#==================================================================================================
# get_weather_hourly.R
#
# Day/night partitioned weather from the NASA POWER *hourly* endpoint.
#
# Sibling collector to get_weather(). Downloads the POWER hourly product,
# partitions each 24 h cycle into day / night / pre-dawn, and returns ONE ROW
# PER ENVIRONMENT-DAY. The 24x data volume stays inside the function.
#
# Why a separate function rather than temporal_api= on get_weather():
#   * POWER caps hourly requests at 15 parameters (daily allows 20); the
#     get_weather() default list has 18 and would fail outright.
#   * T2M_MAX / T2M_MIN / FROST_DAYS are daily aggregates and do not exist
#     hourly. Returning them conditionally would break param_temperature(),
#     param_atmospheric() and processWTH(), which read them by name.
#   * Hourly rate limits are stricter.
#
# LIMITATION: POWER hourly is MERRA-2 reanalysis, not observation. The diurnal
# cycle at 0.5 x 0.625 deg is modelled; sub-daily fidelity is weaker than the
# daily aggregate.
#==================================================================================================


# ============================================================================
# Internal helpers
# ============================================================================

#' POWER fill values are large negative round numbers (-999, -9999).
#' Must be stripped BEFORE aggregation: a single fill destroys a night mean.
#' @keywords internal
#' @noRd
.gwh_strip_na <- function(x) {
  if (!is.numeric(x)) return(x)
  x[x < -900] <- NA_real_
  x
}

#' Mean that returns NA (not NaN) for empty or all-NA input.
#' Needed for polar nights where the window has zero length.
#' @keywords internal
#' @noRd
.gwh_mean <- function(v) {
  v <- v[!is.na(v)]
  if (!length(v)) return(NA_real_)
  mean(v)
}
#' @keywords internal
#' @noRd
.gwh_sd <- function(v) {
  v <- v[!is.na(v)]
  if (length(v) < 2L) return(NA_real_)
  stats::sd(v)
}
#' @keywords internal
#' @noRd
.gwh_min <- function(v) {
  v <- v[!is.na(v)]
  if (!length(v)) return(NA_real_)
  min(v)
}
#' @keywords internal
#' @noRd
.gwh_max <- function(v) {
  v <- v[!is.na(v)]
  if (!length(v)) return(NA_real_)
  max(v)
}

#' Radiation-type POWER parameters. These are ~0 at night by construction, so
#' night / pre-dawn statistics carry no information and are kept for DAYTIME
#' only.
#' @keywords internal
#' @noRd
.gwh_is_radiation <- function(v) {
  grepl("ALLSKY|CLRSKY|_SW_|_LW_|_PAR|SRAD|TOA_", v)
}

#' Build an hour window that wraps correctly around midnight.
#'
#' A naive `from:to` does NOT error for 22:5 -- it silently returns an
#' 18-element DESCENDING sequence (22,21,...,5), i.e. the daytime hours in
#' reverse. That is the dangerous failure mode this guards against.
#' @keywords internal
#' @noRd
.gwh_hour_window <- function(from, to) {
  stopifnot(from %in% 0:23, to %in% 0:23)
  if (from <= to) from:to else c(from:23L, 0L:to)
}

#' Teten / FAO-56 saturation vapour pressure (kPa) for temperature in C.
#' @keywords internal
#' @noRd
.gwh_teten <- function(tc) 0.6108 * exp((17.27 * tc) / (tc + 237.3))

#' Forsythe et al. (1995) daylength, in hours.
#'
#' The sunset hour angle is clamped so polar day returns 24 and polar night
#' returns 0 instead of NaN.
#' @keywords internal
#' @noRd
.gwh_daylength <- function(lat, doy, p = 0.8333) {
  theta <- 0.2163108 + 2 * atan(0.9671396 * tan(0.00860 * (doy - 186)))
  phi   <- asin(0.39795 * cos(theta))
  arg   <- (sin(p * pi / 180) + sin(lat * pi / 180) * sin(phi)) /
           (cos(lat * pi / 180) * cos(phi))
  arg   <- pmax(-1, pmin(1, arg))            # clamp -> polar day/night safe
  24 - (24 / pi) * acos(arg)
}

#' Assign each hour to the night that BEGINS on a given evening.
#'
#' Evening hours (>= 12) keep their own date; morning hours (< 12) are
#' attributed to the previous date. Without this a single night is split
#' across two calendar days and every night mean is computed from two halves
#' of two different nights.
#' @keywords internal
#' @noRd
.gwh_night_id <- function(dates, hours, night.hours) {
  inw <- hours %in% night.hours
  nid <- rep(NA_character_, length(hours))
  eve <- inw & hours >= 12L
  mor <- inw & hours <  12L
  nid[eve] <- as.character(dates[eve])
  nid[mor] <- as.character(dates[mor] - 1L)
  nid
}

#' Retry wrapper. Mirrors .retry() in env_data_collection_src.R.
#' @keywords internal
#' @noRd
.gwh_retry <- function(expr, tries = 3, quiet = FALSE) {
  for (i in seq_len(tries)) {
    out <- try(force(expr), silent = TRUE)
    if (!inherits(out, "try-error")) return(out)
    if (!quiet) message("    attempt ", i, "/", tries, " failed; retrying ...")
    Sys.sleep(min(2^i, 10))
  }
  stop("POWER request failed after ", tries, " attempts: ",
       attr(out, "condition")$message, call. = FALSE)
}


# ============================================================================
# Main function
# ============================================================================

#' @title Day/Night Partitioned Hourly Weather from NASA POWER
#'
#' @description
#' Sibling collector to \code{\link{get_weather}} that downloads the NASA POWER
#' \emph{hourly} endpoint, partitions each 24 h cycle into day / night /
#' pre-dawn, and returns one row per environment-day. It is a PARALLEL producer
#' of a covariate block; it does NOT feed \code{processWTH()} and must not be
#' wired into the daily pipeline by analogy.
#'
#' @param env.id character. Environment identifiers.
#' @param lat,lon numeric. Coordinates, recycled if length 1.
#' @param start.day,end.day character/Date. Collection window, recycled if
#'   length 1.
#' @param variables.names character. POWER parameters. \strong{Maximum 15} on
#'   the hourly endpoint; validated before the request is sent.
#' @param night.def character. How to define night:
#'   \describe{
#'     \item{\code{"solar"}}{Sunset/sunrise from Forsythe daylength. Night
#'       length varies with season and latitude.}
#'     \item{\code{"fixed"}}{User window via \code{night.hours}. Constant and
#'       directly comparable across sites and dates.}
#'     \item{\code{"radiation"}}{Hours with \code{ALLSKY_SFC_SW_DWN <= 0}.
#'       Empirical; costs one parameter slot.}
#'   }
#' @param night.hours integer. Hours 0-23 for \code{night.def = "fixed"}.
#'   Wrap-around (e.g. \code{c(22,23,0,1,2,3,4,5)}) is handled.
#' @param predawn.hours integer. Window for pre-dawn statistics. Default 2:5.
#' @param time_standard character. \code{"LST"} (default) or \code{"UTC"}.
#'   Keep LST: hour 3 is then ~3 a.m. solar time at every site, so pre-dawn
#'   windows are comparable across longitudes.
#' @param summarise boolean. \code{TRUE} (default) returns one row per
#'   environment-day. \code{FALSE} returns raw hourly rows plus \code{period}
#'   and \code{night_id} -- the 24x volume path.
#' @param thresholds named numeric. Night-hour counts above each threshold,
#'   e.g. \code{c(T2M = 22)} gives \code{hours_above_T2M}.
#' @param tries integer. Retry attempts per site.
#' @param save boolean. Write one per-site CSV.
#' @param dir.path character. Output directory when \code{save = TRUE}.
#' @param verbose boolean. Print progress messages.
#'
#' @return
#' A data.frame of class \code{c("weather_hourly", "data.frame")}. With
#' \code{summarise = TRUE}, one row per environment-day, carrying day/night/
#' pre-dawn statistics per variable. \strong{Radiation-type parameters
#' (shortwave, longwave, PAR) are summarised for DAYTIME only}: they are ~0 at
#' night by construction, so night statistics would add no information.
#'
#' @details
#' \strong{POWER hourly is MERRA-2 reanalysis, not observation.} The diurnal
#' cycle at 0.5 x 0.625 degrees is modelled and sub-daily fidelity is weaker
#' than the daily aggregate.
#'
#' A daily \code{T2M_MIN} is not a substitute for a night mean: the two differ
#' materially, and on a \code{Q10 ~ 2} process a few degrees is a large
#' difference in implied rate (Peng et al. 2004).
#'
#' @examples
#' \dontrun{
#' ## Reuse the coordinates of a shipped daily dataset
#' sites <- unique(maizeWTH[, c("env", "LON", "LAT")])
#' h <- get_weather_hourly(env.id = sites$env, lat = sites$LAT, lon = sites$LON,
#'                         start.day = "2016-01-01", end.day = "2016-01-31",
#'                         night.def = "solar")
#'
#' ## Night temperature is not the daily minimum
#' with(h, cor(T2M_night_mean, T2M_day_min, use = "complete.obs"))
#'
#' ## Fixed night window, wrap-around handled correctly
#' h2 <- get_weather_hourly(sites$env, sites$LAT, sites$LON,
#'                          "2016-01-01", "2016-01-31",
#'                          night.def = "fixed", night.hours = c(22:23, 0:5))
#' }
#'
#' @seealso \code{\link{get_weather}}, \code{\link{processWTH}},
#'   \code{\link{summaryWTH}}
#'
#' @references
#' Peng, S. et al. (2004) Rice yields decline with higher night temperature
#' from global warming. \emph{PNAS} 101(27), 9971-9975.
#'
#' Forsythe, W.C. et al. (1995) A model comparison for daylength as a function
#' of latitude and day of the year. \emph{Ecological Modelling} 80, 87-95.
#'
#' Allen, R.G. et al. (1998) Crop evapotranspiration. \emph{FAO Irrigation and
#' Drainage Paper} 56.
#'
#' @export
get_weather_hourly <- function(env.id, lat, lon, start.day, end.day,
                               variables.names = NULL,
                               night.def     = c("solar", "fixed", "radiation"),
                               night.hours   = NULL,
                               predawn.hours = 2:5,
                               time_standard = c("LST", "UTC"),
                               summarise     = TRUE,
                               thresholds    = NULL,
                               tries = 3, save = FALSE, dir.path = NULL,
                               verbose = TRUE) {

  .et_banner("get_weather_hourly", "day/night weather from NASA POWER", verbose)
  .need_pkg("nasapower", "get_weather_hourly")

  night.def     <- match.arg(night.def)
  time_standard <- match.arg(time_standard)

  # ---- validation (mirrors get_weather) ------------------------------------
  env.id <- as.character(env.id)
  n_env  <- length(env.id)
  if (length(lat) == 1L) lat <- rep(lat, n_env)
  if (length(lon) == 1L) lon <- rep(lon, n_env)
  if (length(start.day) == 1L) start.day <- rep(start.day, n_env)
  if (length(end.day)   == 1L) end.day   <- rep(end.day,   n_env)

  if (length(lat) != n_env || length(lon) != n_env)
    stop("'lat' and 'lon' must be length 1 or length(env.id).", call. = FALSE)
  if (length(start.day) != n_env || length(end.day) != n_env)
    stop("'start.day' and 'end.day' must be length 1 or length(env.id).",
         call. = FALSE)
  if (any(as.Date(end.day) < as.Date(start.day)))
    stop("'end.day' precedes 'start.day' for: ",
         paste(env.id[as.Date(end.day) < as.Date(start.day)], collapse = ", "),
         call. = FALSE)
  if (any(lat < -90 | lat > 90, na.rm = TRUE))
    stop("'lat' outside [-90, 90]. Check lat/lon are not swapped.", call. = FALSE)
  if (any(lon < -180 | lon > 180, na.rm = TRUE))
    stop("'lon' outside [-180, 180]. Check lat/lon are not swapped.", call. = FALSE)

  # ---- parameters: enforce the 15-par hourly cap BEFORE the request --------
  if (is.null(variables.names)) {
    variables.names <- c("T2M", "T2MDEW", "RH2M", "QV2M", "PS", "WS2M",
                         "PRECTOTCORR", "ALLSKY_SFC_SW_DWN",
                         "ALLSKY_SFC_LW_DWN", "ALLSKY_SFC_PAR_TOT",
                         "TS", "T2MWET")
  }
  variables.names[grepl("PRECTOT", variables.names)] <- "PRECTOTCORR"
  variables.names <- unique(variables.names)

  if (night.def == "radiation" &&
      !"ALLSKY_SFC_SW_DWN" %in% variables.names) {
    variables.names <- c(variables.names, "ALLSKY_SFC_SW_DWN")
    if (verbose) message("night.def = 'radiation': added ALLSKY_SFC_SW_DWN.")
  }

  if (length(variables.names) > 15L)
    stop("The POWER hourly endpoint accepts at most 15 parameters; ",
         length(variables.names), " supplied:\n  ",
         paste(variables.names, collapse = ", "),
         "\nDrop ", length(variables.names) - 15L, " parameter(s).",
         call. = FALSE)

  if (time_standard == "UTC" && night.def != "radiation")
    warning("time_standard = 'UTC' with night.def = '", night.def,
            "': night windows will NOT be comparable across longitudes. ",
            "Use 'LST' unless you have a specific reason.", call. = FALSE)

  if (night.def == "fixed") {
    if (is.null(night.hours))
      stop("night.def = 'fixed' requires 'night.hours', e.g. c(22:23, 0:5).",
           call. = FALSE)
    if (!all(night.hours %in% 0:23))
      stop("'night.hours' must be integers in [0, 23].", call. = FALSE)
    if (length(night.hours) >= 18L)
      warning("'night.hours' has ", length(night.hours), " hours. If you wrote ",
              "`22:5`, R produced a descending 18-hour DAYTIME sequence. ",
              "Use c(22:23, 0:5) instead.", call. = FALSE)
  }

  if (isTRUE(save)) {
    if (is.null(dir.path)) dir.path <- getwd()
    dir.create(dir.path, showWarnings = FALSE, recursive = TRUE)
  }

  if (verbose) {
    .et_step(paste0("environmental units ......... ", n_env), verbose)
    .et_step(paste0("parameters .................. ", length(variables.names),
                    "/15"), verbose)
    .et_step(paste0("night definition ............ ", night.def), verbose)
    .et_step(paste0("time standard ............... ", time_standard), verbose)
    .et_step(paste0("return ...................... ",
                    if (summarise) "one row per environment-day"
                    else "raw hourly"), verbose)
  }

  # ---- per-site download + partition ---------------------------------------
  out <- vector("list", n_env)

  for (k in seq_len(n_env)) {

    CL <- .gwh_retry(
      as.data.frame(
        nasapower::get_power(community     = "ag",
                             lonlat        = c(lon[k], lat[k]),
                             pars          = variables.names,
                             dates         = c(start.day[k], end.day[k]),
                             temporal_api  = "hourly",
                             time_standard = time_standard)),
      tries = tries, quiet = !verbose)

    names(CL)[names(CL) == "PRECTOTCORR"] <- "PRECTOT"
    vars <- intersect(
      c(sub("PRECTOTCORR", "PRECTOT", variables.names)), names(CL))

    # POWER hourly returns YEAR / MO / DY / HR (verify against a live call;
    # the error below is the guard if the schema ever changes).
    if (!all(c("YEAR", "MO", "DY", "HR") %in% names(CL)))
      stop("Unexpected hourly schema for '", env.id[k],
           "': expected YEAR, MO, DY, HR.", call. = FALSE)

    CL$date <- as.Date(sprintf("%04d-%02d-%02d", CL$YEAR, CL$MO, CL$DY))
    CL$hour <- as.integer(CL$HR)
    CL      <- CL[order(CL$date, CL$hour), , drop = FALSE]

    # strip fill values BEFORE any aggregation
    for (v in vars) CL[[v]] <- .gwh_strip_na(CL[[v]])

    # ---- classify each hour as day or night --------------------------------
    if (night.def == "fixed") {
      CL$is_night <- CL$hour %in% night.hours

    } else if (night.def == "radiation") {
      if (!"ALLSKY_SFC_SW_DWN" %in% names(CL))
        stop("night.def = 'radiation' but ALLSKY_SFC_SW_DWN was not returned.",
             call. = FALSE)
      CL$is_night <- !is.na(CL$ALLSKY_SFC_SW_DWN) & CL$ALLSKY_SFC_SW_DWN <= 0

    } else {  # "solar"
      doy <- as.integer(format(CL$date, "%j"))
      dl  <- .gwh_daylength(lat[k], doy)            # hours of daylight
      sunrise <- 12 - dl / 2                        # centred on solar noon
      sunset  <- 12 + dl / 2
      CL$is_night <- !(CL$hour >= sunrise & CL$hour < sunset)
    }

    # polar diagnostics
    n_by_date <- tapply(CL$is_night, CL$date, sum)
    if (any(n_by_date == 0L | n_by_date == 24L)) {
      bad <- names(n_by_date)[n_by_date == 0L | n_by_date == 24L]
      warning("Environment '", env.id[k], "': ", length(bad),
              " date(s) with a zero-length or 24-hour night (polar day/night). ",
              "Night statistics returned as NA. First: ", bad[1], call. = FALSE)
    }

    # ---- night_id: keep a night on ONE date --------------------------------
    nh_used     <- sort(unique(CL$hour[CL$is_night]))
    CL$night_id <- .gwh_night_id(CL$date, CL$hour, nh_used)
    CL$period   <- ifelse(CL$is_night, "night", "day")

    CL$env <- env.id[k]

    if (!summarise) {
      keep <- c("env", "date", "hour", "period", "night_id", vars)
      out[[k]] <- CL[, keep, drop = FALSE]
      if (verbose)
        .et_step(paste0("[", k, "/", n_env, "] ", env.id[k], " | ",
                        nrow(CL), " hourly rows"), verbose)
      next
    }

    # ---- summarise to one row per environment-day --------------------------
    dates  <- sort(unique(CL$date))
    expect <- length(nh_used)

    rows <- lapply(dates, function(dd) {

      day_rows   <- CL[CL$date == dd & !CL$is_night, , drop = FALSE]
      night_rows <- CL[!is.na(CL$night_id) &
                         CL$night_id == as.character(dd), , drop = FALSE]
      pre_rows   <- night_rows[night_rows$hour %in% predawn.hours, , drop = FALSE]

      base <- data.frame(
        env            = env.id[k],
        date           = dd,
        daysFromStart  = as.numeric(dd - as.Date(start.day[k])) + 1,
        night_hours    = nrow(night_rows),
        day_hours      = nrow(day_rows),
        night_complete = nrow(night_rows) == expect && expect > 0L,
        stringsAsFactors = FALSE
      )

      stats <- list()
      for (v in vars) {
        stats[[paste0(v, "_day_mean")]] <- .gwh_mean(day_rows[[v]])
        stats[[paste0(v, "_day_min")]]  <- .gwh_min(day_rows[[v]])
        stats[[paste0(v, "_day_max")]]  <- .gwh_max(day_rows[[v]])
        stats[[paste0(v, "_day_sd")]]   <- .gwh_sd(day_rows[[v]])
        # Radiation is ~0 at night; night/pre-dawn stats add no information.
        if (!.gwh_is_radiation(v)) {
          stats[[paste0(v, "_night_mean")]]   <- .gwh_mean(night_rows[[v]])
          stats[[paste0(v, "_night_min")]]    <- .gwh_min(night_rows[[v]])
          stats[[paste0(v, "_night_max")]]    <- .gwh_max(night_rows[[v]])
          stats[[paste0(v, "_night_sd")]]     <- .gwh_sd(night_rows[[v]])
          stats[[paste0(v, "_predawn_mean")]] <- .gwh_mean(pre_rows[[v]])
        }
      }

      # --- derived: hourly VPD, then averaged (the correct construction) ----
      if (all(c("T2M", "T2MDEW") %in% vars)) {
        vpd_n <- .gwh_teten(night_rows$T2M) - .gwh_teten(night_rows$T2MDEW)
        vpd_d <- .gwh_teten(day_rows$T2M)   - .gwh_teten(day_rows$T2MDEW)
        vpd_p <- .gwh_teten(pre_rows$T2M)   - .gwh_teten(pre_rows$T2MDEW)
        stats$VPD_night_mean   <- .gwh_mean(pmax(vpd_n, 0))
        stats$VPD_day_mean     <- .gwh_mean(pmax(vpd_d, 0))
        stats$VPD_predawn_mean <- .gwh_mean(pmax(vpd_p, 0))
      }

      # --- diurnal temperature range from hourly extremes -------------------
      if ("T2M" %in% vars) {
        allday <- CL[CL$date == dd, , drop = FALSE]
        stats$DTR <- .gwh_max(allday$T2M) - .gwh_min(allday$T2M)
      }

      # --- threshold hour counts (duration, not peak) -----------------------
      if (!is.null(thresholds)) {
        for (tv in names(thresholds)) {
          if (!tv %in% vars) next
          val <- night_rows[[tv]]
          stats[[paste0("hours_above_", tv)]] <-
            if (!length(val)) NA_integer_
            else sum(val > thresholds[[tv]], na.rm = TRUE)
        }
      }

      cbind(base, as.data.frame(stats, stringsAsFactors = FALSE))
    })

    res <- do.call(rbind, rows)

    if (isTRUE(save))
      utils::write.csv(res,
                       file.path(dir.path, paste0(env.id[k], "_hourly.csv")),
                       row.names = FALSE)

    out[[k]] <- res

    if (verbose)
      .et_step(paste0("[", k, "/", n_env, "] ", env.id[k], " | ", nrow(res),
                      " days | night ", min(res$night_hours), "-",
                      max(res$night_hours), " h | complete ",
                      sum(res$night_complete), "/", nrow(res)), verbose)
  }

  final <- do.call(rbind, out)
  rownames(final) <- NULL
  class(final) <- c("weather_hourly", "data.frame")
  attr(final, "night.def")     <- night.def
  attr(final, "time_standard") <- time_standard
  attr(final, "predawn.hours") <- predawn.hours
  attr(final, "variables")     <- variables.names
  final
}
