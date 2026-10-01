#==================================================================================================
# env_data_collection.R
#
# Data-collection and derived-parameter tools for the 'et' package.
#
# Consolidated from: get_weather.R, get_soil.R, get_spatial.R, extract_GIS.R,
#                    param_radiation.R, param_atmospheric.R, param_temperature.R,
#                    processWTH.R, summaryWTH.R
#
# Crash-resilient retrieval (Section 10):
#   get_weather_resumable() - checkpointed get_weather()
#   get_soil_resumable()    - checkpointed get_soil()
#   restart_from_log()      - resume an interrupted run
#   read_progress_log()     - inspect what succeeded / failed
#
# Author: G Costa-Neto and contributors; refactor 2026
#==================================================================================================


# =========================================================================
# SECTION 0 - internal helpers (not exported)
# =========================================================================

#' Degrees to radians
#' @keywords internal
#' @noRd
.deg2rad <- function(deg) deg * pi / 180

#' Abort when an optional package is absent
#'
#' Replaces the previous `install.packages()` calls, which modified the user's
#' library without consent (a CRAN policy violation).
#'
#' @keywords internal
#' @noRd
.need_pkg <- function(pkg, fun) {
  missing <- pkg[!vapply(pkg, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing))
    stop(fun, "() requires the package(s): ", paste(missing, collapse = ", "),
         ".\n  Install with: install.packages(c(",
         paste0('"', missing, '"', collapse = ", "), "))",
         call. = FALSE)
  invisible(TRUE)
}

#' NASA POWER missing-value sentinels
#'
#' POWER encodes missing data with several magic numbers depending on the
#' parameter and the product version. Confirmed in the wild:
#' \itemize{
#'   \item \code{-999} - the standard fill for most parameters
#'   \item \code{-99}  - used by some radiation parameters
#'   \item \code{-9999}, \code{-99999} - seen in some gridded products
#'   \item \code{-9000} - observed in \code{FROST_DAYS}
#' }
#'
#' @keywords internal
#' @noRd
.POWER_NA <- c(-99, -999, -9000, -9999, -99999)

#' Lower bound below which any value is treated as a fill code
#'
#' Enumerating sentinels alone is brittle: \code{FROST_DAYS = -9000} slipped
#' through a \code{c(-999, -99)} list and produced a mean of exactly -9000.
#' Every POWER fill is a large-magnitude negative round number, while the most
#' extreme physically plausible value in this data (Antarctic \code{T2M_MIN},
#' about -90 C) is an order of magnitude smaller. A threshold of -900 therefore
#' catches unknown future fills with roughly a 9x safety margin.
#'
#' @keywords internal
#' @noRd
.POWER_NA_FLOOR <- -900

#' Replace NASA POWER fill values with NA
#'
#' Applies two complementary rules: exact matches against known sentinels (which
#' catches \code{-99}, a value above the floor), and a magnitude threshold (which
#' catches \code{-9000} and any fill code POWER introduces later).
#'
#' @param x numeric vector.
#' @param sentinels numeric. Exact fill values to remove.
#' @param floor numeric. Values at or below this are treated as fills.
#' @return \code{x} with fill values replaced by \code{NA}.
#' @keywords internal
#' @noRd
.strip_power_na <- function(x, sentinels = .POWER_NA, floor = .POWER_NA_FLOOR) {
  if (!is.numeric(x)) return(x)
  x[!is.na(x) & (x %in% sentinels | x <= floor)] <- NA
  x
}

#' Warn about columns that are entirely missing after fill-stripping
#'
#' A parameter can be unavailable at a given site (\code{FROST_DAYS} in the
#' tropics, snow variables at low latitude). Previously these surfaced as
#' plausible-looking numbers; now they become all-\code{NA}, and the user is told
#' rather than left to discover it in a summary table.
#'
#' @keywords internal
#' @noRd
.warn_empty_cols <- function(df, cols, verbose = TRUE) {
  if (!length(cols)) return(invisible(character(0)))
  cols <- intersect(cols, names(df))
  empty <- cols[vapply(df[cols], function(z) all(is.na(z)), logical(1))]
  if (length(empty) && verbose) {
    message("Note: no valid data returned for: ", paste(empty, collapse = ", "))
    message("      These columns are entirely NA (parameter unavailable at this ",
            "site/period).")
  }
  invisible(empty)
}

#' Extraterrestrial radiation and daylength
#'
#' FAO-56 Eqs. 21-25 for extraterrestrial radiation (Ra), and Forsythe et al.
#' (1995) for daylength (N).
#'
#' The sunset hour angle \eqn{\omega_s = \arccos(-\tan\varphi\tan\delta)} is
#' undefined during polar day/night; the argument is clamped to [-1, 1] so Ra
#' tends to 0 (polar night) or a full-day value (polar day) instead of NaN.
#'
#' Two different daylength formulas previously coexisted in the package
#' (`(24/pi)*ws` in param_radiation, Forsythe in get_weather), disagreeing by up
#' to ~0.5 h at 60 deg. Forsythe is used throughout as it accounts for
#' atmospheric refraction.
#'
#' @param J integer day of year.
#' @param lat numeric latitude in decimal degrees.
#' @return data.frame with Ra (MJ/m2/day) and N (hours).
#' @keywords internal
#' @noRd
.ra_daylength <- function(J, lat) {
  rlat  <- .deg2rad(lat)
  delta <- 0.409 * sin(2 * pi / 365 * J - 1.39)
  dr    <- 1 + 0.033 * cos(2 * pi / 365 * J)

  x  <- pmin(pmax(-tan(rlat) * tan(delta), -1), 1)
  ws <- acos(x)
  Ra <- (1440 / pi) * 0.0820 * dr *
    (ws * sin(rlat) * sin(delta) + cos(rlat) * cos(delta) * sin(ws))
  Ra <- pmax(Ra, 0)

  P <- asin(0.39795 * cos(0.2163108 + 2 * atan(0.9671396 * tan(0.00860 * (J - 186)))))
  z <- (sin(0.8333 * pi / 180) + sin(.deg2rad(lat)) * sin(P)) /
    (cos(.deg2rad(lat)) * cos(P))
  z  <- pmin(pmax(z, -1), 1)
  N  <- 24 - (24 / pi) * acos(z)

  data.frame(Ra = Ra, N = N)
}

#' Tetens saturation vapour pressure (kPa)
#' @keywords internal
#' @noRd
.teten <- function(Temp) 0.61078 * exp((17.27 * Temp) / (Temp + 237.3))

#' Atmospheric pressure from elevation (kPa), FAO-56 Eq. 7
#' @keywords internal
#' @noRd
.atm_pressure <- function(elevation) 101.3 * ((293 - 0.0065 * elevation) / 293)^5.26

#' Retry an HTTP-bound expression with exponential backoff
#'
#' Both remote endpoints (NASA POWER, SoilGrids) return transient 5xx errors under
#' load. Previously a single failure aborted the whole download.
#'
#' @keywords internal
#' @noRd
.retry <- function(expr, tries = 3L, wait = 5, quiet = FALSE) {
  for (k in seq_len(tries)) {
    out <- try(force(expr), silent = TRUE)
    if (!inherits(out, "try-error")) return(out)
    if (k == tries) stop(attr(out, "condition")$message, call. = FALSE)
    if (!quiet)
      message("  request failed (attempt ", k, "/", tries, "); retrying in ",
              wait * k, "s ...")
    Sys.sleep(wait * k)
  }
}

#' Split a vector into chunks of a given length
#' @keywords internal
#' @noRd
.split_chunk <- function(vec, length) split(vec, ceiling(seq_along(vec) / length))

#' Row-bind a list of data.frames, filling absent columns with NA
#'
#' Base \code{rbind()} aborts with \dQuote{names do not match previous names}
#' whenever the elements differ in their columns. That happens in practice: NASA
#' POWER occasionally omits a parameter for one site but not another, and a
#' partially failed download leaves \code{NULL} entries in the results list.
#'
#' This is a dependency-free replacement for \code{plyr::ldply()} and the
#' \code{rbind.fill()} behaviour it relies on. Column order follows first
#' appearance across the list; \code{NULL} and zero-row elements are dropped;
#' factors are demoted to character so differing level sets cannot produce
#' mismatched integer codes.
#'
#' @param x list of data.frames.
#' @return A single data.frame, or an empty data.frame when nothing is bindable.
#' @keywords internal
#' @noRd
.rbind_fill <- function(x) {
  if (!length(x)) return(data.frame())
  x <- Filter(function(d) !is.null(d) && NROW(d) > 0L, lapply(x, as.data.frame))
  if (!length(x)) return(data.frame())

  cols <- unique(unlist(lapply(x, names), use.names = FALSE))

  x <- lapply(x, function(d) {
    fac <- vapply(d, is.factor, logical(1))
    if (any(fac)) d[fac] <- lapply(d[fac], as.character)
    d
  })

  filled <- lapply(x, function(d) {
    miss <- setdiff(cols, names(d))
    if (length(miss)) d[miss] <- NA
    d[cols]
  })

  out <- do.call(rbind, filled)
  rownames(out) <- NULL
  out
}

#' Reshape a long data.frame to wide format
#'
#' Dependency-free replacement for \code{reshape2::dcast()} covering the single
#' case the package needs: one value column spread across a key column, keyed by
#' one or more id columns.
#'
#' Absent id/key combinations become \code{fill} (\code{NA} by default), matching
#' \code{dcast}. This matters for SoilGrids, which can return a property for one
#' site but not another. New columns are sorted by name; id columns retain their
#' original type, so numeric coordinates are not coerced to character.
#'
#' Unlike \code{dcast}, duplicate id/key pairs are not silently aggregated -- the
#' last value wins and a warning is raised, since duplicates indicate a repeated
#' query rather than data needing summarisation.
#'
#' @param data data.frame in long format.
#' @param id.cols character. Columns identifying each output row.
#' @param key.col character. Column whose values become new column names.
#' @param value.col character. Column supplying the cell values.
#' @param fill value used where an id/key combination is absent.
#' @return A wide data.frame.
#' @keywords internal
#' @noRd
.long_to_wide <- function(data, id.cols, key.col, value.col, fill = NA_real_) {
  id.cols <- intersect(id.cols, names(data))
  if (!length(id.cols))
    stop("None of the id columns are present in the data.", call. = FALSE)
  if (!key.col %in% names(data))
    stop("Key column '", key.col, "' not found.", call. = FALSE)
  if (!value.col %in% names(data))
    stop("Value column '", value.col, "' not found.", call. = FALSE)

  idf  <- data[, id.cols, drop = FALSE]
  key  <- do.call(paste, c(lapply(idf, as.character), sep = "\r"))
  ukey <- unique(key)

  base <- idf[match(ukey, key), , drop = FALSE]
  rownames(base) <- NULL

  kv   <- as.character(data[[key.col]])
  cols <- sort(unique(kv))

  ri <- match(key, ukey)
  ci <- match(kv, cols)

  mat <- matrix(fill, nrow = length(ukey), ncol = length(cols),
                dimnames = list(NULL, cols))
  if (any(duplicated(cbind(ri, ci))))
    warning(sum(duplicated(cbind(ri, ci))), " duplicate id/key combination(s) ",
            "in the wide reshape; the last value was kept.", call. = FALSE)
  mat[cbind(ri, ci)] <- data[[value.col]]

  cbind(base, as.data.frame(mat, stringsAsFactors = FALSE))
}

# =========================================================================
# SECTION 1 - get_weather()
# =========================================================================

#' @title Easily Collect Worldwide Daily Weather Data
#'
#' @description
#' Imports daily weather data from the NASA POWER API and augments it with derived
#' agro-meteorological indices (VPD, photoperiod, extraterrestrial radiation,
#' temperature-humidity indices).
#'
#' @author Germano Costa Neto and Giovanni Galli, modified by Tiago Olivoto
#'
#' @param env.id vector (character). Identifier of the site/environment. If
#'   \code{NULL}, environments are named \code{env1 ... envN}.
#' @param lat vector (numeric). Latitude in decimal degrees (WGS84).
#' @param lon vector (numeric). Longitude in decimal degrees (WGS84).
#' @param start.day,end.day vector (character or Date). First and last date of
#'   collection, e.g. \code{"2015-02-15"}. Recycled if length 1.
#' @param variables.names vector (character). POWER parameters to request. See Details.
#' @param dir.path character. Output directory when \code{save = TRUE}. Defaults to
#'   \code{getwd()}.
#' @param save boolean. Write one .csv per environment.
#' @param parallel boolean. Download chunks in parallel. Default \code{FALSE}
#'   (serial), which is gentler on the NASA POWER rate limit; set \code{TRUE}
#'   to speed up large jobs.
#' @param workers integer. Number of parallel processes. Defaults to 90\% of
#'   available cores, capped at \code{chunk_size}.
#' @param chunk_size integer. Points per chunk. Default 29; raising this may exceed
#'   the API rate limit.
#' @param sleep numeric. Seconds to pause between chunks. Default 60.
#' @param tries integer. Retry attempts per failed request. Default 3.
#' @param verbose boolean. Print progress messages.
#'
#' @return
#' A data.frame with one row per environment-day. Always contains \code{env},
#' \code{LON}, \code{LAT}, \code{YYYYMMDD}, \code{DOY} and \code{daysFromStart},
#' followed by the requested and derived variables.
#'
#' @details
#' \strong{daysFromStart is date-based.} It is computed as
#' \code{as.numeric(YYYYMMDD - start.day) + 1}, not as the row position. The previous
#' implementation used \code{1:nrow()}, which silently mis-numbered every row after a
#' gap whenever the API returned fewer days than requested -- assigning observations
#' to the wrong phenological stage in \code{summaryWTH}.
#'
#' \strong{Missing values.} POWER fill values (-999 and -99) are converted to
#' \code{NA} before any derived variable is computed.
#'
#' Commonly used parameters include \code{T2M}, \code{T2M_MAX}, \code{T2M_MIN},
#' \code{T2MDEW}, \code{PRECTOTCORR} (returned as \code{PRECTOT}), \code{RH2M},
#' \code{WS2M}, \code{ALLSKY_SFC_SW_DWN}, \code{ALLSKY_SFC_SW_DNI},
#' \code{ALLSKY_SFC_PAR_TOT}, \code{GWETROOT}, \code{GWETTOP}, \code{EVPTRNS}.
#'
#' Derived variables, added when their inputs are present:
#' \itemize{
#'   \item \code{VPD}: vapour pressure deficit (kPa), from T2MDEW, T2M_MAX, T2M_MIN
#'   \item \code{N}: photoperiod (h), Forsythe et al. (1995)
#'   \item \code{RTA}: extraterrestrial radiation (MJ/m2/day), FAO-56
#'   \item \code{n}: actual sunshine duration (h)
#'   \item \code{TH1}, \code{TH2}: temperature-humidity indices (NRC 1971; Yousef 1985)
#'   \item \code{PAR_TEMP}: ratio of PAR to mean temperature
#' }
#'
#' @examples
#' \dontrun{
#' ## Single location
#' get_weather(env.id = "NM", lat = -13.05, lon = -56.05,
#'             start.day = "2015-02-15", end.day = "2015-06-15",
#'             variables.names = c("T2M", "PRECTOTCORR"))
#'
#' ## Two locations with different planting dates
#' get_weather(env.id = c("NM", "SO"),
#'             lat = c(-13.05, -12.32), lon = c(-56.05, -55.42),
#'             start.day = c("2015-02-15", "2015-02-13"),
#'             end.day   = rep("2015-06-15", 2))
#' }
#'
#' @seealso \code{summaryWTH}, \code{processWTH}, \code{get_soil},
#'   \code{\link{get_weather_hourly}}
#'
#' @references
#' Sparks, A. (2018). nasapower: NASA-POWER Data from R. \emph{JOSS} 3(30), 1035.
#'
#' Forsythe, W.C., Rykiel, E.J., Stahl, R.S., Wu, H., Schoolfield, R.M. (1995).
#' A model comparison for daylength as a function of latitude and day of the year.
#' \emph{Ecological Modelling} 80, 87-95.
#'
#' @importFrom utils write.csv
#' @export
get_weather <- function(env.id = NULL, lat = NULL, lon = NULL,
                        start.day = NULL, end.day = NULL,
                        variables.names = NULL, dir.path = NULL, save = FALSE,
                        parallel = FALSE, workers = NULL, chunk_size = 29,
                        sleep = 60, tries = 3L, verbose = TRUE) {

  .et_banner("get_weather", "downloads daily weather from NASA POWER", verbose)
  .need_pkg("nasapower", "get_weather")
  if (isTRUE(parallel)) .need_pkg("parallel", "get_weather")

  if (is.null(lat) || is.null(lon))
    stop("'lat' and 'lon' are required.", call. = FALSE)
  if (length(lat) != length(lon))
    stop("'lat' and 'lon' must have the same length.", call. = FALSE)

  n_env <- length(lat)
  if (is.null(env.id)) env.id <- paste0("env", seq_len(n_env))
  # Keep as character: factors leak integer codes through split() and binding.
  env.id <- as.character(env.id)
  if (length(env.id) != n_env)
    stop("'env.id' must have the same length as 'lat' and 'lon'.", call. = FALSE)
  if (anyDuplicated(env.id))
    warning("Duplicated env.id values; downstream joins may be ambiguous.", call. = FALSE)

  if (is.null(dir.path)) dir.path <- getwd()

  if (is.null(start.day)) {
    start.day <- Sys.Date() - 1000
    if (verbose) message("start.day is NULL; matched as ", start.day)
  }
  start.day <- as.Date(start.day)
  if (is.null(end.day)) {
    end.day <- start.day + 30
    if (verbose) message("end.day is NULL; matched as ", paste(unique(end.day), collapse = ", "))
  }
  end.day <- as.Date(end.day)

  if (length(start.day) == 1L) start.day <- rep(start.day, n_env)
  if (length(end.day)   == 1L) end.day   <- rep(end.day,   n_env)
  if (length(start.day) != n_env || length(end.day) != n_env)
    stop("'start.day' and 'end.day' must be length 1 or length(lat).", call. = FALSE)
  if (any(end.day < start.day))
    stop("'end.day' precedes 'start.day' for: ",
         paste(env.id[end.day < start.day], collapse = ", "), call. = FALSE)

  if (any(lat < -90 | lat > 90, na.rm = TRUE))
    stop("'lat' outside [-90, 90]. Check lat/lon are not swapped.", call. = FALSE)
  if (any(lon < -180 | lon > 180, na.rm = TRUE))
    stop("'lon' outside [-180, 180]. Check lat/lon are not swapped.", call. = FALSE)

  if (is.null(variables.names)) {
    variables.names <- c("T2M", "T2M_MAX", "T2M_MIN", "T2MDEW",
                         "ALLSKY_SFC_LW_DWN", "ALLSKY_SFC_SW_DWN",
                         "ALLSKY_SFC_SW_DNI", "ALLSKY_SFC_PAR_TOT",
                         "ALLSKY_SFC_UVA", "ALLSKY_SFC_UVB",
                         "PRECTOTCORR", "EVPTRNS", "QV2M", "RH2M",
                         "GWETROOT", "GWETTOP", "FROST_DAYS", "WS2M")
  }
  # POWER serves corrected precipitation as PRECTOTCORR
  variables.names[grepl("PRECTOT", variables.names)] <- "PRECTOTCORR"
  variables.names <- unique(variables.names)

  # ---- per-site download ---------------------------------------------------
  get_helper <- function(lon, lat, variables.names, start.day, end.day,
                         env.id, save, dir.path, tries, verbose) {
    CL <- .retry(
      data.frame(nasapower::get_power(community = "ag",
                                      lonlat = c(lon, lat),
                                      pars = variables.names,
                                      dates = c(start.day, end.day),
                                      temporal_api = "daily")),
      tries = tries, quiet = !verbose)

    names(CL)[names(CL) == "PRECTOTCORR"] <- "PRECTOT"

    # Date-based, NOT positional: robust to gaps in the returned series.
    if ("YYYYMMDD" %in% names(CL)) {
      CL$daysFromStart <- as.numeric(as.Date(CL$YYYYMMDD) - as.Date(start.day)) + 1
    } else {
      CL$daysFromStart <- seq_len(nrow(CL))
      warning("No YYYYMMDD column returned for '", env.id,
              "'; daysFromStart fell back to row position.", call. = FALSE)
    }

    CL$env <- env.id
    CL <- CL[, c("env", setdiff(names(CL), "env")), drop = FALSE]

    if (isTRUE(save))
      utils::write.csv(CL, file = file.path(dir.path, paste0(env.id, ".csv")),
                       row.names = FALSE)
    CL
  }

  startTime <- Sys.time()
  if (verbose) {
    message("---------------------------------------------------------------")
    message("get_weather() - daily weather from NASA POWER")
    message("https://docs.ropensci.org/nasapower")
    message("---------------------------------------------------------------")
    message("Start .............................. ", format(startTime, "%a %b %d %X %Y"))
    message("Environmental units ................ ", n_env)
    message("Parallel ........................... [ ", if (parallel) "x" else " ", " ]")
  }

  results <- vector("list", n_env)

  if (isFALSE(parallel)) {
    init_time <- Sys.time(); iter <- 0L
    for (i in seq_len(n_env)) {
      iter <- iter + 1L
      if (iter >= 30L &&
          as.numeric(difftime(Sys.time(), init_time, units = "secs")) > 60) {
        if (verbose) message("Waiting ", sleep, "s for a new query to the API.")
        Sys.sleep(sleep); iter <- 0L; init_time <- Sys.time()
      }
      results[[i]] <- get_helper(lon[i], lat[i], variables.names,
                                 start.day[i], end.day[i], env.id[i],
                                 save, dir.path, tries, verbose)
      if (verbose) message("  [", i, "/", n_env, "] ", env.id[i], " downloaded")
    }
  } else {
    idx_chunks <- .split_chunk(seq_len(n_env), chunk_size)
    nworkers <- if (is.null(workers)) max(1L, trunc(parallel::detectCores() * 0.9)) else workers
    nworkers <- max(1L, min(nworkers, chunk_size))
    clust <- parallel::makeCluster(nworkers)
    on.exit(parallel::stopCluster(clust), add = TRUE)
    if (verbose) message("Threads ............................ ", nworkers)

    out <- vector("list", length(idx_chunks))
    for (i in seq_along(idx_chunks)) {
      ii <- idx_chunks[[i]]
      lat_t <- lat[ii]; lon_t <- lon[ii]; env_t <- env.id[ii]
      sd_t <- start.day[ii]; ed_t <- end.day[ii]
      parallel::clusterExport(
        clust,
        varlist = c("get_helper", "lat_t", "lon_t", "env_t", "sd_t", "ed_t",
                    "variables.names", "save", "dir.path", "tries", "verbose",
                    ".retry"),
        envir = environment())
      out[[i]] <- parallel::parLapply(clust, seq_along(ii), function(j) {
        get_helper(lon_t[j], lat_t[j], variables.names, sd_t[j], ed_t[j],
                   env_t[j], save, dir.path, tries, verbose)
      })
      if (verbose) message("  chunk ", i, "/", length(idx_chunks),
                           " (", length(ii), " points) downloaded")
      if (i < length(idx_chunks)) {
        if (verbose) message("Waiting ", sleep, "s for a new query to the API.")
        Sys.sleep(sleep)
      }
    }
    results <- unlist(out, recursive = FALSE)
  }

  .out <- .rbind_fill(results)
  if (!nrow(.out)) stop("No data returned by the API.", call. = FALSE)

  # ---- strip fill values BEFORE deriving anything --------------------------
  meta_cols <- c("env", "LON", "LAT", "YEAR", "MM", "DD", "DOY",
                 "YYYYMMDD", "daysFromStart")
  num_cols  <- setdiff(names(.out)[vapply(.out, is.numeric, logical(1))], meta_cols)
  for (cc in num_cols) .out[[cc]] <- .strip_power_na(.out[[cc]])
  .warn_empty_cols(.out, num_cols, verbose = verbose)

  variables.names[variables.names == "PRECTOTCORR"] <- "PRECTOT"
  have <- function(...) all(c(...) %in% names(.out))
  tick <- function(lbl, ok)
    if (verbose) message(lbl, " [ ", if (ok) "x" else " ", " ]")

  ok <- have("T2MDEW", "T2M_MIN", "T2M_MAX")
  if (ok) { .out$VPD <- (.teten(.out$T2M_MIN) + .teten(.out$T2M_MAX)) / 2 -
    .teten(.out$T2MDEW)
  variables.names <- c(variables.names, "VPD") }
  tick("Computing VPD ......................", ok)

  ok <- have("DOY", "LAT")
  if (ok) {
    rd <- .ra_daylength(.out$DOY, .out$LAT)
    .out$N   <- round(rd$N, 3)
    .out$RTA <- round(rd$Ra, 3)
    variables.names <- c(variables.names, "N", "RTA")
  }
  tick("Computing N, RTA ...................", ok)

  ok <- have("RTA", "ALLSKY_SFC_SW_DNI")
  if (ok) {
    r <- .out$ALLSKY_SFC_SW_DNI / .out$RTA
    r[!is.finite(r)] <- NA
    .out$n <- round(pmin(pmax(.out$N * r, 0), .out$N), 3)
    variables.names <- c(variables.names, "n")
  }
  tick("Computing n ........................", ok)

  ok <- have("T2M", "RH2M")
  if (ok) { .out$TH1 <- (1.8 * .out$T2M + 32) -
    (0.55 - 0.0055 * .out$RH2M) * (1.8 * .out$T2M - 26)
  variables.names <- c(variables.names, "TH1") }
  tick("Computing TH1 ......................", ok)

  ok <- have("T2M", "T2MDEW")
  if (ok) { .out$TH2 <- .out$T2M + 0.36 * .out$T2MDEW + 41.2
  variables.names <- c(variables.names, "TH2") }
  tick("Computing TH2 ......................", ok)

  ok <- have("ALLSKY_SFC_PAR_TOT", "T2M")
  if (ok) { pt <- .out$ALLSKY_SFC_PAR_TOT / .out$T2M
  pt[!is.finite(pt)] <- NA
  .out$PAR_TEMP <- pt
  variables.names <- c(variables.names, "PAR_TEMP") }
  tick("Computing PAR_TEMP .................", ok)

  keep <- c(intersect(meta_cols, names(.out)),
            intersect(unique(variables.names), names(.out)))
  .out <- .out[, keep, drop = FALSE]

  endTime <- Sys.time()
  if (verbose) {
    message("Environmental units downloaded ..... ",
            length(unique(.out$env)), "/", n_env)
    message("Daily weather features ............. ",
            length(intersect(unique(variables.names), names(.out))))
    message("Done ............................... ", format(endTime, "%a %b %d %X %Y"))
    message("Total time ......................... ",
            round(difftime(endTime, startTime, units = "secs"), 2), "s")
  }
  .out
}


# =========================================================================
# SECTION 2 - get_spatial()
# =========================================================================

#' @title Extract Point Estimates from Raster Files
#'
#' @description
#' Extracts point values from a raster at the coordinates of each environment. This
#' function supersedes \code{extract_GIS()}, which relied on the retired \pkg{raster}
#' and \pkg{sp} packages.
#'
#' @author Germano Costa Neto
#'
#' @param digital.raster SpatRaster, or a file path readable by \code{terra::rast()}.
#' @param which.raster.number integer. Layer(s) to keep when the input is multi-layer.
#' @param env.dataframe data.frame containing the environment id and coordinates.
#' @param env.id character. Name of the environment-id column in \code{env.dataframe}.
#' @param lat,lng character. Names of the latitude/longitude columns, or numeric
#'   vectors of coordinates used directly.
#' @param name.feature character. Name(s) for the extracted layer(s).
#' @param merge boolean. Join the extracted values back onto \code{env.dataframe}.
#' @param crs numeric or character. CRS of the input coordinates. Default 4326 (WGS84).
#' @param verbose boolean. Print progress messages.
#'
#' @return
#' A data.frame of extracted values, joined to \code{env.dataframe} when
#' \code{merge = TRUE}.
#'
#' @details
#' Coordinates are validated before extraction: latitude must lie in [-90, 90] and
#' longitude in [-180, 180]. This catches the common lat/lon transposition, which
#' otherwise yields silent \code{NA}s or values from the wrong hemisphere.
#'
#' @examples
#' \dontrun{
#' ref <- data.frame(env = "NM", lat = -13.05, lng = -56.05)
#' get_spatial(env.dataframe = ref, digital.raster = my_elevation_raster,
#'             lat = "lat", lng = "lng", env.id = "env",
#'             name.feature = "Elevation")
#' }
#'
#' @seealso \code{get_weather}, \code{get_soil}
#'
#' @export
get_spatial <- function(digital.raster = NULL, which.raster.number = NULL,
                        env.dataframe = NULL, env.id = NULL,
                        lat = NULL, lng = NULL, name.feature = NULL,
                        merge = TRUE, crs = 4326, verbose = TRUE) {

  .et_banner("get_spatial", "extracts point values from raster files", verbose)
  .need_pkg(c("terra", "sf"), "get_spatial")

  if (is.null(digital.raster)) stop("Please supply 'digital.raster'.", call. = FALSE)
  if (is.null(lat) || is.null(lng))
    stop("Provide 'lat' and 'lng' (column names in env.dataframe, or numeric vectors).",
         call. = FALSE)

  startTime <- Sys.time()
  if (verbose) {
    message("Based on sf::st_as_sf() and terra::extract()")
  }

  if (!inherits(digital.raster, "SpatRaster"))
    digital.raster <- terra::rast(digital.raster)
  if (!is.null(which.raster.number))
    digital.raster <- digital.raster[[which.raster.number]]
  if (!is.null(name.feature)) {
    if (length(name.feature) != terra::nlyr(digital.raster))
      stop("'name.feature' has length ", length(name.feature),
           " but the raster has ", terra::nlyr(digital.raster), " layer(s).",
           call. = FALSE)
    names(digital.raster) <- name.feature
  }
  name.feature <- names(digital.raster)

  if (is.numeric(lat) && is.numeric(lng)) {
    .coords <- data.frame(x = lng, y = lat)
    ids <- if (!is.null(env.id) && length(env.id) == nrow(.coords)) env.id else
      paste0("env", seq_len(nrow(.coords)))
  } else {
    if (is.null(env.dataframe))
      stop("'env.dataframe' is required when 'lat'/'lng' are column names.", call. = FALSE)
    miss <- setdiff(c(lat, lng, env.id), names(env.dataframe))
    if (length(miss))
      stop("Column(s) not found in env.dataframe: ", paste(miss, collapse = ", "),
           call. = FALSE)
    .coords <- data.frame(x = env.dataframe[[lng]], y = env.dataframe[[lat]])
    ids <- if (is.null(env.id)) paste0("env", seq_len(nrow(.coords))) else
      env.dataframe[[env.id]]
  }

  if (anyNA(.coords))
    stop("Coordinates contain NA values.", call. = FALSE)
  if (any(.coords$y < -90 | .coords$y > 90))
    stop("Latitude outside [-90, 90]. Are 'lat' and 'lng' swapped?", call. = FALSE)
  if (any(.coords$x < -180 | .coords$x > 180))
    stop("Longitude outside [-180, 180]. Are 'lat' and 'lng' swapped?", call. = FALSE)

  if (verbose) message("Environmental units ................ ", length(unique(ids)))

  sp_vector <- sf::st_as_sf(.coords, coords = c("x", "y"), crs = crs)
  extracted <- terra::extract(digital.raster, sp_vector)
  extracted$ID <- NULL

  n_na <- sum(!stats::complete.cases(extracted))
  if (n_na && verbose)
    message("Note: ", n_na, " point(s) fell outside the raster extent (NA).")

  id_col <- if (is.null(env.id)) "env" else env.id
  out <- cbind(stats::setNames(data.frame(ids, stringsAsFactors = FALSE), id_col),
               extracted)

  if (isTRUE(merge) && !is.null(env.dataframe) && !is.null(env.id))
    out <- merge(env.dataframe, out, by = env.id)

  endTime <- Sys.time()
  if (verbose) {
    message("Raster features .................... ", length(name.feature))
    message("Done ............................... ", format(endTime, "%a %b %d %X %Y"))
    message("Total time ......................... ",
            round(difftime(endTime, startTime, units = "secs"), 2), "s")
  }
  out
}


# =========================================================================
# SECTION 3 - get_soil()
# =========================================================================

#' @title Collect Soil Features from SoilGrids (250 m)
#'
#' @description
#' Queries the ISRIC SoilGrids REST API for soil properties at given coordinates.
#'
#' @author Brandon Monier, modified by Germano Costa Neto
#'
#' @param env.id vector (character). Environment identifiers.
#' @param lat,lon vector (numeric). Coordinates in decimal degrees (WGS84).
#' @param variables.names vector (character). SoilGrids properties. See Details.
#' @param depths vector (character). Depth intervals to request.
#' @param stat character. Which summary to return: \code{"mean"} (default),
#'   \code{"q_5"}, \code{"q_50"}, \code{"q_95"} or \code{"uncertainty"}.
#' @param wide boolean. Return one row per environment with one column per
#'   property-depth (default), else long format.
#' @param sleep numeric. Minimum seconds between successive API requests.
#'   Minimum 12 (SoilGrids fair-use limit of 5 requests/min). Network time counts
#'   toward this interval, and duplicate coordinates are fetched only once.
#' @param tries integer. Retry attempts per failed request.
#' @param timeout numeric. Maximum seconds to wait for a single SoilGrids
#'   response before the request is aborted and retried. Prevents one stalled
#'   connection from blocking the whole sequential download. Default 60.
#' @param verbose boolean. Print progress messages.
#'
#' @return A data.frame of soil features.
#'
#' @details
#' Accepted properties: \code{bdod}, \code{cec}, \code{cfvo}, \code{clay},
#' \code{nitrogen}, \code{ocd}, \code{ocs}, \code{phh2o}, \code{sand},
#' \code{silt}, \code{soc}, \code{wv0010}, \code{wv0033}, \code{wv1500}.
#'
#' SoilGrids returns values in mapped units (e.g. clay in g/kg x 10). No unit
#' conversion is applied; the \code{unit} column records what was served.
#'
#' @examples
#' \dontrun{
#' get_soil(env.id = "NM", lat = -13.05, lon = -56.05,
#'          variables.names = c("clay", "nitrogen"))
#' }
#'
#' @seealso \code{get_weather}, \code{get_spatial}
#'
#' @references
#' Poggio, L. et al. (2021). SoilGrids 2.0: producing soil information for the globe
#' with quantified spatial uncertainty. \emph{SOIL} 7, 217-240.
#'
#' @export
get_soil <- function(env.id = NULL, lat = NULL, lon = NULL,
                     variables.names = "clay",
                     depths = c("0-5cm", "5-15cm", "15-30cm",
                                "30-60cm", "60-100cm", "100-200cm"),
                     stat = c("mean", "q_5", "q_50", "q_95", "uncertainty"),
                     wide = TRUE, sleep = 12, tries = 3L, timeout = 60,
                     verbose = TRUE) {

  .et_banner("get_soil", "downloads soil properties from SoilGrids", verbose)
  .need_pkg(c("httr", "jsonlite"), "get_soil")
  stat <- match.arg(stat)
  sleep <- max(sleep, 12)  # SoilGrids fair-use: 5 requests/min

  if (is.null(lat) || is.null(lon))
    stop("'lat' and 'lon' are required.", call. = FALSE)
  if (length(lat) != length(lon))
    stop("'lat' and 'lon' must have the same length.", call. = FALSE)
  if (any(lat < -90 | lat > 90)) stop("Latitude outside [-90, 90].", call. = FALSE)
  if (any(lon < -180 | lon > 180)) stop("Longitude outside [-180, 180].", call. = FALSE)

  n_env <- length(lat)
  if (is.null(env.id)) env.id <- paste0("env", seq_len(n_env))
  env.id <- as.character(env.id)

  accepted <- c("bdod", "cec", "clay", "cfvo", "nitrogen", "ocd", "ocs",
                "phh2o", "sand", "silt", "soc", "wv0010", "wv0033", "wv1500")
  bad <- setdiff(variables.names, accepted)
  if (length(bad))
    stop("Unsupported soil propert", if (length(bad) > 1) "ies: " else "y: ",
         paste(bad, collapse = ", "),
         "\n  Accepted: ", paste(accepted, collapse = ", "), call. = FALSE)

  base_url <- "https://rest.isric.org/soilgrids/v2.0/properties/query"

  fetch_one <- function(lon, lat) {
    q <- paste0(base_url, "?lon=", lon, "&lat=", lat,
                paste0("&property=", variables.names, collapse = ""),
                paste0("&depth=", depths, collapse = ""),
                "&value=Q0.05&value=Q0.5&value=Q0.95&value=mean&value=uncertainty")
    txt <- .retry({
      r <- httr::GET(q, httr::timeout(timeout))
      if (httr::status_code(r) != 200)
        stop("SoilGrids returned HTTP ", httr::status_code(r))
      httr::content(r, "text", encoding = "UTF-8")
    }, tries = tries, quiet = !verbose)

    js <- jsonlite::fromJSON(txt)
    if (is.null(js$properties$layers) || !length(js$properties$layers$name))
      return(NULL)

    lay <- js$properties$layers
    do.call(rbind, lapply(seq_along(lay$name), function(i) {
      dep <- lay$depths[[i]]
      if (is.null(dep) || !nrow(dep)) return(NULL)
      vals <- dep$values
      data.frame(
        property    = lay$name[i],
        unit        = lay$unit_measure$mapped_units[i],
        depth       = dep$label,
        q_5         = if (!is.null(vals[["Q0.05"]])) vals[["Q0.05"]] else NA_real_,
        q_50        = if (!is.null(vals[["Q0.5"]]))  vals[["Q0.5"]]  else NA_real_,
        q_95        = if (!is.null(vals[["Q0.95"]])) vals[["Q0.95"]] else NA_real_,
        mean        = if (!is.null(vals[["mean"]]))  vals[["mean"]]  else NA_real_,
        uncertainty = if (!is.null(vals[["uncertainty"]])) vals[["uncertainty"]] else NA_real_,
        lon = lon, lat = lat,
        stringsAsFactors = FALSE)
    }))
  }

  startTime <- Sys.time()

  # SoilGrids depends only on lon/lat, but several env.id often share one site
  # (e.g. trials repeated across years); fetch each distinct coordinate once.
  coord.key <- paste(lon, lat, sep = "|")
  uniq.idx  <- which(!duplicated(coord.key))
  n_uniq    <- length(uniq.idx)

  if (verbose) {
    message("---------------------------------------------------------------")
    message("get_soil() - soil features from SoilGrids (ISRIC)")
    message("---------------------------------------------------------------")
    message("Environmental units ................ ", n_env)
    message("Unique coordinates ................. ", n_uniq)
  }

  # Fetch unique points, spacing request *starts* 'sleep' seconds apart so the
  # network time is absorbed into the interval instead of added on top of it.
  cache <- vector("list", n_uniq)
  names(cache) <- coord.key[uniq.idx]
  for (j in seq_len(n_uniq)) {
    i  <- uniq.idx[j]
    t0 <- Sys.time()
    got <- try(fetch_one(lon[i], lat[i]), silent = TRUE)
    if (!inherits(got, "try-error") && !is.null(got)) cache[[j]] <- got
    if (verbose) message("  [", j, "/", n_uniq, "] ", lon[i], ", ", lat[i])
    if (j < n_uniq) {
      wait <- sleep - as.numeric(difftime(Sys.time(), t0, units = "secs"))
      if (wait > 0) Sys.sleep(wait)
    }
  }

  # Map each cached point back to every environment that shares its coordinate.
  res <- vector("list", n_env)
  for (i in seq_len(n_env)) {
    got <- cache[[coord.key[i]]]
    if (is.null(got)) {
      warning("No soil data for '", env.id[i], "'.", call. = FALSE)
    } else {
      got$env <- env.id[i]
      res[[i]] <- got
    }
  }

  out <- do.call(rbind, res)
  if (is.null(out) || !nrow(out))
    stop("SoilGrids returned no usable data.", call. = FALSE)

  out$feature <- paste0(out$property, "|", gsub("-", "_", out$depth))
  out <- out[, c("env", "lon", "lat", "property", "depth", "feature", "unit",
                 "q_5", "q_50", "q_95", "mean", "uncertainty")]

  if (isTRUE(wide)) {
    out <- .long_to_wide(out, id.cols = c("env", "lon", "lat"),
                         key.col = "feature", value.col = stat)
  }

  if (verbose) {
    message("Done ............................... ",
            format(Sys.time(), "%a %b %d %X %Y"))
    message("Total time ......................... ",
            round(difftime(Sys.time(), startTime, units = "secs"), 2), "s")
  }
  out
}


# =========================================================================
# SECTION 4 - param_radiation()
# =========================================================================

#' @title Estimate Basic Solar-Radiation Parameters
#'
#' @description
#' Computes extraterrestrial radiation, daylength and actual sunshine duration from
#' day of year and latitude, following FAO-56.
#'
#' @author Germano Costa Neto
#'
#' @param env.data data.frame. A \code{get_weather} output or equivalent.
#' @param day.id character. Column holding day of year. Default \code{"DOY"}.
#' @param latitude character. Column holding latitude. Default \code{"LAT"}.
#' @param solarRadiation character. Column holding incoming shortwave radiation.
#'   Default \code{"ALLSKY_SFC_SW_DWN"}.
#' @param merge boolean. Bind results onto \code{env.data}. Default \code{TRUE}.
#' @param verbose boolean. Print a summary of what was computed.
#'
#' @return
#' A data.frame with \code{n} (actual sunshine hours), \code{N} (daylength, h) and
#' \code{RTA} (extraterrestrial radiation, MJ/m2/day).
#'
#' @details
#' Daylength uses Forsythe et al. (1995), consistent with \code{get_weather}.
#' The two functions previously used different formulas that disagreed by up to
#' 0.5 h at high latitude.
#'
#' The sunset hour angle is clamped so polar day/night yields 0 rather than
#' \code{NaN}. Missing radiation is median-imputed within the supplied data, and
#' \code{n} is bounded to \code{[0, N]}.
#'
#' @examples
#' \dontrun{
#' env.data <- get_weather(lat = -13.05, lon = -56.05)
#' param_radiation(env.data)
#' }
#'
#' @seealso \code{param_atmospheric}, \code{param_temperature},
#'   \code{processWTH}
#'
#' @importFrom stats median
#' @export
param_radiation <- function(env.data, day.id = "DOY", latitude = "LAT",
                            solarRadiation = NULL, merge = TRUE, verbose = TRUE) {

  .et_banner("param_radiation", "derives radiation parameters", verbose)
  if (is.null(solarRadiation)) solarRadiation <- "ALLSKY_SFC_SW_DWN"
  need <- c(day.id, latitude, solarRadiation)
  miss <- setdiff(need, names(env.data))
  if (length(miss))
    stop("param_radiation() needs column(s): ", paste(miss, collapse = ", "),
         call. = FALSE)

  DOY <- env.data[[day.id]]
  LAT <- env.data[[latitude]]
  DWN <- .strip_power_na(env.data[[solarRadiation]])

  if (all(is.na(DWN)))
    stop("'", solarRadiation, "' is entirely missing.", call. = FALSE)
  if (anyNA(DWN)) DWN[is.na(DWN)] <- stats::median(DWN, na.rm = TRUE)

  rd <- .ra_daylength(DOY, LAT)
  ratio <- DWN / rd$Ra
  ratio[!is.finite(ratio)] <- NA
  n <- pmin(pmax(rd$N * ratio, 0), rd$N)

  if (verbose) {
    message("----------------------------------------------------------------------")
    message("Extraterrestrial radiation (RTA, MJ/m^2/day)")
    message("Daylight hours (N, hours)")
    message("Actual duration of sunshine (n, hours)")
    message("----------------------------------------------------------------------")
  }

  res <- data.frame(n = n, N = rd$N, RTA = rd$Ra)
  if (isTRUE(merge)) cbind(env.data, res) else res
}


# =========================================================================
# SECTION 5 - param_atmospheric()
# =========================================================================

#' @title Estimate Atmospheric Parameters Related to Evapotranspiration
#'
#' @description
#' Computes vapour pressure deficit, the slope of the saturation vapour pressure
#' curve, Priestley-Taylor potential evapotranspiration and the precipitation
#' balance, following FAO-56.
#'
#' @author Germano Costa Neto
#'
#' @param env.data data.frame. A \code{get_weather} output or equivalent.
#' @param PREC,Tdew,Tmin,Tmax,RH,Rad character. Column names for precipitation,
#'   dew point, minimum and maximum temperature, relative humidity and solar
#'   radiation. Defaults follow NASA POWER naming.
#' @param G numeric. Soil heat flux (MJ/m2/day). Default 0.
#' @param Alt numeric scalar, or character naming an elevation column in
#'   \code{env.data}. Default 600 m.
#' @param alpha numeric. Priestley-Taylor coefficient. Default 1.26.
#' @param merge boolean. Bind results onto \code{env.data}.
#' @param verbose boolean. Print a summary of what was computed.
#'
#' @return
#' A data.frame with \code{VPD} (kPa), \code{SPV} (kPa/C), \code{ETP} (mm/day) and
#' \code{PETP} (mm/day).
#'
#' @details
#' \strong{Elevation now works.} \code{Alt} accepts either a numeric elevation or the
#' name of a column. Previously the default \code{Alt = 600} was never \code{NULL},
#' so the elevation column was never read and every site was treated as 600 m --
#' a 26\% error in atmospheric pressure at 2500 m, propagating into the
#' psychrometric constant and ETP.
#'
#' \code{VPD} is the difference between mean saturation vapour pressure (from Tmin
#' and Tmax) and actual vapour pressure (from dew point). Values are floored at 0;
#' negatives can only arise from a dew point above air temperature, which is
#' unphysical and indicates a data problem. A warning is raised when this occurs.
#'
#' @examples
#' \dontrun{
#' env.data <- get_weather(lat = -13.05, lon = -56.05)
#' param_atmospheric(env.data, Alt = 850)
#' }
#'
#' @seealso \code{param_radiation}, \code{param_temperature},
#'   \code{processWTH}
#'
#' @export
param_atmospheric <- function(env.data, PREC = NULL, Tdew = NULL,
                              Tmin = NULL, Tmax = NULL, RH = NULL, Rad = NULL,
                              G = 0, Alt = 600, alpha = 1.26,
                              merge = FALSE, verbose = TRUE) {

  .et_banner("param_atmospheric", "derives atmospheric-demand parameters", verbose)
  if (is.null(PREC)) PREC <- "PRECTOT"
  if (is.null(Tdew)) Tdew <- "T2MDEW"
  if (is.null(Tmin)) Tmin <- "T2M_MIN"
  if (is.null(Tmax)) Tmax <- "T2M_MAX"
  if (is.null(RH))   RH   <- "RH2M"
  if (is.null(Rad))  Rad  <- "ALLSKY_SFC_SW_DWN"

  need <- c(PREC, Tdew, Tmin, Tmax, Rad)
  miss <- setdiff(need, names(env.data))
  if (length(miss))
    stop("param_atmospheric() needs column(s): ", paste(miss, collapse = ", "),
         call. = FALSE)

  prec <- .strip_power_na(env.data[[PREC]])
  tdew <- .strip_power_na(env.data[[Tdew]])
  tmin <- .strip_power_na(env.data[[Tmin]])
  tmax <- .strip_power_na(env.data[[Tmax]])
  rad  <- .strip_power_na(env.data[[Rad]])

  # Elevation: numeric scalar/vector, or a column name.
  if (is.character(Alt)) {
    if (!Alt %in% names(env.data))
      stop("Elevation column '", Alt, "' not found in env.data.", call. = FALSE)
    elev <- .strip_power_na(env.data[[Alt]])
  } else {
    elev <- Alt
  }
  if (is.null(G)) G <- 0

  tmed  <- (tmin + tmax) / 2
  atp   <- .atm_pressure(elev)
  psy   <- 0.665e-3 * atp
  slope <- 4098 * (0.6108 * exp((17.27 * tmed) / (tmed + 237.3))) / (tmed + 237.3)^2

  vpd_raw <- (.teten(tmin) + .teten(tmax)) / 2 - .teten(tdew)
  n_neg <- sum(vpd_raw < 0, na.rm = TRUE)
  if (n_neg)
    warning(n_neg, " negative VPD value(s) floored at 0: dew point exceeds air ",
            "temperature, which indicates a data-quality problem.", call. = FALSE)
  vpd <- pmax(vpd_raw, 0)

  w   <- slope / (slope + psy)
  eto <- alpha * w * (rad - G) * 0.408
  eto <- pmax(eto, 0)
  petp <- prec - eto

  if (verbose) {
    message("----------------------------------------------------------------------")
    message("Slope of saturation vapour pressure curve (SPV, kPa/C)")
    message("Vapour pressure deficit (VPD, kPa)")
    message("Potential evapotranspiration (ETP, mm/day)")
    message("Precipitation balance P - ETP (PETP, mm/day)")
    message("----------------------------------------------------------------------")
  }

  res <- data.frame(VPD = vpd, SPV = slope, ETP = eto, PETP = petp)
  if (isTRUE(merge)) cbind(env.data, res) else res
}


# =========================================================================
# SECTION 6 - param_temperature()
# =========================================================================

#' @title Estimate Thermal Parameters
#'
#' @description
#' Computes growing degree days, the temperature-limited radiation use efficiency
#' factor, and daily temperature range.
#'
#' @author Germano Costa Neto
#'
#' @param env.data data.frame. A \code{get_weather} output or equivalent.
#' @param Tmax,Tmin character. Column names for maximum and minimum air temperature
#'   (Celsius).
#' @param Tbase1 numeric. Base temperature below which development stops.
#' @param Tbase2 numeric. Upper temperature at which development stops.
#' @param Topt1,Topt2 numeric. Lower and upper bounds of the optimal range.
#' @param method character. GDD formulation: \code{"standard"} (default) or
#'   \code{"legacy"}. See Details.
#' @param merge boolean. Bind results onto \code{env.data}.
#' @param verbose boolean. Print a summary of what was computed.
#'
#' @return
#' A data.frame with \code{GDD} (C/day), \code{FRUE} (0-1) and \code{T2M_RANGE} (C).
#'
#' @details
#' \strong{GDD methods.}
#' \describe{
#'   \item{\code{"standard"}}{\eqn{GDD = \max(\frac{\min(T_{max}, T_{base2}) + T_{min}}{2} - T_{base1},\ 0)}.
#'     Tmax is capped at the upper threshold; the result is floored at zero. This is
#'     the conventional agronomic definition.}
#'   \item{\code{"legacy"}}{Clamps \emph{both} Tmax and Tmin into
#'     \eqn{[T_{base1}, T_{base2}]} before averaging. This is the behaviour of
#'     earlier EnvRtype versions and corresponds to the "Method 2" variant. It
#'     credits heat units on days colder than \eqn{T_{base1}} -- e.g. with
#'     Tmax 12, Tmin 5, Tbase1 9 it returns 1.5 where the standard method returns 0.}
#' }
#' The default changed from legacy to standard; results will differ on cold days.
#' Pass \code{method = "legacy"} to reproduce earlier output.
#'
#' \code{FRUE} rises linearly from 0 at \code{Tbase1} to 1 at \code{Topt1}, holds at
#' 1 across the optimum, and falls to 0 at \code{Tbase2}.
#'
#' @examples
#' \dontrun{
#' env.data <- get_weather(lat = -13.05, lon = -56.05)
#' param_temperature(env.data)
#' param_temperature(env.data, method = "legacy")  # pre-refactor behaviour
#' }
#'
#' @seealso \code{param_radiation}, \code{param_atmospheric},
#'   \code{processWTH}
#'
#' @export
param_temperature <- function(env.data, Tmax = NULL, Tmin = NULL,
                              Tbase1 = 9, Tbase2 = 45, Topt1 = 26, Topt2 = 32,
                              method = c("standard", "legacy"),
                              merge = FALSE, verbose = TRUE) {

  .et_banner("param_temperature", "derives temperature-response parameters", verbose)
  method <- match.arg(method)
  if (is.null(Tmin)) Tmin <- "T2M_MIN"
  if (is.null(Tmax)) Tmax <- "T2M_MAX"

  miss <- setdiff(c(Tmin, Tmax), names(env.data))
  if (length(miss))
    stop("param_temperature() needs column(s): ", paste(miss, collapse = ", "),
         call. = FALSE)

  if (!(Tbase1 < Topt1 && Topt1 <= Topt2 && Topt2 < Tbase2))
    stop("Cardinal temperatures must satisfy Tbase1 < Topt1 <= Topt2 < Tbase2.",
         call. = FALSE)

  tmin <- .strip_power_na(env.data[[Tmin]])
  tmax <- .strip_power_na(env.data[[Tmax]])

  n_inv <- sum(tmax < tmin, na.rm = TRUE)
  if (n_inv)
    warning(n_inv, " day(s) with Tmax < Tmin; check the input columns.", call. = FALSE)

  tmed <- (tmin + tmax) / 2

  # FRUE: vectorised trapezoid over the cardinal temperatures
  frue <- rep(1, length(tmed))
  lo <- !is.na(tmed) & tmed < Topt1
  hi <- !is.na(tmed) & tmed > Topt2
  frue[lo] <- (tmed[lo] - Tbase1) / (Topt1 - Tbase1)
  frue[hi] <- (Tbase2 - tmed[hi]) / (Tbase2 - Topt2)
  frue <- pmin(pmax(frue, 0), 1)
  frue[is.na(tmed)] <- NA

  gdd <- if (method == "standard") {
    pmax((pmin(tmax, Tbase2) + tmin) / 2 - Tbase1, 0)
  } else {
    clamp <- function(x) pmin(pmax(x, Tbase1), Tbase2)
    (clamp(tmax) + clamp(tmin)) / 2 - Tbase1
  }

  if (verbose) {
    message("----------------------------------------------------------------------")
    message("Growing Degree Day (GDD, C/day) - method: ", method)
    message("Effect of temperature on radiation use efficiency (FRUE, 0-1)")
    message("Daily temperature range (T2M_RANGE, C)")
    message("----------------------------------------------------------------------")
  }

  res <- data.frame(GDD = gdd, FRUE = frue, T2M_RANGE = tmax - tmin)
  if (isTRUE(merge)) cbind(env.data, res) else res
}


# =========================================================================
# SECTION 7 - processWTH()
# =========================================================================

#' @title Enrich a get_weather() Output with Derived Parameters
#'
#' @description
#' Convenience wrapper chaining \code{param_radiation},
#' \code{param_atmospheric} and \code{param_temperature}.
#'
#' @author Germano Costa Neto
#'
#' @param env.data data.frame. A \code{get_weather} output.
#' @param PREC,Tdew,Tmin,Tmax,RH,Rad character. Column names passed through to the
#'   underlying functions.
#' @param G numeric. Soil heat flux. Default 0.
#' @param Alt numeric or character. Elevation in metres, or the name of an elevation
#'   column.
#' @param Tbase1,Tbase2,Topt1,Topt2 numeric. Cardinal temperatures.
#' @param method character. GDD method, see \code{param_temperature}.
#' @param alpha numeric. Priestley-Taylor coefficient.
#' @param verbose boolean. Print progress messages.
#'
#' @return
#' The input data.frame with radiation, atmospheric and thermal parameters appended.
#'
#' @details
#' If \code{env.data} already carries \code{n}, \code{N} or \code{RTA} from
#' \code{get_weather}, those columns are dropped before recomputation to avoid
#' duplicated names. Note that \code{get_weather()} derives \code{n} from
#' \code{ALLSKY_SFC_SW_DNI} whereas \code{param_radiation()} uses
#' \code{ALLSKY_SFC_SW_DWN}, so the recomputed values may differ.
#'
#' @examples
#' \dontrun{
#' env.data <- get_weather(lat = -13.05, lon = -56.05)
#' processWTH(env.data)
#' }
#'
#' @seealso \code{get_weather}, \code{summaryWTH}
#'
#' @export
processWTH <- function(env.data, PREC = NULL, Tdew = NULL, Tmin = NULL, Tmax = NULL,
                       RH = NULL, Rad = NULL, G = 0, Alt = 600,
                       Tbase1 = 9, Tbase2 = 45, Topt1 = 26, Topt2 = 32,
                       method = c("standard", "legacy"), alpha = 1.26,
                       verbose = TRUE) {

  .et_banner("processWTH", "derives all agro-meteorological parameters", verbose)
  method <- match.arg(method)

  # Every column these three functions produce. get_weather() already emits
  # VPD, n, N and RTA, so without this the cbind()s below append duplicate
  # names -- and df$VPD then silently resolves to whichever came first, even
  # though param_atmospheric() floors VPD at 0 and get_weather() does not.
  produced <- c("n", "N", "RTA", "VPD", "SPV", "ETP", "PETP",
                "GDD", "FRUE", "T2M_RANGE")
  dup <- intersect(produced, names(env.data))
  if (length(dup)) {
    if (verbose)
      message("Dropping existing column(s) before recomputation: ",
              paste(dup, collapse = ", "))
    env.data <- env.data[, setdiff(names(env.data), dup), drop = FALSE]
  }

  env.data <- param_radiation(env.data = env.data, merge = TRUE, verbose = FALSE)
  .et_step("radiation parameters computed", verbose)
  env.data <- param_atmospheric(env.data = env.data, PREC = PREC, Tdew = Tdew,
                                Tmin = Tmin, Tmax = Tmax, RH = RH, Rad = Rad,
                                G = G, Alt = Alt, alpha = alpha,
                                merge = TRUE, verbose = FALSE)
  .et_step("atmospheric-demand parameters computed", verbose)
  env.data <- param_temperature(env.data = env.data, Tmax = Tmax, Tmin = Tmin,
                                Tbase1 = Tbase1, Tbase2 = Tbase2,
                                Topt1 = Topt1, Topt2 = Topt2, method = method,
                                merge = TRUE, verbose = FALSE)
  .et_step("temperature-response parameters computed", verbose)
  env.data
}


# =========================================================================
# SECTION 8 - summaryWTH()
# =========================================================================

#' @title Basic Summary Statistics for Environmental Data
#'
#' @description
#' Summarises a \code{get_weather} output by environment and, optionally, by
#' user-defined time intervals such as phenological stages.
#'
#' @author Germano Costa Neto
#'
#' @param env.data data.frame. A \code{get_weather} output.
#' @param id.names vector (character). Columns treated as identifiers.
#' @param env.id character. Environment-id column. Default \code{"env"}.
#' @param days.id character. Column holding days from start. Default
#'   \code{"daysFromStart"}.
#' @param var.id vector (character). Variables to summarise. Defaults to every
#'   non-identifier column.
#' @param statistic character. One of \code{"all"}, \code{"sum"}, \code{"mean"},
#'   \code{"quantile"}.
#' @param probs vector (numeric). Quantile probabilities. Default
#'   \code{c(0.25, 0.50, 0.75)}.
#' @param by.interval boolean. Summarise within time intervals.
#' @param time.window vector (numeric). Interval breakpoints in days.
#' @param names.window vector (character). Interval names.
#' @param verbose boolean. If \code{TRUE} (default) prints a progress banner.
#'
#' @return
#' A data.frame in long format with one row per environment (and interval) and
#' variable.
#'
#' @details
#' Quantiles are computed with \code{stats::aggregate} rather than the previous
#' nested \code{foreach} loops, which registered no parallel backend and therefore
#' ran sequentially while emitting a warning on every call.
#'
#' @examples
#' \dontrun{
#' env.data <- get_weather(lat = -13.05, lon = -56.05)
#' summaryWTH(env.data, env.id = "env", statistic = "mean")
#'
#' summaryWTH(env.data, env.id = "env", by.interval = TRUE,
#'            time.window  = c(0, 14, 35, 60, 90, 120),
#'            names.window = c("P-E", "E-V1", "V1-V4", "V4-VT", "VT-GF", "GF-PM"))
#' }
#'
#' @seealso \code{get_weather}, \code{processWTH}
#'
#' @importFrom stats aggregate quantile
#' @export
summaryWTH <- function(env.data, id.names = NULL, env.id = NULL, days.id = NULL,
                       var.id = NULL, statistic = NULL, probs = NULL,
                       by.interval = FALSE, time.window = NULL,
                       names.window = NULL, verbose = TRUE) {
  .et_banner("summaryWTH", "summarises weather by environment and interval", verbose)
  env.data <- as.data.frame(env.data)

  if (is.null(statistic)) statistic <- "all"
  statistic <- statistic[1]
  if (!statistic %in% c("all", "sum", "mean", "quantile")) statistic <- "all"
  if (is.null(probs)) probs <- c(0.25, 0.50, 0.75)

  if (is.null(env.id))  env.id  <- "env"
  if (is.null(days.id)) days.id <- "daysFromStart"
  if (!env.id %in% names(env.data))
    stop("Environment column '", env.id, "' not found.", call. = FALSE)
  names(env.data)[names(env.data) == env.id]  <- "env"
  if (days.id %in% names(env.data))
    names(env.data)[names(env.data) == days.id] <- "daysFromStart"

  if (is.null(id.names))
    id.names <- c("env", "LON", "LAT", "YEAR", "MM", "DD", "DOY",
                  "YYYYMMDD", "daysFromStart")
  id.names <- intersect(id.names, names(env.data))
  if (is.null(var.id)) var.id <- setdiff(names(env.data), id.names)
  var.id <- intersect(var.id, names(env.data))
  if (!length(var.id)) stop("No variables left to summarise.", call. = FALSE)

  if (isTRUE(by.interval)) {
    if (!"daysFromStart" %in% names(env.data))
      stop("by.interval = TRUE requires a '", days.id, "' column.", call. = FALSE)
    env.data$interval <- .env_stage_by_dae(.dae = env.data[["daysFromStart"]],
                                           .breaks = time.window,
                                           .names = names.window)
    keep <- !is.na(env.data$interval)
    if (any(!keep)) env.data <- env.data[keep, , drop = FALSE]
    env.data$interval <- droplevels(env.data$interval)
    id.names <- c(id.names, "interval")
  }
  # Wide -> long without reshape2: keeps the package's most-used summary function
  # dependency-free.
  ts <- do.call(rbind, lapply(var.id, function(v) {
    d <- env.data[, id.names, drop = FALSE]
    d$variable <- factor(v, levels = var.id)
    d$value <- suppressWarnings(as.numeric(env.data[[v]]))
    d
  }))
  rownames(ts) <- NULL

  # Variables with no valid data anywhere yield NaN means and 0 sums, which read
  # as real numbers in the output table. Tell the user instead of letting the
  # zeros pass as measurements.
  allna <- var.id[vapply(var.id, function(v)
    all(is.na(ts$value[ts$variable == v])), logical(1))]
  if (length(allna))
    warning("No valid data for: ", paste(allna, collapse = ", "),
            ". Their summaries are NA/NaN; consider dropping them via 'var.id'.",
            call. = FALSE)

  grp <- if (isTRUE(by.interval)) c("env", "interval", "variable") else
    c("env", "variable")
  gl <- ts[grp]

  agg_mean <- function() {
    # mean(NA, na.rm = TRUE) is NaN; return NA for consistency with the other
    # aggregates so all-missing variables read uniformly across the table.
    o <- stats::aggregate(ts$value, by = gl, FUN = function(z)
      if (all(is.na(z))) NA_real_ else mean(z, na.rm = TRUE))
    names(o) <- c(grp, "mean"); o
  }
  agg_sum <- function() {
    # sum(NA, na.rm = TRUE) is 0, which reads as a genuine zero total (e.g. "no
    # rainfall") when in fact nothing was measured. Return NA in that case.
    o <- stats::aggregate(ts$value, by = gl, FUN = function(z)
      if (all(is.na(z))) NA_real_ else sum(z, na.rm = TRUE))
    names(o) <- c(grp, "sum"); o
  }
  agg_quant <- function() {
    q <- stats::aggregate(ts$value, by = gl, FUN = function(z)
      stats::quantile(z, probs = probs, na.rm = TRUE))
    qm <- as.data.frame(q$x)
    names(qm) <- paste0("prob_", probs)
    cbind(q[grp], qm)
  }

  switch(statistic,
         mean     = agg_mean(),
         sum      = agg_sum(),
         quantile = agg_quant(),
         all      = merge(merge(agg_mean(), agg_sum(), by = grp),
                          agg_quant(), by = grp))
}


#==================================================================================================
# SECTION 10 - CRASH-RESILIENT RETRIEVAL
#
#
#
# --------------------------------------------------------------------------
# PROBLEM
# --------------------------------------------------------------------------
# Downloading hundreds of environments from NASA POWER or SoilGrids takes
# hours. If the API times out, the network drops, or R is killed at
# environment 480 of 500, in-memory results are lost and the job restarts
# from zero.
#
# --------------------------------------------------------------------------
# DESIGN
# --------------------------------------------------------------------------
# Each *_resumable() call owns a run directory:
#
#   <dir.path>/<run.id>/
#     manifest.json     run parameters + fingerprint of the request
#     progress.log      one append-only JSON line per environment
#     parts/<key>.rds   one result file per successful environment
#
# Every environment is written to disk the moment it succeeds; nothing is
# held only in memory. restart_from_log() reads the log, skips what already
# succeeded, retries what failed, and assembles the final object.
#
# The log is APPEND-ONLY and flushed per line. A process killed mid-write
# loses at most the final line, and a truncated trailing line is tolerated on
# read. Rewriting a whole state file each iteration would risk losing
# everything to a single bad kill.
#
# Resumption is keyed on a FINGERPRINT of the request. Changing coordinates,
# dates or variables changes the fingerprint, and restart_from_log() refuses
# to mix incompatible results rather than silently returning a blend of two
# different queries.
#
# --------------------------------------------------------------------------
# DEPENDENCIES
# --------------------------------------------------------------------------
# Base R only. JSON is read and written by small internal helpers
# (.ck_to_json / .ck_from_json) so that a crash-recovery tool never fails to
# load because an optional package is missing. jsonlite is used when present
# but is not required.
#
# Inherited at call time: whatever get_weather()/get_soil() need
# (nasapower, httr, jsonlite).
#==================================================================================================


# =========================================================================
# SECTION 10.0 - internal helpers
# =========================================================================

#' Null-coalescing helper
#' @keywords internal
#' @noRd
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

#' Escape a string for JSON
#' @keywords internal
#' @noRd
.ck_esc <- function(x) {
  x <- gsub("\\", "\\\\", x, fixed = TRUE)
  x <- gsub("\"", "\\\"", x, fixed = TRUE)
  x <- gsub("\n", " ", x, fixed = TRUE)
  x <- gsub("\r", " ", x, fixed = TRUE)
  x <- gsub("\t", " ", x, fixed = TRUE)
  x
}

#' Minimal JSON writer for flat lists
#'
#' Handles the only shapes this module stores: scalars and atomic vectors of
#' character, numeric, integer or logical. Deliberately tiny -- using
#' jsonlite here would make crash recovery depend on a package that may not
#' be installed in the session doing the recovery.
#'
#' @param x named list.
#' @return character(1) JSON object.
#' @keywords internal
#' @noRd
.ck_to_json <- function(x) {
  if (requireNamespace("jsonlite", quietly = TRUE))
    return(as.character(jsonlite::toJSON(x, auto_unbox = TRUE, null = "null")))
  one <- function(v) {
    if (is.null(v) || length(v) == 0) return("null")
    if (is.logical(v))  s <- ifelse(is.na(v), "null", ifelse(v, "true", "false"))
    else if (is.numeric(v)) s <- ifelse(is.na(v), "null", format(v, scientific = FALSE,
                                                                 trim = TRUE))
    else s <- ifelse(is.na(v), "null", paste0("\"", .ck_esc(as.character(v)), "\""))
    if (length(v) == 1L) s else paste0("[", paste(s, collapse = ","), "]")
  }
  paste0("{", paste(sprintf("\"%s\":%s", .ck_esc(names(x)),
                            vapply(x, one, character(1))), collapse = ","), "}")
}

#' Minimal JSON reader for flat objects
#'
#' Parses the flat objects written by \code{.ck_to_json}. Nested structures
#' are not supported and are not produced by this module.
#'
#' @param txt character(1) JSON object.
#' @return named list.
#' @keywords internal
#' @noRd
.ck_from_json <- function(txt) {
  if (requireNamespace("jsonlite", quietly = TRUE))
    return(jsonlite::fromJSON(txt))
  txt <- trimws(txt)
  if (!nzchar(txt) || !grepl("^\\{", txt) || !grepl("\\}$", txt))
    stop("not a JSON object", call. = FALSE)
  body <- substr(txt, 2, nchar(txt) - 1L)

  # Split on commas that sit outside quotes and outside brackets.
  ch <- strsplit(body, "")[[1]]
  inq <- FALSE; esc <- FALSE; dep <- 0L; cuts <- integer(0)
  for (i in seq_along(ch)) {
    c0 <- ch[i]
    if (esc) { esc <- FALSE; next }
    if (c0 == "\\") { esc <- TRUE; next }
    if (c0 == "\"") { inq <- !inq; next }
    if (!inq) {
      if (c0 == "[") dep <- dep + 1L
      else if (c0 == "]") dep <- dep - 1L
      else if (c0 == "," && dep == 0L) cuts <- c(cuts, i)
    }
  }
  starts <- c(1L, cuts + 1L); ends <- c(cuts - 1L, length(ch))
  pieces <- mapply(function(a, b) if (b >= a) paste(ch[a:b], collapse = "") else "",
                   starts, ends, USE.NAMES = FALSE)
  pieces <- pieces[nzchar(trimws(pieces))]
  if (!length(pieces)) return(list())

  out <- list()
  for (p in pieces) {
    m <- regexpr(":", p, fixed = TRUE)
    if (m < 0) next
    k <- trimws(substr(p, 1, m - 1L))
    v <- trimws(substr(p, m + 1L, nchar(p)))
    k <- gsub("^\"|\"$", "", k)
    parse_scalar <- function(s) {
      s <- trimws(s)
      if (identical(s, "null")) return(NA)
      if (identical(s, "true")) return(TRUE)
      if (identical(s, "false")) return(FALSE)
      if (grepl("^\".*\"$", s)) {
        s <- substr(s, 2, nchar(s) - 1L)
        s <- gsub("\\\"", "\"", s, fixed = TRUE)
        s <- gsub("\\\\", "\\", s, fixed = TRUE)
        return(s)
      }
      n <- suppressWarnings(as.numeric(s))
      if (!is.na(n)) n else s
    }
    if (grepl("^\\[", v) && grepl("\\]$", v)) {
      inner <- substr(v, 2, nchar(v) - 1L)
      if (!nzchar(trimws(inner))) { out[[k]] <- character(0); next }
      ic <- strsplit(inner, "")[[1]]
      inq2 <- FALSE; esc2 <- FALSE; cc <- integer(0)
      for (i in seq_along(ic)) {
        c0 <- ic[i]
        if (esc2) { esc2 <- FALSE; next }
        if (c0 == "\\") { esc2 <- TRUE; next }
        if (c0 == "\"") { inq2 <- !inq2; next }
        if (!inq2 && c0 == ",") cc <- c(cc, i)
      }
      st <- c(1L, cc + 1L); en <- c(cc - 1L, length(ic))
      el <- mapply(function(a, b) if (b >= a) paste(ic[a:b], collapse = "") else "",
                   st, en, USE.NAMES = FALSE)
      vals <- lapply(el, parse_scalar)
      out[[k]] <- if (all(vapply(vals, is.character, logical(1))))
        unlist(vals) else unlist(vals)
    } else {
      out[[k]] <- parse_scalar(v)
    }
  }
  out
}

#' Stable fingerprint of a retrieval request
#'
#' Two runs may be resumed into one another only if they describe the same
#' query. The fingerprint folds every argument that changes what the API
#' returns -- coordinates, dates, variables, depths -- into a short string.
#'
#' Without this, resuming after editing the variable list would silently
#' concatenate results from two different requests, producing a table whose
#' columns mean different things in different rows. That is worse than an
#' error, because nothing about the output looks wrong.
#'
#' @param ... named arguments describing the request.
#' @return character(1) fingerprint.
#' @keywords internal
#' @noRd
.ck_fingerprint <- function(...) {
  a <- list(...)
  a <- a[order(names(a))]
  flat <- vapply(a, function(z) paste(as.character(z), collapse = "\u001f"),
                 character(1))
  s <- paste(names(a), flat, sep = "=", collapse = "\u001e")
  iv <- utf8ToInt(s)
  if (length(iv) > 200000L) iv <- iv[seq_len(200000L)]
  h1 <- 5381; h2 <- 52711
  for (k in iv) {
    h1 <- (h1 * 33 + k) %% 2147483647
    h2 <- (h2 * 31 + k) %% 2147483647
  }
  sprintf("%08x%08x", as.integer(h1), as.integer(h2))
}

#' Sanitise an environment id for use as a file name
#'
#' Environment ids come from users and routinely contain slashes, spaces or
#' accents, none of which belong in a path.
#'
#' @param x character vector.
#' @return character vector safe for file names, de-duplicated.
#' @keywords internal
#' @noRd
.ck_safe_key <- function(x) {
  y <- gsub("[^A-Za-z0-9._-]", "_", as.character(x))
  y[!nzchar(y)] <- "_"
  if (anyDuplicated(y)) y <- paste0(y, "__", seq_along(y))
  y
}

#' Append one JSON record to the progress log
#'
#' Opens in append mode, writes one line, flushes, closes. Deliberately not a
#' held-open connection: a buffered connection is exactly what gets lost when
#' the process dies.
#'
#' @param path character. Log file.
#' @param rec named list. Record to serialise.
#' @keywords internal
#' @noRd
.ck_log_append <- function(path, rec) {
  con <- file(path, open = "at")
  on.exit(close(con), add = TRUE)
  writeLines(.ck_to_json(rec), con)
  flush(con)
  invisible(TRUE)
}

#' Read a progress log, tolerating a truncated final line
#'
#' A hard kill can leave the last line half-written. That line is dropped
#' with a warning rather than aborting the read -- surviving an ungraceful
#' death is the entire point of the log.
#'
#' @param path character. Log file.
#' @return data.frame with \code{env}, \code{key}, \code{status}, \code{time},
#'   \code{message}; zero rows when the log is absent.
#' @keywords internal
#' @noRd
.ck_log_read <- function(path) {
  empty <- data.frame(env = character(0), key = character(0),
                      status = character(0), time = character(0),
                      message = character(0), stringsAsFactors = FALSE)
  if (!file.exists(path)) return(empty)
  ln <- readLines(path, warn = FALSE)
  ln <- ln[nzchar(trimws(ln))]
  if (!length(ln)) return(empty)

  recs <- lapply(seq_along(ln), function(i) {
    r <- tryCatch(.ck_from_json(ln[i]), error = function(e) NULL)
    if (is.null(r) || is.null(r$env)) {
      if (i == length(ln))
        warning("Progress log ends in a truncated record; ignoring it. ",
                "This is expected after a hard crash.", call. = FALSE)
      else
        warning("Unparseable record at line ", i, " of the progress log; ",
                "ignoring it.", call. = FALSE)
      return(NULL)
    }
    data.frame(env = as.character(r$env %||% NA),
               key = as.character(r$key %||% NA),
               status = as.character(r$status %||% NA),
               time = as.character(r$time %||% NA),
               message = as.character(r$message %||% ""),
               stringsAsFactors = FALSE)
  })
  recs <- recs[!vapply(recs, is.null, logical(1))]
  if (!length(recs)) return(empty)
  do.call(rbind, recs)
}

#' Latest status per environment
#'
#' An environment may appear repeatedly (failed, retried, succeeded). Only
#' the last record counts.
#'
#' @param log data.frame from \code{.ck_log_read}.
#' @return data.frame, one row per environment.
#' @keywords internal
#' @noRd
.ck_latest <- function(log) {
  if (!nrow(log)) return(log)
  log$.i <- seq_len(nrow(log))
  sp <- split(log, log$env)
  out <- do.call(rbind, lapply(sp, function(d) d[which.max(d$.i), , drop = FALSE]))
  out$.i <- NULL
  rownames(out) <- NULL
  out
}

#' Prepare or validate a run directory
#'
#' @param dir.path character. Parent directory.
#' @param run.id character. Run name.
#' @param fingerprint character. Request fingerprint.
#' @param meta named list. Extra manifest fields.
#' @param resume logical. Whether an existing run may be reused.
#' @return list of paths.
#' @keywords internal
#' @noRd
.ck_init_run <- function(dir.path, run.id, fingerprint, meta, resume) {
  root  <- file.path(dir.path, run.id)
  parts <- file.path(root, "parts")
  logf  <- file.path(root, "progress.log")
  manf  <- file.path(root, "manifest.json")

  if (!dir.exists(parts))
    dir.create(parts, recursive = TRUE, showWarnings = FALSE)

  if (file.exists(manf)) {
    old <- tryCatch(.ck_from_json(paste(readLines(manf, warn = FALSE),
                                        collapse = "")),
                    error = function(e) NULL)
    if (!is.null(old) && !identical(as.character(old$fingerprint), fingerprint))
      stop("Run '", run.id, "' already exists in\n  ", root,
           "\nbut describes a DIFFERENT request.\n",
           "  stored fingerprint : ", as.character(old$fingerprint), "\n",
           "  current fingerprint: ", fingerprint, "\n",
           "  Arguments that change the query (coordinates, dates, variables) ",
           "must match exactly.\n",
           "  Choose another run.id, or delete that directory.", call. = FALSE)
  } else {
    man <- c(list(run.id = run.id, fingerprint = fingerprint,
                  created = format(Sys.time(), "%Y-%m-%dT%H:%M:%S")), meta)
    writeLines(.ck_to_json(man), manf)
  }

  list(root = root, parts = parts, log = logf, manifest = manf)
}

#' Per-environment retrieval loop with checkpointing
#'
#' Shared by both resumable wrappers so crash semantics are identical for
#' weather and soil.
#'
#' @param keys character. Safe file keys.
#' @param envs character. Environment ids.
#' @param fetch function(i) returning a data.frame.
#' @param paths list from \code{.ck_init_run}.
#' @param sleep numeric. Seconds between calls.
#' @param verbose logical.
#' @param skip character. Ids already completed.
#' @return named integer vector of successes and failures.
#' @keywords internal
#' @noRd
.ck_run_loop <- function(keys, envs, fetch, paths, sleep, verbose, skip,
                         abort.after = 5L) {
  n <- length(envs); nok <- 0L; nbad <- 0L; streak <- 0L
  for (i in seq_len(n)) {
    if (envs[i] %in% skip) {
      if (verbose) message("  [", i, "/", n, "] ", envs[i], " - cached, skipped")
      next
    }
    t0 <- Sys.time()
    # An interrupt (Ctrl-C, SIGINT) means the operator is stopping the job:
    # propagate immediately rather than logging the remaining sites as
    # failures. Errors are caught per-environment -- see the streak check.
    got <- withCallingHandlers(
      tryCatch(fetch(i), error = function(e) e),
      interrupt = function(e) {
        message("\nInterrupted. Progress is saved; resume with ",
                "restart_from_log().")
        invokeRestart("abort")
      })

    if (inherits(got, "error") || is.null(got) ||
        !is.data.frame(got) || !nrow(got)) {
      msg <- if (inherits(got, "error")) conditionMessage(got) else
        "empty or NULL result"
      nbad <- nbad + 1L
      streak <- streak + 1L
      .ck_log_append(paths$log,
                     list(env = envs[i], key = keys[i], status = "failed",
                          time = format(t0, "%Y-%m-%dT%H:%M:%S"),
                          message = substr(msg, 1, 500)))
      warning("Environment '", envs[i], "' failed: ", substr(msg, 1, 200),
              call. = FALSE)
      if (verbose) message("  [", i, "/", n, "] ", envs[i], " - FAILED")

      # A long unbroken run of failures is not bad luck at individual sites --
      # it means the API is down, credentials expired, or the network dropped.
      # Continuing would burn through the remaining sites and write a log in
      # which everything looks permanently failed, defeating the resume. Stop
      # while the log still distinguishes real failures from the outage.
      if (!is.null(abort.after) && streak >= abort.after) {
        stop("Aborted after ", streak, " consecutive failures ",
             "(last: ", substr(msg, 1, 120), ").\n",
             "  This usually means the service is unreachable rather than ",
             "a problem with individual sites.\n",
             "  Progress so far is saved. Resume with ",
             "restart_from_log(run.id = \"", basename(paths$root), "\").",
             call. = FALSE)
      }
    } else {
      streak <- 0L
      saveRDS(got, file.path(paths$parts, paste0(keys[i], ".rds")))
      nok <- nok + 1L
      .ck_log_append(paths$log,
                     list(env = envs[i], key = keys[i], status = "ok",
                          time = format(t0, "%Y-%m-%dT%H:%M:%S"),
                          message = paste0(nrow(got), " rows")))
      if (verbose) message("  [", i, "/", n, "] ", envs[i], " - ok (",
                           nrow(got), " rows)")
    }
    if (i < n && sleep > 0) Sys.sleep(sleep)
  }
  c(ok = nok, failed = nbad)
}

#' Assemble saved parts into one data.frame
#'
#' @param paths list from \code{.ck_init_run}.
#' @param log data.frame of latest statuses.
#' @return data.frame, or NULL when nothing succeeded.
#' @keywords internal
#' @noRd
.ck_assemble <- function(paths, log) {
  done <- log[log$status == "ok", , drop = FALSE]
  if (!nrow(done)) return(NULL)
  fs <- file.path(paths$parts, paste0(done$key, ".rds"))
  keep <- file.exists(fs)
  if (!all(keep))
    warning(sum(!keep), " part file(s) recorded as complete are missing from\n  ",
            paths$parts, "\n  They will be re-fetched on the next ",
            "restart_from_log() call.", call. = FALSE)
  fs <- fs[keep]
  if (!length(fs)) return(NULL)
  lst <- lapply(fs, readRDS)
  binder <- if (exists(".rbind_fill", mode = "function"))
    get(".rbind_fill") else function(z) do.call(rbind, z)
  binder(lst)
}


# =========================================================================
# SECTION 10.1 - get_weather_resumable()
# =========================================================================

#' @title Crash-Resilient Wrapper for get_weather()
#'
#' @description
#' Downloads daily weather one environment at a time, writing each result to
#' disk and appending a line to a progress log as soon as it succeeds. If the
#' session dies part-way through, \code{restart_from_log} resumes from
#' the last completed environment instead of restarting the whole job.
#'
#' Use this instead of \code{get_weather} whenever losing the run
#' would hurt -- in practice, more than a few dozen environments.
#'
#' @param env.id,lat,lon,start.day,end.day,variables.names Passed to
#'   \code{get_weather}. \code{lat} and \code{lon} are required.
#' @param dir.path character. Parent directory for run directories. Default
#'   \code{getwd()}.
#' @param run.id character. Name of this run's directory. Default
#'   \code{"weather_run"}. Use a distinct name per job.
#' @param sleep numeric. Seconds between environments. Default 1.
#' @param tries integer. Retry attempts per environment before it is logged
#'   as failed. Default 3.
#' @param resume logical. If \code{TRUE} (default), environments already
#'   marked \code{ok} in an existing log are skipped.
#' @param verbose logical. Print progress. Default \code{TRUE}.
#'
#' @details
#' \strong{Why per-environment files.} Results are written as one \code{.rds}
#' per environment under \code{parts/}, never accumulated only in memory. A
#' process killed at any point loses at most the environment in flight.
#'
#' \strong{Why an append-only log.} Each line is appended and flushed
#' individually. Rewriting one state file every iteration would risk the whole
#' file to a kill mid-write; appending risks only the last line, and the
#' reader tolerates a truncated final record.
#'
#' \strong{Why a fingerprint.} The manifest stores a hash of coordinates,
#' dates and variables. Resuming with different arguments is refused, because
#' silently concatenating two different queries produces a table that looks
#' fine and means nothing.
#'
#' \strong{Failures do not abort the run.} A failed environment is logged and
#' the loop continues. Inspect with \code{read_progress_log}, then
#' retry with \code{restart_from_log}.
#'
#' \strong{Serial by design.} Unlike \code{get_weather}, this wrapper
#' does not parallelise. Interleaved writes from several workers make progress
#' ordering ambiguous, and NASA POWER rate-limits aggressive clients anyway.
#' The cost is wall-clock time; the benefit is a log you can trust. For speed
#' on a job you can afford to lose, use \code{get_weather(parallel = TRUE)}.
#'
#' @return
#' A \code{data.frame} of daily weather for every environment that succeeded,
#' with attributes \code{run.dir}, \code{n_ok}, \code{n_failed} and
#' \code{failed}.
#'
#' @seealso \code{restart_from_log}, \code{read_progress_log},
#'   \code{get_weather}
#'
#' @examples
#' \dontrun{
#' sites <- data.frame(env = paste0("E", 1:200),
#'                     lat = runif(200, -30, -5),
#'                     lon = runif(200, -60, -40))
#'
#' wth <- get_weather_resumable(
#'   env.id = sites$env, lat = sites$lat, lon = sites$lon,
#'   start.day = "2023-10-01", end.day = "2024-02-01",
#'   run.id = "trial2024")
#'
#' ## --- session dies at environment 137 ---
#' ## In a NEW session, no arguments needed:
#' wth <- restart_from_log(run.id = "trial2024")
#'
#' read_progress_log(run.id = "trial2024")
#' attr(wth, "failed")
#'
#' ## Straight into the characterization pipeline
#' wth <- processWTH(wth)
#' W   <- W_matrix(wth, env.id = "env", var.id = c("T2M", "PRECTOT"))
#' }
#'
#' @export
get_weather_resumable <- function(env.id = NULL, lat = NULL, lon = NULL,
                                  start.day = NULL, end.day = NULL,
                                  variables.names = NULL,
                                  dir.path = NULL, run.id = "weather_run",
                                  sleep = 1, tries = 3L,
                                  resume = TRUE, verbose = TRUE) {

  .et_banner("get_weather_resumable",
             "downloads weather with resumable progress logging", verbose)
  if (!exists("get_weather", mode = "function"))
    stop("get_weather() not found. Source env_data_collection.R first.",
         call. = FALSE)
  if (is.null(lat) || is.null(lon))
    stop("'lat' and 'lon' are required.", call. = FALSE)
  if (length(lat) != length(lon))
    stop("'lat' and 'lon' must have the same length.", call. = FALSE)

  n <- length(lat)
  if (is.null(env.id)) env.id <- paste0("env", seq_len(n))
  env.id <- as.character(env.id)
  if (length(env.id) != n)
    stop("'env.id' must have the same length as 'lat' and 'lon'.", call. = FALSE)
  if (anyDuplicated(env.id))
    stop("Duplicated env.id values. Checkpointing keys results on env.id, ",
         "so ids must be unique.", call. = FALSE)
  if (is.null(dir.path)) dir.path <- getwd()

  sd_v <- if (is.null(start.day)) rep(NA_character_, n) else
    as.character(rep(start.day, length.out = n))
  ed_v <- if (is.null(end.day)) rep(NA_character_, n) else
    as.character(rep(end.day, length.out = n))

  # The fingerprint describes the QUERY SHAPE, not which sites are being
  # fetched on this call. A resume legitimately passes only the missing
  # subset, so env/lat/lon must NOT be hashed -- including them would make
  # every partial resume look like a different request. Per-site identity is
  # still protected: results are keyed by env.id, and the manifest stores the
  # full site list for inspection.
  fp <- .ck_fingerprint(kind = "weather",
                        start = unique(sd_v[!is.na(sd_v)]),
                        end = unique(ed_v[!is.na(ed_v)]),
                        vars = if (is.null(variables.names)) "default"
                        else sort(variables.names))

  paths <- .ck_init_run(dir.path, run.id, fp,
                        meta = list(kind = "weather", n = n, env.id = env.id,
                                    lat = as.numeric(lat), lon = as.numeric(lon),
                                    start.day = sd_v, end.day = ed_v,
                                    variables.names = variables.names %||% "default",
                                    sleep = sleep, tries = tries),
                        resume = resume)

  log0 <- .ck_latest(.ck_log_read(paths$log))
  skip <- if (isTRUE(resume) && nrow(log0))
    log0$env[log0$status == "ok"] else character(0)
  keys <- .ck_safe_key(env.id)

  if (verbose) {
    message(strrep("-", 63))
    message("get_weather_resumable() - checkpointed NASA POWER retrieval")
    message(strrep("-", 63))
    message("Run directory ...................... ", paths$root)
    message("Environments ....................... ", n)
    if (length(skip))
      message("Already complete (skipped) ......... ", length(skip))
  }

  fetch <- function(i) {
    a <- list(env.id = env.id[i], lat = lat[i], lon = lon[i],
              parallel = FALSE, verbose = FALSE, tries = tries)
    if (!is.na(sd_v[i])) a$start.day <- sd_v[i]
    if (!is.na(ed_v[i])) a$end.day   <- ed_v[i]
    if (!is.null(variables.names)) a$variables.names <- variables.names
    do.call(get_weather, a)
  }

  .ck_run_loop(keys, env.id, fetch, paths, sleep, verbose, skip)

  log1 <- .ck_latest(.ck_log_read(paths$log))
  out  <- .ck_assemble(paths, log1)
  if (is.null(out))
    stop("No environment completed successfully. Inspect the log:\n  ",
         paths$log, call. = FALSE)

  failed <- log1$env[log1$status != "ok"]
  if (verbose) {
    message("Completed .......................... ",
            sum(log1$status == "ok"), "/", n)
    if (length(failed)) {
      message("Failed ............................. ", length(failed))
      message("Resume with: restart_from_log(run.id = \"", run.id, "\")")
    }
    message("Done.")
  }

  attr(out, "run.dir")  <- paths$root
  attr(out, "n_ok")     <- sum(log1$status == "ok")
  attr(out, "n_failed") <- length(failed)
  attr(out, "failed")   <- failed
  out
}


# =========================================================================
# SECTION 10.2 - get_soil_resumable()
# =========================================================================

#' @title Crash-Resilient Wrapper for get_soil()
#'
#' @description
#' Downloads SoilGrids properties one environment at a time with the same
#' checkpointing scheme as \code{get_weather_resumable}: each result is
#' written to disk immediately and recorded in an append-only log, so a
#' crashed run resumes with \code{restart_from_log}.
#'
#' SoilGrids is rate-limited and slower than NASA POWER, which makes long soil
#' jobs the most frequent victims of a mid-run failure.
#'
#' @param env.id,lat,lon,variables.names,depths,stat Passed to
#'   \code{get_soil}.
#' @param dir.path character. Parent directory. Default \code{getwd()}.
#' @param run.id character. Run directory name. Default \code{"soil_run"}.
#' @param sleep numeric. Seconds between environments. Default 20, matching
#'   \code{get_soil}; lowering it much invites throttling by ISRIC.
#' @param tries integer. Retry attempts per environment. Default 3.
#' @param resume logical. Skip environments already logged \code{ok}. Default
#'   \code{TRUE}.
#' @param verbose logical. Print progress. Default \code{TRUE}.
#'
#' @details
#' Results are stored \strong{long} on disk, one row per property-depth, and
#' reshaped to wide only at assembly. Storing long is what makes resumption
#' safe: parts fetched at different times still stack correctly, whereas
#' pre-widened parts could have mismatched columns.
#'
#' See \code{get_weather_resumable} for the full checkpoint design.
#'
#' @return
#' A \code{data.frame} of soil properties (wide by default), with attributes
#' \code{run.dir}, \code{n_ok}, \code{n_failed} and \code{failed}.
#'
#' @seealso \code{restart_from_log}, \code{read_progress_log},
#'   \code{get_soil}
#'
#' @examples
#' \dontrun{
#' grid <- data.frame(env = paste0("P", 1:200),
#'                    lat = runif(200, -13, -12),
#'                    lon = runif(200, -56, -55))
#'
#' soil <- get_soil_resumable(
#'   env.id = grid$env, lat = grid$lat, lon = grid$lon,
#'   variables.names = c("clay", "sand", "silt", "soc", "phh2o"),
#'   run.id = "zoning_grid")
#'
#' ## After a crash, in a new session:
#' soil <- restart_from_log(run.id = "zoning_grid")
#'
#' ## Points in the ocean never resolve -- keep the successes
#' soil <- restart_from_log(run.id = "zoning_grid", retry.failed = FALSE)
#' }
#'
#' @export
get_soil_resumable <- function(env.id = NULL, lat = NULL, lon = NULL,
                               variables.names = "clay",
                               depths = c("0-5cm", "5-15cm", "15-30cm",
                                          "30-60cm", "60-100cm", "100-200cm"),
                               stat = c("mean", "q_5", "q_50", "q_95",
                                        "uncertainty"),
                               dir.path = NULL, run.id = "soil_run",
                               sleep = 20, tries = 3L,
                               resume = TRUE, verbose = TRUE) {

  .et_banner("get_soil_resumable",
             "downloads soil data with resumable progress logging", verbose)
  if (!exists("get_soil", mode = "function"))
    stop("get_soil() not found. Source env_data_collection.R first.",
         call. = FALSE)
  stat <- match.arg(stat)
  if (is.null(lat) || is.null(lon))
    stop("'lat' and 'lon' are required.", call. = FALSE)
  if (length(lat) != length(lon))
    stop("'lat' and 'lon' must have the same length.", call. = FALSE)

  n <- length(lat)
  if (is.null(env.id)) env.id <- paste0("env", seq_len(n))
  env.id <- as.character(env.id)
  if (length(env.id) != n)
    stop("'env.id' must have the same length as 'lat' and 'lon'.", call. = FALSE)
  if (anyDuplicated(env.id))
    stop("Duplicated env.id values. Checkpointing keys results on env.id, ",
         "so ids must be unique.", call. = FALSE)
  if (is.null(dir.path)) dir.path <- getwd()

  # Query shape only -- see the note in get_weather_resumable(). The site
  # list is deliberately excluded so partial resumes match.
  fp <- .ck_fingerprint(kind = "soil",
                        vars = sort(variables.names),
                        depths = sort(depths),
                        stat = stat)

  paths <- .ck_init_run(dir.path, run.id, fp,
                        meta = list(kind = "soil", n = n, env.id = env.id,
                                    lat = as.numeric(lat), lon = as.numeric(lon),
                                    variables.names = variables.names,
                                    depths = depths, stat = stat,
                                    sleep = sleep, tries = tries),
                        resume = resume)

  log0 <- .ck_latest(.ck_log_read(paths$log))
  skip <- if (isTRUE(resume) && nrow(log0))
    log0$env[log0$status == "ok"] else character(0)
  keys <- .ck_safe_key(env.id)

  if (verbose) {
    message(strrep("-", 63))
    message("get_soil_resumable() - checkpointed SoilGrids retrieval")
    message(strrep("-", 63))
    message("Run directory ...................... ", paths$root)
    message("Environments ....................... ", n)
    if (length(skip))
      message("Already complete (skipped) ......... ", length(skip))
  }

  # Long on disk: safe to stack across separate sessions.
  fetch <- function(i) {
    get_soil(env.id = env.id[i], lat = lat[i], lon = lon[i],
             variables.names = variables.names, depths = depths,
             stat = stat, wide = FALSE, sleep = 0, tries = tries,
             verbose = FALSE)
  }

  .ck_run_loop(keys, env.id, fetch, paths, sleep, verbose, skip)

  log1 <- .ck_latest(.ck_log_read(paths$log))
  long <- .ck_assemble(paths, log1)
  if (is.null(long))
    stop("No environment completed successfully. Inspect the log:\n  ",
         paths$log, call. = FALSE)

  out <- if (exists(".long_to_wide", mode = "function"))
    .long_to_wide(long, id.cols = c("env", "lon", "lat"),
                  key.col = "feature", value.col = stat) else long

  failed <- log1$env[log1$status != "ok"]
  if (verbose) {
    message("Completed .......................... ",
            sum(log1$status == "ok"), "/", n)
    if (length(failed)) {
      message("Failed ............................. ", length(failed))
      message("Resume with: restart_from_log(run.id = \"", run.id, "\")")
    }
    message("Done.")
  }

  attr(out, "run.dir")  <- paths$root
  attr(out, "n_ok")     <- sum(log1$status == "ok")
  attr(out, "n_failed") <- length(failed)
  attr(out, "failed")   <- failed
  out
}


# =========================================================================
# SECTION 10.3 - restart_from_log()
# =========================================================================

#' @title Resume an Interrupted get_weather() or get_soil() Run
#'
#' @description
#' Reads the manifest and progress log written by
#' \code{get_weather_resumable} or \code{get_soil_resumable},
#' skips every environment already retrieved, fetches those still missing or
#' failed, and returns the assembled result.
#'
#' Safe to call repeatedly. If everything already succeeded it fetches nothing
#' and simply assembles what is on disk.
#'
#' @param run.id character. The run directory name used originally.
#' @param dir.path character. Parent directory. Default \code{getwd()}.
#' @param retry.failed logical. Retry environments previously logged as
#'   failed. Default \code{TRUE}. Set \code{FALSE} to assemble only the
#'   successes -- useful when failures are permanent, such as coordinates in
#'   the ocean.
#' @param sleep numeric or \code{NULL}. Override the stored delay.
#'   \code{NULL} (default) reuses the manifest value.
#' @param tries integer or \code{NULL}. Override the stored retry count.
#' @param verbose logical. Print progress. Default \code{TRUE}.
#'
#' @details
#' All retrieval parameters come from \code{manifest.json}, so a resumed run
#' necessarily uses the same coordinates, dates and variables as the original.
#' You cannot accidentally resume with different settings -- that is what the
#' stored fingerprint prevents.
#'
#' \strong{Permanently failing environments.} Some failures never resolve: a
#' point in the ocean has no SoilGrids profile. Repeated calls will retry them
#' every time. Once you have confirmed a failure is permanent, pass
#' \code{retry.failed = FALSE}.
#'
#' @return
#' A \code{data.frame} as returned by the original function, with attributes
#' \code{run.dir}, \code{n_ok}, \code{n_failed} and \code{failed}.
#'
#' @seealso \code{get_weather_resumable},
#'   \code{get_soil_resumable}, \code{read_progress_log}
#'
#' @examples
#' \dontrun{
#' ## A long soil job dies at environment 137 of 200.
#' ## In a brand-new session -- no need to remember the arguments:
#' soil <- restart_from_log(run.id = "zoning_grid")
#'
#' ## Inspect what failed and why
#' lg <- read_progress_log(run.id = "zoning_grid")
#' subset(lg, status == "failed")[, c("env", "message")]
#'
#' ## Give up on permanent failures
#' soil <- restart_from_log(run.id = "zoning_grid", retry.failed = FALSE)
#'
#' ## Weather resumes identically
#' wth <- restart_from_log(run.id = "trial2024")
#' }
#'
#' @export
restart_from_log <- function(run.id, dir.path = NULL, retry.failed = TRUE,
                             sleep = NULL, tries = NULL, verbose = TRUE) {

  .et_banner("restart_from_log", "resumes a run from its progress log", verbose)
  if (missing(run.id) || !nzchar(run.id))
    stop("'run.id' is required.", call. = FALSE)
  if (is.null(dir.path)) dir.path <- getwd()

  root <- file.path(dir.path, run.id)
  manf <- file.path(root, "manifest.json")
  if (!dir.exists(root))
    stop("No run directory at\n  ", root,
         "\n  Check 'run.id' and 'dir.path'.", call. = FALSE)
  if (!file.exists(manf))
    stop("No manifest.json in\n  ", root,
         "\n  This directory was not created by a *_resumable() function.",
         call. = FALSE)

  man <- .ck_from_json(paste(readLines(manf, warn = FALSE), collapse = ""))
  kind <- as.character(man$kind)
  if (!kind %in% c("weather", "soil"))
    stop("Unrecognised run kind '", kind, "' in the manifest.", call. = FALSE)

  if (is.null(sleep)) sleep <- as.numeric(man$sleep)
  if (is.null(tries)) tries <- as.integer(man$tries)

  env.id <- as.character(man$env.id)
  lat <- as.numeric(man$lat); lon <- as.numeric(man$lon)

  log <- .ck_latest(.ck_log_read(file.path(root, "progress.log")))
  done <- if (nrow(log)) log$env[log$status == "ok"] else character(0)
  bad  <- if (nrow(log)) log$env[log$status != "ok"] else character(0)

  todo <- setdiff(env.id, done)
  if (!isTRUE(retry.failed)) todo <- setdiff(todo, bad)

  if (verbose) {
    message(strrep("-", 63))
    message("restart_from_log() - resuming '", run.id, "'")
    message(strrep("-", 63))
    message("Kind ............................... ", kind)
    message("Total environments ................. ", length(env.id))
    message("Already complete ................... ", length(done))
    message("Previously failed .................. ", length(bad),
            if (!isTRUE(retry.failed)) " (not retried)" else "")
    message("To fetch now ....................... ", length(todo))
  }

  if (!length(todo)) {
    if (verbose) message("Nothing to fetch; assembling saved parts.")
  } else {
    sel <- match(todo, env.id)
    args <- list(env.id = env.id[sel], lat = lat[sel], lon = lon[sel],
                 dir.path = dir.path, run.id = run.id, sleep = sleep,
                 tries = tries, resume = TRUE, verbose = verbose)
    if (kind == "weather") {
      sd_v <- as.character(man$start.day)[sel]
      ed_v <- as.character(man$end.day)[sel]
      if (!all(is.na(sd_v))) args$start.day <- sd_v
      if (!all(is.na(ed_v))) args$end.day   <- ed_v
      vn <- man$variables.names
      if (!identical(as.character(vn), "default")) args$variables.names <- vn
      suppressWarnings(invisible(do.call(get_weather_resumable, args)))
    } else {
      args$variables.names <- as.character(man$variables.names)
      args$depths <- as.character(man$depths)
      args$stat   <- as.character(man$stat)
      suppressWarnings(invisible(do.call(get_soil_resumable, args)))
    }
  }

  paths <- list(root = root, parts = file.path(root, "parts"),
                log = file.path(root, "progress.log"))
  log2 <- .ck_latest(.ck_log_read(paths$log))
  res  <- .ck_assemble(paths, log2)
  if (is.null(res))
    stop("Nothing could be assembled: no successful parts in\n  ",
         paths$parts, call. = FALSE)

  if (kind == "soil" && exists(".long_to_wide", mode = "function") &&
      "feature" %in% names(res))
    res <- .long_to_wide(res, id.cols = c("env", "lon", "lat"),
                         key.col = "feature",
                         value.col = as.character(man$stat))

  failed <- log2$env[log2$status != "ok"]
  if (verbose) {
    message("Assembled .......................... ",
            sum(log2$status == "ok"), "/", length(env.id), " environments")
    if (length(failed))
      message("Still failing ...................... ",
              paste(utils::head(failed, 5), collapse = ", "),
              if (length(failed) > 5) ", ..." else "")
    message("Done.")
  }

  attr(res, "run.dir")  <- root
  attr(res, "n_ok")     <- sum(log2$status == "ok")
  attr(res, "n_failed") <- length(failed)
  attr(res, "failed")   <- failed
  res
}


# =========================================================================
# SECTION 10.4 - read_progress_log()
# =========================================================================

#' @title Inspect the Progress Log of a Checkpointed Run
#'
#' @description
#' Returns the per-environment retrieval history written by
#' \code{get_weather_resumable} or \code{get_soil_resumable}:
#' what succeeded, what failed, when, and with which error message.
#'
#' @param run.id character. Run directory name.
#' @param dir.path character. Parent directory. Default \code{getwd()}.
#' @param all logical. If \code{TRUE}, return every record including
#'   superseded retries. If \code{FALSE} (default), the latest record per
#'   environment.
#' @param verbose boolean. If \code{TRUE} (default) prints a progress banner.
#'
#' @return
#' A \code{data.frame} with \code{env}, \code{key}, \code{status}
#' (\code{"ok"} or \code{"failed"}), \code{time} and \code{message}.
#'
#' @seealso \code{restart_from_log}
#'
#' @examples
#' \dontrun{
#' lg <- read_progress_log(run.id = "zoning_grid")
#' table(lg$status)
#'
#' ## Why did they fail?
#' subset(lg, status == "failed")[, c("env", "message")]
#'
#' ## Full history including retries
#' read_progress_log(run.id = "zoning_grid", all = TRUE)
#' }
#'
#' @export
read_progress_log <- function(run.id, dir.path = NULL, all = FALSE, verbose = TRUE) {
  .et_banner("read_progress_log", "reads a resumable-run progress log", verbose)
  if (is.null(dir.path)) dir.path <- getwd()
  p <- file.path(dir.path, run.id, "progress.log")
  if (!file.exists(p))
    stop("No progress.log at\n  ", p, call. = FALSE)
  lg <- .ck_log_read(p)
  if (isTRUE(all)) lg else .ck_latest(lg)
}


# =========================================================================
# SECTION - WorldClim static layers: elevation and bioclim
# =========================================================================

#' Shared downloader for a zipped WorldClim 2.1 base grid
#'
#' Base grids ship as .zip containing one .tif per variable/month/index.
#' Cached like the CMIP6 tiles: a present file is never re-downloaded.
#' @keywords internal
#' @noRd
.wc_base_grid <- function(var, resolution, dir.path, verbose = TRUE) {
  stem <- sprintf("wc2.1_%s_%s", resolution, var)
  zipf <- file.path(dir.path, paste0(stem, ".zip"))
  exdir <- file.path(dir.path, stem)

  if (!dir.exists(exdir) || !length(list.files(exdir, "\\.tif$"))) {
    if (!file.exists(zipf)) {
      url <- sprintf("https://geodata.ucdavis.edu/climate/worldclim/2_1/base/%s.zip",
                     stem)
      avail <- .cs_url_ok(url)
      if (identical(avail, FALSE))
        stop("Not published on the WorldClim server:\n  ", url,
             "\n  Check 'var' and 'resolution'. The server responded and the ",
             "file does not exist -- this is NOT a connectivity problem.",
             call. = FALSE)
      .et_step(paste0("downloading ", basename(url),
                      " (this is a large file)"), verbose)
      if (!.cs_fetch(url, zipf, tries = 2L, quiet = !verbose))
        stop("Could not download:\n  ", url,
             "\n  If the host is unreachable (proxy, firewall, offline), ",
             "place the .zip in 'dir.path' manually and re-run.",
             call. = FALSE)
    }
    .et_step("extracting", verbose)
    dir.create(exdir, showWarnings = FALSE, recursive = TRUE)
    utils::unzip(zipf, exdir = exdir)
  }
  tifs <- list.files(exdir, "\\.tif$", full.names = TRUE)
  if (!length(tifs))
    stop("No .tif found after extracting ", zipf, ".", call. = FALSE)
  tifs
}

#' @title Elevation at Trial Sites from WorldClim 2.1
#'
#' @description
#' Returns elevation (m) for each site. Elevation is a strong covariate for
#' temperature and radiation, so it is often worth carrying alongside the
#' weather-derived envirotypes -- but note it is STATIC: it adds
#' between-site structure, never within-season dynamics.
#'
#' @param env.id character. Environment identifiers.
#' @param lat,lon numeric. Site coordinates (decimal degrees).
#' @param resolution character. \code{"10m"}, \code{"5m"}, \code{"2.5m"} or
#'   \code{"30s"}. Default \code{"2.5m"}. 30s is ~9.7 GB -- use only if you
#'   genuinely need ~1 km cells.
#' @param dir.path character. Cache directory. Default is the same cache used
#'   by \code{get_climate_scenario()}.
#' @param verbose boolean.
#'
#' @return data.frame with \code{env}, \code{LAT}, \code{LON},
#'   \code{ALT} (m).
#'
#' @details
#' Sites in the sea or outside the land mask return \code{NA} with a warning
#' naming them, rather than silently propagating NA into downstream kernels.
#'
#' @examples
#' \dontrun{
#' sites <- data.frame(env = c("PIRA", "SETE"),
#'                     lat = c(-22.7, -19.4), lon = c(-47.6, -44.2))
#' alt <- get_elevation(env.id = sites$env, lat = sites$lat, lon = sites$lon)
#' }
#'
#' @references
#' Fick, S.E. & Hijmans, R.J. (2017) WorldClim 2: new 1-km spatial resolution
#' climate surfaces for global land areas. \emph{International Journal of
#' Climatology} 37(12), 4302-4315. \doi{10.1002/joc.5086}
#'
#' @export
get_elevation <- function(env.id = NULL, lat = NULL, lon = NULL,
                          resolution = c("2.5m", "5m", "10m", "30s"),
                          dir.path = NULL, verbose = TRUE) {
  .et_banner("get_elevation", "WorldClim 2.1 elevation", verbose)
  .need_pkg(c("terra"), "get_elevation")
  resolution <- match.arg(resolution)
  if (is.null(lat) || is.null(lon))
    stop("'lat' and 'lon' are required.", call. = FALSE)
  if (length(lat) != length(lon))
    stop("'lat' and 'lon' must be the same length.", call. = FALSE)
  if (is.null(env.id)) env.id <- paste0("env_", seq_along(lat))
  if (length(env.id) != length(lat))
    stop("'env.id' must match the length of 'lat'/'lon'.", call. = FALSE)
  if (is.null(dir.path))
    dir.path <- file.path(tools::R_user_dir("et", "cache"))
  dir.create(dir.path, showWarnings = FALSE, recursive = TRUE)

  tif <- .wc_base_grid("elev", resolution, dir.path, verbose)[1]
  r <- terra::rast(tif)
  v <- .cs_extract(r, lon, lat, expect = 1L)

  out <- data.frame(env = as.character(env.id), LAT = lat, LON = lon,
                    ALT = as.numeric(v[, 1]), stringsAsFactors = FALSE)
  bad <- out$env[is.na(out$ALT)]
  if (length(bad))
    warning(length(bad), " site(s) fell outside the land mask and have ",
            "ALT = NA: ", paste(utils::head(bad, 8), collapse = ", "),
            if (length(bad) > 8) ", ..." else "", ".", call. = FALSE)
  .et_step(paste0("extracted elevation for ", nrow(out), " site(s)"), verbose)
  out
}

#' @title WorldClim Bioclimatic Variables at Trial Sites
#'
#' @description
#' Returns the 19 standard bioclimatic indices (BIO1-BIO19) for each site,
#' from the 1970-2000 WorldClim 2.1 normals.
#'
#' @section What these are, and are not:
#' Bioclim variables are LONG-TERM CLIMATE NORMALS, not the weather of your
#' trial year. They describe the climate a site is drawn from, so they are
#' appropriate for characterising a TPE, clustering locations, or as static
#' site covariates. They are NOT a substitute for
#' \code{get_weather()} + \code{env_phenology()} when you want to explain
#' what happened in a particular season: a site's BIO1 is identical in a
#' drought year and a wet year. Mixing normals and realised weather in one
#' kernel without saying which is which is a common and serious error.
#'
#' @param env.id character. Environment identifiers.
#' @param lat,lon numeric. Site coordinates (decimal degrees).
#' @param vars integer or character. Which indices, e.g. \code{c(1, 12)} or
#'   \code{c("BIO1", "BIO12")}. Default all 19.
#' @param resolution character. \code{"10m"}, \code{"5m"}, \code{"2.5m"} or
#'   \code{"30s"}. Default \code{"2.5m"}.
#' @param dir.path character. Cache directory.
#' @param verbose boolean.
#'
#' @return data.frame with \code{env}, \code{LAT}, \code{LON} and one column
#'   per requested index. Units follow WorldClim: temperature indices in C,
#'   precipitation in mm, BIO3/BIO4/BIO15 dimensionless or percent.
#'
#' @examples
#' \dontrun{
#' sites <- data.frame(env = c("PIRA", "SETE"),
#'                     lat = c(-22.7, -19.4), lon = c(-47.6, -44.2))
#' bc  <- get_bioclim(env.id = sites$env, lat = sites$lat, lon = sites$lon)
#' bc2 <- get_bioclim(sites$env, sites$lat, sites$lon, vars = c(1, 12))
#' }
#'
#' @references
#' Fick, S.E. & Hijmans, R.J. (2017) WorldClim 2: new 1-km spatial resolution
#' climate surfaces for global land areas. \emph{International Journal of
#' Climatology} 37(12), 4302-4315. \doi{10.1002/joc.5086}
#'
#' O'Donnell, M.S. & Ignizio, D.A. (2012) Bioclimatic predictors for
#' supporting ecological applications in the conterminous United States.
#' \emph{U.S. Geological Survey Data Series} 691.
#'
#' @export
get_bioclim <- function(env.id = NULL, lat = NULL, lon = NULL, vars = 1:19,
                        resolution = c("2.5m", "5m", "10m", "30s"),
                        dir.path = NULL, verbose = TRUE) {
  .et_banner("get_bioclim", "WorldClim 2.1 bioclimatic variables", verbose)
  .need_pkg(c("terra"), "get_bioclim")
  resolution <- match.arg(resolution)
  if (is.null(lat) || is.null(lon))
    stop("'lat' and 'lon' are required.", call. = FALSE)
  if (length(lat) != length(lon))
    stop("'lat' and 'lon' must be the same length.", call. = FALSE)
  if (is.null(env.id)) env.id <- paste0("env_", seq_along(lat))
  if (is.character(vars)) vars <- as.integer(sub("^[Bb][Ii][Oo]_?", "", vars))
  if (anyNA(vars) || any(vars < 1 | vars > 19))
    stop("'vars' must be indices 1-19 (or names like 'BIO1').", call. = FALSE)
  vars <- sort(unique(as.integer(vars)))
  if (is.null(dir.path))
    dir.path <- file.path(tools::R_user_dir("et", "cache"))
  dir.create(dir.path, showWarnings = FALSE, recursive = TRUE)

  tifs <- .wc_base_grid("bio", resolution, dir.path, verbose)
  ## Sort by the trailing index so BIO2 precedes BIO10 (lexical sort would not)
  idx <- as.integer(sub(".*_(\\d+)\\.tif$", "\\1", basename(tifs)))
  tifs <- tifs[order(idx)]; idx <- sort(idx)
  keep <- match(vars, idx)
  if (anyNA(keep))
    stop("Missing bioclim layer(s) in the archive: ",
         paste(vars[is.na(keep)], collapse = ", "), ".", call. = FALSE)

  r <- terra::rast(tifs[keep])
  v <- .cs_extract(r, lon, lat, expect = length(vars))
  colnames(v) <- paste0("BIO", vars)

  out <- data.frame(env = as.character(env.id), LAT = lat, LON = lon,
                    v, stringsAsFactors = FALSE, check.names = FALSE)
  bad <- out$env[apply(is.na(v), 1, all)]
  if (length(bad))
    warning(length(bad), " site(s) fell outside the land mask and are all ",
            "NA: ", paste(utils::head(bad, 8), collapse = ", "),
            if (length(bad) > 8) ", ..." else "", ".", call. = FALSE)
  .et_step(paste0("extracted ", length(vars), " index/indices for ",
                  nrow(out), " site(s)"), verbose)
  out
}


# =========================================================================
# SECTION - get_climate_scenario() (CMIP6 downscaled future weather)
# =========================================================================

#' Fetch one WorldClim CMIP6 GeoTIFF into the cache
#'
#' Isolated in its own helper for two reasons: it is the ONLY point in this
#' module that touches the network, and keeping it separate makes
#' get_climate_scenario() testable offline -- a test can replace this
#' function with a stub, or pre-populate the cache so it is never called.
#'
#' @param url character. Remote GeoTIFF URL.
#' @param dest character. Destination path in the cache directory.
#' @param tries integer. Retry attempts, passed to .retry().
#' @param quiet boolean. Suppress download progress.
#' @return \code{TRUE} on success, \code{FALSE} on any failure.
#' Known gaps in the WorldClim CMIP6 downscaled archive
#'
#' Source: worldclim.org/data/cmip6/ coverage tables, verified 2026-09-28.
#' Identical across all four resolutions and all four periods, because these
#' are gaps in the downscaled source rather than tiling artifacts.
#'
#' Each entry lists SSPs that are unusable for delta-change downscaling, which
#' needs tmin AND tmax AND prec. "pr only" counts as unusable: GFDL-ESM4
#' ssp585 ships precipitation but no temperature, so it would fail partway
#' through after a wasted download.
#' @keywords internal
#' @noRd
.CS_GAPS <- list(
  "GFDL-ESM4"       = c("ssp245", "ssp585"),  # ssp245 absent; ssp585 is pr only
  "FIO-ESM-2-0"     = c("ssp370"),
  "HadGEM3-GC31-LL" = c("ssp370")
)

#' Models with complete tn/tx/pr coverage for a given SSP
#' @keywords internal
#' @noRd
.cs_models_for <- function(scenario, pool = .CS_GCMS) {
  bad <- names(.CS_GAPS)[vapply(.CS_GAPS, function(s) scenario %in% s, logical(1))]
  setdiff(pool, bad)
}

#' Is a remote URL actually present? HEAD request; no body transfer.
#'
#' Three-way return, deliberately: TRUE = present, FALSE = definitively absent
#' (4xx), NA = could not tell (timeout, DNS, refused). "Not published" and
#' "cannot reach the host" require different remedies and must never be
#' collapsed into one message.
#' @keywords internal
#' @noRd
.cs_url_ok <- function(url, timeout = 20) {
  if (!requireNamespace("curl", quietly = TRUE)) return(NA)
  tryCatch({
    h <- curl::new_handle(nobody = TRUE, timeout = timeout,
                          followlocation = TRUE)
    code <- curl::curl_fetch_memory(url, handle = h)$status_code
    if (code >= 200 && code < 300) TRUE
    else if (code >= 400 && code < 500) FALSE
    else NA
  }, error = function(e) NA)
}

#' Scenarios published for a resolution x model, from the Apache index.
#' Returns NA_character_ if the listing could not be read.
#' @keywords internal
#' @noRd
.cs_list_scenarios <- function(resolution, model) {
  u <- sprintf("https://geodata.ucdavis.edu/cmip6/%s/%s/", resolution, model)
  txt <- tryCatch(paste(readLines(u, warn = FALSE), collapse = "\n"),
                  error = function(e) NA_character_)
  if (is.na(txt)) return(NA_character_)
  sort(unique(regmatches(txt, gregexpr("ssp[0-9]{3}", txt))[[1]]))
}

#' @keywords internal
#' @noRd
.cs_fetch <- function(url, dest, tries = 2L, quiet = TRUE) {
  ok <- tryCatch({
    .retry(utils::download.file(url, dest, mode = "wb", quiet = quiet),
           tries = tries, quiet = quiet)
    TRUE
  }, error = function(e) FALSE)
  if (!isTRUE(ok) || !file.exists(dest)) {
    if (file.exists(dest)) unlink(dest)
    return(FALSE)
  }
  TRUE
}

#' Bilinear-free nearest-cell extraction from a terra raster
#' @keywords internal
#' @noRd
.cs_extract <- function(r, lon, lat, expect = NULL) {
  pts <- cbind(lon, lat)
  v <- terra::extract(r, pts)

  ## terra returns different shapes depending on input type and version:
  ##   extract(r, matrix)      -> nlyr columns, NO ID
  ##   extract(r, SpatVector)  -> 1 + nlyr columns, first named "ID"
  ## Stripping column 1 unconditionally deletes a real data column (this bug
  ## silently dropped January and reported "got 11"). Detect ID by NAME.
  if (is.data.frame(v)) {
    if (!is.null(names(v)) && identical(names(v)[1], "ID"))
      v <- v[, -1, drop = FALSE]
    v <- as.matrix(v)
  } else {
    v <- as.matrix(v)
  }

  ## Some WorldClim tiles carry a stray extra band. If we know how many
  ## layers we need and there is exactly one spare, drop it from the END
  ## rather than guessing at the front -- month order is Jan..Dec.
  if (!is.null(expect) && ncol(v) == expect + 1L &&
      all(is.na(v[, ncol(v)]))) {
    v <- v[, seq_len(expect), drop = FALSE]
  }
  v
}

#' Valid CMIP6 GCM / SSP / period combinations served by WorldClim 2.1
#' @keywords internal
#' @noRd
.CS_GCMS <- c("ACCESS-CM2", "BCC-CSM2-MR", "CMCC-ESM2", "EC-Earth3-Veg",
              "FIO-ESM-2-0", "GFDL-ESM4", "GISS-E2-1-G", "HadGEM3-GC31-LL",
              "INM-CM5-0", "IPSL-CM6A-LR", "MIROC6", "MPI-ESM1-2-HR",
              "MRI-ESM2-0", "UKESM1-0-LL")
.CS_SSPS    <- c("ssp126", "ssp245", "ssp370", "ssp585")
.CS_PERIODS <- c("2021-2040", "2041-2060", "2061-2080", "2081-2100")

#' @title Future-Climate Envirotyping from CMIP6 Downscaled Scenarios
#'
#' @description
#' Retrieves downscaled CMIP6 monthly climate (WorldClim 2.1) for a set of
#' sites and converts it into a DAILY weather series compatible with
#' \code{processWTH()}, \code{env_phenology()} and \code{water_balance()}, by
#' applying the monthly change signal to an observed baseline series
#' (delta-change downscaling).
#'
#' @author Germano Costa Neto
#'
#' @param env.id vector (character). Site identifiers.
#' @param lat,lon vector (numeric). Site coordinates.
#' @param baseline data.frame. Observed daily weather for the same sites,
#'   from \code{get_weather()}. Supplies the day-to-day variability that
#'   monthly projections cannot. REQUIRED for \code{method = "delta"}.
#' @param model character. CMIP6 GCM name, e.g. \code{"MPI-ESM1-2-HR"}. One of
#'   \code{.CS_GCMS}. Use \code{"ensemble"} to average several.
#' @param models vector (character). Used when \code{model = "ensemble"}.
#'   Default is a four-model spread across climate sensitivities.
#' @param scenario character. \code{"ssp126"}, \code{"ssp245"},
#'   \code{"ssp370"} or \code{"ssp585"}.
#' @param period character. \code{"2021-2040"}, \code{"2041-2060"},
#'   \code{"2061-2080"} or \code{"2081-2100"}.
#' @param resolution character. WorldClim tile resolution: \code{"2.5m"}
#'   (default), \code{"5m"}, \code{"10m"} or \code{"30s"}. Finer is a much
#'   larger download. Resolution does NOT affect which model x scenario
#'   combinations exist: archive gaps are per GCM x SSP and are identical at
#'   every resolution and period. See \code{.CS_GAPS}.
#' @param method character. \code{"delta"} (default) perturbs the baseline;
#'   \code{"raw"} returns the monthly projected values without downscaling.
#' @param dir.path character. Cache directory for downloaded rasters.
#'   Defaults to \code{tools::R_user_dir("et", "cache")}.
#' @param baseline.id,baseline.date character. Columns in \code{baseline}.
#' @param precip.scale character. \code{"multiplicative"} (default) or
#'   \code{"additive"} application of the precipitation delta.
#' @param verbose boolean. Print progress.
#'
#' @return
#' For \code{method = "delta"}, a data.frame in \code{get_weather()} layout
#' (\code{env}, \code{YYYYMMDD}, \code{daysFromStart}, \code{T2M},
#' \code{T2M_MAX}, \code{T2M_MIN}, \code{PRECTOT}, ...) with the scenario
#' signal applied, plus \code{scenario}, \code{model} and \code{period}
#' columns. For \code{method = "raw"}, a monthly data.frame of projected
#' \code{tmin}, \code{tmax} and \code{prec} per site.
#' The attribute \code{"deltas"} holds the monthly change factors that were
#' applied, so the perturbation is auditable.
#'
#' @details
#' \strong{What delta-change does and does not do.} The monthly mean change
#' from the GCM is added to (temperature) or multiplied onto (precipitation)
#' the observed daily series. The result therefore keeps the OBSERVED
#' day-to-day sequence, wet-day frequency and spell structure, and shifts only
#' the monthly means. This is the standard, defensible way to drive a crop
#' model from monthly projections (Hawkins et al. 2013), but note the
#' consequence: \strong{changes in variability, in dry-spell length, and in
#' extreme-event frequency are NOT represented}. A scenario that in reality
#' brings the same total rain in fewer, heavier events will look identical to
#' one that brings it evenly. Heat-stress day counts under delta-change are
#' driven purely by the mean shift, so they are a lower bound.
#'
#' If extremes matter to your question, you need daily GCM output and a
#' quantile-mapping bias correction, not this function.
#'
#' \strong{Model choice is not innocuous.} Across CMIP6, equilibrium climate
#' sensitivity ranges roughly 1.8-5.6 K, and several high-sensitivity models
#' are now regarded as implausibly hot. A single-model result is one draw from
#' a wide distribution. Prefer \code{model = "ensemble"} and report the spread;
#' a projection whose sign flips across models is not a finding.
#'
#' \strong{Network access.} WorldClim is an external host. In restricted
#' environments the download will fail; pre-populate \code{dir.path} with the
#' GeoTIFFs and the function will use the cache.
#'
#' @examples
#' \dontrun{
#' sites <- data.frame(env = c("PIRA", "SETE"),
#'                     lat = c(-22.7, -19.4), lon = c(-47.6, -44.2))
#'
#' base <- get_weather(env.id = sites$env, lat = sites$lat, lon = sites$lon,
#'                     start.day = "2015-10-01", end.day = "2016-03-31")
#'
#' ## Mid-century, intermediate scenario, 4-model ensemble.
#' ## The default 'models' set covers all four SSPs.
#' fut <- get_climate_scenario(env.id = sites$env, lat = sites$lat,
#'                             lon = sites$lon, baseline = base,
#'                             model = "ensemble", scenario = "ssp245",
#'                             period = "2041-2060")
#'
#' ## Which GCMs fully cover a scenario? (GFDL-ESM4 has no ssp245.)
#' .cs_models_for("ssp245")
#'
#' ## GFDL-ESM4 is fine for ssp126 or ssp370, just not ssp245/ssp585:
#' fut126 <- get_climate_scenario(env.id = sites$env, lat = sites$lat,
#'                                lon = sites$lon, baseline = base,
#'                                model = "GFDL-ESM4", scenario = "ssp126",
#'                                period = "2041-2060")
#'
#' ## Drop straight into the existing pipeline
#' fut <- processWTH(fut)
#' ph  <- env_phenology(fut, crop = "maize")
#'
#' ## Compare thermal time now vs mid-century
#' tapply(ph$dTT, ph$env, sum)
#' }
#'
#' @seealso \code{\link{get_weather}}, \code{\link{project_risk}},
#'   \code{\link{env_phenology}}
#'
#' @references
#' Fick, S.E. & Hijmans, R.J. (2017) WorldClim 2. \emph{Int. J. Climatol.}
#' 37(12), 4302-4315.
#'
#' Hawkins, E. et al. (2013) Calibration and bias correction of climate
#' projections for crop modelling. \emph{Agric. For. Meteorol.} 170, 19-31.
#'
#' O'Neill, B.C. et al. (2016) The Scenario Model Intercomparison Project
#' (ScenarioMIP) for CMIP6. \emph{Geosci. Model Dev.} 9, 3461-3482.
#'
#' @export
get_climate_scenario <- function(env.id = NULL, lat = NULL, lon = NULL,
                                 baseline = NULL,
                                 model = "MPI-ESM1-2-HR",
                                 ## Spread across climate sensitivities using
                                 ## only GCMs with complete tn/tx/pr coverage
                                 ## for ALL four SSPs. GFDL-ESM4 is deliberately
                                 ## excluded: it has no ssp245 and ssp585 is
                                 ## precipitation-only (verified 2026-09-28).
                                 models = c("IPSL-CM6A-LR", "MPI-ESM1-2-HR",
                                            "MRI-ESM2-0", "UKESM1-0-LL"),
                                 scenario = "ssp245",
                                 period = "2041-2060",
                                 ## Coverage gaps are per GCM x SSP and do NOT
                                 ## vary with resolution, so this is purely a
                                 ## precision/download-size choice.
                                 resolution = c("2.5m", "5m", "10m", "30s"),
                                 method = c("delta", "raw"),
                                 dir.path = NULL,
                                 baseline.id = "env", baseline.date = "YYYYMMDD",
                                 precip.scale = c("multiplicative", "additive"),
                                 verbose = TRUE) {

  .et_banner("get_climate_scenario", "CMIP6 downscaled scenario weather", verbose)
  .need_pkg(c("terra"), "get_climate_scenario")
  method       <- match.arg(method)
  resolution   <- match.arg(resolution)
  precip.scale <- match.arg(precip.scale)

  if (is.null(lat) || is.null(lon))
    stop("'lat' and 'lon' are required.", call. = FALSE)
  if (length(lat) != length(lon))
    stop("'lat' and 'lon' must have the same length.", call. = FALSE)
  if (any(lat < -90 | lat > 90, na.rm = TRUE) ||
      any(lon < -180 | lon > 180, na.rm = TRUE))
    stop("Coordinates out of range; check lat/lon are not swapped.", call. = FALSE)
  n_env <- length(lat)
  if (is.null(env.id)) env.id <- paste0("env", seq_len(n_env))
  env.id <- as.character(env.id)

  scenario <- match.arg(scenario, .CS_SSPS)
  period   <- match.arg(period, .CS_PERIODS)
  ## Ensemble membership must be fixed and stated. An unpublished member
  ## raises an error above rather than being silently dropped -- otherwise the
  ## same call at two resolutions would average different model sets and
  ## report neither, which is not reproducible.
  use_models <- if (identical(model, "ensemble")) models else model

  ## ---- fail fast on known archive gaps -----------------------------------
  ## Checked BEFORE any download. Coverage holes are per GCM x SSP and are
  ## identical at every resolution, so switching resolution cannot fix them.
  gaps <- use_models[vapply(use_models, function(g)
    scenario %in% .CS_GAPS[[g]], logical(1))]
  if (length(gaps)) {
    alt <- .cs_models_for(scenario)
    stop("The WorldClim CMIP6 archive has no usable ", scenario,
         " data for: ", paste(gaps, collapse = ", "), ".\n",
         "  These are gaps in the downscaled source -- identical at every ",
         "resolution and period, so changing 'resolution' will not help.\n",
         "  (GFDL-ESM4 ssp585 exists but is precipitation-only, with no ",
         "tmin/tmax, so it cannot drive delta-change downscaling.)\n",
         "  GCMs with complete coverage for ", scenario, ": ",
         paste(alt, collapse = ", "), ".",
         call. = FALSE)
  }
  bad <- setdiff(use_models, .CS_GCMS)
  if (length(bad))
    stop("Unknown GCM(s): ", paste(bad, collapse = ", "),
         ".\n  Available: ", paste(.CS_GCMS, collapse = ", "), call. = FALSE)

  if (is.null(dir.path)) {
    dir.path <- tryCatch(tools::R_user_dir("et", "cache"),
                         error = function(e) file.path(tempdir(), "et_cache"))
  }
  dir.create(dir.path, recursive = TRUE, showWarnings = FALSE)

  if (verbose) {
    message(strrep("-", 63))
    message("CMIP6 downscaled scenario (WorldClim 2.1)")
    message("Scenario ........... ", scenario)
    message("Period ............. ", period)
    message("Model(s) ........... ", paste(use_models, collapse = ", "))
    message("Resolution ......... ", resolution)
    message("Cache .............. ", dir.path)
    message(strrep("-", 63))
  }

  ## ---- fetch monthly future climate per model ----------------------------
  base_url <- "https://geodata.ucdavis.edu/cmip6"
  vars <- c("tmin", "tmax", "prec")
  fut <- array(NA_real_, dim = c(n_env, 12L, length(vars), length(use_models)),
               dimnames = list(env.id, month.abb, vars, use_models))

  for (m in seq_along(use_models)) {
    gcm <- use_models[m]
    for (v in vars) {
      fn <- sprintf("wc2.1_%s_%s_%s_%s_%s.tif",
                    resolution, v, gcm, scenario, period)
      dest <- file.path(dir.path, fn)
      if (!file.exists(dest)) {
        url <- sprintf("%s/%s/%s/%s/%s", base_url, resolution, gcm, scenario, fn)

        ## Preflight. Coverage is ragged: verified 2026-09-28,
        ## /cmip6/2.5m/GFDL-ESM4/ publishes ssp126, ssp370 and ssp585 but NOT
        ## ssp245, while /cmip6/10m/MIROC6/ssp245/ is fully populated. There is
        ## no complete matrix to hardcode, so ask the server and fail with a
        ## message that names the real alternatives.
        avail <- .cs_url_ok(url)
        if (identical(avail, FALSE)) {
          sc <- .cs_list_scenarios(resolution, gcm)
          hint <- if (!identical(sc, NA_character_) && length(sc))
            paste0("\n  Published for ", gcm, " at ", resolution, ": ",
                   paste(sc, collapse = ", "), ".")
          else ""
          stop("Not published on the WorldClim server:\n  ",
               "model = '", gcm, "', scenario = '", scenario,
               "', period = '", period, "', resolution = '", resolution, "'.",
               hint,
               "\n  The server responded and the file genuinely does not ",
               "exist -- this is NOT a connectivity problem, so retrying or ",
               "changing network will not help.",
               "\n  Pick a published scenario, another model, or a different ",
               "resolution; coverage varies across all three.",
               call. = FALSE)
        }

        .et_step(paste0("downloading ", fn), verbose)
        if (!.cs_fetch(url, dest, tries = 2L, quiet = !verbose))
          stop("Could not download:\n  ", url,
               "\n  The file is published but the transfer failed. If the host ",
               "is unreachable from this machine (proxy, firewall, offline), ",
               "pre-populate 'dir.path' with the WorldClim CMIP6 GeoTIFFs and ",
               "re-run; the cache is used whenever the file already exists.",
               call. = FALSE)
      }
      r <- terra::rast(dest)
      vals <- .cs_extract(r, lon, lat, expect = 12L)
      if (ncol(vals) != 12L)
        stop("Expected 12 monthly layers in ", fn, "; got ", ncol(vals),
             ".\n  The raster on disk has ", terra::nlyr(r), " layer(s). ",
             "If that is not 12 the cached file is truncated or corrupt -- ",
             "delete it and re-run:\n    unlink(\"", dest, "\")",
             call. = FALSE)
      fut[, , v, m] <- as.matrix(vals)
    }
  }

  # ensemble mean across models
  futm <- apply(fut, c(1, 2, 3), mean, na.rm = TRUE)
  if (anyNA(futm)) {
    nas <- unique(env.id[apply(is.na(futm), 1, any)])
    warning("No scenario data at: ", paste(nas, collapse = ", "),
            ". Sites over water or outside the raster return NA.", call. = FALSE)
  }

  if (method == "raw") {
    out <- do.call(rbind, lapply(seq_len(n_env), function(i)
      data.frame(env = env.id[i], month = 1:12,
                 tmin = futm[i, , "tmin"], tmax = futm[i, , "tmax"],
                 prec = futm[i, , "prec"],
                 scenario = scenario, model = paste(use_models, collapse = "+"),
                 period = period, stringsAsFactors = FALSE)))
    rownames(out) <- NULL
    attr(out, "models") <- use_models
    return(out)
  }

  ## ---- delta-change downscaling -----------------------------------------
  if (is.null(baseline))
    stop("method = 'delta' requires 'baseline', an observed daily series from ",
         "get_weather() for the same sites.", call. = FALSE)
  baseline <- as.data.frame(baseline)
  if (!baseline.id %in% names(baseline))
    stop("Column '", baseline.id, "' not found in 'baseline'.", call. = FALSE)
  if (!baseline.date %in% names(baseline))
    stop("Column '", baseline.date, "' not found in 'baseline'.", call. = FALSE)

  bdate <- as.Date(baseline[[baseline.date]])
  bmon  <- as.numeric(format(bdate, "%m"))
  bsite <- as.character(baseline[[baseline.id]])

  need <- c("T2M_MAX", "T2M_MIN", "PRECTOT")
  miss <- setdiff(need, names(baseline))
  if (length(miss))
    stop("'baseline' is missing: ", paste(miss, collapse = ", "), call. = FALSE)

  unknown <- setdiff(unique(bsite), env.id)
  if (length(unknown))
    stop("Baseline contains site(s) with no coordinates supplied: ",
         paste(unknown, collapse = ", "), call. = FALSE)

  .et_step("computing monthly baseline climatology", verbose)
  btmax <- .strip_power_na(baseline$T2M_MAX)
  btmin <- .strip_power_na(baseline$T2M_MIN)
  bprec <- .strip_power_na(baseline$PRECTOT)

  key <- paste(bsite, bmon, sep = "\r")
  cl_tmax <- tapply(btmax, key, mean, na.rm = TRUE)
  cl_tmin <- tapply(btmin, key, mean, na.rm = TRUE)
  # monthly TOTAL precipitation averaged over years, to match WorldClim units
  yrs <- as.numeric(format(bdate, "%Y"))
  key_y <- paste(bsite, bmon, yrs, sep = "\r")
  mtot <- tapply(bprec, key_y, sum, na.rm = TRUE)
  pk <- do.call(rbind, strsplit(names(mtot), "\r", fixed = TRUE))
  cl_prec <- tapply(as.numeric(mtot), paste(pk[, 1], pk[, 2], sep = "\r"),
                    mean, na.rm = TRUE)

  .et_step("applying monthly delta-change signal", verbose)
  d_tmax <- d_tmin <- d_prec <- rep(NA_real_, nrow(baseline))
  delta_tab <- list()

  for (i in seq_len(n_env)) {
    s <- env.id[i]
    for (mo in 1:12) {
      k <- paste(s, mo, sep = "\r")
      idx <- which(bsite == s & bmon == mo)
      if (!length(idx)) next
      obs_tmax <- cl_tmax[[k]]; obs_tmin <- cl_tmin[[k]]
      obs_prec <- if (k %in% names(cl_prec)) cl_prec[[k]] else NA_real_

      dt_max <- futm[i, mo, "tmax"] - obs_tmax
      dt_min <- futm[i, mo, "tmin"] - obs_tmin
      fp <- futm[i, mo, "prec"]
      dp <- if (precip.scale == "multiplicative") {
        if (is.finite(obs_prec) && obs_prec > 1) fp / obs_prec else 1
      } else {
        if (is.finite(obs_prec)) (fp - obs_prec) / max(length(idx), 1) else 0
      }
      # Guard against runaway multiplicative factors from a near-dry baseline
      if (precip.scale == "multiplicative" && is.finite(dp))
        dp <- min(max(dp, 0), 5)

      d_tmax[idx] <- dt_max; d_tmin[idx] <- dt_min; d_prec[idx] <- dp
      delta_tab[[length(delta_tab) + 1L]] <-
        data.frame(env = s, month = mo, d_tmax = dt_max, d_tmin = dt_min,
                   d_prec = dp, obs_prec = obs_prec, fut_prec = fp,
                   stringsAsFactors = FALSE)
    }
  }

  out <- baseline
  out$T2M_MAX <- btmax + d_tmax
  out$T2M_MIN <- btmin + d_tmin
  out$T2M <- (out$T2M_MAX + out$T2M_MIN) / 2
  out$PRECTOT <- if (precip.scale == "multiplicative") {
    bprec * d_prec
  } else pmax(bprec + d_prec, 0)

  # a warmer atmosphere holds more water: keep dewpoint physically consistent
  # by preserving relative humidity rather than absolute dewpoint.
  if ("T2MDEW" %in% names(out)) {
    bdew <- .strip_power_na(baseline$T2MDEW)
    out$T2MDEW <- bdew + (d_tmax + d_tmin) / 2
  }

  inv <- out$T2M_MAX < out$T2M_MIN
  if (any(inv, na.rm = TRUE))
    warning(sum(inv, na.rm = TRUE), " day(s) ended with Tmax < Tmin after the ",
            "delta shift; tmin and tmax deltas differ in sign there. ",
            "Inspect attr(x, 'deltas').", call. = FALSE)

  out$scenario <- scenario
  out$model    <- paste(use_models, collapse = "+")
  out$period   <- period

  # drop stale derived columns: they describe the baseline, not the scenario
  stale <- intersect(c("n", "N", "RTA", "VPD", "SPV", "ETP", "PETP",
                       "GDD", "FRUE", "T2M_RANGE", "dTT", "cumTT",
                       "stage", "dstage", "stage_day"), names(out))
  if (length(stale)) {
    if (verbose)
      message("Dropping baseline-derived column(s) invalidated by the shift: ",
              paste(stale, collapse = ", "), ". Re-run processWTH().")
    out <- out[, setdiff(names(out), stale), drop = FALSE]
  }

  attr(out, "deltas") <- do.call(rbind, delta_tab)
  attr(out, "scenario") <- list(models = use_models, scenario = scenario,
                                period = period, method = method,
                                precip.scale = precip.scale)

  if (verbose) {
    dd <- attr(out, "deltas")
    message(strrep("-", 63))
    message("Mean applied signal across sites and months")
    message(sprintf("  Tmax  %+5.2f C", mean(dd$d_tmax, na.rm = TRUE)))
    message(sprintf("  Tmin  %+5.2f C", mean(dd$d_tmin, na.rm = TRUE)))
    if (precip.scale == "multiplicative")
      message(sprintf("  Prec  x%4.2f", mean(dd$d_prec, na.rm = TRUE)))
    message(strrep("-", 63))
  }
  out
}


#==================================================================================================
# GAEZ AGRO-ECOLOGICAL ZONES  - retrieve FAO/IIASA GAEZ v4 data
#
# Companion to get_weather() (NASA POWER) and get_soil() (SoilGrids).
#
# GAEZ v4 has NO documented public REST API for point queries: it distributes
# whole-globe GeoTIFF rasters through a JavaScript data portal at
# https://gaez.fao.org/pages/data-access-download. So get_AEZ() is a
# raster-extract tool. It resolves a raster in this order:
#   1. an explicit `file =` path supplied by the caller (recommended);
#   2. otherwise the GAEZ layer BUNDLED with the package under inst/extdata,
#      so the function works offline out of the box;
#   3. only if that is missing, a download from `url` (an UNVERIFIED template
#      -- see aez_default_url()).
#
# The bundled layer is the GAEZ v4 57-class AEZ historical baseline
# (aez_v9v2red_5m_CRUTS32_Hist_8110_100_avg.tif, ~5 arc-min). Supply your own
# via `file =` for a different release or scenario.
#
# Dependencies: terra (raster IO / point extraction), utils (download.file).
#==================================================================================================

#' Abort when an optional package is absent
#' @keywords internal
#' @noRd
.aez_need <- function(pkg, fun) {
  miss <- pkg[!vapply(pkg, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss))
    stop(fun, "() requires the package(s): ", paste(miss, collapse = ", "),
         ".\n  Install with: install.packages(c(",
         paste0('"', miss, '"', collapse = ", "), "))", call. = FALSE)
  invisible(TRUE)
}

# Filename of the GAEZ layer shipped under inst/extdata.
.AEZ_BUNDLED_FILE <- "aez_v9v2red_5m_CRUTS32_Hist_8110_100_avg.tif"

#' Path to the GAEZ raster bundled with the package, or "" if not installed
#' @keywords internal
#' @noRd
.aez_bundled <- function() {
  system.file("extdata", .AEZ_BUNDLED_FILE, package = "etTest")
}

#' The 57-class GAEZ agro-ecological zone legend
#'
#' GAEZ v4's "AEZ classification (57 classes)" layer. Codes are the raster
#' cell values; labels follow the GAEZ v4 Model Documentation. The grouping
#' column collapses the 57 classes into the broad thermal/moisture regimes
#' that are usually what a breeding programme actually wants.
#'
#' Verify against the GAEZ v4 documentation for your specific layer before
#' relying on the labels: FAO has revised class definitions between releases,
#' and a mislabelled zone is worse than an unlabelled one.
#'
#' @keywords internal
#' @noRd
.aez_legend_57 <- function() {
  d <- data.frame(
    code = 1:57,
    label = c(
      "Tropics, lowland; semi-arid",           "Tropics, lowland; sub-humid",
      "Tropics, lowland; humid",               "Tropics, lowland; per-humid",
      "Tropics, highland; semi-arid",          "Tropics, highland; sub-humid",
      "Tropics, highland; humid",              "Tropics, highland; per-humid",
      "Subtropics, warm; semi-arid",           "Subtropics, warm; sub-humid",
      "Subtropics, warm; humid",               "Subtropics, warm; per-humid",
      "Subtropics, moderately cool; semi-arid","Subtropics, moderately cool; sub-humid",
      "Subtropics, moderately cool; humid",    "Subtropics, moderately cool; per-humid",
      "Subtropics, cool; semi-arid",           "Subtropics, cool; sub-humid",
      "Subtropics, cool; humid",               "Subtropics, cool; per-humid",
      "Subtropics, cold; semi-arid",           "Subtropics, cold; sub-humid",
      "Subtropics, cold; humid",               "Subtropics, cold; per-humid",
      "Temperate, moderate; semi-arid",        "Temperate, moderate; sub-humid",
      "Temperate, moderate; humid",            "Temperate, moderate; per-humid",
      "Temperate, cool; semi-arid",            "Temperate, cool; sub-humid",
      "Temperate, cool; humid",                "Temperate, cool; per-humid",
      "Temperate, cold; semi-arid",            "Temperate, cold; sub-humid",
      "Temperate, cold; humid",                "Temperate, cold; per-humid",
      "Boreal, cold; semi-arid",               "Boreal, cold; sub-humid",
      "Boreal, cold; humid",                   "Boreal, cold; per-humid",
      "Boreal, very cold; semi-arid",          "Boreal, very cold; sub-humid",
      "Boreal, very cold; humid",              "Boreal, very cold; per-humid",
      "Arctic / very cold",                    "Desert, tropical",
      "Desert, subtropical",                   "Desert, temperate",
      "Desert, boreal / arctic",               "Land with severe soil constraints",
      "Land with steep terrain",               "Land with poor drainage",
      "Water bodies",                          "Urban / built-up",
      "Protected area",                        "Barren / very sparse vegetation",
      "No data / unclassified"),
    stringsAsFactors = FALSE)

  # Broad regime, parsed from the label rather than re-typed, so the two
  # columns cannot drift out of sync if a label is corrected.
  d$thermal <- sub("[,;].*$", "", d$label)
  d$moisture <- ifelse(grepl("per-humid", d$label), "per-humid",
                       ifelse(grepl("sub-humid", d$label), "sub-humid",
                              ifelse(grepl("humid",     d$label), "humid",
                                     ifelse(grepl("semi-arid", d$label), "semi-arid", NA_character_))))
  d$arable <- !d$code %in% c(50:57)
  d
}

#' @title FAO GAEZ Agro-Ecological Zone Legend
#'
#' @description
#' Returns the class table used to decode GAEZ agro-ecological zone raster values
#' into human-readable zone names.
#'
#' @param classes integer. Which legend to return. Currently only \code{57}
#'   (the GAEZ v4 57-class AEZ layer) is built in.
#'
#' @return
#' A data.frame with \code{code}, \code{label}, \code{thermal}, \code{moisture}
#' and \code{arable} (logical; \code{FALSE} for water, urban, barren and the
#' constraint classes).
#'
#' @details
#' \strong{Verify before relying on these labels.} FAO has revised AEZ class
#' definitions between GAEZ releases. The table here follows the GAEZ v4
#' documentation, but if your raster came from a different release the codes may
#' not line up. \code{\link{get_AEZ}} will warn when a raster contains codes that
#' are absent from the legend -- that warning is the signal to check.
#'
#' @examples
#' leg <- aez_legend()
#' head(leg)
#' table(leg$thermal)
#'
#' @seealso \code{\link{get_AEZ}}
#' @export
aez_legend <- function(classes = 57) {
  if (!identical(as.integer(classes), 57L))
    stop("Only the 57-class GAEZ v4 legend is built in (classes = 57). ",
         "Supply your own lookup via the 'legend' argument of get_AEZ().",
         call. = FALSE)
  .aez_legend_57()
}

# Internal. Downloads the GAEZ GeoTIFF and returns a local path.
#
# THIS IS THE UNVERIFIED PART OF THE FILE. It has never been run against a
# live FAO server: the sandbox used to author this cannot reach the internet,
# and the default URL is a reconstruction (see aez_default_url()). It is only
# reached when no local file is supplied AND the bundled raster is missing.
#
# Swap it out for a proxy-aware or offline version with:
#   assignInNamespace(".aez_fetch", my_fetch, ns = "etTest")
.aez_fetch <- function(url, dir.path = NULL, overwrite = FALSE,
                       timeout = 1800L, verbose = TRUE) {

  if (!is.character(url) || length(url) != 1L || !nzchar(url))
    stop("'url' must be a single non-empty string.", call. = FALSE)

  dir.path <- dir.path %||% file.path(tempdir(), "gaez")
  dir.create(dir.path, recursive = TRUE, showWarnings = FALSE)

  dest <- file.path(dir.path, basename(sub("[?#].*$", "", url)))
  if (!grepl("\\.tif{1,2}$", dest, ignore.case = TRUE))
    dest <- paste0(dest, ".tif")

  # Cache: these layers are static, so never re-download by default.
  if (file.exists(dest) && !overwrite) {
    if (verbose) message("Using cached raster: ", dest)
    return(dest)
  }

  if (verbose)
    message("Downloading GAEZ raster (this is large and may take a while) ...\n  ",
            url)

  old <- options(timeout = max(timeout, getOption("timeout", 60L)))
  on.exit(options(old), add = TRUE)

  tmp <- paste0(dest, ".part")
  res <- try(utils::download.file(url, destfile = tmp, mode = "wb",
                                  quiet = !verbose), silent = TRUE)

  if (inherits(res, "try-error") || !file.exists(tmp)) {
    unlink(tmp)
    stop("Download failed: ", url, "\n  ",
         if (inherits(res, "try-error"))
           sub("\n.*", "", conditionMessage(attr(res, "condition"))) else "",
         "\n  The default URL is a RECONSTRUCTION and may be wrong -- see ",
         "?aez_default_url.\n  Download the layer manually from ",
         "https://gaez.fao.org and pass it with file=.", call. = FALSE)
  }

  # A failed fetch often yields an HTML error page with a .tif name. Catch it
  # here rather than letting terra::rast() fail confusingly later.
  sz <- file.info(tmp)$size
  if (is.na(sz) || sz < 1024) {
    unlink(tmp)
    stop("Downloaded file is only ", sz, " bytes -- almost certainly an error ",
         "page, not a GeoTIFF.\n  Check the URL, or download manually from ",
         "https://gaez.fao.org and pass it with file=.", call. = FALSE)
  }
  hdr <- readBin(tmp, "raw", n = 4L)
  is.tif <- (hdr[1] == as.raw(0x49) && hdr[2] == as.raw(0x49)) ||
    (hdr[1] == as.raw(0x4D) && hdr[2] == as.raw(0x4D))
  if (!is.tif) {
    unlink(tmp)
    stop("Downloaded file is not a TIFF (bad magic bytes). The URL probably ",
         "returned an HTML error page.\n  See ?aez_default_url -- the default ",
         "path is unverified.", call. = FALSE)
  }

  if (!file.rename(tmp, dest)) {
    file.copy(tmp, dest, overwrite = TRUE); unlink(tmp)
  }
  if (verbose) message("Saved: ", dest, " (",
                       round(sz / 1024^2, 1), " MB)")
  dest
}

#' @title Retrieve FAO GAEZ Agro-Ecological Zone Classification for Sites
#'
#' @description
#' Extracts the FAO/IIASA Global Agro-Ecological Zone (GAEZ v4) class at a set of
#' coordinates. By default it reads the GAEZ 57-class raster \strong{bundled with
#' the package}, so it works offline with no arguments beyond the coordinates.
#' Companion to \code{\link{get_weather}} and \code{\link{get_soil}}.
#'
#' @param env.id character vector. Environment/site identifiers.
#' @param lat numeric vector. Latitude in decimal degrees (WGS84), same length as
#'   \code{env.id}.
#' @param lon numeric vector. Longitude in decimal degrees (WGS84), same length as
#'   \code{env.id}.
#' @param file character. Path to a GAEZ GeoTIFF already on disk. If \code{NULL}
#'   (default), \code{get_AEZ} uses the raster \strong{shipped with the package}
#'   (\code{system.file("extdata", "aez_v9v2red_5m_CRUTS32_Hist_8110_100_avg.tif",
#'   package = "etTest")}). Only if that bundled file is unavailable does it fall
#'   back to downloading from \code{url}. Supply \code{file} to use a different
#'   GAEZ release or scenario.
#' @param url character. Download URL for a GAEZ raster, used only when \code{file}
#'   is \code{NULL} \emph{and} the bundled raster is missing. Defaults to
#'   \code{aez_default_url()}, a best-effort template that has \emph{not} been
#'   verified against a live FAO server -- see Details.
#' @param dir.path character. Where to cache a downloaded raster. Default
#'   \code{tempdir()}. Use a persistent directory to avoid re-downloading.
#' @param legend data.frame or NULL. Lookup table with at least \code{code} and
#'   \code{label} columns. \code{NULL} (default) uses \code{\link{aez_legend}()}.
#'   Pass \code{NA} to skip labelling and return raw codes only.
#' @param buffer numeric. If greater than 0, the radius in metres of a circular
#'   neighbourhood summarised around each point (majority class and its share),
#'   rather than a single-cell lookup. Default 0.
#' @param overwrite boolean. Re-download even if the cached file exists. Default
#'   \code{FALSE}.
#' @param timeout integer. Seconds allowed for the download. Default 1800; GAEZ
#'   rasters are large and R's 60-second default will abort them.
#' @param verbose boolean. Progress messages. Default \code{TRUE}.
#'
#' @return
#' A data.frame with one row per environment:
#' \describe{
#'   \item{\code{env}}{environment id}
#'   \item{\code{lat}, \code{lon}}{coordinates as supplied}
#'   \item{\code{aez_code}}{raster cell value}
#'   \item{\code{aez_label}, \code{thermal}, \code{moisture}, \code{arable}}{
#'     from the legend, when one is used}
#'   \item{\code{majority_share}}{when \code{buffer > 0}, the proportion of cells
#'     in the neighbourhood holding the majority class -- a purity measure}
#' }
#' The raster path used is attached as the \code{"source"} attribute.
#'
#' @details
#' \strong{Where the raster comes from.} GAEZ v4 has no documented per-point API:
#' it publishes whole-globe GeoTIFFs through a JavaScript portal. \code{get_AEZ}
#' is therefore a raster-extract tool and resolves its raster in three steps --
#' an explicit \code{file}, else the layer bundled under \code{inst/extdata}, else
#' a download from \code{url}. The bundled layer is the GAEZ v4 57-class AEZ
#' historical baseline (\code{CRUTS32_Hist_8110}, ~5 arc-min), which makes the
#' default call fully offline and reproducible.
#'
#' \strong{The default URL is unverified.} It is only used when both \code{file}
#' and the bundled raster are unavailable. FAO has reorganised GAEZ download paths
#' between releases and the portal is JavaScript-rendered, so no stable path could
#' be confirmed. If a download is triggered and fails:
#' \enumerate{
#'   \item open \url{https://gaez.fao.org/pages/data-access-download},
#'   \item download the AEZ classification layer manually,
#'   \item pass the file with \code{file = "path/to/layer.tif"}.
#' }
#'
#' \strong{Coordinates must be WGS84 lon/lat.} They are reprojected to the raster's
#' CRS automatically, but only if the raster declares one. A raster with an
#' undefined CRS triggers an error rather than a silent mis-extraction.
#'
#' \strong{Interpreting \code{buffer}.} A single-cell lookup at ~9 km resolution can
#' be misleading near a zone boundary. With \code{buffer > 0} the returned
#' \code{majority_share} shows how homogeneous the neighbourhood is; a value near
#' 1 means the site sits well inside one zone, while 0.4 means the classification
#' is effectively arbitrary at that location.
#'
#' @examples
#' \dontrun{
#' sites <- data.frame(
#'   env = c("Ames_IA", "Lincoln_NE", "Piracicaba_BR"),
#'   lat = c(42.03, 40.81, -22.71),
#'   lon = c(-93.62, -96.68, -47.63))
#'
#' ## Default: uses the raster bundled with the package (offline)
#' aez <- get_AEZ(sites$env, sites$lat, sites$lon)
#' aez
#'
#' ## Check zone purity in a 25 km neighbourhood
#' aez25 <- get_AEZ(sites$env, sites$lat, sites$lon, buffer = 25000)
#' aez25[, c("env", "aez_label", "majority_share")]
#'
#' ## Raw codes, no labelling
#' get_AEZ(sites$env, sites$lat, sites$lon, legend = NA)
#'
#' ## A different GAEZ layer you downloaded yourself
#' get_AEZ(sites$env, sites$lat, sites$lon, file = "GAEZ_AEZ_scenario.tif")
#' }
#'
#' @seealso \code{\link{aez_legend}}, \code{\link{get_weather}}, \code{\link{get_soil}}
#'
#' @references
#' FAO & IIASA (2021). \emph{Global Agro-Ecological Zones (GAEZ v4)}.
#' Rome, FAO. \url{https://gaez.fao.org}
#'
#' @importFrom utils download.file
#' @export
get_AEZ <- function(env.id, lat, lon, file = NULL, url = aez_default_url(),
                    dir.path = NULL, legend = NULL, buffer = 0,
                    overwrite = FALSE, timeout = 1800L, verbose = TRUE) {

  # ---- validate ----------------------------------------------------------
  n <- length(env.id)
  if (!n) stop("'env.id' is empty.", call. = FALSE)
  if (length(lat) != n || length(lon) != n)
    stop("'lat' and 'lon' must be the same length as 'env.id' (", n, ").",
         call. = FALSE)
  lat <- as.numeric(lat); lon <- as.numeric(lon)
  if (anyNA(lat) || anyNA(lon))
    stop("'lat'/'lon' contain missing values.", call. = FALSE)
  if (any(abs(lat) > 90))
    stop("'lat' must lie in [-90, 90]. Are lat and lon swapped?", call. = FALSE)
  if (any(abs(lon) > 180))
    stop("'lon' must lie in [-180, 180]. Are lat and lon swapped?", call. = FALSE)
  if (anyDuplicated(env.id))
    warning("Duplicated env.id values; rows will be returned as supplied.",
            call. = FALSE)
  if (!is.numeric(buffer) || length(buffer) != 1 || buffer < 0)
    stop("'buffer' must be a single non-negative number of metres.", call. = FALSE)

  .aez_need("terra", "get_AEZ")

  # ---- obtain the raster -------------------------------------------------
  # Resolution order: caller's file -> raster bundled in inst/extdata ->
  # download. The bundled layer keeps the default call offline; the download
  # (NETWORK-DEPENDENT, UNTESTED) is isolated in .aez_fetch() and reached only
  # as a last resort.
  if (is.null(file)) {
    bundled <- .aez_bundled()
    if (nzchar(bundled)) {
      file <- bundled
      if (verbose) message("Using bundled GAEZ raster: ", basename(file))
    } else {
      file <- .aez_fetch(url = url, dir.path = dir.path, overwrite = overwrite,
                         timeout = timeout, verbose = verbose)
    }
  }

  if (!file.exists(file))
    stop("Raster not found: ", file, call. = FALSE)

  r <- try(terra::rast(file), silent = TRUE)
  if (inherits(r, "try-error"))
    stop("Could not open '", file, "' as a raster.\n  ",
         sub("\n.*", "", conditionMessage(attr(r, "condition"))),
         "\n  Is this really a GeoTIFF? A failed download often leaves an ",
         "HTML error page with a .tif name.", call. = FALSE)

  if (terra::nlyr(r) > 1L) {
    warning("Raster has ", terra::nlyr(r), " layers; using the first.",
            call. = FALSE)
    r <- r[[1]]
  }

  crs.txt <- terra::crs(r)
  if (is.na(crs.txt) || !nzchar(crs.txt))
    stop("The raster has no coordinate reference system defined, so lon/lat ",
         "points cannot be placed on it reliably.\n  Set one explicitly, e.g. ",
         "terra::crs(r) <- \"EPSG:4326\", and pass the corrected raster.",
         call. = FALSE)

  # ---- extract -----------------------------------------------------------
  pts <- terra::vect(data.frame(lon = lon, lat = lat),
                     geom = c("lon", "lat"), crs = "EPSG:4326")
  pts <- terra::project(pts, crs.txt)

  if (buffer > 0) {
    if (verbose) message("Extracting majority class within ", buffer, " m ...")
    # terra::buffer works in metres for projected CRS and handles lon/lat
    # internally via geodesic buffering.
    bf  <- terra::buffer(pts, width = buffer)
    ex  <- terra::extract(r, bf)
    names(ex)[2] <- "value"
    agg <- lapply(split(ex$value, ex$ID), function(v) {
      v <- v[!is.na(v)]
      if (!length(v)) return(c(code = NA_real_, share = NA_real_))
      tb <- sort(table(v), decreasing = TRUE)
      c(code = as.numeric(names(tb)[1]), share = as.numeric(tb[1]) / length(v))
    })
    idx  <- as.integer(names(agg))
    code <- rep(NA_real_, n); share <- rep(NA_real_, n)
    code[idx]  <- vapply(agg, function(z) z[["code"]],  numeric(1))
    share[idx] <- vapply(agg, function(z) z[["share"]], numeric(1))
  } else {
    if (verbose) message("Extracting cell values at ", n, " site(s) ...")
    ex   <- terra::extract(r, pts)
    code <- as.numeric(ex[[2]])
    share <- rep(NA_real_, n)
  }

  out <- data.frame(env = as.character(env.id), lat = lat, lon = lon,
                    aez_code = code, stringsAsFactors = FALSE)
  if (buffer > 0) out$majority_share <- round(share, 3)

  # ---- label -------------------------------------------------------------
  use.legend <- !(length(legend) == 1L && is.na(legend))
  if (use.legend) {
    leg <- if (is.null(legend)) aez_legend() else as.data.frame(legend)
    if (!all(c("code", "label") %in% names(leg)))
      stop("'legend' must have at least 'code' and 'label' columns.", call. = FALSE)
    m <- match(out$aez_code, leg$code)
    out$aez_label <- leg$label[m]
    for (extra in intersect(c("thermal", "moisture", "arable"), names(leg)))
      out[[extra]] <- leg[[extra]][m]

    unknown <- unique(out$aez_code[!is.na(out$aez_code) & is.na(m)])
    if (length(unknown))
      warning("Raster contains code(s) absent from the legend: ",
              paste(sort(unknown), collapse = ", "),
              ".\n  The legend may not match this GAEZ release -- check the ",
              "layer documentation before using the labels.", call. = FALSE)
  }

  n.miss <- sum(is.na(out$aez_code))
  if (n.miss && verbose)
    message("  ", n.miss, " site(s) returned no value (outside the raster, ",
            "or over nodata such as ocean).")

  attr(out, "source") <- file
  if (verbose) {
    message("Done. ", n - n.miss, "/", n, " site(s) classified.")
  }
  out
}

#' @title Default GAEZ Download URL
#'
#' @description
#' Returns the URL template \code{\link{get_AEZ}} uses when no local file is given
#' and the bundled raster is unavailable.
#'
#' @param layer character. Which GAEZ layer. Currently \code{"aez57"}.
#'
#' @return A single character string.
#'
#' @details
#' \strong{This URL is unverified.} It was constructed from the documented GAEZ v4
#' portal structure but could not be tested against a live FAO server, and FAO has
#' reorganised these paths between releases. It is exposed as a function, rather
#' than hidden inside \code{get_AEZ()}, so it can be inspected and overridden
#' without editing the package:
#'
#' \preformatted{
#'   get_AEZ(env, lat, lon, url = "https://.../correct/path.tif")
#' }
#'
#' For reproducible analyses, rely on the bundled raster or download the layer
#' once and pass \code{file=}.
#'
#' @examples
#' aez_default_url()
#'
#' @seealso \code{\link{get_AEZ}}
#' @export
aez_default_url <- function(layer = "aez57") {
  layer <- match.arg(layer, c("aez57"))
  switch(layer,
         aez57 = paste0("https://s3.eu-west-1.amazonaws.com/data.gaezdev.aws.fao.org/",
                        "LR/aez/57_class/aez_v9v2red_5m_ENSEMBLE_rcp2p6_2020s.tif"))
}
