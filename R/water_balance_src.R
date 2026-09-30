#==================================================================================================
# water_balance.R
#
# Daily soil water balance and crop water-stress indices, combining daily
# weather from get_weather() with soil hydraulic properties from get_soil().
#
# Implements the FAO-56 single crop coefficient method:
#
#   Allen, R.G., Pereira, L.S., Raes, D. & Smith, M. (1998)
#   "Crop evapotranspiration: guidelines for computing crop water
#   requirements." FAO Irrigation and Drainage Paper 56. FAO, Rome.
#   ISBN 92-5-104219-5
#
# Companion to env_data_collection.R and soil_classification.R.
#
# --------------------------------------------------------------------------
# SCOPE AND LIMITATIONS -- READ BEFORE USE
# --------------------------------------------------------------------------
# This is the SINGLE coefficient method (FAO-56 Ch. 6-8). It does NOT:
#   * separate soil evaporation from transpiration (no dual Kc, Ch. 7)
#   * model surface runoff -- ALL rainfall is assumed to infiltrate
#   * model capillary rise from a shallow water table
#   * model lateral flow, or any horizontal redistribution
#   * schedule irrigation for you
# On steep, crusted or low-infiltration soils the no-runoff assumption
# overestimates stored water and therefore UNDERESTIMATES stress. Subtract
# runoff from the precipitation column beforehand if that matters.
#
# The Kc and p defaults describe a generic 120-day maize crop. They are
# FAO-56 tabulated starting points, not calibrated values for your cultivar,
# planting density or region. Calibrate before drawing agronomic conclusions.
#
# --------------------------------------------------------------------------
# METHOD REFERENCES
# --------------------------------------------------------------------------
# Water balance, crop coefficients, stress coefficient
#   Allen, R.G., Pereira, L.S., Raes, D. & Smith, M. (1998) Crop
#     evapotranspiration. FAO Irrigation and Drainage Paper 56. FAO, Rome.
#     -- Ch. 6 (Kc), Table 12 (Kc values), Table 22 (p values),
#        Ch. 8 and Eq. 84-85 (root-zone depletion, Ks)
#   Allen, R.G. et al. (2005) FAO-56 dual crop coefficient method for
#     estimating evaporation from soil and application extensions.
#     J. Irrig. Drain. Eng. 131(1), 2-13.
#     doi:10.1061/(ASCE)0733-9437(2005)131:1(2)
#   Pereira, L.S. et al. (2021) Standard single and basal crop coefficients
#     for field crops. Agricultural Water Management 243, 106466.
#     doi:10.1016/j.agwat.2020.106466
#
# Yield response to water deficit (why ETa/ETc is reported)
#   Doorenbos, J. & Kassam, A.H. (1979) Yield response to water. FAO
#     Irrigation and Drainage Paper 33. FAO, Rome.
#   Steduto, P., Hsiao, T.C., Fereres, E. & Raes, D. (2012) Crop yield
#     response to water. FAO Irrigation and Drainage Paper 66. FAO, Rome.
#
# Soil water retention and pedotransfer
#   Poggio, L. et al. (2021) SoilGrids 2.0: producing soil information for the
#     globe with quantified spatial uncertainty. SOIL 7, 217-240.
#     doi:10.5194/soil-7-217-2021
#   Turek, M.E. et al. (2023) Global mapping of volumetric water retention at
#     100, 330 and 15000 cm suction using the WoSIS database. International
#     Soil and Water Conservation Research 11(2), 225-239.
#     doi:10.1016/j.iswcr.2022.08.001
#   Saxton, K.E. & Rawls, W.J. (2006) Soil water characteristic estimates by
#     texture and organic matter for hydrologic solutions. Soil Sci. Soc. Am.
#     J. 70(5), 1569-1578. doi:10.2136/sssaj2005.0117
#
# Root growth
#   Borg, H. & Grimes, D.W. (1986) Depth development of roots with time: an
#     empirical description. Transactions of the ASAE 29(1), 194-197.
#     doi:10.13031/2013.30125
#
# Environmental covariates in GxE analysis (downstream use of these indices)
#   Costa-Neto, G., Fritsche-Neto, R. & Crossa, J. (2021) Nonlinear kernels,
#     dominance, and envirotyping data increase the accuracy of genome-based
#     prediction in multi-environment trials. Heredity 126, 92-106.
#     doi:10.1038/s41437-020-00353-1
#   Resende, R.T. et al. (2021) Enviromics in breeding: applications and
#     perspectives on envirotypic-assisted selection. Theor. Appl. Genet. 134,
#     95-112. doi:10.1007/s00122-020-03684-z
#==================================================================================================


# =========================================================================
# SECTION 0 - internal helpers (not exported)
# =========================================================================

#' SoilGrids mapped-unit conversion factors
#'
#' SoilGrids stores every property as an integer to save space, so the served
#' values are not in conventional units. Dividing by these factors recovers the
#' conventional units documented by ISRIC.
#'
#' Getting this wrong is silent and severe: \code{wv0033 = 314} is a volumetric
#' water content of 0.314 cm3/cm3, not 314. Treating it literally inflates total
#' available water by three orders of magnitude, and every stress index then
#' collapses to "no stress" -- a result that looks clean and is entirely wrong.
#'
#' @format Named numeric vector; names are SoilGrids property codes.
#' @references Poggio, L. et al. (2021) SOIL 7, 217-240.
#'   \doi{10.5194/soil-7-217-2021}
#' @keywords internal
#' @noRd
.SG_FACTORS <- c(
  bdod = 100,    # cg/cm3            -> kg/dm3
  cec  = 10,     # mmol(c)/kg        -> cmol(c)/kg
  cfvo = 10,     # cm3/dm3           -> vol%
  clay = 10,     # g/kg              -> %
  sand = 10,     # g/kg              -> %
  silt = 10,     # g/kg              -> %
  nitrogen = 100,# cg/kg             -> g/kg
  soc  = 10,     # dg/kg             -> g/kg
  ocd  = 10,     # hg/m3             -> kg/m3
  ocs  = 10,     # t/ha              -> kg/m2
  phh2o = 10,    # pH*10             -> pH
  wv0010 = 1000, # 10^-3 cm3/cm3     -> cm3/cm3
  wv0033 = 1000, # 10^-3 cm3/cm3     -> cm3/cm3
  wv1500 = 1000  # 10^-3 cm3/cm3     -> cm3/cm3
)

#' Standard SoilGrids depth intervals, in mm
#'
#' SoilGrids reports six fixed depth intervals. Depths are stored in mm so they
#' combine directly with rooting depth, which is also in mm.
#'
#' @format data.frame with columns \code{layer}, \code{top}, \code{bot}.
#' @keywords internal
#' @noRd
.SG_DEPTHS <- data.frame(
  layer = c("0_5cm", "5_15cm", "15_30cm", "30_60cm", "60_100cm", "100_200cm"),
  top   = c(0, 50, 150, 300, 600, 1000),
  bot   = c(50, 150, 300, 600, 1000, 2000),
  stringsAsFactors = FALSE
)

#' Parse a wide get_soil() row into long property/depth records
#'
#' Column names produced by \code{get_soil(wide = TRUE)} look like
#' \code{"wv0033|30_60cm"}. This splits them on the pipe, applies the ISRIC
#' conversion factor, and attaches the depth interval in mm.
#'
#' @param soil.row one-row data.frame (or list) of wide-format soil values.
#' @param convert logical. Apply \code{.SG_FACTORS}. Default \code{TRUE}.
#' @return data.frame with \code{layer}, \code{property}, \code{value},
#'   \code{top}, \code{bot}.
#' @keywords internal
#' @noRd
.parse_soil_wide <- function(soil.row, convert = TRUE) {
  nm  <- names(soil.row)
  key <- grep("\\|", nm, value = TRUE)
  if (!length(key))
    stop("No 'property|depth' columns found. Was this produced by ",
         "get_soil(wide = TRUE)?", call. = FALSE)

  parts <- do.call(rbind, strsplit(key, "|", fixed = TRUE))
  out <- data.frame(
    property = parts[, 1],
    layer    = parts[, 2],
    value    = as.numeric(unlist(soil.row[key], use.names = FALSE)),
    stringsAsFactors = FALSE
  )

  if (isTRUE(convert)) {
    f <- .SG_FACTORS[out$property]
    unknown <- unique(out$property[is.na(f)])
    if (length(unknown))
      warning("No SoilGrids conversion factor for: ",
              paste(unknown, collapse = ", "), ". Values left unconverted.",
              call. = FALSE)
    f[is.na(f)] <- 1
    out$value <- out$value / as.numeric(f)
  }

  merge(out, .SG_DEPTHS, by = "layer", all.x = TRUE)
}

#' Total available water over a root zone
#'
#' Integrates (field capacity - permanent wilting point) across the SoilGrids
#' layers intersecting the root zone, weighting each layer by the depth it
#' actually contributes:
#' \deqn{TAW = \sum_j (\theta_{FC,j} - \theta_{WP,j}) \, \Delta z_j}
#' where \eqn{\Delta z_j} is the overlap between layer \eqn{j} and the root
#' zone.
#'
#' Taking a single layer instead would misstate TAW badly: SoilGrids layers
#' thicken with depth (50 mm at the surface, 400 mm at 60-100 cm), so an
#' unweighted mean over-weights the thin topsoil layers.
#'
#' @param soil.long output of \code{.parse_soil_wide()}.
#' @param root.depth numeric. Rooting depth in mm.
#' @param fc.var,pwp.var character. Property names for field capacity and
#'   wilting point.
#' @return list with \code{taw} (mm) and \code{layers}, the per-layer
#'   contributions.
#' @keywords internal
#' @noRd
.compute_taw <- function(soil.long, root.depth, fc.var = "wv0033",
                         pwp.var = "wv1500") {
  fc  <- soil.long[soil.long$property == fc.var, ]
  pwp <- soil.long[soil.long$property == pwp.var, ]
  if (!nrow(fc) || !nrow(pwp))
    stop("Soil data must contain '", fc.var, "' and '", pwp.var,
         "'. Request them in get_soil(variables.names = ...).", call. = FALSE)

  lay <- merge(fc[, c("layer", "top", "bot", "value")],
               pwp[, c("layer", "value")],
               by = "layer", suffixes = c(".fc", ".pwp"))
  lay <- lay[order(lay$top), ]

  lay$awc <- lay$value.fc - lay$value.pwp
  bad <- !is.na(lay$awc) & lay$awc <= 0
  if (any(bad))
    warning(sum(bad), " layer(s) with field capacity <= wilting point; ",
            "treated as zero available water.", call. = FALSE)
  lay$awc <- pmax(lay$awc, 0)

  # Depth of each layer that falls inside the root zone.
  lay$eff_mm <- pmax(0, pmin(lay$bot, root.depth) - pmin(lay$top, root.depth))
  lay$taw_mm <- lay$awc * lay$eff_mm

  if (root.depth > max(lay$bot, na.rm = TRUE))
    warning("Root depth (", root.depth, " mm) exceeds the deepest soil layer (",
            max(lay$bot, na.rm = TRUE), " mm). TAW covers the mapped profile only.",
            call. = FALSE)

  list(taw = sum(lay$taw_mm, na.rm = TRUE), layers = lay)
}

#' Linear interpolation of a crop coefficient curve over days after planting
#'
#' FAO-56 defines Kc as a piecewise-linear curve through the initial, crop
#' development, mid-season and late-season stages (Ch. 6, Fig. 25). Values
#' outside the supplied range are held flat (\code{rule = 2}) rather than
#' extrapolated, since extrapolating a Kc curve produces nonsense.
#'
#' @param dap numeric. Days after planting.
#' @param stages numeric. Increasing breakpoints in days after planting.
#' @param kc numeric. Crop coefficients at those breakpoints.
#' @return numeric vector of Kc, same length as \code{dap}.
#' @references Allen et al. (1998), FAO-56 Ch. 6 and Table 12.
#' @keywords internal
#' @noRd
.kc_curve <- function(dap, stages, kc) {
  if (length(stages) != length(kc))
    stop("'kc.stages' and 'kc.values' must have the same length.", call. = FALSE)
  if (is.unsorted(stages))
    stop("'kc.stages' must be increasing.", call. = FALSE)
  stats::approx(x = stages, y = kc, xout = dap, method = "linear",
                rule = 2)$y
}

#' Linear root-growth curve, capped at maximum rooting depth
#'
#' Rooting depth increases linearly from \code{root.init} at planting to
#' \code{root.max} at \code{dap.max}, then stays constant.
#'
#' A linear ramp is a simplification. Observed root descent is closer to
#' sigmoidal (Borg & Grimes 1986), and real depth is limited by compaction,
#' acidity and water-table depth that this function knows nothing about. The
#' consequence is that early-season TAW is somewhat approximate. Since stress
#' during establishment is rarely the binding constraint, the simplification is
#' usually acceptable -- but supply a shallower \code{root.depth} if you know
#' the profile is restricted.
#'
#' @param dap numeric. Days after planting.
#' @param root.init,root.max numeric. Initial and maximum rooting depth (mm).
#' @param dap.max numeric. Days after planting at which \code{root.max} is met.
#' @return numeric vector of rooting depth in mm.
#' @references Borg, H. & Grimes, D.W. (1986) Trans. ASAE 29(1), 194-197.
#' @keywords internal
#' @noRd
.root_curve <- function(dap, root.init, root.max, dap.max) {
  frac <- pmin(1, pmax(0, dap / dap.max))
  root.init + (root.max - root.init) * frac
}


# =========================================================================
# SECTION 1 - water_balance()
# =========================================================================

#' @title Daily Soil Water Balance and Crop Water-Stress Indices
#'
#' @description
#' Runs a FAO-56 single crop coefficient soil water balance by combining daily
#' weather from \code{get_weather} with soil hydraulic properties from
#' \code{get_soil}, returning daily root-zone depletion, a water-stress
#' coefficient, actual evapotranspiration, and a categorical stress level.
#'
#' The stress coefficient \code{Ks} is the main output: bounded in [0, 1],
#' dimensionless, and directly proportional to the reduction in transpiration.
#'
#' @param env.data data.frame. Daily weather, typically from
#'   \code{get_weather} or \code{processWTH}. Must contain the
#'   environment id, a day counter, precipitation, and a reference
#'   evapotranspiration column.
#' @param soil.data data.frame. Soil properties from
#'   \code{get_soil(wide = TRUE)}. Must include \code{wv0033} (field capacity)
#'   and \code{wv1500} (wilting point) for the layers spanning the root zone.
#' @param env.id character. Environment-id column present in both inputs.
#'   Default \code{"env"}.
#' @param days.id character. Column giving days from planting. Default
#'   \code{"daysFromStart"}.
#' @param PREC character. Precipitation column (mm/day). Default
#'   \code{"PRECTOT"}.
#' @param ETo character. Reference evapotranspiration column (mm/day). If
#'   \code{NULL} (default), the function looks for \code{"ETP"} (added by
#'   \code{param_atmospheric}), then \code{"EVPTRNS"}.
#' @param irrigation numeric or character. Irrigation in mm/day: a single value
#'   applied daily, a vector as long as \code{env.data}, or the name of a
#'   column. Default 0.
#' @param root.depth numeric. Maximum rooting depth in mm. Default 1000.
#' @param root.init numeric. Rooting depth at emergence in mm. Default 300.
#' @param dap.root.max numeric. Days after planting at which maximum rooting
#'   depth is reached. Default 60.
#' @param kc.stages,kc.values numeric. Days after planting and the corresponding
#'   crop coefficients, linearly interpolated between. Defaults describe a
#'   generic 120-day maize crop (FAO-56 Table 12). \strong{Calibrate these.}
#' @param p numeric. Soil-water depletion fraction for no stress (FAO-56 Table
#'   22). Default 0.55 for maize.
#' @param p.adjust logical. Adjust \code{p} daily for evaporative demand using
#'   FAO-56 Eq. 84. Default \code{TRUE}.
#' @param initial.depletion numeric. Root-zone depletion at planting in mm, or
#'   \code{NULL} (default) to start at field capacity.
#' @param fc.var,pwp.var character. SoilGrids property names for field capacity
#'   and wilting point. Defaults \code{"wv0033"} and \code{"wv1500"}.
#' @param convert.units logical. Apply ISRIC conversion factors to the soil
#'   data. Set \code{FALSE} only if already converted. Default \code{TRUE}.
#' @param verbose logical. Print progress messages. Default \code{TRUE}.
#'
#' @details
#' \strong{Method.} FAO-56 single crop coefficient approach (Allen et al. 1998,
#' Ch. 8). Root-zone depletion is stepped forward daily:
#' \deqn{D_{r,i} = D_{r,i-1} - P_i - I_i + ET_{a,i} + DP_i}
#' bounded to \eqn{[0, TAW]}. Water in excess of field capacity becomes deep
#' percolation; depletion cannot exceed total available water.
#'
#' The stress coefficient follows FAO-56 Eq. 84:
#' \deqn{K_s = 1 \quad \mathrm{if}\ D_r \le RAW}
#' \deqn{K_s = \frac{TAW - D_r}{TAW - RAW} \quad \mathrm{if}\ D_r > RAW}
#' so transpiration proceeds at the potential rate until readily available
#' water is exhausted, then declines linearly to zero at wilting point.
#'
#' \strong{Causal ordering.} \code{Ks} on day \eqn{i} is computed from
#' \emph{yesterday's} depletion against \emph{today's} capacity. Using the same
#' day's closing depletion would let today's rainfall relieve stress that
#' already reduced today's transpiration -- a subtle look-ahead that
#' systematically understates stress.
#'
#' \strong{Units.} SoilGrids serves integers in mapped units:
#' \code{wv0033 = 314} means 0.314 cm3/cm3. Conversion is automatic; see
#' \code{convert.units}.
#'
#' \strong{Total available water} is integrated across the SoilGrids layers
#' intersecting the root zone, each weighted by the depth it contributes, not
#' taken from a single layer. Because rooting depth grows over the season, TAW
#' is recomputed daily (cached by unique depth for speed).
#'
#' \strong{What this does not model.} No dual crop coefficient, so soil
#' evaporation is not separated from transpiration. No capillary rise, no
#' lateral flow, and \strong{no surface runoff} -- all rainfall is assumed to
#' infiltrate, which overestimates stored water and understates stress on
#' steep, crusted or low-infiltration soils. Subtract runoff from the
#' precipitation column beforehand if that matters. Irrigation is applied as
#' given; the balance does not schedule it.
#'
#' \strong{Interpreting the output.} Regress \code{Ks} against yield: it is
#' bounded, dimensionless, and scales transpiration directly.
#' \code{depletion_frac} is easier to explain to agronomists.
#' \code{stress_level} is a convenience for grouping and tabulation and should
#' not be used as a modelling covariate -- it discards magnitude.
#'
#' @return
#' A data.frame with one row per environment-day, containing the input columns
#' plus:
#' \describe{
#'   \item{\code{TAW}}{total available water in the root zone (mm)}
#'   \item{\code{RAW}}{readily available water, \eqn{p \times TAW} (mm)}
#'   \item{\code{Dr}}{root-zone depletion at end of day (mm)}
#'   \item{\code{ASW}}{available soil water remaining, \eqn{TAW - D_r} (mm)}
#'   \item{\code{depletion_frac}}{\eqn{D_r / TAW}; 0 at field capacity, 1 at
#'     wilting point}
#'   \item{\code{Ks}}{water-stress coefficient; 1 unstressed, 0 fully stressed}
#'   \item{\code{ETc}}{crop evapotranspiration, standard conditions (mm/day)}
#'   \item{\code{ETa}}{actual evapotranspiration, \eqn{K_s \times ET_c}}
#'   \item{\code{deficit}}{\eqn{ET_c - ET_a} (mm/day)}
#'   \item{\code{drainage}}{deep percolation below the root zone (mm/day)}
#'   \item{\code{irrigation}}{irrigation applied (mm/day)}
#'   \item{\code{Kc}}{interpolated crop coefficient}
#'   \item{\code{root_depth}}{rooting depth that day (mm)}
#'   \item{\code{stress_level}}{ordered factor: none, mild, moderate, severe}
#' }
#'
#' @references
#' Allen, R.G., Pereira, L.S., Raes, D. & Smith, M. (1998) Crop
#' evapotranspiration: guidelines for computing crop water requirements.
#' \emph{FAO Irrigation and Drainage Paper 56}. FAO, Rome.
#'
#' Pereira, L.S. et al. (2021) Standard single and basal crop coefficients for
#' field crops. \emph{Agricultural Water Management} 243, 106466.
#' \doi{10.1016/j.agwat.2020.106466}
#'
#' Doorenbos, J. & Kassam, A.H. (1979) Yield response to water. \emph{FAO
#' Irrigation and Drainage Paper 33}. FAO, Rome.
#'
#' Steduto, P., Hsiao, T.C., Fereres, E. & Raes, D. (2012) Crop yield response
#' to water. \emph{FAO Irrigation and Drainage Paper 66}. FAO, Rome.
#'
#' Poggio, L. et al. (2021) SoilGrids 2.0. \emph{SOIL} 7, 217-240.
#' \doi{10.5194/soil-7-217-2021}
#'
#' Turek, M.E. et al. (2023) Global mapping of volumetric water retention at
#' 100, 330 and 15000 cm suction using the WoSIS database. \emph{International
#' Soil and Water Conservation Research} 11(2), 225-239.
#' \doi{10.1016/j.iswcr.2022.08.001}
#'
#' Borg, H. & Grimes, D.W. (1986) Depth development of roots with time.
#' \emph{Transactions of the ASAE} 29(1), 194-197. \doi{10.13031/2013.30125}
#'
#' @seealso \code{summary_water_balance} to condense the daily output,
#'   \code{get_weather} and \code{get_soil} for the inputs,
#'   \code{processWTH} to derive ETo,
#'   \code{soil_classification} for cross-site soil zoning.
#'
#' @examples
#' \dontrun{
#' ## ---------------------------------------------------------------
#' ## 1. Minimal use
#' ## ---------------------------------------------------------------
#' wth  <- get_weather(env.id = "NM", lat = -13.05, lon = -56.05,
#'                     start.day = "2015-02-15", end.day = "2015-06-15")
#' wth  <- processWTH(wth)                 # adds ETP
#' soil <- get_soil(env.id = "NM", lat = -13.05, lon = -56.05,
#'                  variables.names = c("wv0033", "wv1500", "clay"))
#'
#' wb <- water_balance(wth, soil)
#' head(wb[, c("env", "daysFromStart", "TAW", "Dr", "Ks",
#'             "deficit", "stress_level")])
#'
#'
#' ## ---------------------------------------------------------------
#' ## 2. Several environments at once
#' ## ---------------------------------------------------------------
#' sites <- data.frame(
#'   env = c("SOR", "LON", "BAR"),
#'   lat = c(-12.5453, -23.3045, -12.1530),
#'   lon = c(-55.7113, -51.1696, -44.9900))
#'
#' wth  <- processWTH(get_weather(env.id = sites$env, lat = sites$lat,
#'                                lon = sites$lon,
#'                                start.day = "2023-10-15",
#'                                end.day   = "2024-02-12"))
#' soil <- get_soil(env.id = sites$env, lat = sites$lat, lon = sites$lon,
#'                  variables.names = c("wv0033", "wv1500"))
#'
#' wb <- water_balance(wth, soil, root.depth = 1000)
#' summary_water_balance(wb)
#'
#'
#' ## ---------------------------------------------------------------
#' ## 3. Irrigation: constant, per-day, or from a column
#' ## ---------------------------------------------------------------
#' water_balance(wth, soil, irrigation = 5)              # 5 mm every day
#' water_balance(wth, soil, irrigation = rep(0, nrow(wth)))
#'
#' wth$irr <- ifelse(wth$daysFromStart %in% 40:70, 8, 0)  # top-up at flowering
#' wb_irr  <- water_balance(wth, soil, irrigation = "irr")
#'
#' ## How much stress did irrigation remove?
#' c(rainfed   = sum(water_balance(wth, soil, verbose = FALSE)$deficit),
#'   irrigated = sum(wb_irr$deficit))
#'
#'
#' ## ---------------------------------------------------------------
#' ## 4. A different crop: soybean
#' ## ---------------------------------------------------------------
#' ## Kc and p come from FAO-56 Tables 12 and 22. Always check them
#' ## against your own crop, cycle length and region.
#' wb_soy <- water_balance(
#'   wth, soil,
#'   kc.stages = c(1, 25, 50, 95, 120),
#'   kc.values = c(0.40, 0.80, 1.15, 1.15, 0.50),
#'   p         = 0.50,          # soybean, FAO-56 Table 22
#'   root.depth = 800)
#'
#'
#' ## ---------------------------------------------------------------
#' ## 5. Restricted root zone (compaction, shallow profile)
#' ## ---------------------------------------------------------------
#' shallow <- water_balance(wth, soil, root.depth = 400, verbose = FALSE)
#' deep    <- water_balance(wth, soil, root.depth = 1200, verbose = FALSE)
#' c(shallow_TAW = max(shallow$TAW), deep_TAW = max(deep$TAW))
#' ## A shallower profile stores less water and stresses sooner.
#'
#'
#' ## ---------------------------------------------------------------
#' ## 6. Dry planting: start part-depleted rather than at field capacity
#' ## ---------------------------------------------------------------
#' wb_dry <- water_balance(wth, soil, initial.depletion = 60)
#'
#'
#' ## ---------------------------------------------------------------
#' ## 7. Season and growth-stage summaries
#' ## ---------------------------------------------------------------
#' summary_water_balance(wb)
#'
#' summary_water_balance(wb, by.interval = TRUE,
#'                       time.window  = c(0, 14, 35, 60, 90, 120),
#'                       names.window = c("P-E", "E-V1", "V1-V4",
#'                                        "V4-VT", "VT-GF", "GF-PM"))
#' ## Stress at VT-GF (flowering to grain fill) costs far more yield
#' ## than the same stress during vegetative growth.
#'
#'
#' ## ---------------------------------------------------------------
#' ## 8. Using the indices as envirotype covariates
#' ## ---------------------------------------------------------------
#' W <- summary_water_balance(wb, by.interval = TRUE,
#'                            time.window = c(0, 30, 60, 90, 120))
#' wide <- reshape(W[, c("env", "interval", "mean_Ks")],
#'                 idvar = "env", timevar = "interval", direction = "wide")
#' head(wide)   # one row per environment, ready to join to trial data
#' }
#'
#' @importFrom stats approx
#' @export
water_balance <- function(env.data, soil.data,
                          env.id = "env", days.id = "daysFromStart",
                          PREC = "PRECTOT", ETo = NULL,
                          irrigation = 0,
                          root.depth = 1000, root.init = 300,
                          dap.root.max = 60,
                          kc.stages = c(1, 20, 45, 90, 120),
                          kc.values = c(0.40, 0.75, 1.15, 1.15, 0.60),
                          p = 0.55, p.adjust = TRUE,
                          initial.depletion = NULL,
                          fc.var = "wv0033", pwp.var = "wv1500",
                          convert.units = TRUE, verbose = TRUE) {

  .et_banner("water_balance", "runs a FAO-56 daily soil water balance", verbose)
  env.data  <- as.data.frame(env.data)
  soil.data <- as.data.frame(soil.data)

  # ---- resolve columns -----------------------------------------------------
  if (!env.id %in% names(env.data))
    stop("Environment column '", env.id, "' not found in env.data.", call. = FALSE)
  if (!env.id %in% names(soil.data))
    stop("Environment column '", env.id, "' not found in soil.data.", call. = FALSE)
  if (!days.id %in% names(env.data))
    stop("Day column '", days.id, "' not found in env.data.", call. = FALSE)
  if (!PREC %in% names(env.data))
    stop("Precipitation column '", PREC, "' not found in env.data.", call. = FALSE)

  if (is.null(ETo)) {
    ETo <- intersect(c("ETP", "EVPTRNS"), names(env.data))[1]
    if (is.na(ETo))
      stop("No reference evapotranspiration column found. Run processWTH() ",
           "first, or name the column via 'ETo'.", call. = FALSE)
    if (verbose) message("Using '", ETo, "' as reference evapotranspiration.")
  }
  if (!ETo %in% names(env.data))
    stop("ETo column '", ETo, "' not found in env.data.", call. = FALSE)

  if (!(p > 0 && p < 1)) stop("'p' must lie strictly between 0 and 1.", call. = FALSE)
  if (root.init <= 0 || root.depth <= 0)
    stop("Rooting depths must be positive.", call. = FALSE)
  if (root.init > root.depth)
    stop("'root.init' cannot exceed 'root.depth'.", call. = FALSE)
  if (dap.root.max <= 0) stop("'dap.root.max' must be positive.", call. = FALSE)

  # ---- irrigation ----------------------------------------------------------
  if (is.character(irrigation)) {
    if (!irrigation %in% names(env.data))
      stop("Irrigation column '", irrigation, "' not found in env.data.", call. = FALSE)
    irr_all <- env.data[[irrigation]]
  } else if (length(irrigation) == 1L) {
    irr_all <- rep(irrigation, nrow(env.data))
  } else if (length(irrigation) == nrow(env.data)) {
    irr_all <- irrigation
  } else {
    stop("'irrigation' must be length 1, length nrow(env.data), or a column name.",
         call. = FALSE)
  }
  irr_all[is.na(irr_all)] <- 0

  envs <- unique(as.character(env.data[[env.id]]))
  missing_soil <- setdiff(envs, as.character(soil.data[[env.id]]))
  if (length(missing_soil))
    stop("No soil data for environment(s): ", paste(missing_soil, collapse = ", "),
         call. = FALSE)

  if (verbose) {
    message("---------------------------------------------------------------")
    message("water_balance() - FAO-56 soil water balance")
    message("Allen et al. (1998), FAO Irrigation and Drainage Paper 56")
    message("---------------------------------------------------------------")
    message("Environments ....................... ", length(envs))
    message("Rooting depth ...................... ", root.init, " -> ",
            root.depth, " mm by ", dap.root.max, " DAP")
    message("Depletion fraction p ............... ", p,
            if (p.adjust) " (adjusted daily for ETc)" else "")
  }

  out <- vector("list", length(envs))

  for (k in seq_along(envs)) {
    e   <- envs[k]
    idx <- which(as.character(env.data[[env.id]]) == e)
    d   <- env.data[idx, , drop = FALSE]
    ord <- order(d[[days.id]])
    d   <- d[ord, , drop = FALSE]
    irr <- irr_all[idx][ord]

    srow <- soil.data[as.character(soil.data[[env.id]]) == e, , drop = FALSE][1, ]
    slong <- .parse_soil_wide(srow, convert = convert.units)

    dap <- as.numeric(d[[days.id]])
    Kc  <- .kc_curve(dap, kc.stages, kc.values)
    rd  <- .root_curve(dap, root.init, root.depth, dap.root.max)

    # TAW depends on rooting depth, so it changes daily. Cache by unique depth
    # to avoid repeating the layer integration for every day of the season.
    urd <- unique(rd)
    taw_lookup <- vapply(urd, function(z)
      .compute_taw(slong, z, fc.var, pwp.var)$taw, numeric(1))
    TAW <- taw_lookup[match(rd, urd)]

    if (all(TAW <= 0))
      stop("Total available water is zero for '", e,
           "'. Check that ", fc.var, " and ", pwp.var, " were returned by get_soil().",
           call. = FALSE)

    P   <- d[[PREC]]; P[is.na(P)] <- 0
    ET0 <- d[[ETo]];  ET0[is.na(ET0)] <- 0
    ET0 <- pmax(ET0, 0)
    ETc <- Kc * ET0

    # FAO-56 Eq. 84: p rises in low-demand conditions, falls in high demand,
    # bounded to [0.1, 0.8] as the paper recommends.
    p_day <- if (isTRUE(p.adjust)) pmin(0.8, pmax(0.1, p + 0.04 * (5 - ETc))) else
      rep(p, length(ETc))
    RAW <- p_day * TAW

    n   <- nrow(d)
    Dr  <- numeric(n); Ks <- numeric(n); ETa <- numeric(n)
    DP  <- numeric(n)

    prev <- if (is.null(initial.depletion)) 0 else
      min(max(initial.depletion, 0), TAW[1])

    for (i in seq_len(n)) {
      # Stress is governed by YESTERDAY's depletion against TODAY's capacity.
      # Using today's closing depletion would let today's rain relieve stress
      # that already suppressed today's transpiration -- a look-ahead that
      # systematically understates stress.
      Ks[i]  <- if (TAW[i] <= 0) 0 else
        if (prev <= RAW[i]) 1 else
          max(0, (TAW[i] - prev) / (TAW[i] - RAW[i]))
      ETa[i] <- ETc[i] * Ks[i]

      bal   <- prev - P[i] - irr[i] + ETa[i]
      DP[i] <- max(0, -bal)                    # surplus drains away
      Dr[i] <- min(max(bal, 0), TAW[i])        # bounded [0, TAW]
      prev  <- Dr[i]
    }

    res <- d
    res$Kc             <- round(Kc, 4)
    res$root_depth     <- round(rd, 1)
    res$TAW            <- round(TAW, 3)
    res$RAW            <- round(RAW, 3)
    res$Dr             <- round(Dr, 3)
    res$ASW            <- round(pmax(TAW - Dr, 0), 3)
    res$depletion_frac <- round(ifelse(TAW > 0, Dr / TAW, NA_real_), 4)
    res$Ks             <- round(Ks, 4)
    res$ETc            <- round(ETc, 3)
    res$ETa            <- round(ETa, 3)
    res$deficit        <- round(ETc - ETa, 3)
    res$drainage       <- round(DP, 3)
    res$irrigation     <- irr

    out[[k]] <- res

    if (verbose)
      message("  [", k, "/", length(envs), "] ", e,
              " | TAW ", round(min(TAW)), "-", round(max(TAW)), " mm",
              " | stressed days ", sum(Ks < 1), "/", n)
  }

  final <- do.call(rbind, out)
  rownames(final) <- NULL

  final$stress_level <- cut(
    1 - final$Ks,
    breaks = c(-Inf, 1e-6, 0.2, 0.5, Inf),
    labels = c("none", "mild", "moderate", "severe"),
    ordered_result = TRUE
  )

  if (verbose) {
    tb <- table(final$stress_level)
    message("Stress days: ",
            paste(sprintf("%s=%d", names(tb), as.integer(tb)), collapse = "  "))
    message("Done.")
  }

  final
}


# =========================================================================
# SECTION 2 - summary_water_balance()
# =========================================================================

#' @title Summarise Water-Balance Output by Environment or Growth Stage
#'
#' @description
#' Condenses the daily output of \code{water_balance} into one row per
#' environment, or per environment and growth stage, producing the season-level
#' indices most often used as environmental covariates in
#' genotype-by-environment models.
#'
#' @param wb data.frame. Output of \code{water_balance}.
#' @param env.id character. Environment-id column. Default \code{"env"}.
#' @param days.id character. Day counter column. Default
#'   \code{"daysFromStart"}.
#' @param by.interval logical. Summarise within growth stages as well as by
#'   environment. Default \code{FALSE}.
#' @param time.window numeric. Stage breakpoints in days after planting. When
#'   \code{NULL}, 10-day windows are used.
#' @param names.window character. Stage names, one per interval. Ignored if the
#'   length does not match the number of intervals.
#' @param verbose boolean. If \code{TRUE} (default) prints a progress banner.
#'
#' @details
#' \code{ETa_ETc_ratio} is the relative transpiration deficit, the quantity
#' entering most water-limited yield models -- notably the FAO Ky framework of
#' Doorenbos & Kassam (1979), where relative yield loss is proportional to
#' \eqn{1 - ET_a/ET_c}.
#'
#' \code{max_consecutive_stress} distinguishes one damaging spell from the same
#' number of scattered days, which a mean cannot. Ten consecutive stressed days
#' at flowering is a different event from ten isolated days across the season,
#' and yield responds very differently to the two.
#'
#' \strong{Stage timing matters more than season totals.} Water deficit during
#' flowering and grain fill costs disproportionately more yield than the same
#' deficit during vegetative growth (Steduto et al. 2012). Prefer
#' \code{by.interval = TRUE} with breakpoints matched to your crop's phenology
#' over a single season-long mean.
#'
#' @return
#' A data.frame with, per group: \code{n_days}, \code{mean_Ks}, \code{min_Ks},
#' \code{total_deficit}, \code{n_stress_days},
#' \code{max_consecutive_stress}, \code{total_ETc}, \code{total_ETa},
#' \code{ETa_ETc_ratio}, \code{total_rain}, \code{total_irrigation},
#' \code{total_drainage}, \code{mean_ASW}, \code{mean_depletion}.
#'
#' @references
#' Doorenbos, J. & Kassam, A.H. (1979) Yield response to water. \emph{FAO
#' Irrigation and Drainage Paper 33}. FAO, Rome.
#'
#' Steduto, P., Hsiao, T.C., Fereres, E. & Raes, D. (2012) Crop yield response
#' to water. \emph{FAO Irrigation and Drainage Paper 66}. FAO, Rome.
#'
#' Costa-Neto, G., Fritsche-Neto, R. & Crossa, J. (2021) Nonlinear kernels,
#' dominance, and envirotyping data increase the accuracy of genome-based
#' prediction in multi-environment trials. \emph{Heredity} 126, 92-106.
#' \doi{10.1038/s41437-020-00353-1}
#'
#' @seealso \code{water_balance}, \code{summaryWTH}
#'
#' @examples
#' \dontrun{
#' wb <- water_balance(wth, soil)
#'
#' ## Season-level, one row per environment
#' summary_water_balance(wb)
#'
#' ## By maize growth stage
#' summary_water_balance(wb, by.interval = TRUE,
#'                       time.window  = c(0, 14, 35, 60, 90, 120),
#'                       names.window = c("P-E", "E-V1", "V1-V4",
#'                                        "V4-VT", "VT-GF", "GF-PM"))
#'
#' ## Default 10-day windows when no breakpoints are given
#' summary_water_balance(wb, by.interval = TRUE)
#'
#' ## Rank environments by flowering-period stress
#' s <- summary_water_balance(wb, by.interval = TRUE,
#'                            time.window  = c(0, 60, 90, 120),
#'                            names.window = c("veg", "flower", "fill"))
#' s[s$interval == "flower", ][order(s$mean_Ks[s$interval == "flower"]), ]
#' }
#'
#' @export
summary_water_balance <- function(wb, env.id = "env", days.id = "daysFromStart",
                                  by.interval = FALSE, time.window = NULL,
                                  names.window = NULL, verbose = TRUE) {

  .et_banner("summary_water_balance",
             "summarises a water balance by environment/stage", verbose)
  wb <- as.data.frame(wb)
  need <- c(env.id, "Ks", "deficit", "ETc", "ETa", "ASW", "drainage")
  miss <- setdiff(need, names(wb))
  if (length(miss))
    stop("Not a water_balance() output; missing: ", paste(miss, collapse = ", "),
         call. = FALSE)

  grp_cols <- env.id
  if (isTRUE(by.interval)) {
    if (!days.id %in% names(wb))
      stop("by.interval = TRUE requires column '", days.id, "'.", call. = FALSE)
    brk <- if (is.null(time.window))
      seq(0, max(wb[[days.id]], na.rm = TRUE) + 10, by = 10) else time.window
    cutv <- cut(wb[[days.id]], breaks = c(brk, Inf), right = FALSE)
    nm <- if (is.null(names.window)) paste0("Interval_", seq_len(nlevels(cutv))) else
      names.window
    if (length(nm) != nlevels(cutv))
      nm <- paste0("Interval_", seq_len(nlevels(cutv)))
    levels(cutv) <- nm
    wb$interval <- cutv
    grp_cols <- c(env.id, "interval")
  }

  # Longest run of consecutive stressed days.
  longest_run <- function(flag) {
    flag <- flag & !is.na(flag)
    if (!any(flag)) return(0L)
    r <- rle(flag)
    max(r$lengths[r$values])
  }

  key  <- do.call(paste, c(lapply(grp_cols, function(z) as.character(wb[[z]])),
                           sep = "\r"))
  ukey <- unique(key)

  rows <- lapply(ukey, function(k) {
    s <- wb[key == k, , drop = FALSE]
    etc <- sum(s$ETc, na.rm = TRUE)
    eta <- sum(s$ETa, na.rm = TRUE)
    base <- s[1, grp_cols, drop = FALSE]
    cbind(base, data.frame(
      n_days                 = nrow(s),
      mean_Ks                = mean(s$Ks, na.rm = TRUE),
      min_Ks                 = suppressWarnings(min(s$Ks, na.rm = TRUE)),
      total_deficit          = sum(s$deficit, na.rm = TRUE),
      n_stress_days          = sum(s$Ks < 1, na.rm = TRUE),
      max_consecutive_stress = longest_run(s$Ks < 1),
      total_ETc              = etc,
      total_ETa              = eta,
      ETa_ETc_ratio          = if (etc > 0) eta / etc else NA_real_,
      total_rain             = if ("PRECTOT" %in% names(s))
        sum(s$PRECTOT, na.rm = TRUE) else NA_real_,
      total_irrigation       = if ("irrigation" %in% names(s))
        sum(s$irrigation, na.rm = TRUE) else NA_real_,
      total_drainage         = sum(s$drainage, na.rm = TRUE),
      mean_ASW               = mean(s$ASW, na.rm = TRUE),
      mean_depletion         = mean(s$depletion_frac, na.rm = TRUE),
      stringsAsFactors = FALSE))
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  num <- vapply(out, is.numeric, logical(1))
  out[num] <- lapply(out[num], function(z) round(z, 4))
  out[do.call(order, out[grp_cols]), , drop = FALSE]
}
