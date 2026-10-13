require(EnvRtype)

data(G2F_2014_17)

# first 5 environments, as in your data_G2F.R example
d     <- G2F_2014_17[1:5, ]
env   <- d$env
lat   <- d$lat
lon   <- d$lon
start <- d$start
end   <- d$end


#==================================================================================================
# env_remote_sensing.R
#
# Satellite land-surface covariates for enviromic analysis, retrieved via the
# 'modisfast' package (OPeNDAP subsetting of NASA MODIS / VIIRS / GPM).
#
# Module contents
#   get_remote_sensing()  - point time series of MODIS/VIIRS/GPM variables,
#                           returned in the same long format as get_weather()
#   rs_to_daily()         - composite (8-/16-day) -> daily grid, with an
#                           explicit refusal to interpolate across long gaps
#   rs_qa_mask()          - decode and apply per-collection QA bitfields
#   rs_collections()      - the collections this adapter knows how to handle
#
# Internal helpers (not exported) are prefixed with '.rs_'.
#
# --------------------------------------------------------------------------
# OPTIONAL MODULE -- READ BEFORE USE
# --------------------------------------------------------------------------
# 'modisfast' is in Suggests, not Imports. Unlike get_weather(), this module
# CANNOT run without credentials:
#
#   Sys.setenv(EARTHDATA_TOKEN = "<bearer token>")
#   # obtain at urs.earthdata.nasa.gov
#
# LIMITATIONS
#   * LST is SKIN temperature of the land surface, NOT 2 m air temperature.
#     Over a dry canopy at midday it can exceed T2M by >10 K. It is not a
#     drop-in replacement for POWER T2M and must not be fed to
#     param_temperature() as if it were.
#   * Optical products (LST, reflectance, LAI, VI) are MISSING UNDER CLOUD,
#     and that missingness is NOT random -- it is worst in exactly the humid
#     environments where water is not limiting. Every row carries a 'gap'
#     flag and an 'n_pixel' count. Do not average them away.
#   * 8-day and 16-day composites are NOT daily observations. rs_to_daily()
#     will interpolate, but refuses gaps longer than max_gap and flags every
#     interpolated row.
#   * GPM IMERG is microwave/IR and is gap-free; it is the only collection
#     here that is safe to treat as a complete daily series.
#   * NASA endpoint stability is outside our control. modisfast was off CRAN
#     for roughly a year (2024-2025) following a server migration. Do not
#     put this module on the critical path of a reproducible pipeline
#     without caching the downloads.
#
# --------------------------------------------------------------------------
# PHENOTYPE CONTAMINATION -- WHY CANOPY VARIABLES ARE GATED
# --------------------------------------------------------------------------
# NDVI, EVI, LAI, FPAR and surface reflectance over a trial plot are PARTLY A
# FUNCTION OF THE GENOTYPE GROWING THERE. Using them as environmental
# covariables in env_kernel() -> kernel_model() leaks the response into the
# predictor: cross-validation accuracy inflates and will not replicate at
# untested sites, which is the entire purpose of the reaction-norm model.
#
# Collections flagged 'canopy_state' therefore require
# allow_canopy_state = TRUE. Legitimate uses exist (env_dissection(),
# descriptive characterization, off-plot landscape context) -- the gate asks
# for acknowledgement, it does not forbid.
#
# Forcing variables (LST, ET/PET, precipitation) carry no such risk.
#
# --------------------------------------------------------------------------
# METHOD REFERENCES
# --------------------------------------------------------------------------
# modisfast
#   Taconet, P., Degenne, P., Moiroux, N., Dennis, R. (2024) modisfast: An R
#     package for fast and efficient access to MODIS, VIIRS and GPM Earth
#     Observation data. Journal of Open Source Software.
#
# MODIS LST and its QA
#   Wan, Z. (2014) New refinements and validation of the collection-6
#     MODIS land-surface temperature/emissivity product. Remote Sensing of
#     Environment 140, 36-45. doi:10.1016/j.rse.2013.08.027
#   Wan, Z. (2013) Collection-6 MODIS Land Surface Temperature Products
#     Users' Guide. ERI, University of California, Santa Barbara.
#
# Skin vs air temperature (why LST != T2M)
#   Jin, M. & Dickinson, R.E. (2010) Land surface skin temperature climatology.
#     Environmental Research Letters 5, 044004.
#     doi:10.1088/1748-9326/5/4/044004
#
# MODIS evapotranspiration
#   Mu, Q., Zhao, M. & Running, S.W. (2011) Improvements to a MODIS global
#     terrestrial evapotranspiration algorithm. Remote Sensing of Environment
#     115(8), 1781-1800. doi:10.1016/j.rse.2011.02.019
#
# GPM IMERG
#   Huffman, G.J. et al. (2020) Integrated Multi-satellite Retrievals for the
#     Global Precipitation Measurement (GPM) Mission (IMERG). In: Satellite
#     Precipitation Measurement, Springer. doi:10.1007/978-3-030-24568-9_19
#
# Author: G Costa-Neto and contributors
#==================================================================================================


# =========================================================================
# SECTION 0 - collection registry
# =========================================================================

#' Registry of collections this adapter understands
#'
#' Each entry records what the adapter must know that modisfast does not tell
#' it: the scale factor, the physical unit, the compositing period, whether
#' the variable describes canopy state (and is therefore gated), and which QA
#' decoder applies.
#'
#' \code{period} is the nominal compositing interval in days; 1 means daily.
#' \code{valid} is the physically plausible range AFTER scaling, used to catch
#' undocumented fill values that slipped past the QA mask.
#'
#' @keywords internal
#' @noRd
.rs_registry <- function() {
  list(
    "MOD11A1.061" = list(
      label = "MODIS/Terra LST daily 1 km",
      period = 1L, canopy_state = FALSE, qa = "lst",
      vars = list(
        LST_Day_1km   = list(scale = 0.02, offset = 0, unit = "K",
                             valid = c(150, 350), qa_layer = "QC_Day"),
        LST_Night_1km = list(scale = 0.02, offset = 0, unit = "K",
                             valid = c(150, 350), qa_layer = "QC_Night"))),
    "MYD11A1.061" = list(
      label = "MODIS/Aqua LST daily 1 km",
      period = 1L, canopy_state = FALSE, qa = "lst",
      vars = list(
        LST_Day_1km   = list(scale = 0.02, offset = 0, unit = "K",
                             valid = c(150, 350), qa_layer = "QC_Day"),
        LST_Night_1km = list(scale = 0.02, offset = 0, unit = "K",
                             valid = c(150, 350), qa_layer = "QC_Night"))),
    "MOD11A2.061" = list(
      label = "MODIS/Terra LST 8-day 1 km",
      period = 8L, canopy_state = FALSE, qa = "lst",
      vars = list(
        LST_Day_1km   = list(scale = 0.02, offset = 0, unit = "K",
                             valid = c(150, 350), qa_layer = "QC_Day"),
        LST_Night_1km = list(scale = 0.02, offset = 0, unit = "K",
                             valid = c(150, 350), qa_layer = "QC_Night"))),
    "MOD16A2GF.061" = list(
      label = "MODIS/Terra ET gap-filled 8-day 500 m",
      period = 8L, canopy_state = FALSE, qa = "et",
      vars = list(
        ET_500m  = list(scale = 0.1, offset = 0, unit = "mm/8day",
                        valid = c(0, 1000), qa_layer = "ET_QC_500m"),
        PET_500m = list(scale = 0.1, offset = 0, unit = "mm/8day",
                        valid = c(0, 1000), qa_layer = "ET_QC_500m"))),
    "MYD16A2GF.061" = list(
      label = "MODIS/Aqua ET gap-filled 8-day 500 m",
      period = 8L, canopy_state = FALSE, qa = "et",
      vars = list(
        ET_500m  = list(scale = 0.1, offset = 0, unit = "mm/8day",
                        valid = c(0, 1000), qa_layer = "ET_QC_500m"),
        PET_500m = list(scale = 0.1, offset = 0, unit = "mm/8day",
                        valid = c(0, 1000), qa_layer = "ET_QC_500m"))),
    "GPM_3IMERGDF.07" = list(
      label = "GPM IMERG Final daily 0.1 deg",
      period = 1L, canopy_state = FALSE, qa = "none",
      vars = list(
        precipitation    = list(scale = 1, offset = 0, unit = "mm/day",
                                valid = c(0, 2000), qa_layer = NA_character_),
        precipitationCal = list(scale = 1, offset = 0, unit = "mm/day",
                                valid = c(0, 2000), qa_layer = NA_character_))),
    "GPM_3IMERGDL.06" = list(
      label = "GPM IMERG Late daily 0.1 deg",
      period = 1L, canopy_state = FALSE, qa = "none",
      vars = list(
        precipitationCal = list(scale = 1, offset = 0, unit = "mm/day",
                                valid = c(0, 2000), qa_layer = NA_character_))),
    "MCD43A4.061" = list(
      label = "MODIS NBAR daily 500 m (CANOPY STATE)",
      period = 1L, canopy_state = TRUE, qa = "nbar",
      vars = stats::setNames(
        lapply(1:7, function(i)
          list(scale = 1e-4, offset = 0, unit = "reflectance",
               valid = c(-0.01, 1.6),
               qa_layer = paste0("BRDF_Albedo_Band_Mandatory_Quality_Band", i))),
        paste0("Nadir_Reflectance_Band", 1:7))),
    "VNP15A2H.002" = list(
      label = "VIIRS/NPP LAI-FPAR 8-day 500 m (CANOPY STATE)",
      period = 8L, canopy_state = TRUE, qa = "none",
      vars = list(
        LAI  = list(scale = 0.1,  offset = 0, unit = "m2/m2",
                    valid = c(0, 10), qa_layer = NA_character_),
        FPAR = list(scale = 0.01, offset = 0, unit = "fraction",
                    valid = c(0, 1),  qa_layer = NA_character_)))
  )
}


#' @title Collections Supported by get_remote_sensing()
#'
#' @description
#' Lists the MODIS / VIIRS / GPM collections for which this adapter knows the
#' scale factors, units, compositing period and QA decoding. \pkg{modisfast}
#' supports many more; those listed here are the ones wired into the
#' EnvRtype pipeline.
#'
#' @author Germano Costa Neto
#'
#' @return
#' A data.frame with one row per collection-variable pair: \code{collection},
#' \code{variable}, \code{unit}, \code{period_days}, \code{canopy_state}.
#'
#' @details
#' \code{canopy_state = TRUE} marks variables that describe the state of the
#' vegetation rather than the environmental forcing. Those require
#' \code{allow_canopy_state = TRUE} in \code{get_remote_sensing()}; see the
#' contamination note in that function's documentation.
#'
#' Run \code{modisfast::mf_list_collections()} for the full catalogue, and
#' \code{modisfast::mf_list_variables(<collection>)} for the bands of any one
#' collection.
#'
#' @examples
#' rs_collections()
#'
#' @seealso \code{\link{get_remote_sensing}}
#'
#' @export
rs_collections <- function() {
  reg <- .rs_registry()
  out <- do.call(rbind, lapply(names(reg), function(cn) {
    cc <- reg[[cn]]
    data.frame(collection   = cn,
               variable     = names(cc$vars),
               unit         = vapply(cc$vars, function(v) v$unit, character(1)),
               period_days  = cc$period,
               canopy_state = cc$canopy_state,
               stringsAsFactors = FALSE)
  }))
  rownames(out) <- NULL
  out
}


# =========================================================================
# SECTION 1 - internal helpers
# =========================================================================

#' Build an sf polygon ROI from point sites plus a buffer in kilometres
#'
#' A fixed buffer in DEGREES is wrong everywhere except the equator: one
#' degree of longitude is 111 km at the equator and 56 km at 60 deg N. The
#' longitude half-width is therefore divided by cos(latitude).
#'
#' Sites are returned as ONE polygon per site rather than a single bounding
#' box over all of them. A bounding box over continental-scale trial networks
#' would request a data cube orders of magnitude larger than needed, which is
#' precisely the waste modisfast exists to avoid.
#'
#' @param env.id character. Site identifiers.
#' @param lat,lon numeric. Coordinates in WGS84 decimal degrees.
#' @param buffer_km numeric. Half-width of the square window, in km.
#' @return An sf object with columns \code{id} and \code{geom}.
#' @keywords internal
#' @noRd
.rs_sites_to_roi <- function(env.id, lat, lon, buffer_km) {
  .need_pkg("sf", "get_remote_sensing")

  if (!is.numeric(buffer_km) || length(buffer_km) != 1L ||
      is.na(buffer_km) || buffer_km <= 0)
    stop("'buffer_km' must be a single positive number.", call. = FALSE)

  km_per_deg_lat <- 111.32
  dlat <- buffer_km / km_per_deg_lat
  # cos() guarded: at |lat| -> 90 the scaling diverges.
  coslat <- pmax(cos(lat * pi / 180), 0.01)
  dlon <- buffer_km / (km_per_deg_lat * coslat)

  ymin <- pmax(lat - dlat, -90);  ymax <- pmin(lat + dlat, 90)
  xmin <- pmax(lon - dlon, -180); xmax <- pmin(lon + dlon, 180)

  polys <- lapply(seq_along(lat), function(i) {
    sf::st_polygon(list(cbind(
      c(xmin[i], xmax[i], xmax[i], xmin[i], xmin[i]),
      c(ymin[i], ymin[i], ymax[i], ymax[i], ymin[i]))))
  })

  sf::st_sf(id = as.character(env.id),
            geom = sf::st_sfc(polys, crs = 4326))
}


#' Extract bits [hi:lo] from an integer QA word
#'
#' MODIS QA layers pack several independent flags into one 8- or 16-bit
#' integer. bitwAnd/bitwShiftR are base R and integer-safe for the widths
#' used here.
#'
#' @param x integer vector. Packed QA values.
#' @param lo,hi integer. Inclusive bit positions, 0-indexed from the LSB.
#' @return Integer vector of the extracted field.
#' @keywords internal
#' @noRd
.rs_bits <- function(x, lo, hi) {
  if (hi < lo) stop("'hi' must be >= 'lo'.", call. = FALSE)
  width <- hi - lo + 1L
  mask  <- bitwShiftL(1L, width) - 1L
  xi <- suppressWarnings(as.integer(x))
  out <- bitwAnd(bitwShiftR(xi, lo), mask)
  out[is.na(xi)] <- NA_integer_
  out
}


#' Decode the MOD11A1 / MOD11A2 QC_Day and QC_Night bitfields
#'
#' Collection-6 LST QC layout (Wan 2013, Users' Guide Table 5):
#'   bits 0-1  mandatory QA : 0 good, 1 other-quality, 2 not produced (cloud),
#'                            3 not produced (other)
#'   bits 2-3  data quality : 0 good, 1 other quality
#'   bits 4-5  emissivity error
#'   bits 6-7  LST error     : 0 <=1 K, 1 <=2 K, 2 <=3 K, 3 >3 K
#'
#' Default keeps pixels with mandatory QA 0 and LST error <= 2 K. Accepting
#' mandatory QA 1 ("other quality") is the single most common way MODIS LST
#' time series acquire a warm bias, because partially cloudy pixels are
#' retained; it is allowed only by raising max_qa explicitly.
#'
#' @param qa integer vector. Packed QC values.
#' @param max_qa integer. Largest acceptable mandatory-QA code (0 or 1).
#' @param max_lst_err integer. Largest acceptable LST-error class (0-3).
#' @return Logical vector, TRUE where the pixel is acceptable.
#' @keywords internal
#' @noRd
.rs_qa_lst <- function(qa, max_qa = 0L, max_lst_err = 1L) {
  mandatory <- .rs_bits(qa, 0L, 1L)
  lst_err   <- .rs_bits(qa, 6L, 7L)
  ok <- !is.na(mandatory) & mandatory <= max_qa &
    !is.na(lst_err)   & lst_err   <= max_lst_err
  ok[is.na(qa)] <- FALSE
  ok
}


#' Decode the MOD16A2 ET_QC_500m bitfield
#'
#' MOD16 QC layout (Mu et al. 2011; MOD16 Users' Guide):
#'   bit 0     MODLAND : 0 good, 1 other
#'   bits 3-4  cloud state : 0 clear, 1 cloudy, 2 mixed, 3 unknown
#'   bits 5-7  SCF_QC   : 0-1 main method, higher = degraded retrieval
#'
#' Note that MOD16A2GF is the GAP-FILLED product: cloudy intervals have
#' already been filled by the algorithm using the annual LAI/FPAR cycle. The
#' QC layer still records that they were filled, so masking strictly here
#' removes the gap-filled values and reintroduces the gaps. That is usually
#' what you want when validating a water balance against observations; set
#' max_scf = 3 to keep the filled values.
#'
#' @param qa integer vector. Packed QC values.
#' @param max_scf integer. Largest acceptable SCF_QC class (0-7).
#' @param allow_cloudy logical. Keep pixels flagged cloudy or mixed.
#' @return Logical vector, TRUE where the pixel is acceptable.
#' @keywords internal
#' @noRd
.rs_qa_et <- function(qa, max_scf = 1L, allow_cloudy = FALSE) {
  modland <- .rs_bits(qa, 0L, 0L)
  cloud   <- .rs_bits(qa, 3L, 4L)
  scf     <- .rs_bits(qa, 5L, 7L)
  ok <- !is.na(modland) & modland == 0L & !is.na(scf) & scf <= max_scf
  if (!isTRUE(allow_cloudy)) ok <- ok & !is.na(cloud) & cloud == 0L
  ok[is.na(qa)] <- FALSE
  ok
}


#' Decode an MCD43A4 mandatory-quality layer
#'
#' BRDF/Albedo mandatory quality: 0 = full BRDF inversion (best),
#' 1 = magnitude inversion (lower quality). Values above 1 are fill.
#'
#' @param qa integer vector.
#' @param max_qa integer. Largest acceptable code.
#' @return Logical vector.
#' @keywords internal
#' @noRd
.rs_qa_nbar <- function(qa, max_qa = 0L) {
  v <- suppressWarnings(as.integer(qa))
  ok <- !is.na(v) & v <= max_qa
  ok
}


#' Dispatch to the QA decoder named in the registry
#' @keywords internal
#' @noRd
.rs_qa_dispatch <- function(kind, qa, ...) {
  switch(kind,
         lst  = .rs_qa_lst(qa, ...),
         et   = .rs_qa_et(qa, ...),
         nbar = .rs_qa_nbar(qa, ...),
         none = rep(TRUE, length(qa)),
         stop("Unknown QA decoder: '", kind, "'.", call. = FALSE))
}


#' Apply scale/offset and clamp to the physically valid range
#'
#' Fill values that the QA mask did not catch are set to NA rather than
#' silently clamped: a clamped fill value is indistinguishable from a real
#' extreme and would corrupt any statistic computed downstream.
#'
#' @param x numeric vector of raw DN values.
#' @param spec list. Registry entry for the variable.
#' @return Numeric vector in physical units.
#' @keywords internal
#' @noRd
.rs_scale <- function(x, spec) {
  v <- as.numeric(x) * spec$scale + spec$offset
  bad <- !is.na(v) & (v < spec$valid[1] | v > spec$valid[2])
  v[bad] <- NA_real_
  v
}


#' Parse acquisition dates from modisfast SpatRaster layer names
#'
#' modisfast names layers with the variable and the date, but the exact
#' pattern varies by collection and version. Rather than depend on one
#' pattern, this scans for the first recognisable date token in each name:
#' YYYY-MM-DD, YYYYMMDD, or a MODIS AYYYYDDD julian stamp.
#'
#' @param nms character vector of layer names.
#' @return Date vector, NA where no date could be parsed.
#' @keywords internal
#' @noRd
.rs_parse_dates <- function(nms) {
  out <- rep(as.Date(NA), length(nms))

  m <- regmatches(nms, regexpr("[0-9]{4}-[0-9]{2}-[0-9]{2}", nms))
  hit <- regexpr("[0-9]{4}-[0-9]{2}-[0-9]{2}", nms) > 0
  if (any(hit)) out[hit] <- as.Date(m, format = "%Y-%m-%d")

  todo <- is.na(out)
  if (any(todo)) {
    j <- regexpr("A[0-9]{7}", nms)
    hit2 <- todo & j > 0
    if (any(hit2)) {
      tok <- regmatches(nms[hit2], regexpr("A[0-9]{7}", nms[hit2]))
      out[hit2] <- as.Date(substr(tok, 2, 8), format = "%Y%j")
    }
  }

  todo <- is.na(out)
  if (any(todo)) {
    k <- regexpr("(?<![0-9])[0-9]{8}(?![0-9])", nms, perl = TRUE)
    hit3 <- todo & k > 0
    if (any(hit3)) {
      tok <- regmatches(nms[hit3],
                        regexpr("(?<![0-9])[0-9]{8}(?![0-9])", nms[hit3],
                                perl = TRUE))
      out[hit3] <- as.Date(tok, format = "%Y%m%d")
    }
  }
  out
}


#' Summarise the pixels of one site-date window into one row
#'
#' @param vals numeric vector of QA-accepted, scaled pixel values.
#' @param n_total integer. Pixels in the window before masking.
#' @param fun function. Aggregator applied across pixels.
#' @return A one-row list with value, n_pixel, n_total and gap.
#' @keywords internal
#' @noRd
.rs_aggregate_window <- function(vals, n_total, fun = mean) {
  keep <- vals[!is.na(vals)]
  if (!length(keep))
    return(list(value = NA_real_, n_pixel = 0L,
                n_total = as.integer(n_total), gap = TRUE))
  list(value = as.numeric(fun(keep)), n_pixel = length(keep),
       n_total = as.integer(n_total), gap = FALSE)
}


# =========================================================================
# SECTION 2 - rs_qa_mask()
# =========================================================================

#' @title Apply MODIS / VIIRS Quality Bitfields to a Raster Stack
#'
#' @description
#' Masks a \code{SpatRaster} of science values using the collection's quality
#' layers, decoding the packed bitfields according to the product's users'
#' guide. Pixels failing the quality test are set to \code{NA}.
#'
#' @author Germano Costa Neto
#'
#' @param r SpatRaster containing both the science layers and their QA layers,
#'   as returned by \code{modisfast::mf_import_data()}.
#' @param collection character. Collection id, e.g. \code{"MOD11A1.061"}.
#' @param variables character. Science variables to mask. Default: all
#'   registry variables present in \code{r}.
#' @param strict boolean. If \code{TRUE} (default) only the highest quality
#'   class is retained. If \code{FALSE} the thresholds are relaxed by one
#'   class, which roughly doubles data volume at the cost of a warm bias in
#'   LST and a wet bias in ET.
#' @param verbose boolean. Report the rejection rate per layer.
#'
#' @return
#' A \code{SpatRaster} with the same science layers, quality-masked. QA layers
#' are dropped from the result.
#'
#' @details
#' Decoders implemented: MODIS LST (\code{QC_Day} / \code{QC_Night}), MOD16
#' evapotranspiration (\code{ET_QC_500m}), and MCD43A4 BRDF mandatory quality.
#' GPM IMERG carries no per-pixel QA and is passed through unchanged.
#'
#' When a QA layer named in the registry is absent from \code{r} the science
#' layer is passed through unmasked and a warning is raised. Silently
#' returning unmasked data would be worse: cloud-contaminated LST looks
#' entirely plausible.
#'
#' @examples
#' \dontrun{
#' r  <- modisfast::mf_import_data(path = d, collection = "MOD11A1.061")
#' rm <- rs_qa_mask(r, collection = "MOD11A1.061")
#' }
#'
#' @seealso \code{\link{get_remote_sensing}}
#'
#' @export
rs_qa_mask <- function(r, collection, variables = NULL, strict = TRUE,
                       verbose = TRUE) {
  .need_pkg("terra", "rs_qa_mask")

  reg <- .rs_registry()
  if (!collection %in% names(reg))
    stop("Collection '", collection, "' is not in the adapter registry.\n",
         "  See rs_collections() for supported collections.", call. = FALSE)
  cc <- reg[[collection]]
  if (identical(cc$qa, "none")) {
    if (verbose) .et_step(paste0(collection, " carries no per-pixel QA"), verbose)
    return(r)
  }

  nms <- names(r)
  if (is.null(variables)) variables <- names(cc$vars)
  variables <- intersect(variables, names(cc$vars))

  keep_layers <- list()
  for (v in variables) {
    sci_idx <- grep(v, nms, fixed = TRUE)
    # Exclude QA layers whose names contain the variable name.
    sci_idx <- sci_idx[!grepl("QC|Quality", nms[sci_idx])]
    if (!length(sci_idx)) next

    qa_name <- cc$vars[[v]]$qa_layer
    qa_idx  <- if (!is.na(qa_name)) grep(qa_name, nms, fixed = TRUE) else integer(0)

    sci <- r[[sci_idx]]
    if (!length(qa_idx)) {
      warning("QA layer '", qa_name, "' not found for variable '", v,
              "'; returning UNMASKED values. Cloud-contaminated pixels are ",
              "retained.", call. = FALSE)
      keep_layers[[v]] <- sci
      next
    }
    qa <- r[[qa_idx]]
    if (terra::nlyr(qa) != terra::nlyr(sci)) {
      warning("QA layer count (", terra::nlyr(qa), ") does not match science ",
              "layer count (", terra::nlyr(sci), ") for '", v,
              "'; returning UNMASKED values.", call. = FALSE)
      keep_layers[[v]] <- sci
      next
    }

    args <- switch(cc$qa,
                   lst  = list(max_qa = if (strict) 0L else 1L,
                               max_lst_err = if (strict) 1L else 2L),
                   et   = list(max_scf = if (strict) 1L else 3L,
                               allow_cloudy = !strict),
                   nbar = list(max_qa = if (strict) 0L else 1L),
                   list())

    ok <- terra::app(qa, function(x)
      as.numeric(do.call(.rs_qa_dispatch, c(list(cc$qa, x), args))))
    masked <- terra::mask(sci, ok, maskvalues = c(0, NA))

    if (verbose) {
      tot <- terra::global(!is.na(sci), "sum", na.rm = TRUE)[1, 1]
      kep <- terra::global(!is.na(masked), "sum", na.rm = TRUE)[1, 1]
      pct <- if (is.finite(tot) && tot > 0) 100 * (1 - kep / tot) else NA_real_
      .et_step(sprintf("%-22s QA rejected %5.1f%% of pixels", v, pct), verbose)
    }
    keep_layers[[v]] <- masked
  }

  if (!length(keep_layers))
    stop("No science layers matched the registry for collection '",
         collection, "'.", call. = FALSE)
  out <- terra::rast(keep_layers)
  out
}


# =========================================================================
# SECTION 3 - get_remote_sensing()
# =========================================================================

#' @title Collect Satellite Land-Surface Covariates for Trial Sites
#'
#' @description
#' Downloads MODIS, VIIRS or GPM data for a set of trial sites via
#' \pkg{modisfast} OPeNDAP subsetting, applies the product quality mask,
#' aggregates the pixels within a buffer around each site, and returns a long
#' data.frame in the same shape as \code{\link{get_weather}} so that the
#' result flows directly into \code{\link{W_matrix}},
#' \code{\link{env_typing}} and the kernel pipeline.
#'
#' @author Germano Costa Neto
#'
#' @param env.id vector (character). Identification of the site/environment.
#' @param lat vector (numeric). Latitude in WGS84 decimal degrees.
#' @param lon vector (numeric). Longitude in WGS84 decimal degrees.
#' @param start.day vector (character or Date). First day, \code{YYYY-MM-DD}.
#'   Length 1 or \code{length(lat)}.
#' @param end.day vector (character or Date). Last day. Length 1 or
#'   \code{length(lat)}.
#' @param collection character. Collection id; see \code{\link{rs_collections}}.
#' @param variables character. Bands to retrieve. Default: all registry
#'   variables for the collection.
#' @param buffer_km numeric. Half-width in km of the square window centred on
#'   each site. Must exceed the pixel size of the product, otherwise the
#'   window may contain a single pixel and \code{n_pixel} loses its meaning.
#' @param aggregate function. Applied across the valid pixels of a window.
#'   Default \code{mean}. \code{median} is more robust where the window
#'   straddles a land-cover boundary.
#' @param qa_mask boolean. Apply the product quality bitfield. Default
#'   \code{TRUE}. Setting \code{FALSE} retains cloud-contaminated pixels.
#' @param strict_qa boolean. Passed to \code{\link{rs_qa_mask}}.
#' @param allow_canopy_state boolean. Required to be \code{TRUE} for
#'   collections describing vegetation state. See Details.
#' @param parallel boolean. Parallel downloads inside \pkg{modisfast}.
#' @param num_workers integer. Download workers when \code{parallel = TRUE}.
#' @param save_dir character. Directory for the downloaded subsets. Defaults
#'   to a session temporary directory; supply a persistent path to cache.
#' @param verbose boolean. Print progress.
#'
#' @return
#' A data.frame with one row per environment-date:
#' \describe{
#'   \item{env}{environment id}
#'   \item{LAT, LON}{site coordinates}
#'   \item{YYYYMMDD}{acquisition date (Date)}
#'   \item{daysFromStart}{days since \code{start.day} for that site, 1-based}
#'   \item{<variables>}{one column per requested band, in physical units}
#'   \item{n_pixel}{valid pixels after QA masking}
#'   \item{n_total}{pixels in the window before masking}
#'   \item{gap}{logical; \code{TRUE} when no valid pixel was available}
#' }
#'
#' @details
#' \strong{Credentials.} An Earthdata bearer token is mandatory:
#' \code{Sys.setenv(EARTHDATA_TOKEN = "...")}, obtainable at
#' urs.earthdata.nasa.gov. Unlike \code{get_weather()}, this
#' function cannot run anonymously.
#'
#' \strong{LST is not air temperature.} Land surface temperature is the
#' radiometric skin temperature of the canopy and soil. Over a dry surface at
#' midday it routinely exceeds 2 m air temperature by more than 10 K. Do not
#' substitute it for \code{T2M} in \code{param_temperature()} or in any
#' cardinal-temperature envirotyping without a local calibration.
#'
#' \strong{Cloud gaps are not random.} Optical retrievals fail under cloud,
#' and cloud is most persistent in humid environments. A season mean computed
#' over clear days only is biased toward dry, hot conditions. The \code{gap}
#' and \code{n_pixel} columns exist so that this bias stays visible; inspect
#' \code{tapply(x$gap, x$env, mean)} before using the data.
#'
#' \strong{Phenotype contamination.} Vegetation indices, LAI, FPAR and surface
#' reflectance measured over a trial plot are partly determined by the
#' genotype growing there. Using them as environmental covariables in
#' \code{\link{env_kernel}} leaks the response into the predictor and inflates
#' cross-validation accuracy that will not replicate at untested sites. Such
#' collections are gated behind \code{allow_canopy_state = TRUE}. Forcing
#' variables -- LST, ET, precipitation -- are not gated.
#'
#' @examples
#' \dontrun{
#' Sys.setenv(EARTHDATA_TOKEN = "<your token>")
#'
#' data(G2F_2014_17)
#' d <- G2F_2014_17[1:5, ]
#'
#' lst <- get_remote_sensing(env.id = d$env, lat = d$lat, lon = d$lon,
#'                           start.day = d$start, end.day = d$end,
#'                           collection = "MOD11A1.061",
#'                           variables = c("LST_Day_1km", "LST_Night_1km"),
#'                           buffer_km = 2)
#'
#' # always inspect the gap fraction before using optical data
#' tapply(lst$gap, lst$env, mean)
#' }
#'
#' @seealso \code{\link{rs_collections}}, \code{\link{rs_to_daily}},
#'   \code{\link{rs_qa_mask}}, \code{\link{get_weather}}
#'
#' @export
get_remote_sensing <- function(env.id = NULL, lat = NULL, lon = NULL,
                               start.day = NULL, end.day = NULL,
                               collection = "MOD11A1.061",
                               variables = NULL,
                               buffer_km = 2,
                               aggregate = mean,
                               qa_mask = TRUE,
                               strict_qa = TRUE,
                               allow_canopy_state = FALSE,
                               parallel = TRUE, num_workers = 3L,
                               save_dir = NULL,
                               verbose = TRUE) {

 # .et_banner("get_remote_sensing",
 #            "satellite land-surface covariates via modisfast", verbose)
#  .need_pkg(c("modisfast", "terra", "sf"), "get_remote_sensing")

  # ---- credentials ---------------------------------------------------------
  if (!nzchar(Sys.getenv("EARTHDATA_TOKEN")))
    stop("An Earthdata bearer token is required.\n",
         "  1. Register / sign in at urs.earthdata.nasa.gov\n",
         "  2. Generate a token, then:\n",
         '     Sys.setenv(EARTHDATA_TOKEN = "<your token>")',
         call. = FALSE)

  # ---- registry and gating -------------------------------------------------
  reg <- .rs_registry()
  if (!collection %in% names(reg))
    stop("Collection '", collection, "' is not supported by this adapter.\n",
         "  Supported: ", paste(names(reg), collapse = ", "), "\n",
         "  See rs_collections() for the variable-level listing.", call. = FALSE)
  cc <- reg[[collection]]

  if (isTRUE(cc$canopy_state) && !isTRUE(allow_canopy_state))
    stop("Collection '", collection, "' (", cc$label, ") yields CANOPY-STATE\n",
         "  variables. Over a trial plot these are partly a function of the\n",
         "  GENOTYPE being grown, so using them as environmental covariables\n",
         "  in env_kernel() leaks the response into the predictor and inflates\n",
         "  cross-validation accuracy that will not replicate at untested\n",
         "  environments.\n\n",
         "  Legitimate uses exist (env_dissection(), descriptive\n",
         "  characterization, off-plot landscape context). To proceed:\n",
         "    get_remote_sensing(..., allow_canopy_state = TRUE)",
         call. = FALSE)

  if (is.null(variables)) variables <- names(cc$vars)
  unknown <- setdiff(variables, names(cc$vars))
  if (length(unknown))
    stop("Variable(s) not in the registry for '", collection, "': ",
         paste(unknown, collapse = ", "), "\n",
         "  Known: ", paste(names(cc$vars), collapse = ", "), "\n",
         "  (modisfast::mf_list_variables('", collection,
         "') lists every band the product carries; only those with a known\n",
         "  scale factor and QA decoder are wired into this adapter.)",
         call. = FALSE)

  # ---- argument validation, mirroring get_weather() ------------------------
  if (is.null(lat) || is.null(lon))
    stop("'lat' and 'lon' are required.", call. = FALSE)
  if (length(lat) != length(lon))
    stop("'lat' and 'lon' must have the same length.", call. = FALSE)

  n_env <- length(lat)
  if (is.null(env.id)) env.id <- paste0("env", seq_len(n_env))
  env.id <- as.character(env.id)
  if (length(env.id) != n_env)
    stop("'env.id' must have the same length as 'lat' and 'lon'.", call. = FALSE)
  if (anyDuplicated(env.id))
    warning("Duplicated env.id values; downstream joins may be ambiguous.",
            call. = FALSE)

  if (any(lat < -90 | lat > 90, na.rm = TRUE))
    stop("'lat' outside [-90, 90]. Check lat/lon are not swapped.", call. = FALSE)
  if (any(lon < -180 | lon > 180, na.rm = TRUE))
    stop("'lon' outside [-180, 180]. Check lat/lon are not swapped.", call. = FALSE)
  if (anyNA(lat) || anyNA(lon))
    stop("Coordinates contain NA values.", call. = FALSE)

  if (is.null(start.day) || is.null(end.day))
    stop("'start.day' and 'end.day' are required.", call. = FALSE)
  start.day <- as.Date(start.day); end.day <- as.Date(end.day)
  if (length(start.day) == 1L) start.day <- rep(start.day, n_env)
  if (length(end.day)   == 1L) end.day   <- rep(end.day,   n_env)
  if (length(start.day) != n_env || length(end.day) != n_env)
    stop("'start.day' and 'end.day' must be length 1 or length(lat).",
         call. = FALSE)
  if (any(end.day < start.day))
    stop("'end.day' precedes 'start.day' for: ",
         paste(env.id[end.day < start.day], collapse = ", "), call. = FALSE)

  if (is.null(save_dir)) save_dir <- file.path(tempdir(), "EnvRtype_rs")
  dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)

  startTime <- Sys.time()
  if (verbose) {
    message("---------------------------------------------------------------")
    message("get_remote_sensing() - satellite covariates via modisfast")
    message("github.com/ptaconet/modisfast")
    message("---------------------------------------------------------------")
    message("Start .............................. ",
            format(startTime, "%a %b %d %X %Y"))
    message("Collection ......................... ", collection)
    message("                                     ", cc$label)
    message("Variables .......................... ",
            paste(variables, collapse = ", "))
    message("Environmental units ................ ", n_env)
    message("Buffer (km) ........................ ", buffer_km)
    message("QA mask ............................ [ ",
            if (qa_mask) "x" else " ", " ]")
    message("Parallel ........................... [ ",
            if (parallel) "x" else " ", " ]")
  }

  # ---- QA layers must be downloaded alongside the science bands ------------
  qa_needed <- unique(stats::na.omit(
    vapply(cc$vars[variables], function(v) v$qa_layer, character(1))))
  dl_vars <- unique(c(variables, if (isTRUE(qa_mask)) qa_needed))

  # ---- per-site retrieval --------------------------------------------------
  # One modisfast call per site. Sites in a trial network are far apart; a
  # single ROI spanning all of them would request the convex hull of the
  # network, defeating the purpose of OPeNDAP subsetting.
  results <- vector("list", n_env)

  for (i in seq_len(n_env)) {
    roi_i <- .rs_sites_to_roi(env.id[i], lat[i], lon[i], buffer_km)

    res_i <- tryCatch({
      urls <- modisfast::mf_get_url(
        collection = collection,
        variables  = dl_vars,
        roi        = roi_i,
        time_range = c(start.day[i], end.day[i]))

      dl <- modisfast::mf_download_data(
        urls, path = save_dir,
        parallel = parallel, num_workers = num_workers)

      r <- modisfast::mf_import_data(
        path       = dirname(dl$destfile[1]),
        collection = collection,
        proj_epsg  = 4326)

      if (isTRUE(qa_mask))
        r <- rs_qa_mask(r, collection = collection, variables = variables,
                        strict = strict_qa, verbose = FALSE)

      .rs_extract_site(r, env.id[i], lat[i], lon[i], start.day[i],
                       collection, variables, aggregate)
    },
    error = function(e) {
      warning("Site '", env.id[i], "' failed: ", conditionMessage(e),
              call. = FALSE)
      NULL
    })

    results[[i]] <- res_i
    if (verbose) {
      ok <- !is.null(res_i)
      message("  [", i, "/", n_env, "] ", env.id[i], " ",
              if (ok) paste0("-> ", nrow(res_i), " dates") else "FAILED")
    }
  }

  results <- results[!vapply(results, is.null, logical(1))]
  if (!length(results))
    stop("No data returned for any site. Check the token, the collection id ",
         "and the date range against the product's temporal extent ",
         "(modisfast::mf_list_collections()).", call. = FALSE)

  out <- .rbind_fill(results)
  out <- out[order(out$env, out$YYYYMMDD), , drop = FALSE]
  rownames(out) <- NULL

  endTime <- Sys.time()
  if (verbose) {
    gapfrac <- mean(out$gap, na.rm = TRUE)
    message("Environmental units retrieved ...... ",
            length(unique(out$env)), "/", n_env)
    message("Rows ............................... ", nrow(out))
    message("Gap fraction (no valid pixel) ...... ",
            sprintf("%.1f%%", 100 * gapfrac))
    if (is.finite(gapfrac) && gapfrac > 0.5)
      message("  NOTE: more than half the series is missing. Cloud loss is ",
              "non-random;\n        season means over the remaining days are ",
              "biased toward clear, dry weather.")
    message("Done ............................... ",
            format(endTime, "%a %b %d %X %Y"))
    message("Total time ......................... ",
            round(difftime(endTime, startTime, units = "secs"), 2), "s")
  }
  out
}


#' Extract and aggregate one site's raster stack into long format
#'
#' Split out of get_remote_sensing() so it can be unit-tested against a
#' synthetic SpatRaster with no network access.
#'
#' @param r SpatRaster, QA-masked, covering one site.
#' @param env_id character(1), lat_i/lon_i numeric(1), start_i Date(1).
#' @param collection character(1). Registry key.
#' @param variables character. Science variables to extract.
#' @param aggregate function applied across pixels.
#' @return data.frame in the documented long format.
#' @keywords internal
#' @noRd
.rs_extract_site <- function(r, env_id, lat_i, lon_i, start_i,
                             collection, variables, aggregate = mean) {
  .need_pkg("terra", "get_remote_sensing")
  cc <- .rs_registry()[[collection]]
  nms <- names(r)

  per_var <- list()
  for (v in variables) {
    idx <- grep(v, nms, fixed = TRUE)
    idx <- idx[!grepl("QC|Quality", nms[idx])]
    if (!length(idx)) next

    dates <- .rs_parse_dates(nms[idx])
    spec  <- cc$vars[[v]]

    # terra::values() -> matrix [cell, layer]; whole window is the site.
    vals <- terra::values(r[[idx]], mat = TRUE)
    n_total <- nrow(vals)

    rows <- lapply(seq_along(idx), function(j) {
      scaled <- .rs_scale(vals[, j], spec)
      .rs_aggregate_window(scaled, n_total, fun = aggregate)
    })

    per_var[[v]] <- data.frame(
      YYYYMMDD = dates,
      value    = vapply(rows, function(z) z$value,   numeric(1)),
      n_pixel  = vapply(rows, function(z) z$n_pixel, integer(1)),
      n_total  = vapply(rows, function(z) z$n_total, integer(1)),
      gap      = vapply(rows, function(z) z$gap,     logical(1)),
      stringsAsFactors = FALSE)
    names(per_var[[v]])[2] <- v
  }

  if (!length(per_var))
    stop("No requested variable found among layers: ",
         paste(utils::head(nms, 10), collapse = ", "), call. = FALSE)

  # Merge variables on date; n_pixel/gap are taken from the first variable
  # and recomputed as the across-variable worst case.
  base <- per_var[[1]][, c("YYYYMMDD", names(per_var)[1], "n_pixel",
                           "n_total", "gap"), drop = FALSE]
  if (length(per_var) > 1L) {
    for (k in 2:length(per_var)) {
      nm <- names(per_var)[k]
      add <- per_var[[k]][, c("YYYYMMDD", nm, "n_pixel", "gap"), drop = FALSE]
      names(add)[names(add) == "n_pixel"] <- ".np2"
      names(add)[names(add) == "gap"]     <- ".gap2"
      base <- merge(base, add, by = "YYYYMMDD", all = TRUE)
      base$n_pixel <- pmin(base$n_pixel, base$.np2, na.rm = TRUE)
      base$gap     <- base$gap | base$.gap2
      base$.np2 <- NULL; base$.gap2 <- NULL
    }
  }

  base <- base[!is.na(base$YYYYMMDD), , drop = FALSE]
  base <- base[order(base$YYYYMMDD), , drop = FALSE]

  out <- data.frame(
    env           = as.character(env_id),
    LAT           = as.numeric(lat_i),
    LON           = as.numeric(lon_i),
    YYYYMMDD      = base$YYYYMMDD,
    daysFromStart = as.numeric(base$YYYYMMDD - as.Date(start_i)) + 1,
    stringsAsFactors = FALSE)

  sci <- setdiff(names(base), c("YYYYMMDD", "n_pixel", "n_total", "gap"))
  for (v in sci) out[[v]] <- base[[v]]
  out$n_pixel <- base$n_pixel
  out$n_total <- base$n_total
  out$gap     <- base$gap
  rownames(out) <- NULL
  out
}


# =========================================================================
# SECTION 4 - rs_to_daily()
# =========================================================================

#' @title Expand Composite Satellite Series onto a Daily Grid
#'
#' @description
#' Converts an 8-day or 16-day composite series returned by
#' \code{\link{get_remote_sensing}} into one row per day, so that it can be
#' merged with \code{\link{get_weather}} output. Interpolation is explicit,
#' bounded by \code{max_gap}, and every synthesised row is flagged.
#'
#' @author Germano Costa Neto
#'
#' @param rs data.frame. Output of \code{\link{get_remote_sensing}}.
#' @param variables character. Columns to interpolate. Default: all numeric
#'   science columns.
#' @param method character. \code{"linear"} (default), \code{"constant"}
#'   (last observation carried forward), \code{"spline"}, or \code{"none"}
#'   (expand the grid, leave gaps as \code{NA}).
#' @param max_gap integer. Longest gap, in days, across which interpolation is
#'   permitted. Gaps longer than this are left \code{NA}. Default 16.
#' @param divide_by_interval boolean. If \code{TRUE}, divide flux variables by
#'   the composite length so an 8-day total becomes a daily rate. Applies only
#'   to collections whose unit ends in \code{"/8day"}.
#' @param period integer. Composite length in days. Inferred from the median
#'   spacing of the observations when \code{NULL}.
#' @param verbose boolean. Print progress.
#'
#' @return
#' A data.frame with one row per environment-day, carrying an
#' \code{interpolated} logical column and the original \code{gap} flag on
#' observed dates.
#'
#' @details
#' \strong{Interpolation is not observation.} An 8-day MODIS composite is a
#' single retrieval representing the best pixel in that window, not a daily
#' mean. Expanding it to daily creates seven rows per observation that
#' contain no new information. The \code{interpolated} column must be carried
#' into any downstream diagnostic: validating a model against interpolated
#' values partly validates the interpolation.
#'
#' \code{max_gap} exists because linear interpolation across a 40-day cloud
#' gap produces a smooth, plausible and entirely fictional curve. The default
#' of 16 days permits bridging one missing 8-day composite and refuses more.
#'
#' @examples
#' # synthetic 8-day series with a long gap
#' d <- data.frame(env = "E1",
#'                 YYYYMMDD = as.Date("2020-05-01") + c(0, 8, 16, 48, 56),
#'                 ET_500m = c(10, 14, 18, 22, 20),
#'                 gap = FALSE)
#' out <- rs_to_daily(d, method = "linear", max_gap = 16, verbose = FALSE)
#' # the 32-day gap is NOT bridged
#' sum(is.na(out$ET_500m))
#'
#' @seealso \code{\link{get_remote_sensing}}
#'
#' @export
rs_to_daily <- function(rs, variables = NULL, method = c("linear", "constant",
                                                         "spline", "none"),
                        max_gap = 16L, divide_by_interval = FALSE,
                        period = NULL, verbose = TRUE) {

  method <- match.arg(method)
  .et_banner("rs_to_daily", "composite -> daily grid, gap-bounded", verbose)

  if (!is.data.frame(rs) || !nrow(rs))
    stop("'rs' must be a non-empty data.frame.", call. = FALSE)
  need <- c("env", "YYYYMMDD")
  miss <- setdiff(need, names(rs))
  if (length(miss))
    stop("Column(s) missing from 'rs': ", paste(miss, collapse = ", "),
         call. = FALSE)
  rs$YYYYMMDD <- as.Date(rs$YYYYMMDD)

  meta <- c("env", "LAT", "LON", "YYYYMMDD", "daysFromStart",
            "n_pixel", "n_total", "gap", "interpolated")
  if (is.null(variables))
    variables <- setdiff(names(rs)[vapply(rs, is.numeric, logical(1))], meta)
  variables <- intersect(variables, names(rs))
  if (!length(variables))
    stop("No numeric science columns found to interpolate.", call. = FALSE)

  if (!is.numeric(max_gap) || length(max_gap) != 1L || max_gap < 1)
    stop("'max_gap' must be a single positive number of days.", call. = FALSE)

  out <- lapply(split(rs, rs$env), function(d) {
    d <- d[order(d$YYYYMMDD), , drop = FALSE]
    d <- d[!duplicated(d$YYYYMMDD), , drop = FALSE]

    per <- period
    if (is.null(per)) {
      dd <- as.numeric(diff(d$YYYYMMDD))
      per <- if (length(dd)) max(1L, round(stats::median(dd))) else 1L
    }

    grid <- seq(min(d$YYYYMMDD), max(d$YYYYMMDD), by = "day")
    res <- data.frame(env = d$env[1], YYYYMMDD = grid, stringsAsFactors = FALSE)
    if ("LAT" %in% names(d)) res$LAT <- d$LAT[1]
    if ("LON" %in% names(d)) res$LON <- d$LON[1]

    observed <- grid %in% d$YYYYMMDD

    # Distance in days to the nearest observation that is not NA for this
    # variable; computed per variable because gaps differ between bands.
    for (v in variables) {
      obs_ok <- !is.na(d[[v]])
      xs <- as.numeric(d$YYYYMMDD[obs_ok])
      ys <- d[[v]][obs_ok]

      yv <- rep(NA_real_, length(grid))
      if (length(xs) >= 1L) {
        xg <- as.numeric(grid)
        if (identical(method, "none") || length(xs) == 1L) {
          hit <- match(xg, xs)
          yv <- ifelse(is.na(hit), NA_real_, ys[hit])
        } else if (identical(method, "spline") && length(xs) >= 4L) {
          yv <- stats::spline(xs, ys, xout = xg)$y
        } else {
          rule_method <- if (identical(method, "constant")) "constant" else "linear"
          yv <- stats::approx(xs, ys, xout = xg, method = rule_method,
                              rule = 1)$y
        }

        # Enforce max_gap: blank any grid point whose bracketing observations
        # are further apart than max_gap.
        if (length(xs) >= 2L && !identical(method, "none")) {
          lo <- findInterval(xg, xs)
          span <- rep(NA_real_, length(xg))
          inside <- lo >= 1L & lo < length(xs)
          span[inside] <- xs[lo[inside] + 1L] - xs[lo[inside]]
          too_far <- !is.na(span) & span > max_gap
          exact   <- xg %in% xs
          yv[too_far & !exact] <- NA_real_
        }
      }

      if (isTRUE(divide_by_interval)) yv <- yv / per
      res[[v]] <- yv
    }

    res$observed <- observed
    res$interpolated <- !observed & !is.na(res[[variables[1]]])

    if ("gap" %in% names(d)) {
      g <- rep(NA, length(grid))
      g[observed] <- d$gap[match(grid[observed], d$YYYYMMDD)]
      res$gap <- g
    }
    if ("n_pixel" %in% names(d)) {
      np <- rep(NA_integer_, length(grid))
      np[observed] <- d$n_pixel[match(grid[observed], d$YYYYMMDD)]
      res$n_pixel <- np
    }
    if ("daysFromStart" %in% names(d)) {
      origin <- d$YYYYMMDD[1] - (d$daysFromStart[1] - 1)
      res$daysFromStart <- as.numeric(grid - origin) + 1
    }
    res
  })

  out <- .rbind_fill(out)
  out <- out[order(out$env, out$YYYYMMDD), , drop = FALSE]
  rownames(out) <- NULL

  if (verbose) {
    fi <- mean(out$interpolated, na.rm = TRUE)
    nn <- sum(is.na(out[[variables[1]]]))
    .et_step(sprintf("rows ................ %d", nrow(out)), verbose)
    .et_step(sprintf("interpolated ........ %.1f%%", 100 * fi), verbose)
    .et_step(sprintf("left NA (> max_gap) . %d", nn), verbose)
  }
  out
}


require(modisfast)


Sys.setenv(EARTHDATA_TOKEN = "eyJ0eXAiOiJKV1Qi...")   # urs.earthdata.nasa.gov

data(G2F_2014_17)
d <- G2F_2014_17[1:5, ]



lst <- get_remote_sensing(
  env.id     = d$env,
  lat        = d$lat,
  lon        = d$lon,
  start.day  = d$start,
  end.day    = d$end,
  collection = "MOD11A1.061",
  variables  = c("LST_Day_1km", "LST_Night_1km"),
  buffer_km  = 2,
  aggregate  = mean,
  qa_mask    = TRUE,
  strict_qa  = TRUE
)



rain <- get_remote_sensing(
  env.id = d$env, lat = d$lat, lon = d$lon,
  start.day = d$start, end.day = d$end,
  collection = "GPM_3IMERGDF.07",
  variables  = "precipitationCal",
  buffer_km  = 15        # IMERG is 0.1 deg (~10 km) — buffer MUST exceed pixel size
)

stopifnot(!any(rain$gap))   # should hold for IMERG
