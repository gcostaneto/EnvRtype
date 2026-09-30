#==================================================================================================
# env_dissection.R
#
# Dimensionality reduction of interval-resolved environmental data, following
# Della Coletta et al. (2023) GENETICS 224(4):iyad103,
# "Linking genetic and environmental factors through marker effect networks to
#  understand trait plasticity".
#
# Module contents
#   env_indices()        - build the environment x (factor x interval) index matrix
#   env_pca()            - PCA of the index matrix (prcomp, as in the paper)
#   env_pca_biplot()     - environments in PC space (paper Fig. 1c)
#   env_pca_scree()      - variance explained, for choosing how many PCs to keep
#   env_cor_heatmap()    - correlation among indices (paper Fig. 1b)
#   env_loading_curve()  - loadings through the season, loess-smoothed (paper Fig. 6c)
#   env_pc_associate()   - correlate PCs with an external per-environment variable
#
# --------------------------------------------------------------------------
# WHAT THE PAPER DID (Methods, "Obtaining weather factors and indices")
# --------------------------------------------------------------------------
# 17 weather factors were extracted with EnvRtype across 9 environments in
# 3-day intervals from planting to 151 days after planting. Fixing the last
# interval at 151 DAP gives every environment the same number of intervals:
#   51 intervals x 17 factors = 867 environmental indices.
# Pearson correlation among all indices was computed with cor(); PCA was
# performed with prcomp() to reduce the dimensionality of this highly
# correlated set. PCs 1-9 were retained on the basis of a scree plot. The
# loadings of an associated PC were then plotted through the growing season
# with a loess curve, one line per weather factor (Fig. 6c), to show which
# factors at which part of the season drive that PC.
#
# --------------------------------------------------------------------------
# DESIGN NOTES / DEVIATIONS
# --------------------------------------------------------------------------
# * The paper does not state whether prcomp() was called with scale. = TRUE.
#   The 17 factors are on wildly different scales (radiation in MJ/m2/day vs
#   precipitation in mm vs a unitless FRUE in [0,1]), so an unscaled PCA is
#   dominated by whichever factor has the largest variance. This module
#   therefore defaults to scale. = TRUE and exposes the choice. Set
#   scale. = FALSE to reproduce a covariance-matrix PCA.
#
# * Zero-variance indices are dropped before scaling, not silently passed to
#   prcomp() where they produce NaN. Which ones were dropped is returned.
#
# * The number of PCs to keep is NOT hardcoded to 9. That was a scree-plot
#   judgement on one dataset; env_pca_scree() is provided so the same
#   judgement can be made on yours.
#
# --------------------------------------------------------------------------
# DEPENDENCIES
# --------------------------------------------------------------------------
# Base R (stats, graphics, grDevices) only. ggplot2 is used for plots when
# installed; every plotting function falls back to base graphics otherwise.
#==================================================================================================


# =========================================================================
# SECTION 0 - internal helpers
# =========================================================================

#' Null-coalescing helper
#' @keywords internal
#' @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a

#' Coerce to a numeric matrix with informative errors
#'
#' @param x matrix or data.frame.
#' @param arg character. Argument name used in messages.
#' @return A numeric matrix.
#' @keywords internal
#' @noRd
.rel_as_matrix <- function(x, arg = "x") {
  if (is.null(x)) stop("'", arg, "' is missing.", call. = FALSE)
  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.matrix(x)) x <- as.matrix(x)
  if (!is.numeric(x)) {
    storage.mode(x) <- "numeric"
    if (!is.numeric(x))
      stop("'", arg, "' must be numeric.", call. = FALSE)
  }
  x
}

#' Split "<factor>_<interval>" index names into their two parts
#'
#' Index names are built by \code{env_indices()} as \code{factor_sep_interval},
#' e.g. \code{"T2M_MAX_i017"}. Splitting on the LAST separator is deliberate:
#' factor names themselves contain underscores (T2M_MAX, T2M_MIN), so splitting
#' on the first would truncate them.
#'
#' @param nm character vector of index names.
#' @param sep character. Separator used by \code{env_indices()}.
#' @return A data.frame with \code{index}, \code{factor} and \code{interval}.
#' @keywords internal
#' @noRd
.rel_split_names <- function(nm, sep = "_") {
  pos <- regexpr(paste0(sep, "[^", sep, "]*$"), nm)
  ok  <- pos > 0
  fac <- ifelse(ok, substr(nm, 1, pos - 1), nm)
  itv <- ifelse(ok, substr(nm, pos + nchar(sep), nchar(nm)), NA_character_)
  data.frame(index = nm, factor = fac, interval = itv,
             stringsAsFactors = FALSE)
}

#' Map an interval label back to a numeric mid-point day
#'
#' @param itv character vector of interval labels, e.g. \code{"i017"}.
#' @param map named numeric vector from \code{env_indices()}
#'   (\code{"interval.mid"} attribute); \code{NULL} falls back to the
#'   integer embedded in the label.
#' @return Numeric vector of days.
#' @keywords internal
#' @noRd
.rel_interval_day <- function(itv, map = NULL) {
  if (!is.null(map) && all(itv %in% names(map))) return(unname(map[itv]))
  suppressWarnings(as.numeric(gsub("[^0-9.]", "", itv)))
}

#' Is ggplot2 usable?
#' @keywords internal
#' @noRd
.rel_has_ggplot <- function() requireNamespace("ggplot2", quietly = TRUE)


# =========================================================================
# SECTION 1 - env_indices()
# =========================================================================

#' @title Build Interval-Resolved Environmental Indices
#'
#' @description
#' Converts daily environmental data into the environment x index matrix used by
#' Della Coletta et al. (2023): each weather factor is summarised within fixed-width
#' time intervals from planting, producing one column per
#' \emph{factor x interval} combination.
#'
#' @author Implementation following Della Coletta et al. (2023)
#'
#' @param env.data data.frame of daily weather, e.g. a \code{\link{get_weather}}
#'   or \code{\link{processWTH}} output.
#' @param env.id character. Column identifying the environment. Default \code{"env"}.
#' @param day.id character. Column giving days from planting. Default \code{"daysFromStart"}.
#' @param var.id character vector. Which weather factors to use. \code{NULL} (default)
#'   uses every numeric column that is not an id column.
#' @param interval integer. Width of each time window in days. Default 3, as in the paper.
#' @param end.day integer. Last day retained. Default 151, as in the paper. Days beyond
#'   this are discarded so that every environment contributes the same intervals.
#' @param statistic character. Within-interval summary: \code{"mean"} (default),
#'   \code{"sum"}, \code{"min"}, \code{"max"} or \code{"median"}. Precipitation is
#'   arguably better summed than averaged -- see \code{statistic.by}.
#' @param statistic.by named character vector. Per-factor override of \code{statistic},
#'   e.g. \code{c(PRECTOT = "sum")}. Names are factor names, values are statistics.
#' @param id.names character vector. Additional columns to treat as identifiers
#'   (never summarised). Default \code{NULL}.
#' @param sep character. Separator between factor and interval in the output column
#'   names. Default \code{"_"}.
#' @param require.complete boolean. If \code{TRUE} (default) an error is raised when
#'   any environment x interval x factor cell is missing, because an incomplete matrix
#'   silently biases the PCA. Set \code{FALSE} to allow \code{NA}s through.
#' @param verbose boolean. Print a summary. Default \code{TRUE}.
#'
#' @return
#' A numeric matrix with environments in rows and \eqn{n_{factor} \times n_{interval}}
#' indices in columns. Attributes:
#' \describe{
#'   \item{\code{"factors"}}{the factor names used}
#'   \item{\code{"intervals"}}{the interval labels}
#'   \item{\code{"interval.mid"}}{named numeric vector of interval mid-point days}
#'   \item{\code{"interval.range"}}{data.frame of interval start/end days}
#'   \item{\code{"n.missing"}}{number of missing cells}
#' }
#'
#' @details
#' \strong{Why a fixed end day.} Environments differ in season length. Without a common
#' cut-off, environments contribute different numbers of intervals and the matrix is
#' ragged. The paper fixed the last interval at 151 days after planting; the same device
#' is used here via \code{end.day}.
#'
#' \strong{Interval labelling.} With \code{interval = 3} and \code{end.day = 151}, days
#' 1-3 form interval 1, days 4-6 interval 2, and so on; the final partial window
#' (day 151) forms interval 51. This reproduces the 51 intervals of the paper, and with
#' 17 factors gives the 867 indices reported.
#'
#' \strong{Mean vs sum.} The paper does not state the within-interval statistic.
#' A mean is used by default because it is defined for every factor, but for
#' accumulating quantities such as precipitation a sum is more interpretable;
#' \code{statistic.by} allows a per-factor choice.
#'
#' @examples
#' \donttest{
#' data("maizeWTH")
#'
#' ## 3-day intervals across the shared season. The bundled sample does not
#' ## reach the paper's 151 days in every environment, so end.day is trimmed to
#' ## a fully observed window (the default 151 would leave gaps and error out).
#' W <- env_indices(maizeWTH, env.id = "env", day.id = "daysFromStart",
#'                  end.day = 130)
#' dim(W)
#' attr(W, "intervals")[1:5]
#'
#' ## Sum precipitation within each window, average everything else
#' W2 <- env_indices(maizeWTH, statistic = "mean",
#'                   statistic.by = c(PRECTOT = "sum"), end.day = 130)
#'
#' ## Weekly windows instead
#' W3 <- env_indices(maizeWTH, interval = 7, end.day = 140)
#' }
#'
#' @seealso \code{\link{env_pca}}, \code{\link{env_cor_heatmap}}
#'
#' @references
#' Della Coletta R., Liese S.E., Fernandes S.B., Mikel M.A., Bohn M.O., Lipka A.E.,
#' Hirsch C.N. (2023). Linking genetic and environmental factors through marker effect
#' networks to understand trait plasticity. \emph{GENETICS} 224(4), iyad103.
#'
#' @importFrom stats aggregate median
#' @export
env_indices <- function(env.data, env.id = "env", day.id = "daysFromStart",
                        var.id = NULL, interval = 3L, end.day = 151L,
                        statistic = c("mean", "sum", "min", "max", "median"),
                        statistic.by = NULL, id.names = NULL, sep = "_",
                        require.complete = TRUE, verbose = TRUE) {

  .et_banner("env_indices", "builds interval-resolved environmental indices", verbose)
  statistic <- match.arg(statistic)
  env.data  <- as.data.frame(env.data)

  if (!env.id %in% names(env.data))
    stop("Environment column '", env.id, "' not found in env.data.", call. = FALSE)
  if (!day.id %in% names(env.data))
    stop("Day column '", day.id, "' not found in env.data.", call. = FALSE)
  if (!is.numeric(interval) || interval < 1)
    stop("'interval' must be a positive number of days.", call. = FALSE)

  # ---- choose the weather factors ----------------------------------------
  drop.cols <- unique(c(env.id, day.id, id.names,
                        "LON", "LAT", "YEAR", "MM", "DD", "DOY", "YYYYMMDD",
                        "daysFromStart"))
  if (is.null(var.id)) {
    num <- vapply(env.data, is.numeric, logical(1))
    var.id <- setdiff(names(env.data)[num], drop.cols)
  } else {
    miss <- setdiff(var.id, names(env.data))
    if (length(miss))
      stop("Weather factors not found in env.data: ",
           paste(miss, collapse = ", "), call. = FALSE)
  }
  if (!length(var.id))
    stop("No weather factors to summarise. Supply 'var.id' explicitly.", call. = FALSE)

  # ---- restrict to the common season window ------------------------------
  d <- suppressWarnings(as.numeric(env.data[[day.id]]))
  if (all(is.na(d)))
    stop("Column '", day.id, "' contains no usable numeric days.", call. = FALSE)

  keep <- !is.na(d) & d >= 1 & d <= end.day
  n.drop <- sum(!keep)
  env.data <- env.data[keep, , drop = FALSE]
  d <- d[keep]
  if (!nrow(env.data))
    stop("No records fall within days 1..", end.day, ".", call. = FALSE)

  # ---- assign each day to an interval ------------------------------------
  # Day 1 must land in interval 1, so index from zero before dividing.
  idx <- as.integer(floor((d - 1) / interval) + 1L)
  n.int <- as.integer(ceiling(end.day / interval))
  # Pad to at least 3 digits so interval labels (and therefore column names)
  # have the same form whether a run has 51 intervals or 151. Data-dependent
  # width would silently rename columns when end.day or interval changes.
  width <- max(3L, nchar(as.character(n.int)))
  lab <- sprintf(paste0("i%0", width, "d"), seq_len(n.int))

  start.day <- (seq_len(n.int) - 1L) * interval + 1L
  stop.day  <- pmin(seq_len(n.int) * interval, end.day)
  mid.day   <- (start.day + stop.day) / 2
  names(mid.day) <- lab

  env.data$.interval <- factor(lab[idx], levels = lab)
  env.data$.env      <- as.character(env.data[[env.id]])

  envs <- sort(unique(env.data$.env))

  # ---- summarise each factor within each environment x interval ----------
  fun.for <- function(v) {
    s <- if (!is.null(statistic.by) && v %in% names(statistic.by))
      statistic.by[[v]] else statistic
    switch(s,
           mean   = function(z) if (all(is.na(z))) NA_real_ else mean(z, na.rm = TRUE),
           # sum(NA, na.rm = TRUE) is 0, which reads as a real zero total
           # (e.g. "no rain") when nothing was actually measured.
           sum    = function(z) if (all(is.na(z))) NA_real_ else sum(z, na.rm = TRUE),
           min    = function(z) if (all(is.na(z))) NA_real_ else min(z, na.rm = TRUE),
           max    = function(z) if (all(is.na(z))) NA_real_ else max(z, na.rm = TRUE),
           median = function(z) if (all(is.na(z))) NA_real_ else stats::median(z, na.rm = TRUE),
           stop("Unknown statistic '", s, "' for factor '", v, "'.", call. = FALSE))
  }

  out <- matrix(NA_real_, nrow = length(envs), ncol = length(var.id) * n.int,
                dimnames = list(envs,
                                paste0(rep(var.id, each = n.int), sep,
                                       rep(lab, times = length(var.id)))))

  ei <- match(env.data$.env, envs)
  ii <- as.integer(env.data$.interval)

  for (k in seq_along(var.id)) {
    v  <- var.id[k]
    fk <- fun.for(v)
    vals <- suppressWarnings(as.numeric(env.data[[v]]))
    agg  <- tapply(vals, list(ei, ii), fk)
    # tapply drops empty combinations; place values by name, not position.
    ri <- as.integer(rownames(agg))
    ci <- as.integer(colnames(agg))
    block <- (k - 1L) * n.int
    out[ri, block + ci] <- agg
  }

  n.missing <- sum(is.na(out))
  if (n.missing && isTRUE(require.complete)) {
    bad <- which(is.na(out), arr.ind = TRUE)
    eg  <- utils::head(paste0(rownames(out)[bad[, 1]], " / ",
                              colnames(out)[bad[, 2]]), 5)
    stop(n.missing, " of ", length(out), " environment x index cells are missing.\n",
         "  e.g. ", paste(eg, collapse = "; "), "\n",
         "  An incomplete matrix biases the PCA. Shorten 'end.day', widen ",
         "'interval',\n  drop sparse factors, or set require.complete = FALSE ",
         "to proceed anyway.", call. = FALSE)
  }

  attr(out, "factors")        <- var.id
  attr(out, "intervals")      <- lab
  attr(out, "interval.mid")   <- mid.day
  attr(out, "interval.range") <- data.frame(interval = lab, start = start.day,
                                            end = stop.day, mid = mid.day,
                                            row.names = NULL)
  attr(out, "n.missing")      <- n.missing
  attr(out, "sep")            <- sep

  if (isTRUE(verbose)) {
    message(strrep("-", 60))
    message("env_indices(): ", length(envs), " environments x ",
            length(var.id), " factors x ", n.int, " intervals = ",
            ncol(out), " indices")
    message("  interval width : ", interval, " day(s), days 1..", end.day)
    if (n.drop)
      message("  discarded      : ", n.drop, " record(s) outside 1..", end.day)
    if (n.missing)
      message("  missing cells  : ", n.missing)
    message(strrep("-", 60))
  }

  out
}


# =========================================================================
# SECTION 2 - env_pca()
# =========================================================================

#' @title PCA of Interval-Resolved Environmental Indices
#'
#' @description
#' Runs \code{\link[stats]{prcomp}} on the environment x index matrix, as in
#' Della Coletta et al. (2023), and returns scores, loadings and variance
#' explained in a form the plotting functions in this module can consume.
#'
#' @param W matrix. Environment x index matrix, typically from \code{\link{env_indices}}.
#' @param center boolean. Centre the indices. Default \code{TRUE}.
#' @param scale. boolean. Scale the indices to unit variance. Default \code{TRUE};
#'   see Details for why this differs from a bare \code{prcomp} call.
#' @param rank. integer or NULL. Maximum number of PCs to compute. \code{NULL} (default)
#'   computes all available, which is \code{min(nrow - 1, ncol)}.
#' @param drop.constant boolean. If \code{TRUE} (default) indices with zero variance
#'   are removed before the PCA rather than producing \code{NaN} under scaling.
#' @param verbose boolean. Print a summary. Default \code{TRUE}.
#'
#' @return
#' A list of class \code{"env_pca"}:
#' \describe{
#'   \item{\code{scores}}{environment x PC matrix (\code{prcomp$x})}
#'   \item{\code{loadings}}{index x PC matrix (\code{prcomp$rotation})}
#'   \item{\code{variance}}{data.frame with \code{PC}, \code{sd}, \code{prop}, \code{cum}}
#'   \item{\code{prcomp}}{the underlying \code{prcomp} object}
#'   \item{\code{dropped}}{indices removed as constant}
#'   \item{\code{meta}}{data.frame mapping each index to its factor, interval and day}
#' }
#'
#' @details
#' \strong{Scaling.} The paper reports only that \code{prcomp} was used. Because the
#' 17 factors are on incomparable scales -- radiation in MJ/m2/day, precipitation in mm,
#' FRUE unitless in [0,1] -- an unscaled PCA is dominated by the highest-variance
#' factor and the leading PCs largely describe that one variable. \code{scale. = TRUE}
#' is therefore the default here. Set \code{scale. = FALSE} for a covariance-matrix PCA.
#'
#' \strong{Rank.} With \eqn{q} environments and \eqn{k \gg q} indices, at most
#' \eqn{q - 1} PCs carry variance. With the 9 environments of the paper that is 8, so
#' the reported "ninth PC" explaining \eqn{<0.01\%} is effectively numerical noise.
#' Expect very few informative PCs unless you have many environments.
#'
#' @examples
#' \donttest{
#' data("maizeWTH")
#' W  <- env_indices(maizeWTH, end.day = 130)
#' pc <- env_pca(W)
#' head(pc$variance)
#' env_pca_scree(pc)
#' env_pca_biplot(pc)
#' }
#'
#' @seealso \code{\link{env_indices}}, \code{\link{env_pca_biplot}},
#'   \code{\link{env_loading_curve}}
#'
#' @importFrom stats prcomp sd
#' @export
env_pca <- function(W, center = TRUE, scale. = TRUE, rank. = NULL,
                    drop.constant = TRUE, verbose = TRUE) {

  .et_banner("env_pca", "runs PCA on the environmental index matrix", verbose)
  W <- .rel_as_matrix(W, "W")
  if (nrow(W) < 3)
    stop("At least three environments are required for a PCA (got ", nrow(W), ").",
         call. = FALSE)
  if (anyNA(W))
    stop("W contains missing values. Re-run env_indices() with a shorter 'end.day' ",
         "or wider 'interval', or impute before calling env_pca().", call. = FALSE)

  sep <- attr(W, "sep") %||% "_"
  mid <- attr(W, "interval.mid")

  dropped <- character(0)
  if (isTRUE(drop.constant)) {
    sdv <- apply(W, 2, stats::sd)
    bad <- !is.finite(sdv) | sdv < .Machine$double.eps
    if (any(bad)) {
      dropped <- colnames(W)[bad]
      W <- W[, !bad, drop = FALSE]
      if (!ncol(W))
        stop("Every index is constant across environments; nothing to decompose.",
             call. = FALSE)
    }
  }

  max.rank <- min(nrow(W) - 1L, ncol(W))
  if (is.null(rank.)) rank. <- max.rank
  rank. <- min(rank., max.rank)

  pc <- stats::prcomp(W, center = center, scale. = scale., rank. = rank.)

  sdev <- pc$sdev
  prop <- sdev^2 / sum(sdev^2)
  keep <- seq_len(min(length(sdev), rank.))
  variance <- data.frame(PC = paste0("PC", keep),
                         sd = sdev[keep],
                         prop = prop[keep],
                         cum = cumsum(prop)[keep],
                         row.names = NULL, stringsAsFactors = FALSE)

  meta <- .rel_split_names(rownames(pc$rotation), sep = sep)
  meta$day <- .rel_interval_day(meta$interval, mid)

  out <- list(scores = pc$x, loadings = pc$rotation, variance = variance,
              prcomp = pc, dropped = dropped, meta = meta)
  class(out) <- "env_pca"

  if (isTRUE(verbose)) {
    message(strrep("-", 60))
    message("env_pca(): ", nrow(W), " environments x ", ncol(W), " indices")
    if (length(dropped))
      message("  dropped ", length(dropped), " constant index/indices")
    message("  informative PCs: ", max.rank, " (= n.env - 1)")
    message("  PC1 ", sprintf("%.1f%%", 100 * variance$prop[1]),
            if (nrow(variance) > 1)
              paste0(", PC2 ", sprintf("%.1f%%", 100 * variance$prop[2])) else "",
            if (nrow(variance) > 1)
              paste0(", cumulative ", sprintf("%.1f%%", 100 * variance$cum[2])) else "")
    message(strrep("-", 60))
  }

  out
}

#' @export
print.env_pca <- function(x, ...) {
  cat("env_pca:", nrow(x$scores), "environments,",
      nrow(x$loadings), "indices\n")
  cat("Variance explained (first", min(6, nrow(x$variance)), "PCs):\n")
  v <- utils::head(x$variance, 6)
  v$prop <- sprintf("%.2f%%", 100 * v$prop)
  v$cum  <- sprintf("%.2f%%", 100 * v$cum)
  print(v, row.names = FALSE)
  invisible(x)
}


# =========================================================================
# SECTION 3 - scree and biplot
# =========================================================================

#' @title Scree Plot of Environmental PCA
#'
#' @description
#' Plots the proportion of variance explained per PC, the visualisation the paper
#' used to decide how many PCs to retain.
#'
#' @param x an \code{\link{env_pca}} object.
#' @param n.pc integer. How many PCs to show. Default \code{min(15, available)}.
#' @param cumulative boolean. Overlay the cumulative proportion. Default \code{TRUE}.
#' @param title character. Plot title.
#' @param verbose boolean. If \code{TRUE} (default) prints a progress banner.
#'
#' @return A ggplot object when \pkg{ggplot2} is available (invisibly), otherwise
#'   \code{NULL} after drawing a base plot.
#'
#' @examples
#' \donttest{
#' pc <- env_pca(env_indices(maizeWTH, end.day = 130))
#' env_pca_scree(pc)
#' }
#'
#' @seealso \code{\link{env_pca}}
#' @importFrom graphics barplot lines points axis legend par text
#' @export
env_pca_scree <- function(x, n.pc = NULL, cumulative = TRUE,
                          title = "Variance explained by environmental PCs",
                          verbose = TRUE) {

  .et_banner("env_pca_scree", "plots variance explained per PC", verbose)
  if (!inherits(x, "env_pca")) stop("'x' must be an env_pca object.", call. = FALSE)
  v <- x$variance
  if (is.null(n.pc)) n.pc <- min(15L, nrow(v))
  v <- v[seq_len(min(n.pc, nrow(v))), , drop = FALSE]
  v$PC <- factor(v$PC, levels = v$PC)

  if (.rel_has_ggplot()) {
    p <- ggplot2::ggplot(v, ggplot2::aes(x = .data[["PC"]], y = .data[["prop"]])) +
      ggplot2::geom_col(fill = "grey35") +
      ggplot2::scale_y_continuous(labels = function(z) paste0(round(100 * z, 1), "%")) +
      ggplot2::labs(x = NULL, y = "Variance explained", title = title) +
      ggplot2::theme_bw()
    if (isTRUE(cumulative))
      p <- p +
        ggplot2::geom_line(ggplot2::aes(y = .data[["cum"]], group = 1),
                           colour = "firebrick") +
        ggplot2::geom_point(ggplot2::aes(y = .data[["cum"]]), colour = "firebrick")
    print(p)
    return(invisible(p))
  }

  op <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(op), add = TRUE)
  bp <- graphics::barplot(v$prop, names.arg = as.character(v$PC), las = 2,
                          ylim = c(0, 1), col = "grey35",
                          ylab = "Variance explained", main = title)
  if (isTRUE(cumulative)) {
    graphics::lines(bp, v$cum, col = "firebrick", lwd = 2)
    graphics::points(bp, v$cum, col = "firebrick", pch = 16)
  }
  invisible(NULL)
}

#' @title Environments in Principal Component Space
#'
#' @description
#' Scatter of environments on two PCs, reproducing Fig. 1c of Della Coletta et al.
#' (2023), where each point is one growth environment labelled by its code.
#'
#' @param x an \code{\link{env_pca}} object.
#' @param pc integer vector of length 2. Which PCs to plot. Default \code{c(1, 2)}.
#' @param label boolean. Draw environment names. Default \code{TRUE}.
#' @param group factor or named vector. Optional grouping used to colour points,
#'   e.g. a mega-environment assignment. Names are matched to environment names.
#' @param title character. Plot title.
#' @param verbose boolean. If \code{TRUE} (default) prints a progress banner.
#'
#' @return A ggplot object when \pkg{ggplot2} is available (invisibly), otherwise
#'   \code{NULL} after drawing a base plot.
#'
#' @examples
#' \donttest{
#' pc <- env_pca(env_indices(maizeWTH, end.day = 130))
#' env_pca_biplot(pc)
#' env_pca_biplot(pc, pc = c(1, 3))
#' }
#'
#' @seealso \code{\link{env_pca}}, \code{\link{env_pca_scree}}
#' @importFrom graphics plot text abline
#' @export
env_pca_biplot <- function(x, pc = c(1, 2), label = TRUE, group = NULL,
                           title = "Environments in PC space", verbose = TRUE) {

  .et_banner("env_pca_biplot", "plots environments in PC space", verbose)
  if (!inherits(x, "env_pca")) stop("'x' must be an env_pca object.", call. = FALSE)
  if (length(pc) != 2) stop("'pc' must be two PC numbers, e.g. c(1, 2).", call. = FALSE)
  if (max(pc) > ncol(x$scores))
    stop("Only ", ncol(x$scores), " PCs available; cannot plot PC", max(pc), ".",
         call. = FALSE)

  s <- as.data.frame(x$scores[, pc, drop = FALSE])
  names(s) <- c("x", "y")
  s$env <- rownames(x$scores)

  v <- x$variance
  labs <- sprintf("PC%d (%.1f%%)", pc, 100 * v$prop[pc])

  s$grp <- if (is.null(group)) "all" else {
    g <- if (!is.null(names(group))) group[s$env] else group
    as.character(g)
  }

  if (.rel_has_ggplot()) {
    p <- ggplot2::ggplot(s, ggplot2::aes(x = .data[["x"]], y = .data[["y"]])) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey85") +
      ggplot2::geom_vline(xintercept = 0, colour = "grey85") +
      ggplot2::labs(x = labs[1], y = labs[2], title = title) +
      ggplot2::theme_bw()
    p <- if (is.null(group))
      p + ggplot2::geom_point(size = 3, colour = "steelblue4") else
      p + ggplot2::geom_point(ggplot2::aes(colour = .data[["grp"]]), size = 3) +
          ggplot2::labs(colour = NULL)
    if (isTRUE(label)) {
      p <- if (requireNamespace("ggrepel", quietly = TRUE))
        p + ggrepel::geom_text_repel(ggplot2::aes(label = .data[["env"]]), size = 3) else
        p + ggplot2::geom_text(ggplot2::aes(label = .data[["env"]]),
                               vjust = -0.8, size = 3)
    }
    print(p)
    return(invisible(p))
  }

  op <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(op), add = TRUE)
  cols <- if (is.null(group)) "steelblue4" else
    grDevices::hcl.colors(length(unique(s$grp)), "Dark 3")[factor(s$grp)]
  graphics::plot(s$x, s$y, pch = 16, cex = 1.4, col = cols,
                 xlab = labs[1], ylab = labs[2], main = title)
  graphics::abline(h = 0, v = 0, col = "grey85")
  if (isTRUE(label))
    graphics::text(s$x, s$y, labels = s$env, pos = 3, cex = 0.8)
  invisible(NULL)
}


# =========================================================================
# SECTION 4 - correlation heatmap (paper Fig. 1b)
# =========================================================================

#' @title Correlation Heatmap of Environmental Indices
#'
#' @description
#' Pearson correlation among all environmental indices, as in Fig. 1b of
#' Della Coletta et al. (2023), which motivated the dimensionality reduction:
#' the 867 indices were highly correlated.
#'
#' @param W matrix. Environment x index matrix from \code{\link{env_indices}}.
#' @param order character. Column ordering: \code{"factor"} (default, groups all
#'   intervals of a factor together, giving the block structure of Fig. 1b),
#'   \code{"interval"} (chronological across factors) or \code{"none"}.
#' @param max.index integer. Guard against enormous plots; if the matrix has more
#'   columns than this, a regularly spaced subset is shown. Default 400.
#' @param title character. Plot title.
#' @param verbose boolean. If \code{TRUE} (default) prints a progress banner.
#'
#' @return Invisibly, the correlation matrix that was plotted.
#'
#' @details
#' The mean absolute off-diagonal correlation is printed as a one-number summary of
#' how redundant the index set is -- the quantity that justifies running a PCA at all.
#'
#' @examples
#' \donttest{
#' W <- env_indices(maizeWTH, end.day = 130)
#' R <- env_cor_heatmap(W)
#' }
#'
#' @seealso \code{\link{env_indices}}, \code{\link{env_pca}}
#' @importFrom stats cor
#' @importFrom graphics image axis box
#' @importFrom grDevices hcl.colors
#' @export
env_cor_heatmap <- function(W, order = c("factor", "interval", "none"),
                            max.index = 400L,
                            title = "Correlation among environmental indices",
                            verbose = TRUE) {

  .et_banner("env_cor_heatmap", "plots correlation among environmental indices", verbose)
  order <- match.arg(order)
  W <- .rel_as_matrix(W, "W")
  sep <- attr(W, "sep") %||% "_"
  mid <- attr(W, "interval.mid")

  meta <- .rel_split_names(colnames(W), sep = sep)
  meta$day <- .rel_interval_day(meta$interval, mid)

  ord <- switch(order,
                factor   = order(meta$factor, meta$day),
                interval = order(meta$day, meta$factor),
                none     = seq_len(ncol(W)))
  W <- W[, ord, drop = FALSE]
  meta <- meta[ord, , drop = FALSE]

  if (ncol(W) > max.index) {
    sel <- round(seq(1, ncol(W), length.out = max.index))
    W <- W[, sel, drop = FALSE]
    meta <- meta[sel, , drop = FALSE]
    message("Showing ", max.index, " of ", length(ord),
            " indices (raise 'max.index' to see more).")
  }

  R <- suppressWarnings(stats::cor(W, use = "pairwise.complete.obs"))
  offdiag <- R[upper.tri(R)]
  message("Mean |r| among indices: ",
          sprintf("%.3f", mean(abs(offdiag), na.rm = TRUE)),
          "  (high values justify the PCA)")

  op <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(op), add = TRUE)
  graphics::par(mar = c(5, 5, 4, 2))
  graphics::image(seq_len(ncol(R)), seq_len(nrow(R)), R,
                  zlim = c(-1, 1), col = grDevices::hcl.colors(101, "Blue-Red 3"),
                  xlab = "", ylab = "", axes = FALSE, main = title)
  # One tick per factor block, placed at the block centre.
  bl <- tapply(seq_len(nrow(meta)), meta$factor, mean)
  graphics::axis(1, at = bl, labels = names(bl), las = 2, cex.axis = 0.7, tick = FALSE)
  graphics::axis(2, at = bl, labels = names(bl), las = 2, cex.axis = 0.7, tick = FALSE)
  graphics::box()
  invisible(R)
}


# =========================================================================
# SECTION 5 - loading curves through the season (paper Fig. 6c)
# =========================================================================

#' @title Seasonal Loading Curves of an Environmental PC
#'
#' @description
#' Reproduces Fig. 6c of Della Coletta et al. (2023): for a chosen PC, the loading of
#' every index is plotted against its position in the growing season, with one loess
#' curve per weather factor. This answers "which weather factors, at which point in
#' the season, drive this PC".
#'
#' @param x an \code{\link{env_pca}} object.
#' @param pc integer. Which PC to display. Default 1.
#' @param factors character vector. Restrict to these weather factors. \code{NULL}
#'   (default) shows all.
#' @param span numeric. Loess span. Default 0.4; larger is smoother.
#' @param points boolean. Draw the raw per-interval loadings under the curves.
#'   Default \code{TRUE}.
#' @param facet boolean. One panel per factor instead of overlaying. Default \code{FALSE}.
#' @param title character. Plot title. \code{NULL} builds one automatically.
#' @param verbose boolean. If \code{TRUE} (default) prints a progress banner.
#'
#' @return
#' Invisibly, a data.frame of the plotted loadings (\code{index}, \code{factor},
#' \code{day}, \code{loading}), so the curves can be re-plotted or tabulated.
#' When \pkg{ggplot2} is available the ggplot object is attached as the
#' \code{"plot"} attribute.
#'
#' @details
#' \strong{Interpretation.} A factor whose curve departs strongly from zero over a
#' stretch of the season contributes heavily to that PC during that window. Sign is
#' arbitrary up to the usual PCA reflection: only relative signs within a PC are
#' meaningful.
#'
#' \strong{Loess caveat.} The smoother is cosmetic. With 51 intervals per factor the
#' curve is well determined, but with few intervals it can imply structure that the
#' underlying points do not support -- keep \code{points = TRUE} to see the evidence.
#'
#' @examples
#' \donttest{
#' pc <- env_pca(env_indices(maizeWTH, end.day = 130))
#'
#' ## All factors on one panel
#' env_loading_curve(pc, pc = 1)
#'
#' ## Temperature factors only, one panel each
#' env_loading_curve(pc, pc = 2,
#'                   factors = c("T2M", "T2M_MAX", "T2M_MIN"), facet = TRUE)
#'
#' ## Recover the plotted values
#' L <- env_loading_curve(pc, pc = 1)
#' head(L[order(-abs(L$loading)), ])
#' }
#'
#' @seealso \code{\link{env_pca}}, \code{\link{env_pc_associate}}
#'
#' @importFrom stats loess predict
#' @importFrom graphics plot lines points legend abline par
#' @importFrom grDevices hcl.colors
#' @export
env_loading_curve <- function(x, pc = 1, factors = NULL, span = 0.4,
                              points = TRUE, facet = FALSE, title = NULL,
                              verbose = TRUE) {

  .et_banner("env_loading_curve", "plots seasonal PC loading curves", verbose)
  if (!inherits(x, "env_pca")) stop("'x' must be an env_pca object.", call. = FALSE)
  if (length(pc) != 1) stop("'pc' must be a single PC number.", call. = FALSE)
  if (pc > ncol(x$loadings))
    stop("Only ", ncol(x$loadings), " PCs available; cannot plot PC", pc, ".",
         call. = FALSE)

  L <- data.frame(index   = rownames(x$loadings),
                  factor  = x$meta$factor,
                  day     = x$meta$day,
                  loading = as.numeric(x$loadings[, pc]),
                  stringsAsFactors = FALSE)

  if (!is.null(factors)) {
    miss <- setdiff(factors, unique(L$factor))
    if (length(miss))
      stop("Factors not present in the PCA: ", paste(miss, collapse = ", "),
           call. = FALSE)
    L <- L[L$factor %in% factors, , drop = FALSE]
  }
  if (all(is.na(L$day)))
    stop("Index names carry no interval information; was W built by env_indices()?",
         call. = FALSE)
  L <- L[!is.na(L$day), , drop = FALSE]

  vexp <- 100 * x$variance$prop[pc]
  if (is.null(title))
    title <- sprintf("Loadings of PC%d (%.1f%% of variance) through the season",
                     pc, vexp)

  if (.rel_has_ggplot()) {
    p <- ggplot2::ggplot(L, ggplot2::aes(x = .data[["day"]], y = .data[["loading"]],
                                         colour = .data[["factor"]]))
    if (isTRUE(points))
      p <- p + ggplot2::geom_point(alpha = 0.35, size = 1)
    p <- p +
      ggplot2::geom_hline(yintercept = 0, colour = "grey60", linetype = 2) +
      ggplot2::geom_smooth(method = "loess", span = span, se = FALSE,
                           formula = y ~ x, linewidth = 0.9) +
      ggplot2::labs(x = "Days after planting", y = paste0("PC", pc, " loading"),
                    colour = NULL, title = title) +
      ggplot2::theme_bw()
    if (isTRUE(facet))
      p <- p + ggplot2::facet_wrap(~ factor) + ggplot2::theme(legend.position = "none")
    print(p)
    attr(L, "plot") <- p
    return(invisible(L))
  }

  # ---- base-graphics fallback ----
  fac <- sort(unique(L$factor))
  cols <- grDevices::hcl.colors(max(2, length(fac)), "Dark 3")
  op <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(op), add = TRUE)

  if (isTRUE(facet)) {
    nr <- ceiling(sqrt(length(fac)))
    graphics::par(mfrow = c(nr, ceiling(length(fac) / nr)), mar = c(4, 4, 2, 1))
    for (i in seq_along(fac)) {
      s <- L[L$factor == fac[i], ]
      s <- s[order(s$day), ]
      graphics::plot(s$day, s$loading, type = "n", xlab = "Days after planting",
                     ylab = paste0("PC", pc, " loading"), main = fac[i])
      graphics::abline(h = 0, col = "grey60", lty = 2)
      if (isTRUE(points)) graphics::points(s$day, s$loading, pch = 16, cex = .5,
                                           col = cols[i])
      .rel_add_loess(s, span, cols[i])
    }
    return(invisible(L))
  }

  graphics::par(mar = c(4, 4, 3, 8))
  graphics::plot(L$day, L$loading, type = "n", xlab = "Days after planting",
                 ylab = paste0("PC", pc, " loading"), main = title)
  graphics::abline(h = 0, col = "grey60", lty = 2)
  for (i in seq_along(fac)) {
    s <- L[L$factor == fac[i], ]
    s <- s[order(s$day), ]
    if (isTRUE(points))
      graphics::points(s$day, s$loading, pch = 16, cex = .4, col = cols[i])
    .rel_add_loess(s, span, cols[i])
  }
  graphics::legend("topright", legend = fac, col = cols, lwd = 2, bty = "n",
                   cex = 0.7, inset = c(-0.28, 0), xpd = TRUE)
  invisible(L)
}

#' Add a loess curve, falling back to a plain line when loess cannot fit
#'
#' @keywords internal
#' @noRd
.rel_add_loess <- function(s, span, col) {
  if (nrow(s) < 4) { graphics::lines(s$day, s$loading, col = col, lwd = 2); return(invisible()) }
  fit <- try(stats::loess(loading ~ day, data = s, span = span), silent = TRUE)
  if (inherits(fit, "try-error")) {
    graphics::lines(s$day, s$loading, col = col, lwd = 2)
  } else {
    xx <- seq(min(s$day), max(s$day), length.out = 200)
    graphics::lines(xx, stats::predict(fit, data.frame(day = xx)), col = col, lwd = 2)
  }
  invisible()
}


# =========================================================================
# SECTION 6 - associating PCs with an external variable
# =========================================================================

#' @title Correlate Environmental PCs with a Per-Environment Variable
#'
#' @description
#' Pearson correlation between each environmental PC and an external
#' per-environment quantity -- in the paper, module eigenmarkers; equally, a trait
#' mean, heritability or predictive ability per environment. P-values are
#' FDR-adjusted across PCs, as in the paper.
#'
#' @param x an \code{\link{env_pca}} object.
#' @param y numeric vector or matrix. One value (or column of values) per environment.
#'   Names, or rownames, are matched against the environments in \code{x}.
#' @param n.pc integer. How many PCs to test. Default \code{min(9, available)},
#'   9 being the paper's choice.
#' @param method character. Correlation method passed to \code{\link[stats]{cor.test}}.
#'   Default \code{"pearson"}.
#' @param p.adjust character. Multiple-testing correction. Default \code{"BH"}
#'   (Benjamini-Hochberg), as in the paper.
#' @param verbose boolean. Report the strongest association. Default \code{TRUE}.
#'
#' @return
#' A data.frame with \code{column}, \code{PC}, \code{r}, \code{p}, \code{p.adj} and
#' \code{n}, sorted by adjusted p-value.
#'
#' @details
#' \strong{Power warning.} With \eqn{q} environments each correlation has \eqn{q - 2}
#' degrees of freedom. At the paper's \eqn{q = 9} even \eqn{r = 0.86} gives an
#' unadjusted \eqn{p \approx 0.003}, and a handful of PCs tested against a handful of
#' modules exhausts the evidence quickly. Treat these associations as hypothesis-
#' generating, not confirmatory, and prefer many environments.
#'
#' \strong{FDR does not cover module selection.} The adjustment applied here accounts
#' for testing several PCs. It does \emph{not} account for the eigenmarker being PC1 of
#' a module that was itself chosen by clustering the same effect matrix. That selection
#' inflates correlations and is invisible to any p-value computed from the observed
#' data. See the worked example below for a permutation check, and for a demonstration
#' that the permutation check is itself insufficient.
#'
#' @examples
#' \donttest{
#' pc <- env_pca(env_indices(maizeWTH, end.day = 130))
#'
#' ## Trait mean per environment
#' ybar <- tapply(maizeYield$value, maizeYield$env, mean)
#' env_pc_associate(pc, ybar)
#' }
#'
#' \dontrun{
#' ## ----------------------------------------------------------------------
#' ## WORKED EXAMPLE: marker effect networks -> eigenmarkers -> environmental PCs
#' ##
#' ## The full analysis chain of Della Coletta et al. (2023), on simulated data
#' ## where the answer is known in advance so the method can be checked.
#' ##
#' ## Planted truth:
#' ##   module A markers respond to mid-season TEMPERATURE
#' ##   module B markers respond to late-season PRECIPITATION
#' ##   module C markers have real but environment-INDEPENDENT effects
#' ## ----------------------------------------------------------------------
#'
#' ## --- 1. Marker effects per environment (RR-BLUP) ----------------------
#' ## Solve via the n x n system: with p markers >> n genotypes this is far
#' ## cheaper than the p x p form and gives the same ridge solution.
#' rrblup_effects <- function(y, Z, h2 = 0.5) {
#'   yc <- as.numeric(y) - mean(y)
#'   lambda <- ncol(Z) * (1 - h2) / h2
#'   a <- solve(tcrossprod(Z) + lambda * diag(nrow(Z)), yc)
#'   as.numeric(crossprod(Z, a))
#' }
#' beta.hat <- sapply(colnames(Y), function(e) rrblup_effects(Y[, e], M, h2 = 0.55))
#' dimnames(beta.hat) <- list(colnames(M), colnames(Y))
#'
#' ## --- 2. Marker effect network -> modules ------------------------------
#' ## Markers whose effects do not VARY across environments cannot covary with
#' ## anything; screening them out first avoids NA correlations and spurious
#' ## modules. The paper used WGCNA; this is the same idea in base R.
#' eff.sd <- apply(beta.hat, 1, sd)
#' Bk     <- beta.hat[eff.sd > quantile(eff.sd, 0.35), , drop = FALSE]
#' adj    <- abs(cor(t(Bk)))^6                 # soft-thresholding power
#' mods   <- cutree(hclust(as.dist(1 - adj), method = "average"), h = 0.92)
#' mods   <- mods[mods %in% as.integer(names(table(mods))[table(mods) >= 15])]
#'
#' ## --- 3. Eigenmarkers --------------------------------------------------
#' ## PC1 of a module's effect profiles, the analogue of a WGCNA module
#' ## eigengene. Sign is arbitrary in PCA, so orient it deterministically --
#' ## otherwise a module can flip between runs and reverse every downstream
#' ## correlation.
#' eigenmarker <- function(B.mod) {
#'   X  <- scale(t(B.mod))
#'   em <- prcomp(X, center = FALSE, scale. = FALSE)$x[, 1]
#'   if (mean(cor(em, X)) < 0) em <- -em
#'   em
#' }
#' ids <- sort(unique(mods))
#' EM  <- sapply(ids, function(m) eigenmarker(Bk[names(mods)[mods == m], , drop = FALSE]))
#' dimnames(EM) <- list(colnames(beta.hat), paste0("ME", ids))
#'
#' ## --- 4. Environmental PCA and association -----------------------------
#' W  <- env_indices(wth, interval = 3, end.day = 151,
#'                   statistic.by = c(PRECTOT = "sum", ETP = "sum"))
#' pc <- env_pca(W)
#' assoc <- env_pc_associate(pc, EM, n.pc = 9, p.adjust = "BH")
#' head(assoc)
#' ##  column  PC     r        p    p.adj  n
#' ##     ME1 PC1  0.99 1.91e-17 6.88e-16 20
#' ##     ME9 PC2  0.94 1.27e-09 2.29e-08 20
#' ##    ME12 PC2 -0.73 2.77e-04 3.33e-03 20
#'
#' ## --- 5. Name the drivers of an associated PC (paper Fig. 6c) ----------
#' k <- 1
#' L <- data.frame(factor = pc$meta$factor, day = pc$meta$day,
#'                 loading = pc$loadings[, k])
#' sort(tapply(abs(L$loading), L$factor, mean), decreasing = TRUE)[1:4]
#' ##  PC1 -> T2M, GDD, T2M_MIN, T2M_MAX   (temperature: module A recovered)
#' ##  PC2 -> PRECTOT, PETP, RH2M          (precipitation: module B recovered)
#' env_loading_curve(pc, pc = k)
#'
#' ## --- 6. Permutation null, and why it is NOT enough --------------------
#' ## Permuting environment labels of the eigenmarker breaks the
#' ## environment<->module link while preserving both structures.
#' perm.max <- replicate(500, {
#'   EMp <- EM[sample(nrow(EM)), , drop = FALSE]
#'   rownames(EMp) <- rownames(EM)
#'   max(abs(env_pc_associate(pc, EMp, n.pc = 9, verbose = FALSE)$r))
#' })
#' quantile(perm.max, 0.95)   # ~0.66: any |r| below this is noise-compatible
#'
#' ## IMPORTANT NEGATIVE RESULT. In this simulation module ME12 contained no
#' ## enriched markers at all -- it is a clustering artefact -- yet it reached
#' ## r = -0.73 (FDR q = 0.003) AND passed this permutation test (p = 0.008).
#' ## The permutation preserves the module definitions, so it cannot detect a
#' ## module that was spuriously assembled in the first place. Guard against
#' ## this with module stability across resampled genotypes, or by requiring
#' ## an independent set of environments -- not with permutation alone.
#' }
#'
#' @importFrom stats cor.test p.adjust
#' @export
env_pc_associate <- function(x, y, n.pc = NULL, method = "pearson",
                             p.adjust = "BH", verbose = TRUE) {

  .et_banner("env_pc_associate", "correlates PCs with an external variable", verbose)
  if (!inherits(x, "env_pca")) stop("'x' must be an env_pca object.", call. = FALSE)

  if (is.null(dim(y)))
    y <- matrix(y, ncol = 1, dimnames = list(names(y), "y"))
  y <- .rel_as_matrix(y, "y")
  if (is.null(colnames(y))) colnames(y) <- paste0("y", seq_len(ncol(y)))

  S <- x$scores
  if (!is.null(rownames(y)) && !is.null(rownames(S))) {
    common <- intersect(rownames(S), rownames(y))
    if (length(common) < 3)
      stop("Fewer than three environments shared between the PCA and 'y'.",
           call. = FALSE)
    if (isTRUE(verbose) && length(common) < nrow(S))
      message("Using ", length(common), " environments shared by the PCA and 'y'.")
    S <- S[common, , drop = FALSE]
    y <- y[common, , drop = FALSE]
  } else if (nrow(y) != nrow(S)) {
    stop("'y' has ", nrow(y), " rows but the PCA has ", nrow(S),
         " environments, and names are missing so they cannot be matched.",
         call. = FALSE)
  }

  if (is.null(n.pc)) n.pc <- min(9L, ncol(S))
  n.pc <- min(n.pc, ncol(S))

  res <- do.call(rbind, lapply(seq_len(ncol(y)), function(j) {
    do.call(rbind, lapply(seq_len(n.pc), function(k) {
      ct <- try(stats::cor.test(S[, k], y[, j], method = method), silent = TRUE)
      if (inherits(ct, "try-error"))
        return(data.frame(column = colnames(y)[j], PC = paste0("PC", k),
                          r = NA_real_, p = NA_real_,
                          n = sum(stats::complete.cases(S[, k], y[, j])),
                          stringsAsFactors = FALSE))
      data.frame(column = colnames(y)[j], PC = paste0("PC", k),
                 r = unname(ct$estimate), p = ct$p.value,
                 n = sum(stats::complete.cases(S[, k], y[, j])),
                 stringsAsFactors = FALSE)
    }))
  }))

  res$p.adj <- stats::p.adjust(res$p, method = p.adjust)
  res <- res[order(res$p.adj, -abs(res$r)), c("column", "PC", "r", "p", "p.adj", "n")]
  rownames(res) <- NULL

  if (isTRUE(verbose) && nrow(res)) {
    b <- res[1, ]
    message(sprintf("Strongest: %s ~ %s, r = %.2f, p = %.3g, FDR p = %.3g (n = %d)",
                    b$column, b$PC, b$r, b$p, b$p.adj, b$n))
    if (b$n < 12)
      message("  Note: n = ", b$n,
              " environments. Low power; treat as hypothesis-generating.")
  }

  res
}