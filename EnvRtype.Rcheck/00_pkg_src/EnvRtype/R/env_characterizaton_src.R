#==================================================================================================
# Title.    : Environmental Characterization Module
# Author.   : G Costa-Neto
# Updated at: 2026-09 (EnvRtype 1.2.0)
#
# Module contents
#   W_matrix()              - environmental covariable matrix (q environments x k covariables)
#   env_typing()            - envirotypes by cardinals, quantiles or data-driven mining
#   env_cluster()           - mega-environment discovery by K-means + silhouette
#   env_target_importance() - environmental drivers of a target-based grouping (Random Forest)
#
# Internal helpers (not exported) are prefixed with '.env_'.
#==================================================================================================


# =========================================================================
# SECTION 0 - shared internal helpers
# =========================================================================

#' Coerce an environmental covariable object to a numeric matrix
#'
#' @param x matrix or data.frame of environmental covariables.
#' @param arg character. Argument name used in error messages.
#' @return A numeric matrix.
#' @keywords internal
#' @noRd
.env_as_matrix <- function(x, arg = "W") {
  if (is.null(x)) stop("'", arg, "' is missing.", call. = FALSE)
  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.matrix(x)) x <- as.matrix(x)
  if (!is.numeric(x)) {
    storage.mode(x) <- "numeric"
    if (!is.numeric(x)) stop("'", arg, "' must be numeric.", call. = FALSE)
  }
  x
}

#' Select the number of clusters by maximising the average silhouette width
#'
#' Shared by \code{env_cluster()} and \code{env_target_importance()} so that both
#' functions use exactly the same rule to pick \code{k}.
#'
#' @param X numeric matrix used for clustering (observations in rows).
#' @param k.max integer. Largest number of clusters tested.
#' @param nstart integer. Random starts passed to \code{stats::kmeans}.
#' @param verbose boolean. Print the chosen \code{k}.
#' @param label character. Wording used in the progress message.
#' @return A list with \code{k} and the \code{silhouette} data.frame.
#' @keywords internal
#' @noRd
.env_best_k <- function(X, k.max = 10, nstart = 25, verbose = TRUE,
                        label = "clusters") {
  if (!requireNamespace("cluster", quietly = TRUE))
    stop("Package 'cluster' is required for silhouette optimisation. ",
         "Install it with install.packages('cluster') or supply 'k'.", call. = FALSE)

  k.max <- min(k.max, nrow(X) - 1L)
  if (k.max < 2L) stop("Not enough environments to form clusters.", call. = FALSE)

  d <- stats::dist(X)
  widths <- vapply(2:k.max, function(g) {
    km <- stats::kmeans(X, centers = g, nstart = nstart)
    mean(cluster::silhouette(km$cluster, d)[, 3])
  }, numeric(1))

  sil <- data.frame(k = 2:k.max, avg.silhouette = widths)
  k   <- sil$k[which.max(sil$avg.silhouette)]
  if (isTRUE(verbose))
    message("Optimal number of ", label, " (silhouette): k = ", k)

  list(k = k, silhouette = sil)
}

#' Assign observations to time intervals (development stages)
#'
#' \code{cut()} is called with an open-ended final break, so the user-supplied
#' \code{.breaks} define the LEFT edges and one extra trailing interval
#' \code{[last, Inf)} is always appended.
#'
#' @param .dae numeric vector of days from start.
#' @param .breaks numeric vector of temporal breaks.
#' @param .names character vector of interval names.
#' @param .min.window integer. Minimum number of records required to keep the
#'   open-ended trailing interval. Use 0 or NULL to always keep it.
#' @return A factor of interval membership.
#' @keywords internal
#' @noRd
.env_stage_by_dae <- function(.dae = NULL, .breaks = NULL, .names = NULL,
                              .min.window = 2L) {
  if (is.null(.dae)) stop("'days.id' values are missing.", call. = FALSE)
  if (is.null(.breaks))
    .breaks <- seq(from = 1 - min(.dae, na.rm = TRUE),
                   to = max(.dae, na.rm = TRUE) + 10, by = 10)
  .breaks <- sort(unique(as.numeric(.breaks)))

  .breaks.full <- c(.breaks, Inf)
  n.int <- length(.breaks.full) - 1L

  if (is.null(.names)) {
    .names <- paste0("Interval_", .breaks)
  } else if (length(.names) == n.int - 1L) {
    .names <- c(.names, paste0("Interval_", .breaks[length(.breaks)]))
  } else if (length(.names) != n.int) {
    stop("names.window has ", length(.names), " entries but time.window = c(",
         paste(.breaks, collapse = ","), ") implies ", n.int,
         " intervals (an open-ended final interval [",
         .breaks[length(.breaks)], ",Inf) is always appended). Supply ",
         n.int, " names, or ", n.int - 1L,
         " and the last one will be generated.", call. = FALSE)
  }

  pstage <- cut(x = .dae, breaks = .breaks.full, right = FALSE)
  levels(pstage) <- .names

  # Drop a degenerate trailing window fitted on a handful of points.
  if (!is.null(.min.window) && .min.window > 0) {
    tb <- table(pstage)
    if (length(tb) > 1L && tb[length(tb)] < .min.window)
      pstage <- factor(pstage, levels = .names[-length(.names)])
  }
  pstage
}

#' Melt a wide weather data.frame into long format
#'
#' @param .GeTw data.frame. A \code{get_weather()}-like object.
#' @param var.id character vector of variables to melt.
#' @param by.interval boolean. Compute temporal intervals.
#' @param days character. Name of the days-from-start column.
#' @param time.window,names.window interval definition.
#' @param id.names character vector of id columns.
#' @param env.id character. Environment id column.
#' @param min.window integer. See \code{.env_stage_by_dae}.
#' @return A long data.frame with \code{env}, \code{variable} and \code{value}.
#' @keywords internal
#' @noRd
.env_melt_wth <- function(.GeTw, var.id = NULL, by.interval = FALSE, days = NULL,
                          time.window = NULL, names.window = NULL, id.names = NULL,
                          env.id = NULL, min.window = 2L) {

  .GeTw <- as.data.frame(.GeTw)
  if (is.null(env.id)) env.id <- ".id"
  if (is.null(days))   days   <- "daysFromStart"

  names(.GeTw)[names(.GeTw) %in% env.id] <- "env"
  names(.GeTw)[names(.GeTw) %in% days]   <- "daysFromStart"

  if (is.null(id.names)) {
    # default = NASA-POWER column set from get_weather(); keep only what is
    # present so that arbitrary data.frames can be used as well.
    id.names <- c("env", "LON", "LAT", "YEAR", "MM", "DD", "DOY",
                  "YYYYMMDD", "daysFromStart")
    id.names <- id.names[id.names %in% names(.GeTw)]
  } else {
    .miss <- id.names[!id.names %in% names(.GeTw)]
    if (length(.miss))
      stop("id variables not found in data: ", paste(.miss, collapse = ", "),
           call. = FALSE)
  }

  if (is.null(var.id)) var.id <- names(.GeTw)[!names(.GeTw) %in% id.names]
  .missv <- var.id[!var.id %in% names(.GeTw)]
  if (length(.missv))
    stop("variables not found in data: ", paste(.missv, collapse = ", "),
         call. = FALSE)

  if (isTRUE(by.interval)) {
    if (!"daysFromStart" %in% names(.GeTw))
      stop("by.interval = TRUE requires a days-from-start column ",
           "(see the 'days.id' argument).", call. = FALSE)
    .GeTw$interval <- .env_stage_by_dae(.dae = .GeTw[, "daysFromStart"],
                                        .breaks = time.window,
                                        .names = names.window,
                                        .min.window = min.window)
    # rows falling in a dropped (empty) trailing window become NA -> remove
    .keep <- !is.na(.GeTw$interval)
    if (any(!.keep)) .GeTw <- .GeTw[.keep, , drop = FALSE]
    id.names <- c(id.names, "interval")
  }

  ts <- suppressWarnings(
    reshape2::melt(.GeTw, measure.vars = var.id, id.vars = id.names))
  ts$value <- suppressWarnings(as.numeric(ts$value))
  if (isTRUE(by.interval)) ts$interval <- droplevels(as.factor(ts$interval))
  ts
}


# =========================================================================
# SECTION 1 - W_matrix()
# =========================================================================

#' @title Build the Matrix of Environmental Covariables for Reaction Norm
#'
#' @description
#' Estimates a weather covariable realized matrix \strong{W} with dimensions
#' \eqn{q \times k}, for \eqn{q} environments and \eqn{k} covariables. This matrix is
#' the basic input of the reaction-norm and enviromic-assisted models in
#' \pkg{EnvRtype} (\code{env_kernel}, \code{env_cluster},
#' \code{env_target_importance}, \code{get_kernel}).
#'
#' @author Germano Costa Neto
#'
#' @param env.data data.frame of environmental variables from \code{get_weather},
#'   or an already summarised object when \code{is.processed = TRUE}.
#' @param is.processed boolean. Indicates whether the data.frame was previously processed
#'   with \code{summaryWTH} and already contains means, medians, etc.
#' @param id.names character. Columns used as id for the environmental variables.
#' @param env.id character. Column used as id for environments.
#' @param var.id vector (character). Which variables will be used in the analysis.
#' @param probs vector (numeric). Probability quantiles in \eqn{[0,1]}, used when
#'   \code{statistic = 'quantile'}. If \code{NULL}, \code{c(0.25, 0.50, 0.75)}.
#' @param by.interval boolean. Indicates if temporal intervals must be computed inside
#'   each environment. Default \code{FALSE}.
#' @param time.window vector (numeric). If \code{by.interval = TRUE}, the temporal breaks.
#' @param names.window vector (character). If \code{by.interval = TRUE}, the interval names.
#' @param center boolean. Indicates whether the matrix should be centred. Default \code{TRUE}.
#' @param scale boolean. If \code{TRUE}, variables assume \eqn{x \sim N(0,1)}. Default \code{TRUE}.
#' @param sd.tol numeric. Maximum standard deviation tolerated in quality control. Default 10.
#' @param statistic vector (character). Statistic to be computed,
#'   \code{c('all','sum','mean','quantile')}. Default \code{'mean'}.
#' @param tol numeric. Numerical tolerance; variables whose standard deviation is at or
#'   below \code{tol} are treated as near-constant and flagged for removal. Default 1E-3.
#' @param QC boolean. Indicates whether Quality Control is applied. QC removes variables
#'   with \code{sd(x) > sd.tol} and near-constant variables.
#' @param impute character. How to handle missing values before scaling: \code{"none"}
#'   (default, keep NAs), \code{"mean"} (column mean) or \code{"drop"} (remove columns with any NA).
#' @param verbose boolean. If \code{TRUE} (default) prints quality-control messages.
#'
#' @return
#' An environmental covariable realized matrix with dimensions \eqn{q \times k}.
#' The centring and scaling values, plus the list of removed markers, are attached as
#' attributes (\code{"scaled:center"}, \code{"scaled:scale"}, \code{"removed"}) so that new
#' environments can be projected onto the same space.
#'
#' @details
#' Quality control follows Morais Junior et al. (2018): covariables whose standard deviation
#' across environments exceeds \code{sd.tol} are discarded, as are near-constant covariables
#' (\eqn{sd \le tol}) which carry no information about environmental differences.
#'
#' Unlike earlier versions, the numerical tolerance is \emph{not} added to the data before
#' scaling (which silently shifted unscaled outputs); it is used only to detect
#' near-constant columns.
#'
#' @examples
#' \donttest{
#' data("maizeWTH")
#' env.data <- maizeWTH[maizeWTH$daysFromStart < 100, ]
#'
#' ## Mean-centred and scaled matrix (default statistic = 'mean')
#' W <- W_matrix(env.data = env.data)
#' dim(W)
#'
#' ## Adding time windows (one block of covariables per development stage)
#' W <- W_matrix(env.data = env.data, by.interval = TRUE,
#'               time.window = c(0, 14, 35, 60, 90, 120))
#'
#' ## Selecting the statistic to be used
#' W <- W_matrix(env.data = env.data, by.interval = TRUE, statistic = 'quantile',
#'               time.window = c(0, 14, 35, 60, 90, 120))
#'
#' ## With Quality Control based on the maximum sd tolerated
#' W <- W_matrix(env.data = env.data, QC = TRUE, sd.tol = 3)
#' attr(W, "removed")
#'
#' ## Creating W for specific variables
#' W <- W_matrix(env.data = env.data, var.id = c('T2M_MAX', 'T2M_MIN', 'T2M'))
#'
#' ## Combining with summaryWTH by using is.processed = TRUE
#' data <- summaryWTH(env.data, env.id = 'env', statistic = 'quantile')
#' W <- W_matrix(env.data = data, is.processed = TRUE)
#' }
#'
#' @seealso \code{summaryWTH}, \code{env_kernel}, \code{env_typing},
#'   \code{\link{T_matrix}}
#'
#' @references
#' Costa-Neto G., Galli G., Carvalho H.F., Crossa J., Fritsche-Neto R. (2021).
#' EnvRtype: a software to interplay enviromics and quantitative genomics in agriculture.
#' \emph{G3} 11(4), jkab040.
#'
#' @param copula character or NULL. If not \code{NULL}, the covariable matrix is replaced by
#'   copula-based indices computed by \code{env_copula}. One of \code{"none"} (default,
#'   same as \code{NULL}), \code{"pobs"}, \code{"joint"}, \code{"survival"}, \code{"kendall"}
#'   or \code{"all"}. See \code{env_copula} for the meaning of each index.
#' @param copula.args list. Extra arguments passed to \code{env_copula}, e.g.
#'   \code{list(groups = list(heat = c("T2M_mean","T2M_MAX_mean")), ties.method = "average")}.
#'
#' @importFrom stats sd
#' @importFrom reshape2 acast melt
#' @export
W_matrix <- function(env.data, is.processed = FALSE, id.names = NULL, env.id = NULL,
                     var.id = NULL, probs = NULL, by.interval = NULL,
                     time.window = NULL, names.window = NULL,
                     center = TRUE, scale = TRUE, sd.tol = 10, statistic = NULL,
                     tol = 1E-3, QC = FALSE,
                     impute = c("none", "mean", "drop"), verbose = TRUE,
                     copula = NULL, copula.args = list()) {

  .et_banner("W_matrix", "builds the environmental covariable (W) matrix", verbose)
  impute <- match.arg(impute)
  if (is.null(statistic))   statistic   <- "mean"
  if (is.null(by.interval)) by.interval <- FALSE

  if (isFALSE(is.processed)) {
    .et_step("summarising weather data into environmental covariables", verbose)
    W <- summaryWTH(env.data = env.data, id.names = id.names, env.id = env.id,
                    statistic = statistic, probs = probs, var.id = var.id,
                    by.interval = by.interval, time.window = time.window,
                    names.window = names.window, verbose = FALSE)

    colid <- c("env", "variable", "interval")
    W <- reshape2::melt(W, measure.vars = names(W)[!names(W) %in% colid],
                        variable.name = "stat")
    W$variable <- paste(W$variable, W$stat, sep = "_")
    W <- if (isTRUE(by.interval)) {
      reshape2::acast(W, env ~ variable + interval, value.var = "value")
    } else {
      reshape2::acast(W, env ~ variable, value.var = "value")
    }
  }

  # ---- optional copula reparameterisation of the covariable space ----------
  # Done BEFORE centring/scaling: copula indices are rank-based and already
  # live on [0,1], so scaling them afterwards is opt-in via 'center'/'scale'.
  if (!is.null(copula) && !identical(copula, "none")) {
    .et_step("reparameterising covariables with copula index", verbose)
    W <- do.call(env_copula,
                 c(list(W = W, index = copula, verbose = verbose), copula.args))
    if (isTRUE(verbose))
      message("W replaced by copula index '", copula, "' (", ncol(W), " columns).")
  }

  .et_step("centring, scaling and quality-controlling W", verbose)
  .env_w_scale(env.data = W, center = center, scale = scale, sd.tol = sd.tol,
               tol = tol, QC = QC, impute = impute, verbose = verbose)
}

#' Centre, scale and quality-control an environmental covariable matrix
#'
#' @inheritParams W_matrix
#' @return A scaled matrix carrying the \code{"removed"} attribute.
#' @keywords internal
#' @noRd
.env_w_scale <- function(env.data, center = TRUE, scale = TRUE, sd.tol = 10,
                         tol = 1E-3, QC = FALSE, impute = "none", verbose = TRUE) {

  env.data <- .env_as_matrix(env.data, "env.data")

  # ---- missing-value handling, before any statistic is computed ----
  if (anyNA(env.data)) {
    if (impute == "mean") {
      cm  <- colMeans(env.data, na.rm = TRUE)
      idx <- which(is.na(env.data), arr.ind = TRUE)
      env.data[idx] <- cm[idx[, 2]]
    } else if (impute == "drop") {
      env.data <- env.data[, colSums(is.na(env.data)) == 0, drop = FALSE]
    }
  }
  if (!ncol(env.data))
    stop("No environmental covariables left after missing-value handling.", call. = FALSE)

  # ---- quality control: too-variable and near-constant covariables ----
  sdA <- apply(env.data, 2, stats::sd, na.rm = TRUE)
  t   <- ncol(env.data)
  removed <- unique(c(names(sdA[sdA > sd.tol]),
                      names(sdA[sdA <= tol | is.na(sdA)])))

  # NOTE: the tolerance is NOT added to the data (that shifted unscaled output);
  # it is used only to flag near-constant columns above.
  env.data <- scale(env.data, center = center, scale = scale)

  if (isTRUE(QC)) {
    keep <- !colnames(env.data) %in% removed
    if (!any(keep))
      stop("Quality control removed every covariable; relax 'sd.tol' or set QC = FALSE.",
           call. = FALSE)
    ctr <- attr(env.data, "scaled:center")
    scl <- attr(env.data, "scaled:scale")
    env.data <- env.data[, keep, drop = FALSE]
    if (!is.null(ctr)) attr(env.data, "scaled:center") <- ctr[keep]
    if (!is.null(scl)) attr(env.data, "scaled:scale")  <- scl[keep]

    if (isTRUE(verbose)) {
      message(strrep("-", 48))
      message("Quality Control based on sd.tol = ", sd.tol)
      message("Removed variables: ", length(removed), " from ", t)
      if (length(removed)) message(paste(removed, collapse = "\n"))
      message(strrep("-", 48))
    }
  }

  attr(env.data, "removed") <- removed
  env.data
}

# =========================================================================
# SECTION 1B - env_copula()
# =========================================================================

#' @title Copula-based Dependence Indices for an Environmental Covariable Matrix
#'
#' @description
#' Performs a copula analysis of an environment x covariable matrix (\code{W}) and returns a
#' matrix of \emph{copula indices} that can replace the raw covariables. Sklar's theorem lets
#' the joint distribution of the environmental variables be split into (i) the marginal
#' behaviour of each variable and (ii) the dependence structure linking them (the copula).
#' By moving to the copula scale, each environment is described by where it sits in the
#' \emph{multivariate} environmental distribution rather than by raw physical units.
#'
#' @author Germano Costa Neto
#'
#' @param W matrix or data.frame. Environmental covariables with environments in rows and
#'   variables in columns (typically a \code{W_matrix} output).
#' @param index character. Which copula index to return:
#'   \describe{
#'     \item{\code{"pobs"}}{Pseudo-observations \eqn{u_{ij} = r_{ij}/(q+1)}, the marginal
#'       probability-integral transform. Same dimension as \code{W}; one column per variable.}
#'     \item{\code{"joint"}}{Multivariate empirical copula \eqn{C_q(u_i)}, the joint
#'       non-exceedance probability of environment \eqn{i}: the probability that all
#'       variables are simultaneously at or below their observed levels. One column per block.}
#'     \item{\code{"survival"}}{Survival (joint exceedance) copula: the probability that all
#'       variables are simultaneously at or above their observed levels. Highlights
#'       compound-stress environments. One column per block.}
#'     \item{\code{"kendall"}}{Kendall function \eqn{K_C(z) = P(C(U) \le z)} evaluated at each
#'       environment's joint probability. A scalar multivariate return level ordering
#'       environments from jointly-benign to jointly-extreme. One column per block.}
#'     \item{\code{"all"}}{Column-binds \code{pobs}, \code{joint}, \code{survival}, \code{kendall}.}
#'   }
#' @param groups list or NULL. Optional named list splitting the columns of \code{W} into
#'   variable blocks, e.g. \code{list(heat = c("T2M_mean","T2M_MAX_mean"), water = "PRECTOT_mean")}.
#'   Elements may be column names or column indices. A separate multivariate copula is fitted
#'   per block, so joint/survival/Kendall describe dependence \emph{within} a block. If
#'   \code{NULL} (default) all columns form a single block named \code{"all"}.
#' @param ties.method character. How ranks are computed when there are ties, passed to
#'   \code{rank}. Default \code{"average"}.
#' @param normalize boolean. If \code{TRUE} (default) the Kendall index is rescaled to
#'   \eqn{[0,1]} by its empirical range when non-degenerate; \code{FALSE} keeps raw \eqn{K_C}.
#' @param verbose boolean. If \code{TRUE} (default) prints progress messages.
#'
#' @return
#' A numeric matrix with \eqn{q} environments in rows (rownames preserved from \code{W}) and
#' the requested copula indices in columns. Attributes \code{"copula.index"} and
#' \code{"copula.groups"} record the index produced and the block structure used.
#'
#' @details
#' \strong{Why rank-based.} All indices are computed from ranks, so they are invariant to any
#' monotone transformation of the marginals: storing precipitation in mm or log(mm) gives the
#' same result, and no distributional assumption is imposed. This is the main practical gain
#' over a z-scored \strong{W}, where one extreme site can dominate the Euclidean geometry.
#'
#' \strong{Pseudo-observations.} For \eqn{q} environments and variable \eqn{j},
#' \deqn{u_{ij} = r_{ij}/(q+1)}
#' with \eqn{r_{ij}} the rank of environment \eqn{i}. Dividing by \eqn{q+1} keeps values
#' strictly inside \eqn{(0,1)}, avoiding boundary problems in downstream density evaluation.
#'
#' \strong{Empirical copula.} The joint index is the empirical copula evaluated at each
#' environment's own pseudo-observation, a nonparametric estimator of \eqn{C} that requires no
#' copula family to be chosen.
#'
#' \strong{Kendall function.} \eqn{K_C} is the CDF of \eqn{Z = C(U)}, estimated by the
#' Genest-Rivest empirical construction. It collapses a \eqn{d}-dimensional dependence
#' structure into one interpretable number per environment and underpins multivariate return
#' periods in hydrology.
#'
#' \strong{Sample-size caveat.} These are rank statistics \emph{across environments}, so their
#' resolution is bounded by \eqn{q}. With five environments each pseudo-observation can only
#' take values in 1/6,...,5/6 and the joint copula is very coarse; the indices become
#' informative with dozens of environments. A warning is issued when \eqn{q < 10}.
#'
#' \strong{Scope.} The dependence modelled is that \emph{among environmental variables across
#' the environment panel}. Geographic coordinates are not used and no spatial autocorrelation
#' (variogram / distance decay) is fitted. Use \code{groups} to control which variables share
#' a copula.
#'
#' @examples
#' \donttest{
#' data("maizeWTH")
#' env.data <- maizeWTH[maizeWTH$daysFromStart < 100, ]
#'
#' W <- W_matrix(env.data = env.data, var.id = c("T2M", "T2M_MAX", "PRECTOT"),
#'               statistic = "mean")
#'
#' ## 1. Marginal copula scale (pseudo-observations)
#' round(env_copula(W, index = "pobs"), 3)
#'
#' ## 2. Joint non-exceedance probability per environment
#' env_copula(W, index = "joint")
#'
#' ## 3. Compound-stress (joint exceedance) index
#' env_copula(W, index = "survival")
#'
#' ## 4. Kendall multivariate return level
#' env_copula(W, index = "kendall")
#'
#' ## 5. Blocks: heat and water handled as separate copulas
#' env_copula(W, index = "all",
#'            groups = list(heat  = c("T2M_mean", "T2M_MAX_mean"),
#'                          water = "PRECTOT_mean"))
#'
#' ## 6. Straight from W_matrix() via the 'copula' argument
#' Wc <- W_matrix(env.data = env.data, var.id = c("T2M", "T2M_MAX", "PRECTOT"),
#'                statistic = "mean", copula = "kendall",
#'                center = FALSE, scale = FALSE)
#'
#' ## 7. Copula-scale kernel for reaction-norm models
#' Kc <- env_kernel(env.data = env_copula(W, index = "pobs"), gaussian = TRUE)
#' round(Kc$envCov, 2)
#' }
#'
#' @seealso \code{W_matrix}, \code{env_kernel}, \code{env_cluster}
#'
#' @references
#' Sklar, A. (1959). Fonctions de repartition a n dimensions et leurs marges.
#' \emph{Publ. Inst. Statist. Univ. Paris} 8, 229-231.
#'
#' Genest, C., Rivest, L.-P. (1993). Statistical inference procedures for bivariate
#' Archimedean copulas. \emph{JASA} 88, 1034-1043.
#'
#' Salvadori, G., De Michele, C., Kottegoda, N.T., Rosso, R. (2007).
#' \emph{Extremes in Nature: An Approach Using Copulas}. Springer.
#'
#' @export
env_copula <- function(W, index = c("pobs", "joint", "survival", "kendall", "all"),
                       groups = NULL, ties.method = "average",
                       normalize = TRUE, verbose = TRUE) {

  .et_banner("env_copula", "reparameterises covariables into copula indices", verbose)
  index <- match.arg(index)
  W <- .env_as_matrix(W, "W")

  q <- nrow(W)
  if (q < 3)
    stop("At least three environments are required for a copula analysis (q = ", q, ").",
         call. = FALSE)
  if (q < 10 && isTRUE(verbose))
    warning("Only ", q, " environments: copula indices are rank statistics across ",
            "environments and will be very coarse. Interpret with care.", call. = FALSE)
  if (anyNA(W))
    stop("W contains missing values; impute them (see W_matrix(impute=)) before env_copula().",
         call. = FALSE)

  if (is.null(colnames(W))) colnames(W) <- paste0("V", seq_len(ncol(W)))
  if (is.null(rownames(W))) rownames(W) <- paste0("env", seq_len(q))

  # ---- resolve the variable blocks ---------------------------------------
  if (is.null(groups)) {
    groups <- list(all = colnames(W))
  } else {
    if (!is.list(groups)) stop("'groups' must be a list.", call. = FALSE)
    if (is.null(names(groups)) || any(!nzchar(names(groups))))
      names(groups) <- paste0("G", seq_along(groups))
    groups <- lapply(groups, function(g) {
      if (is.numeric(g)) {
        if (any(g < 1 | g > ncol(W)))
          stop("'groups' contains out-of-range column indices.", call. = FALSE)
        colnames(W)[g]
      } else {
        miss <- setdiff(g, colnames(W))
        if (length(miss))
          stop("'groups' refers to columns absent from W: ",
               paste(miss, collapse = ", "), call. = FALSE)
        g
      }
    })
  }

  if (isTRUE(verbose))
    message("env_copula: ", q, " environments, ", ncol(W), " covariables, ",
            length(groups), " block(s); index = '", index, "'.")

  # ---- marginal probability integral transform ----------------------------
  U <- .env_pobs(W, ties.method = ties.method)

  if (index == "pobs") {
    colnames(U) <- paste0("cop_u_", colnames(W))
    attr(U, "copula.index")  <- index
    attr(U, "copula.groups") <- groups
    return(U)
  }

  # ---- per-block multivariate indices -------------------------------------
  joint <- surv <- kend <- matrix(NA_real_, q, length(groups),
                                  dimnames = list(rownames(W), names(groups)))

  for (g in seq_along(groups)) {
    Ug <- U[, groups[[g]], drop = FALSE]
    Cj <- .env_emp_copula(Ug)
    joint[, g] <- Cj
    surv[, g]  <- .env_emp_survival(Ug)
    kend[, g]  <- .env_kendall_fun(Cj, normalize = normalize)
  }

  colnames(joint) <- paste0("cop_C_",    names(groups))
  colnames(surv)  <- paste0("cop_surv_", names(groups))
  colnames(kend)  <- paste0("cop_K_",    names(groups))

  out <- switch(index,
                joint    = joint,
                survival = surv,
                kendall  = kend,
                all      = cbind(`colnames<-`(U, paste0("cop_u_", colnames(W))),
                                 joint, surv, kend))

  attr(out, "copula.index")  <- index
  attr(out, "copula.groups") <- groups
  out
}

#' Pseudo-observations (marginal probability integral transform)
#'
#' @param W numeric matrix, environments in rows.
#' @param ties.method passed to \code{rank}.
#' @return Matrix of the same dimension with entries in (0,1).
#' @keywords internal
#' @noRd
.env_pobs <- function(W, ties.method = "average") {
  q <- nrow(W)
  U <- apply(W, 2, function(z) rank(z, ties.method = ties.method) / (q + 1))
  matrix(U, nrow = q, dimnames = dimnames(W))
}
#' Multivariate empirical copula evaluated at the sample points
#'
#' Counts are divided by \code{q + 1}, matching the divisor used for the
#' pseudo-observations in \code{.env_pobs}. Using \code{q} here instead would put
#' C and u on marginally different scales and let C exceed its Frechet-Hoeffding
#' upper bound min(u) by up to 1/q.
#'
#' @param U matrix of pseudo-observations.
#' @return Numeric vector of joint non-exceedance probabilities.
#' @keywords internal
#' @noRd
.env_emp_copula <- function(U) {
  q <- nrow(U)
  d <- ncol(U)
  vapply(seq_len(q), function(i) {
    sum(colSums(t(U) <= U[i, ]) == d) / (q + 1)
  }, numeric(1))
}

#' Multivariate empirical survival copula evaluated at the sample points
#'
#' Uses the same \code{q + 1} divisor as \code{.env_emp_copula} for consistency.
#'
#' @param U matrix of pseudo-observations.
#' @return Numeric vector of joint exceedance probabilities.
#' @keywords internal
#' @noRd
.env_emp_survival <- function(U) {
  q <- nrow(U)
  d <- ncol(U)
  vapply(seq_len(q), function(i) {
    sum(colSums(t(U) >= U[i, ]) == d) / (q + 1)
  }, numeric(1))
}

#' Empirical Kendall distribution function at given joint probabilities
#'
#' Genest-Rivest empirical estimator of \eqn{K_C(z) = P(C(U) \le z)}.
#'
#' @param Cz numeric vector of joint probabilities from \code{.env_emp_copula}.
#' @param normalize boolean. Rescale the result to [0,1] by its range.
#' @return Numeric vector of Kendall-function values.
#' @keywords internal
#' @noRd
.env_kendall_fun <- function(Cz, normalize = TRUE) {
  q <- length(Cz)
  K <- vapply(Cz, function(z) sum(Cz <= z) / (q + 1), numeric(1))
  if (isTRUE(normalize)) {
    rg <- range(K, na.rm = TRUE)
    # A constant K carries no information; rescaling would divide by zero.
    if (diff(rg) > .Machine$double.eps) K <- (K - rg[1]) / diff(rg)
  }
  K
}
# =========================================================================
# SECTION 2 - env_typing()
# =========================================================================

#' @title Environmental Typologies based on Cardinal, Quantilic Limits or Data-Driven Mining
#'
#' @description
#' Returns environmental typologies (\emph{envirotypes}) that can be used as envirotype
#' markers. Typologies are given by cardinals (discrete intervals for each variable), by
#' empirical quantiles, or -- when \code{envirotype_mining = TRUE} -- learned from the data:
#' for each variable a K-means clustering is run over a grid of candidate \code{k} values and
#' the optimal number of envirotypes is selected by the Calinski-Harabasz (variance ratio)
#' criterion.
#'
#' @author Germano Costa Neto
#'
#' @param env.data data.frame of environmental variables from \code{get_weather}.
#' @param var.id character. Which variables will be used in the analysis.
#' @param env.id character. Column used as id for environments.
#' @param cardinals list (numeric). Cardinal thresholds per variable. If \code{NULL}, see
#'   \code{quantiles}. Ignored when \code{envirotype_mining = TRUE}.
#' @param days.id character. Name of the column indicating the days from start.
#' @param time.window vector (numeric). If \code{by.interval = TRUE}, the temporal breaks.
#' @param names.window vector (character). If \code{by.interval = TRUE}, the interval names.
#'   Note that \code{cut()} is applied with an open-ended final break, so
#'   \code{time.window = c(0,40,80,120)} defines FOUR intervals; supply either 4 names or 3,
#'   in which case the trailing one is auto-named.
#' @param quantiles vector (numeric). Probability quantiles used when \code{cardinals} is
#'   \code{NULL}. Default \code{c(0.01, .25, .50, .75, .99)}. Ignored when mining.
#' @param id.names vector (character). Columns used as id for the environmental variables.
#' @param by.interval boolean. Compute temporal intervals inside each environment. Default \code{FALSE}.
#' @param scale boolean. If \code{TRUE}, variables assume \eqn{x \sim N(0,1)}. Default \code{FALSE}.
#' @param format character. Output shape, \code{'long'} (default) or \code{'wide'}.
#' @param ratio boolean. If \code{TRUE}, \code{Freq} is rescaled to a relative frequency in
#'   \eqn{[0,1]} within each environment x interval x variable group. Default \code{FALSE}.
#' @param joint boolean. Only used when \code{envirotype_mining = TRUE}. If \code{FALSE}
#'   (default) one univariate K-means is run per variable. If \code{TRUE}, the variables in
#'   \code{var.id} are COMBINED: the vector of all variables observed at the same time point is
#'   clustered jointly, so each envirotype is a recurring multi-variable state of the
#'   environment (analogous to a haplotype across loci).
#' @param joint.by.interval boolean. Only used when \code{joint = TRUE}. If \code{TRUE}
#'   (default) a separate K-means is fitted inside each time interval; if \code{FALSE} a
#'   single K-means is fitted over all intervals pooled.
#' @param combine.features boolean. Only valid when \code{envirotype_mining = TRUE} and
#'   \code{joint = FALSE}. If \code{TRUE}, envirotypes are mined for every unique
#'   COMBINATION of the variables in \code{var.id}: each variable on its own, then every
#'   pair, every triple, and so on. A combination of two or more variables is clustered
#'   jointly (as in \code{joint = TRUE}); a single-variable combination reduces to the
#'   univariate case. Each combination receives a compact generic id (\code{FC001},
#'   \code{FC002}, ...) and the returned \code{feature_combinations} table maps every id
#'   back to its variables and time intervals. Default \code{FALSE}.
#' @param combine.max.order integer or \code{NULL}. Largest combination size to enumerate
#'   when \code{combine.features = TRUE}. Because the number of combinations grows as a
#'   power set (\eqn{2^n - 1}), this defaults to \code{2} (single variables and pairs
#'   only), which stays tractable for dozens of variables. Increase it to include triples
#'   and higher (e.g. \code{3}), or set \code{NULL} to enumerate all
#'   \code{length(var.id)} orders. A guard stops the run if the requested number of
#'   combinations is impractically large.
#' @param envirotype_mining boolean. If \code{TRUE}, envirotype classes are mined from the
#'   data by K-means with \code{k} chosen by the Calinski-Harabasz index. Default \code{FALSE}.
#' @param k.range integer vector. Candidate numbers of envirotypes tested when mining. Default \code{2:10}.
#' @param nstart numeric. Random starts passed to \code{stats::kmeans}. Default 25.
#' @param iter.max numeric. Maximum iterations passed to \code{stats::kmeans}. Default 100.
#' @param seed numeric. Random seed for reproducible K-means. Default 1234; \code{NULL} skips seeding.
#' @param plot_envirotypes boolean. If \code{TRUE}, draws a barplot of the frequency (or ratio)
#'   of each envirotype per environment, faceted by variable (univariate) or by time interval
#'   (joint). Uses \pkg{ggplot2} when available, base graphics otherwise. Default \code{FALSE}.
#' @param min.window integer. Minimum number of observations required in the open-ended
#'   trailing interval. Default 2; use 0 or \code{NULL} to always keep it.
#' @param verbose boolean. If \code{TRUE} (default) prints progress messages.
#'
#' @return
#' If \code{envirotype_mining = FALSE} (default) and \code{plot_envirotypes = FALSE}: a
#' data.frame (or matrix, if \code{format = 'wide'}) of envirotype frequencies.
#'
#' Otherwise a \code{list} with:
#' \describe{
#'   \item{\code{typologies}}{data.frame (long) or matrix (wide) of envirotype frequencies.}
#'   \item{\code{envirotype_description}}{data.frame describing the numeric range captured by
#'     each envirotype: variable, label, cluster index, \code{n}, \code{min}, \code{max},
#'     \code{mean}, \code{sd}, \code{median} and the \code{lower}/\code{upper} cut boundaries.}
#'   \item{\code{mining_summary}}{data.frame with the optimal \code{k}, the CH index at the
#'     optimum and the full CH profile over \code{k.range}.}
#'   \item{\code{feature_combinations}}{only when \code{combine.features = TRUE}: a
#'     data.frame mapping each generic combination id to its number of features
#'     (\code{n_features}), the \code{features} themselves and the time \code{intervals}.}
#'   \item{\code{plot}}{the plot object, or \code{NULL}.}
#' }
#'
#' @details
#' The Calinski-Harabasz index for a partition into \eqn{k} groups of \eqn{n} observations is
#' \deqn{CH(k) = \frac{BGSS/(k-1)}{WGSS/(n-k)}}
#' where BGSS is the between-group and WGSS the within-group sum of squares. The \eqn{k}
#' maximising CH is retained. Clusters are relabelled in increasing order of their centres so
#' that \code{Envirotype_<var>_00001} is always the lowest range of the variable. For joint
#' (multivariate) mining there is no natural ordering, so centroids are ordered along their
#' first principal component and variables are standardised within the window before
#' clustering so that no feature dominates the Euclidean distance.
#'
#' @examples
#' \donttest{
#' data("maizeWTH")
#' env.data <- maizeWTH[maizeWTH$daysFromStart < 100, ]
#'
#' ## 1. Generic time intervals
#' env_typing(env.data = env.data, env.id = 'env', var.id = 'T2M', by.interval = TRUE)
#'
#' ## 2. Specific time intervals with names
#' env_typing(env.data = env.data, env.id = 'env', var.id = 'T2M',
#'            by.interval = TRUE,
#'            time.window  = c(0, 15, 35, 65, 90),
#'            names.window = c('1-initial growing', '2-leaf expansion I',
#'                             '3-leaf expansion II', '4-flowering', '5-grain filling'))
#'
#' ## 3. Cardinal (ecophysiological) thresholds
#' env_typing(env.data = env.data, env.id = 'env',
#'            var.id    = c('T2M', 'PRECTOT'),
#'            cardinals = list(T2M = c(0, 9, 22, 32, 45), PRECTOT = c(0, 5, 10, 25, 100)))
#'
#' ## 4. Data-driven univariate mining (K-means + Calinski-Harabasz)
#' out <- env_typing(env.data = env.data, env.id = 'env',
#'                   var.id = c('T2M', 'PRECTOT'),
#'                   envirotype_mining = TRUE, k.range = 2:8)
#' out$typologies
#' out$envirotype_description
#' out$mining_summary
#'
#' ## 5. Joint (multi-variable) envirotypes per development stage
#' outj <- env_typing(env.data = env.data, env.id = 'env',
#'                    var.id = c('T2M', 'PRECTOT'),
#'                    by.interval = TRUE, time.window = c(0, 30, 60, 90),
#'                    envirotype_mining = TRUE, joint = TRUE)
#' head(outj$envirotype_description)
#'
#' ## 6. Relative frequencies, wide format (ready for env_kernel)
#' ET <- env_typing(env.data = env.data, env.id = 'env', var.id = 'T2M',
#'                  format = 'wide', ratio = TRUE)
#'
#' ## 7. Envirotypes for feature combinations (singles + pairs by default)
#' outc <- env_typing(env.data = env.data, env.id = 'env',
#'                    var.id = c('T2M', 'PRECTOT', 'VPD'),
#'                    envirotype_mining = TRUE, combine.features = TRUE)
#' outc$feature_combinations        # generic id -> which variables
#' head(outc$envirotype_description)
#'
#' ## include triples too (combine.max.order = 3), or NULL for all orders
#' outc3 <- env_typing(env.data = env.data, env.id = 'env',
#'                     var.id = c('T2M', 'PRECTOT', 'VPD'),
#'                     envirotype_mining = TRUE, combine.features = TRUE,
#'                     combine.max.order = 3)
#' }
#'
#' @seealso \code{W_matrix}, \code{env_kernel}, \code{\link{T_matrix}}
#'
#' @importFrom stats quantile kmeans median sd ave as.formula prcomp
#' @importFrom utils combn
#' @importFrom grDevices hcl.colors
#' @importFrom graphics barplot par
#' @importFrom reshape2 acast dcast melt
#' @export
env_typing <- function(env.data, var.id, env.id, cardinals = NULL, days.id = NULL,
                       time.window = NULL, names.window = NULL, quantiles = NULL,
                       id.names = NULL, by.interval = FALSE, scale = FALSE,
                       format = NULL, ratio = FALSE, joint = FALSE,
                       joint.by.interval = TRUE, combine.features = FALSE,
                       combine.max.order = 2L, envirotype_mining = FALSE,
                       k.range = 2:10, nstart = 25, iter.max = 100, seed = 1234,
                       plot_envirotypes = FALSE, min.window = 2L, verbose = TRUE) {

  .et_banner("env_typing", "mines environmental types from weather data", verbose)
  if (is.null(format)) format <- "long"
  format <- match.arg(format, c("long", "wide"))

  if (isTRUE(combine.features)) {
    if (!isTRUE(envirotype_mining))
      stop("combine.features = TRUE requires envirotype_mining = TRUE.",
           call. = FALSE)
    if (isTRUE(joint))
      stop("combine.features = TRUE is only valid when joint = FALSE.",
           call. = FALSE)
  }

  env.data <- as.data.frame(env.data)
  if (isTRUE(scale))
    env.data[, var.id] <- scale(env.data[, var.id], center = TRUE, scale = TRUE)

  .et_step("assigning records to time-interval stages", verbose)
  .GET <- .env_melt_wth(.GeTw = env.data, days = days.id, by.interval = by.interval,
                        time.window = time.window, names.window = names.window,
                        id.names = id.names, var.id = var.id, env.id = env.id,
                        min.window = min.window)

  if (isFALSE(by.interval)) .GET <- data.frame(.GET, interval = "by environment")

  env  <- as.character(unique(.GET$env))
  int  <- as.character(unique(.GET$interval))
  vars <- as.character(unique(.GET$variable))

  ## ================================================================
  ## BRANCH 1 - data-driven envirotype mining
  ## ================================================================
  if (isTRUE(envirotype_mining)) {
    if (!is.null(cardinals))
      message("'cardinals' is ignored when envirotype_mining = TRUE")

    .GET$env.variable <- NA_character_

    if (isTRUE(combine.features))
      return(.env_typing_combine(.GET = .GET, vars = vars, env = env, int = int,
                                 by.interval = by.interval, k.range = k.range,
                                 nstart = nstart, iter.max = iter.max, seed = seed,
                                 ratio = ratio, format = format,
                                 plot_envirotypes = plot_envirotypes,
                                 verbose = verbose, max.order = combine.max.order))

    if (isTRUE(joint))
      return(.env_typing_joint(.GET = .GET, vars = vars, env = env, int = int,
                               k.range = k.range, nstart = nstart, iter.max = iter.max,
                               seed = seed, joint.by.interval = joint.by.interval,
                               ratio = ratio, format = format,
                               plot_envirotypes = plot_envirotypes, verbose = verbose))

    return(.env_typing_univariate(.GET = .GET, vars = vars, env = env, int = int,
                                  by.interval = by.interval, k.range = k.range,
                                  nstart = nstart, iter.max = iter.max, seed = seed,
                                  ratio = ratio, format = format,
                                  plot_envirotypes = plot_envirotypes, verbose = verbose))
  }

  ## ================================================================
  ## BRANCH 2 - cardinals / quantiles
  ## ================================================================
  if (is.null(quantiles)) quantiles <- c(.01, .25, .50, .75, .99)

  if (is.null(cardinals)) {
    cardinals <- vector("list", length(vars))
    names(cardinals) <- vars
  }
  if (is.null(names(cardinals))) names(cardinals) <- vars

  for (i in seq_along(vars)) {
    v <- vars[i]
    if (is.null(cardinals[[v]])) {
      SM <- round(as.numeric(stats::quantile(.GET$value[.GET$variable %in% v],
                                             quantiles, na.rm = TRUE)), 3)
      cardinals[[v]] <- SM[!duplicated(SM)]
    }
  }

  results <- do.call(rbind, lapply(seq_along(vars), function(s) {
    do.call(rbind, lapply(seq_along(env), function(j) {
      do.call(rbind, lapply(seq_along(int), function(t) {
        idx <- which(.GET$env == env[j] & .GET$variable == vars[s] &
                       .GET$interval == int[t])
        data.frame(table(cut(x = .GET$value[idx],
                             breaks = cardinals[[vars[s]]], right = TRUE)),
                   env = env[j], interval = int[t], var = vars[s],
                   stringsAsFactors = FALSE)
      }))
    }))
  }))

  names(results)[1] <- "env.variable"
  results$env.variable <- paste0(results$var, "_", results$env.variable)
  if (isTRUE(by.interval))
    results$env.variable <- paste0(results$env.variable, "_", results$interval)

  if (isTRUE(ratio)) results <- .env_freq2ratio(results)

  .p <- NULL
  if (isTRUE(plot_envirotypes))
    .p <- .env_plot_envirotypes(results, .split = "var", .ratio = ratio)

  if (format == "wide") {
    results <- reshape2::acast(results, env ~ env.variable, stats::median,
                               value.var = "Freq")
    results[is.na(results)] <- 0
  }

  if (isTRUE(plot_envirotypes)) return(list(typologies = results, plot = .p))
  results
}

#' Convert raw counts to relative frequencies within env x interval x variable
#'
#' @param .res long data.frame of envirotype frequencies.
#' @return The same data.frame with \code{Freq} rescaled to \eqn{[0,1]}.
#' @keywords internal
#' @noRd
.env_freq2ratio <- function(.res) {
  grp <- paste(.res$env, .res$interval, .res$var, sep = "\r")
  tot <- stats::ave(.res$Freq, grp, FUN = function(z) sum(z, na.rm = TRUE))
  .res$Freq <- as.numeric(ifelse(tot > 0, .res$Freq / tot, 0))
  .res
}

#' Barplot of envirotype frequency per environment
#'
#' @param .res long data.frame of envirotype frequencies.
#' @param .split character. Faceting column (\code{'var'} or \code{'interval'}).
#' @param .ratio boolean. Whether frequencies are relative.
#' @return A ggplot object (invisibly) or \code{NULL}.
#' @keywords internal
#' @noRd
.env_plot_envirotypes <- function(.res, .split = "var", .ratio = FALSE) {

  .res <- .res[!is.na(.res$env.variable), ]
  if (!nrow(.res)) return(invisible(NULL))
  ylab <- if (isTRUE(.ratio)) "Relative frequency" else "Frequency (days)"

  if (!.split %in% names(.res) || length(unique(.res[[.split]])) < 1) .split <- NULL

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    p <- ggplot2::ggplot(.res,
                         ggplot2::aes(x = .data[["env"]], y = .data[["Freq"]],
                                      fill = .data[["env.variable"]])) +
      ggplot2::geom_bar(stat = "identity", position = "stack", colour = NA) +
      ggplot2::labs(x = "Environment", y = ylab, fill = "Envirotype") +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                     legend.position = "right")
    if (!is.null(.split))
      p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", .split)),
                                   scales = "free_y")
    print(p)
    return(invisible(p))
  }

  # ---- base-graphics fallback ----
  sp <- if (is.null(.split)) rep("all", nrow(.res)) else as.character(.res[[.split]])
  lv <- unique(sp)
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  graphics::par(mfrow = c(1, max(1, length(lv))), mar = c(7, 4, 3, 1))
  for (l in lv) {
    sub <- .res[sp == l, ]
    m <- tapply(sub$Freq, list(sub$env.variable, sub$env), sum)
    m[is.na(m)] <- 0
    graphics::barplot(m, beside = FALSE, las = 2, main = l, ylab = ylab,
                      col = grDevices::hcl.colors(nrow(m), "Zissou 1"))
  }
  invisible(NULL)
}

#' Calinski-Harabasz index for a K-means partition
#'
#' @param .km a \code{kmeans} object.
#' @param .x the clustered data matrix.
#' @return Numeric CH index, or \code{NA}.
#' @keywords internal
#' @noRd
.env_ch_index <- function(.km, .x) {
  n <- nrow(as.matrix(.x))
  k <- length(.km$size)
  if (k < 2 || k >= n) return(NA_real_)
  WGSS <- sum(.km$withinss)
  BGSS <- .km$totss - WGSS
  if (!is.finite(WGSS) || WGSS <= .Machine$double.eps) return(NA_real_)
  (BGSS / (k - 1)) / (WGSS / (n - k))
}

#' Mine envirotypes by K-means with k chosen by Calinski-Harabasz
#'
#' @param .x numeric vector (one variable) or matrix/data.frame (joint variables).
#' @param .k.range integer vector of candidate k.
#' @param .nstart,.iter.max K-means control.
#' @param .seed integer or NULL.
#' @param .scale.joint boolean. Standardise columns before joint clustering.
#' @return A list with cluster assignment, retained rows, optimal k, centres and CH profile.
#' @keywords internal
#' @noRd
.env_mine_envirotypes <- function(.x, .k.range = 2:10, .nstart = 25, .iter.max = 100,
                                  .seed = NULL, .scale.joint = TRUE) {

  .X <- as.matrix(as.data.frame(.x))
  storage.mode(.X) <- "numeric"
  p <- ncol(.X)

  keep <- which(apply(.X, 1, function(r) all(is.finite(r))))
  Xo   <- .X[keep, , drop = FALSE]

  # Put variables on a common scale for multivariate mining, otherwise the
  # Euclidean distance is dominated by the largest-variance feature.
  ctr <- rep(0, p); scl <- rep(1, p)
  if (isTRUE(.scale.joint) && p > 1 && nrow(Xo) > 1) {
    ctr <- colMeans(Xo)
    scl <- apply(Xo, 2, stats::sd)
    scl[!is.finite(scl) | scl < .Machine$double.eps] <- 1
  }
  Xs <- sweep(sweep(Xo, 2, ctr, "-"), 2, scl, "/")

  nuniq <- nrow(unique(as.data.frame(Xs)))

  if (nuniq < 2 || nrow(Xo) < 2) {
    ctrs <- matrix(colMeans(Xo), nrow = 1, dimnames = list(NULL, colnames(.X)))
    return(list(cluster = rep(1L, nrow(Xo)), X = Xo, keep = keep, k = 1L,
                centers = ctrs, CH = NA_real_,
                CH.profile = data.frame(k = NA_integer_, CH = NA_real_)))
  }

  kr <- .k.range[.k.range >= 2 & .k.range < nuniq]
  if (!length(kr)) kr <- 2:max(2, min(nuniq - 1, 3))

  if (!is.null(.seed)) set.seed(.seed)

  CH.prof <- data.frame(k = kr, CH = NA_real_)
  fits <- vector("list", length(kr))
  for (a in seq_along(kr)) {
    fit <- try(stats::kmeans(x = Xs, centers = kr[a], nstart = .nstart,
                             iter.max = .iter.max), silent = TRUE)
    if (inherits(fit, "try-error")) next
    fits[[a]] <- fit
    CH.prof$CH[a] <- .env_ch_index(fit, Xs)
  }

  if (all(is.na(CH.prof$CH))) {
    ctrs <- matrix(colMeans(Xo), nrow = 1, dimnames = list(NULL, colnames(.X)))
    return(list(cluster = rep(1L, nrow(Xo)), X = Xo, keep = keep, k = 1L,
                centers = ctrs, CH = NA_real_, CH.profile = CH.prof))
  }

  best <- which.max(CH.prof$CH)
  kmb  <- fits[[best]]

  # Deterministic, interpretable labels: order clusters by their centre
  # (univariate) or along the dominant gradient (first PC, multivariate).
  Cs <- kmb$centers
  score <- if (p == 1) {
    as.numeric(Cs[, 1])
  } else {
    pc <- try(stats::prcomp(Cs, center = TRUE, scale. = FALSE), silent = TRUE)
    if (inherits(pc, "try-error")) rowSums(Cs) else pc$x[, 1]
  }
  ord   <- order(score)
  remap <- match(kmb$cluster, ord)

  ctrs <- sweep(sweep(Cs[ord, , drop = FALSE], 2, scl, "*"), 2, ctr, "+")
  colnames(ctrs) <- colnames(.X)

  list(cluster = remap, X = Xo, keep = keep, k = CH.prof$k[best],
       centers = ctrs, CH = CH.prof$CH[best], CH.profile = CH.prof)
}

#' Univariate envirotype mining branch of env_typing()
#'
#' @keywords internal
#' @noRd
.env_typing_univariate <- function(.GET, vars, env, int, by.interval, k.range,
                                   nstart, iter.max, seed, ratio, format,
                                   plot_envirotypes, verbose) {

  description  <- vector("list", length(vars))
  summary.mine <- vector("list", length(vars))

  for (s in seq_along(vars)) {
    if (isTRUE(verbose)) message("Mining envirotypes for: ", vars[s])

    idx <- which(.GET$variable == vars[s])
    xv  <- as.numeric(.GET$value[idx])

    mine <- .env_mine_envirotypes(.x = xv, .k.range = k.range, .nstart = nstart,
                                  .iter.max = iter.max, .seed = seed)

    .GET$env.variable[idx[mine$keep]] <-
      paste0("Envirotype_", vars[s], "_", formatC(mine$cluster, width = 5, flag = "0"))

    dsc <- do.call(rbind, lapply(sort(unique(mine$cluster)), function(cc) {
      vals <- mine$X[mine$cluster == cc, 1]
      data.frame(variable   = vars[s],
                 envirotype = paste0("Envirotype_", vars[s], "_",
                                     formatC(cc, width = 5, flag = "0")),
                 cluster = cc, n = length(vals),
                 min = min(vals, na.rm = TRUE), max = max(vals, na.rm = TRUE),
                 mean = mean(vals, na.rm = TRUE), sd = stats::sd(vals),
                 median = stats::median(vals, na.rm = TRUE),
                 center = mine$centers[cc, 1], stringsAsFactors = FALSE)
    }))

    dsc <- dsc[order(dsc$cluster), ]
    dsc$lower <- c(-Inf, (dsc$max[-nrow(dsc)] + dsc$min[-1]) / 2)
    dsc$upper <- c((dsc$max[-nrow(dsc)] + dsc$min[-1]) / 2, Inf)
    rownames(dsc) <- NULL
    description[[s]] <- dsc

    summary.mine[[s]] <- data.frame(
      variable = vars[s], k.optimal = mine$k, CH.index = mine$CH,
      k.evaluated = paste(mine$CH.profile$k, collapse = ","),
      CH.profile  = paste(round(mine$CH.profile$CH, 3), collapse = ","),
      stringsAsFactors = FALSE)
  }

  description  <- do.call(rbind, description)
  summary.mine <- do.call(rbind, summary.mine)

  .use <- .GET[!is.na(.GET$env.variable), ]
  results <- do.call(rbind, lapply(seq_along(vars), function(s) {
    lv <- description$envirotype[description$variable == vars[s]]
    do.call(rbind, lapply(seq_along(env), function(j) {
      do.call(rbind, lapply(seq_along(int), function(t) {
        sel <- .use$env.variable[.use$env == env[j] &
                                   .use$variable == vars[s] &
                                   .use$interval == int[t]]
        data.frame(table(factor(sel, levels = lv)),
                   env = env[j], interval = int[t], var = vars[s],
                   stringsAsFactors = FALSE)
      }))
    }))
  }))

  names(results)[1] <- "env.variable"
  results$env.variable <- as.character(results$env.variable)
  if (isTRUE(by.interval))
    results$env.variable <- paste0(results$env.variable, "_", results$interval)

  if (isTRUE(ratio)) results <- .env_freq2ratio(results)

  .p <- NULL
  if (isTRUE(plot_envirotypes))
    .p <- .env_plot_envirotypes(results, .split = "var", .ratio = ratio)

  if (format == "wide") {
    results <- reshape2::acast(results, env ~ env.variable, stats::median,
                               value.var = "Freq")
    results[is.na(results)] <- 0
  }

  list(typologies = results, envirotype_description = description,
       mining_summary = summary.mine, plot = .p)
}

#' Joint (multi-variable) envirotype mining branch of env_typing()
#'
#' @keywords internal
#' @noRd
.env_typing_joint <- function(.GET, vars, env, int, k.range, nstart, iter.max,
                              seed, joint.by.interval, ratio, format,
                              plot_envirotypes, verbose) {

  .GET$.row <- paste(.GET$env, .GET$interval, .GET$daysFromStart, sep = "\r")
  W <- reshape2::dcast(.GET, .row + env + interval ~ variable,
                       value.var = "value", fun.aggregate = mean)
  jvars <- vars[vars %in% names(W)]
  if (length(jvars) < 2)
    warning("joint = TRUE with fewer than 2 variables; result is univariate.",
            call. = FALSE)

  wins <- if (isTRUE(joint.by.interval)) int else "ALL"
  description  <- vector("list", length(wins))
  summary.mine <- vector("list", length(wins))
  W$env.variable <- NA_character_

  for (w in seq_along(wins)) {
    widx <- if (isTRUE(joint.by.interval)) which(W$interval == wins[w]) else seq_len(nrow(W))
    if (!length(widx)) next
    if (isTRUE(verbose)) message("Mining joint envirotypes for window: ", wins[w])

    mine <- .env_mine_envirotypes(.x = W[widx, jvars, drop = FALSE],
                                  .k.range = k.range, .nstart = nstart,
                                  .iter.max = iter.max, .seed = seed,
                                  .scale.joint = TRUE)

    stem <- if (isTRUE(joint.by.interval)) paste0("Envirotype_", wins[w]) else
      paste0("Envirotype_", paste(jvars, collapse = "."))

    W$env.variable[widx[mine$keep]] <-
      paste0(stem, "_", formatC(mine$cluster, width = 5, flag = "0"))

    description[[w]] <- do.call(rbind, lapply(sort(unique(mine$cluster)), function(cc) {
      rows <- mine$X[mine$cluster == cc, , drop = FALSE]
      do.call(rbind, lapply(seq_along(jvars), function(q) {
        vals <- rows[, q]
        data.frame(window = wins[w],
                   envirotype = paste0(stem, "_", formatC(cc, width = 5, flag = "0")),
                   cluster = cc, variable = jvars[q], n = length(vals),
                   min = min(vals, na.rm = TRUE), max = max(vals, na.rm = TRUE),
                   mean = mean(vals, na.rm = TRUE), sd = stats::sd(vals),
                   median = stats::median(vals, na.rm = TRUE),
                   center = mine$centers[cc, q], stringsAsFactors = FALSE)
      }))
    }))

    summary.mine[[w]] <- data.frame(
      window = wins[w], variables = paste(jvars, collapse = ","),
      k.optimal = mine$k, CH.index = mine$CH,
      k.evaluated = paste(mine$CH.profile$k, collapse = ","),
      CH.profile  = paste(round(mine$CH.profile$CH, 3), collapse = ","),
      stringsAsFactors = FALSE)
  }

  description  <- do.call(rbind, description)
  summary.mine <- do.call(rbind, summary.mine)
  rownames(description) <- NULL

  .use <- W[!is.na(W$env.variable), ]
  results <- do.call(rbind, lapply(seq_along(env), function(j) {
    do.call(rbind, lapply(seq_along(int), function(t) {
      lv <- unique(description$envirotype[
        if (isTRUE(joint.by.interval)) description$window == int[t] else TRUE])
      if (!length(lv)) return(NULL)
      sel <- .use$env.variable[.use$env == env[j] & .use$interval == int[t]]
      data.frame(table(factor(sel, levels = lv)),
                 env = env[j], interval = int[t],
                 var = paste(jvars, collapse = "."), stringsAsFactors = FALSE)
    }))
  }))

  names(results)[1] <- "env.variable"
  results$env.variable <- as.character(results$env.variable)

  if (isTRUE(ratio)) results <- .env_freq2ratio(results)

  # joint envirotypes already span every variable -> split by INTERVAL
  .p <- NULL
  if (isTRUE(plot_envirotypes))
    .p <- .env_plot_envirotypes(results, .split = "interval", .ratio = ratio)

  if (format == "wide") {
    results <- reshape2::acast(results, env ~ env.variable, stats::median,
                               value.var = "Freq")
    results[is.na(results)] <- 0
  }

  list(typologies = results, envirotype_description = description,
       mining_summary = summary.mine, plot = .p)
}

#' Feature-combination envirotype mining branch of env_typing()
#'
#' Mines envirotypes for every unique combination of the variables up to
#' \code{max.order} (each single variable, every pair, every triple, ...). Each
#' combination of size >= 2 is clustered jointly; singletons reduce to univariate
#' clustering. Combinations receive compact generic ids (\code{FC001}, ...) decoded
#' by the returned \code{feature_combinations} table.
#'
#' @keywords internal
#' @noRd
.env_typing_combine <- function(.GET, vars, env, int, by.interval, k.range,
                                nstart, iter.max, seed, ratio, format,
                                plot_envirotypes, verbose, max.order) {

  vars <- unique(vars)                       # never combine a variable with itself
  n <- length(vars)
  max.order <- if (is.null(max.order)) n else min(as.integer(max.order), n)
  if (is.na(max.order) || max.order < 1L)
    stop("combine.max.order must be a positive integer.", call. = FALSE)

  combos <- unlist(lapply(seq_len(max.order),
                          function(m) utils::combn(vars, m, simplify = FALSE)),
                   recursive = FALSE)
  ## defensive: keep only unique feature sets (order-independent, no repeats)
  combos <- combos[!duplicated(lapply(combos, sort))]
  n_combos <- length(combos)
  if (n_combos > 5000L)
    stop(n_combos, " feature combinations requested. Enumerating the power set of ",
         n, " variables is impractical; lower `combine.max.order` (currently ",
         max.order, ") or use fewer variables in `var.id`.", call. = FALSE)

  ids <- paste0("FC", formatC(seq_len(n_combos),
                              width = max(3L, nchar(n_combos)), flag = "0"))
  feature_combinations <- data.frame(
    combo_id   = ids,
    n_features = vapply(combos, length, integer(1)),
    features   = vapply(combos, function(z) paste(z, collapse = ", "),
                        character(1)),
    intervals  = paste(int, collapse = ", "),
    stringsAsFactors = FALSE)

  ## one row per env x interval x timepoint, variables spread into columns.
  ## the within-group sequence keys the pivot without needing a day column.
  .GET$.k <- stats::ave(seq_len(nrow(.GET)),
                        paste(.GET$env, .GET$interval, .GET$variable, sep = "\r"),
                        FUN = seq_along)
  W <- reshape2::dcast(.GET, .k + env + interval ~ variable,
                       value.var = "value", fun.aggregate = mean)

  typ <- dsc <- smy <- vector("list", n_combos)

  for (a in seq_len(n_combos)) {
    id <- ids[a]
    jv <- combos[[a]][combos[[a]] %in% names(W)]
    if (isTRUE(verbose))
      message("Mining combination ", id, " (", paste(jv, collapse = "+"), ")")

    mine <- .env_mine_envirotypes(.x = W[, jv, drop = FALSE], .k.range = k.range,
                                  .nstart = nstart, .iter.max = iter.max,
                                  .seed = seed, .scale.joint = TRUE)

    lab <- rep(NA_character_, nrow(W))
    lab[mine$keep] <- paste0("Envirotype_", id, "_",
                             formatC(mine$cluster, width = 5, flag = "0"))

    dsc[[a]] <- do.call(rbind, lapply(sort(unique(mine$cluster)), function(g) {
      rows <- mine$X[mine$cluster == g, , drop = FALSE]
      do.call(rbind, lapply(seq_along(jv), function(q) {
        vals <- rows[, q]
        data.frame(combo_id = id, features = paste(jv, collapse = ", "),
                   envirotype = paste0("Envirotype_", id, "_",
                                       formatC(g, width = 5, flag = "0")),
                   cluster = g, variable = jv[q], n = length(vals),
                   min = min(vals, na.rm = TRUE), max = max(vals, na.rm = TRUE),
                   mean = mean(vals, na.rm = TRUE), sd = stats::sd(vals),
                   median = stats::median(vals, na.rm = TRUE),
                   center = mine$centers[g, q], stringsAsFactors = FALSE)
      }))
    }))

    smy[[a]] <- data.frame(
      combo_id = id, features = paste(jv, collapse = ", "),
      n_features = length(jv), k.optimal = mine$k, CH.index = mine$CH,
      k.evaluated = paste(mine$CH.profile$k, collapse = ","),
      CH.profile  = paste(round(mine$CH.profile$CH, 3), collapse = ","),
      stringsAsFactors = FALSE)

    lv <- sort(unique(stats::na.omit(lab)))
    Wc <- data.frame(env = W$env, interval = W$interval, lab = lab,
                     stringsAsFactors = FALSE)
    typ[[a]] <- do.call(rbind, lapply(seq_along(env), function(j) {
      do.call(rbind, lapply(seq_along(int), function(t) {
        sel <- Wc$lab[Wc$env == env[j] & Wc$interval == int[t]]
        data.frame(table(factor(sel, levels = lv)),
                   env = env[j], interval = int[t], var = id,
                   stringsAsFactors = FALSE)
      }))
    }))
  }

  results      <- do.call(rbind, typ)
  description  <- do.call(rbind, dsc)
  summary.mine <- do.call(rbind, smy)
  rownames(description) <- NULL

  names(results)[1] <- "env.variable"
  results$env.variable <- as.character(results$env.variable)
  if (isTRUE(by.interval))
    results$env.variable <- paste0(results$env.variable, "_", results$interval)

  if (isTRUE(ratio)) results <- .env_freq2ratio(results)

  .p <- NULL
  if (isTRUE(plot_envirotypes))
    .p <- .env_plot_envirotypes(results, .split = "var", .ratio = ratio)

  if (format == "wide") {
    results <- reshape2::acast(results, env ~ env.variable, stats::median,
                               value.var = "Freq")
    results[is.na(results)] <- 0
  }

  list(typologies = results, envirotype_description = description,
       mining_summary = summary.mine,
       feature_combinations = feature_combinations, plot = .p)
}


# =========================================================================
# SECTION 4 - env_cluster()
# =========================================================================

#' @title Envirotype clustering into mega-environments
#'
#' @description
#' Clusters environments into homogeneous groups ("mega-environments") from an environmental
#' covariable matrix (\code{W}) using K-means. The number of groups can be set by the user or
#' selected automatically by maximising the average silhouette width. An environmental
#' relationship (kinship) matrix is also returned for convenience.
#'
#' @author Germano Costa Neto
#'
#' @param W matrix or data.frame. Environmental covariables (environments in rows, variables in
#'   columns), typically a \code{W_matrix} output. An envirotype frequency matrix from
#'   \code{\link{T_matrix}} works identically.
#' @param k integer or NULL. Number of mega-environments. If \code{NULL} (default) it is chosen
#'   by silhouette.
#' @param k.max integer. Maximum number of clusters tested when \code{k} is \code{NULL}. Default 10.
#' @param scale boolean. If \code{TRUE} (default) \code{W} is centred and scaled before clustering.
#' @param nstart integer. Number of random starts for K-means. Default 25.
#' @param seed integer or NULL. Optional seed for reproducibility.
#' @param verbose boolean. If \code{TRUE} (default) prints progress messages.
#'
#' @return
#' A list with \code{clusters} (named integer vector), \code{k}, \code{silhouette} (data.frame
#' of average silhouette per tested \code{k}, or \code{NULL} when \code{k} was supplied),
#' \code{kmeans} (the fitted object) and \code{env.kinship} (an environment x environment
#' correlation matrix).
#'
#' @details
#' Mega-environments are repeatable groups of locations that rank genotypes similarly. Here they
#' are approximated by grouping environments that share a similar enviromic profile, which is
#' useful to stratify a Target Population of Environments (TPE) before multi-environment
#' prediction. The silhouette width \eqn{s(i)} contrasts within- and between-cluster distances;
#' the \code{k} maximising its average is retained.
#'
#' @examples
#' \donttest{
#' data("maizeWTH")
#'
#' W <- W_matrix(env.data = maizeWTH[maizeWTH$daysFromStart < 100, ],
#'               var.id = c("T2M", "T2M_MAX", "PRECTOT"), statistic = "mean")
#'
#' ## Automatic number of mega-environments (chosen by silhouette)
#' mega <- env_cluster(W, seed = 1)
#' mega$k
#' mega$clusters
#' mega$silhouette
#'
#' ## Force a fixed number of groups
#' mega2 <- env_cluster(W, k = 2, seed = 1)
#' mega2$clusters
#'
#' ## Environmental relationship (kinship) matrix
#' round(mega$env.kinship, 2)
#'
#' ## Use the grouping downstream (e.g. as a fixed effect or to stratify CV)
#' table(mega$clusters)
#' }
#'
#' @seealso \code{W_matrix}, \code{env_kernel}, \code{env_target_importance},
#'   \code{\link{T_matrix}}
#'
#' @importFrom stats kmeans dist cor setNames
#' @export
env_cluster <- function(W, k = NULL, k.max = 10, scale = TRUE,
                        nstart = 25, seed = NULL, verbose = TRUE) {

  .et_banner("env_cluster", "groups environments into mega-environments", verbose)
  if (!is.null(seed)) set.seed(seed)

  W <- .env_as_matrix(W, "W")
  if (anyNA(W))
    stop("W contains missing values; use W_matrix(impute = ) or remove them first.",
         call. = FALSE)
  if (nrow(W) < 3)
    stop("At least three environments are required to build mega-environments.",
         call. = FALSE)

  Ws <- if (isTRUE(scale)) scale(W) else W

  sil_tab <- NULL
  if (is.null(k)) {
    .et_step("selecting the number of clusters by silhouette", verbose)
    best    <- .env_best_k(Ws, k.max = k.max, nstart = nstart, verbose = verbose,
                           label = "mega-environments")
    k       <- best$k
    sil_tab <- best$silhouette
  }

  .et_step("running k-means clustering", verbose)
  km <- stats::kmeans(Ws, centers = k, nstart = nstart)

  list(clusters    = stats::setNames(km$cluster, rownames(W)),
       k           = k,
       silhouette  = sil_tab,
       kmeans      = km,
       env.kinship = stats::cor(t(Ws)))
}


# =========================================================================
# SECTION 5 - env_target_importance()
# =========================================================================

#' @title Environmental drivers of phenotypic/target clustering
#'
#' @description
#' Groups environments by their target metrics (\code{Y}: phenotypic means, heritabilities,
#' selection accuracies, etc.) using K-means, then trains a Random Forest to predict those
#' clusters from the environmental covariable matrix (\code{W}). The Random Forest variable
#' importance ranks which environmental variables best explain the target-based grouping of
#' environments. The number of clusters can be supplied by the user or optimised automatically
#' with the average silhouette width.
#'
#' @author Germano Costa Neto
#'
#' @param Y numeric vector, matrix or data.frame. Target metrics with environments in rows
#'   (e.g. adjusted means, \eqn{H^2}, predictive ability). A named vector is also accepted.
#' @param W matrix or data.frame. Environmental covariables with environments in rows and
#'   variables in columns (e.g. a \code{W_matrix} output).
#' @param k integer or NULL. Number of clusters. If \code{NULL} (default) it is chosen by
#'   maximising the average silhouette width over \code{2:k.max}.
#' @param k.max integer. Maximum number of clusters tested when \code{k} is \code{NULL}. Default 10.
#' @param scale.Y boolean. If \code{TRUE} (default) \code{Y} is centred and scaled before clustering.
#' @param ntree integer. Number of trees in the Random Forest. Default 1000.
#' @param importance.type character. \code{"gini"} (MeanDecreaseGini, default) or
#'   \code{"permutation"} (MeanDecreaseAccuracy).
#' @param nstart integer. Number of random starts for K-means. Default 25.
#' @param seed integer or NULL. Optional seed for reproducibility.
#' @param save boolean. If \code{TRUE}, writes the importance table and cluster assignment as \code{.csv}.
#' @param dir.path character. Output directory when \code{save = TRUE}. Defaults to the working directory.
#' @param verbose boolean. If \code{TRUE} (default) prints progress messages.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{clusters}: named integer vector of the cluster of each environment.
#'   \item \code{k}: the number of clusters used.
#'   \item \code{silhouette}: data.frame of average silhouette width per tested \code{k}
#'     (\code{NULL} if \code{k} was given).
#'   \item \code{kmeans}: the \code{kmeans} object.
#'   \item \code{randomForest}: the fitted \code{randomForest} object.
#'   \item \code{importance}: data.frame of environmental variables ranked by importance (descending).
#' }
#'
#' @details
#' This is a supervised counterpart of \code{env_cluster}: instead of grouping
#' environments by their enviromic profile, environments are grouped by \emph{what the breeder
#' cares about} (yield level, heritability, accuracy) and the enviromic data are then asked to
#' explain that grouping. Variables at the top of \code{importance} are candidate drivers of the
#' observed environmental stratification and good candidates for \code{env_threshold_gxe}.
#'
#' Permutation importance (\code{"permutation"}) is generally preferred when covariables differ
#' in scale or in number of categories, since MeanDecreaseGini is biased towards high-cardinality
#' predictors.
#'
#' @references
#' Costa-Neto, G., da Matta, D., Fernandes, I. K., & Heinemann, A. B. (2023).
#' Environmental clusters defining breeding zones for tropical irrigated rice in Brazil.
#' \emph{Agronomy Journal}, 115(5). \doi{10.1002/agj2.21481}
#'
#' @examples
#' \donttest{
#' if (requireNamespace("randomForest", quietly = TRUE)) {
#'   data("maizeYield"); data("maizeWTH")
#'
#'   ## Environment-level target: mean grain yield per environment
#'   Y <- tapply(maizeYield$value, maizeYield$env, mean)
#'
#'   ## Environmental covariables (environments x covariables)
#'   W <- W_matrix(env.data = maizeWTH[maizeWTH$daysFromStart < 100, ],
#'                 var.id = c("T2M", "T2M_MAX", "PRECTOT"), statistic = "mean")
#'
#'   ## Which environmental variables best explain the target-based grouping?
#'   out <- env_target_importance(Y = Y, W = W, seed = 1)
#'   out$k
#'   out$clusters
#'   head(out$importance)
#'
#'   ## Permutation importance and a fixed number of clusters
#'   out2 <- env_target_importance(Y = Y, W = W, k = 2,
#'                                 importance.type = "permutation", seed = 1)
#'   head(out2$importance)
#'
#'   ## Multi-trait target (e.g. mean and heritability per environment)
#'   Ymat <- cbind(mean = Y, sd = tapply(maizeYield$value, maizeYield$env, stats::sd))
#'   out3 <- env_target_importance(Y = Ymat, W = W, seed = 1)
#'   out3$importance
#' }
#' }
#'
#' @seealso \code{env_cluster}, \code{W_matrix}
#'
#' @importFrom stats kmeans dist setNames
#' @importFrom utils write.csv
#' @export
env_target_importance <- function(Y, W, k = NULL, k.max = 10, scale.Y = TRUE,
                                  ntree = 1000,
                                  importance.type = c("gini", "permutation"),
                                  nstart = 25, seed = NULL, save = FALSE,
                                  dir.path = NULL, verbose = TRUE) {

  .et_banner("env_target_importance",
             "ranks covariables by importance to a target", verbose)
  importance.type <- match.arg(importance.type)
  if (!requireNamespace("randomForest", quietly = TRUE))
    stop("Package 'randomForest' is required. ",
         "Install it with install.packages('randomForest').", call. = FALSE)

  if (!is.null(seed)) set.seed(seed)

  # ---- shape inputs -------------------------------------------------------
  if (is.vector(Y) && is.null(dim(Y)))
    Y <- matrix(Y, ncol = 1, dimnames = list(names(Y), "target"))
  Y <- .env_as_matrix(Y, "Y")
  W <- .env_as_matrix(W, "W")
  if (anyNA(Y)) stop("Y contains missing values; impute or remove them first.", call. = FALSE)

  # ---- align environments between Y and W ---------------------------------
  if (!is.null(rownames(Y)) && !is.null(rownames(W))) {
    common <- intersect(rownames(Y), rownames(W))
    if (length(common) < 2)
      stop("Fewer than two environments shared between Y and W row names.", call. = FALSE)
    if (isTRUE(verbose) && length(common) < nrow(Y))
      message("Using ", length(common), " environments shared by Y and W.")
    Y <- Y[common, , drop = FALSE]
    W <- W[common, , drop = FALSE]
  } else if (nrow(Y) != nrow(W)) {
    stop("Y and W must have the same number of environments (rows) or matching row names.",
         call. = FALSE)
  }
  if (anyNA(W))
    stop("W contains missing values; use W_matrix(impute = ) or remove them first.",
         call. = FALSE)

  Ys <- if (isTRUE(scale.Y)) scale(Y) else Y

  # ---- number of clusters via silhouette ----------------------------------
  sil_tab <- NULL
  if (is.null(k)) {
    best    <- .env_best_k(Ys, k.max = k.max, nstart = nstart, verbose = verbose,
                           label = "target clusters")
    k       <- best$k
    sil_tab <- best$silhouette
  }

  # ---- K-means on the target metrics --------------------------------------
  .et_step("clustering environments on the target metrics", verbose)
  km <- stats::kmeans(Ys, centers = k, nstart = nstart)
  clusters <- stats::setNames(km$cluster, rownames(W))

  # ---- Random Forest: clusters ~ W ----------------------------------------
  .et_step("fitting random forest of clusters on covariables", verbose)
  rf_data <- data.frame(.cluster = factor(clusters), W, check.names = TRUE)
  rf <- randomForest::randomForest(.cluster ~ ., data = rf_data,
                                   ntree = ntree, importance = TRUE)

  imp <- randomForest::importance(rf)
  imp.col <- if (importance.type == "permutation" &&
                 "MeanDecreaseAccuracy" %in% colnames(imp))
    "MeanDecreaseAccuracy" else "MeanDecreaseGini"

  importance <- data.frame(variable = rownames(imp),
                           importance = imp[, imp.col],
                           row.names = NULL, stringsAsFactors = FALSE)
  importance <- importance[order(-importance$importance), ]
  rownames(importance) <- NULL

  if (isTRUE(save)) {
    if (is.null(dir.path)) dir.path <- getwd()
    utils::write.csv(importance, file.path(dir.path, "env_target_importance.csv"),
                     row.names = FALSE)
    utils::write.csv(data.frame(env = names(clusters), cluster = clusters),
                     file.path(dir.path, "env_target_clusters.csv"), row.names = FALSE)
    if (isTRUE(verbose)) message("Results written to ", normalizePath(dir.path))
  }

  list(clusters     = clusters,
       k            = k,
       silhouette   = sil_tab,
       kmeans       = km,
       randomForest = rf,
       importance   = importance)
}


# =========================================================================
# SECTION - thermal-time phenology internal helpers
# =========================================================================

#' Built-in thermal-time stage templates
#'
#' A registry of cardinal temperatures and cumulative thermal-time thresholds
#' (C d) at the END of each named developmental stage, for the field crops most
#' often handled in multi-environment trials.
#'
#' Each entry carries:
#' \describe{
#'   \item{Tbase, Tupper}{cardinal temperatures for thermal time (C)}
#'   \item{stages}{named, strictly increasing cumulative thermal time (C d)}
#'   \item{scale}{the staging system the labels belong to}
#'   \item{day.type}{"short", "long" or "neutral" photoperiod response}
#'   \item{p.opt}{optimum daylength (h) for the photoperiod factor}
#'   \item{label}{common name}
#'   \item{source}{where the thresholds come from}
#'   \item{note}{the main caveat for that crop}
#' }
#'
#' READ THIS BEFORE TRUSTING A TEMPLATE. These are extension-service or
#' textbook tabulations for a generic mid-season cultivar in a well-studied
#' production region. They are STARTING POINTS for a consistent relative
#' scale, not calibrated phenology predictions. Maturity group alone moves
#' maize R6 by 300-500 C d and soybean R8 by more than that; wheat and barley
#' templates ignore vernalisation entirely, which dominates development in
#' winter types. Calibrate `stages` against observed dates before drawing
#' agronomic conclusions. Use show_phenology() to inspect and compare.
#'
#' @keywords internal
#' @noRd
.PHEN_TEMPLATES <- list(

  ## ---------------- C4 cereals ----------------
  maize = list(
    Tbase = 10, Tupper = 30, scale = "Abendroth V/R",
    day.type = "short", p.opt = 12.5, label = "Maize (Zea mays)",
    stages = c(VE = 125, V6 = 475, V12 = 740, VT = 1135,
               R1 = 1250, R3 = 1660, R5 = 2110, R6 = 2700),
    source = "Abendroth et al. (2011) PMR 1009, Iowa State Univ. Extension",
    note = paste("Tuned to a ~110-115 day US Corn Belt hybrid. Short-season",
                 "hybrids reach R6 near 2400 C d, full-season near 3000.")
  ),
  sorghum = list(
    Tbase = 10, Tupper = 34, scale = "Vanderlip GS",
    day.type = "short", p.opt = 12.5, label = "Sorghum (Sorghum bicolor)",
    stages = c(EM = 130, GP3 = 450, FI = 800, BT = 1150,
               ANT = 1350, SD = 1750, PM = 2300),
    source = "Vanderlip (1993) How a Sorghum Plant Develops, KSU S-3",
    note = paste("Photoperiod-sensitive tropical landraces will not flower on",
                 "these thresholds at long daylengths; set photoperiod = TRUE.")
  ),
  pearl_millet = list(
    Tbase = 10, Tupper = 34, scale = "Maiti GS",
    day.type = "short", p.opt = 12.0, label = "Pearl millet (Pennisetum glaucum)",
    stages = c(EM = 100, GP3 = 350, PI = 600, BT = 850,
               ANT = 1000, GF = 1350, PM = 1700),
    source = "Maiti & Bidinger (1981) ICRISAT Research Bulletin 6",
    note = "Short-duration crop; Sahelian landraces are strongly photoperiodic."
  ),
  sugarcane = list(
    Tbase = 12, Tupper = 32, scale = "phenophase",
    day.type = "neutral", p.opt = 12.5, label = "Sugarcane (Saccharum spp.)",
    stages = c(GERM = 400, TILL = 1600, GRAND = 4500, MAT = 6500),
    source = "Inman-Bamber (1994) Field Crops Res. 36, 41-51",
    note = paste("Plant-cane cycle of 12-18 months; ratoon cycles are shorter.",
                 "Thresholds are far more site-dependent than for annuals.")
  ),

  ## ---------------- C3 small grains ----------------
  wheat = list(
    Tbase = 0, Tupper = 26, scale = "Zadoks-anchored",
    day.type = "long", p.opt = 16.0, label = "Wheat (Triticum aestivum)",
    stages = c(EM = 150, TS = 500, SE = 900, HD = 1250,
               ANT = 1400, GF = 2000, PM = 2400),
    source = "Zadoks et al. (1974) Weed Res. 14, 415-421; McMaster & Wilhelm (1997)",
    note = paste("NO VERNALISATION TERM. Winter wheat sown in autumn will be",
                 "staged far too fast. Valid for spring types, or recalibrate.")
  ),
  barley = list(
    Tbase = 0, Tupper = 26, scale = "Zadoks-anchored",
    day.type = "long", p.opt = 16.0, label = "Barley (Hordeum vulgare)",
    stages = c(EM = 130, TS = 450, SE = 800, HD = 1100,
               ANT = 1220, GF = 1750, PM = 2100),
    source = "Zadoks et al. (1974); Bauer et al. (1984) Agron. J. 76, 829-835",
    note = "Earlier than wheat throughout. Same vernalisation caveat."
  ),
  oat = list(
    Tbase = 0, Tupper = 25, scale = "Zadoks-anchored",
    day.type = "long", p.opt = 16.0, label = "Oat (Avena sativa)",
    stages = c(EM = 140, TS = 480, SE = 850, HD = 1180,
               ANT = 1300, GF = 1850, PM = 2250),
    source = "Sonego et al. (2000) NZ J. Crop Hortic. Sci. 28, 177-186",
    note = "Intermediate between barley and wheat; sensitive to heat at anthesis."
  ),
  rice = list(
    Tbase = 10, Tupper = 30, scale = "Counce P/R",
    day.type = "short", p.opt = 12.5, label = "Rice (Oryza sativa)",
    stages = c(EM = 100, TIL = 450, PI = 900, HD = 1350,
               ANT = 1450, GF = 1900, PM = 2300),
    source = "Counce et al. (2000) Crop Sci. 40, 436-443; IRRI RKB",
    note = paste("Indica lowland defaults. Tbase 10 C understates cold injury",
                 "in temperate japonica, where 12-15 C is more appropriate.")
  ),

  ## ---------------- Grain legumes ----------------
  soybean = list(
    Tbase = 10, Tupper = 30, scale = "Fehr-Caviness V/R",
    day.type = "short", p.opt = 12.8, label = "Soybean (Glycine max)",
    stages = c(VE = 120, V3 = 400, R1 = 800, R3 = 1150,
               R5 = 1500, R6 = 1900, R8 = 2400),
    source = "Fehr & Caviness (1977) Iowa State Univ. Spec. Rep. 80",
    note = paste("STRONGLY photoperiodic. Maturity groups 000-X differ by more",
                 "than 1000 C d to R8; set photoperiod = TRUE across latitudes.")
  ),
  common_bean = list(
    Tbase = 10, Tupper = 30, scale = "CIAT V/R",
    day.type = "short", p.opt = 12.5, label = "Common bean (Phaseolus vulgaris)",
    stages = c(V1 = 130, V3 = 350, R5 = 600, R6 = 750,
               R7 = 950, R8 = 1250, R9 = 1500),
    source = "van Schoonhoven & Pastor-Corrales (1987) CIAT standard system",
    note = "Short-cycle crop; determinate and indeterminate habits differ widely."
  ),
  chickpea = list(
    Tbase = 0, Tupper = 30, scale = "phenophase",
    day.type = "long", p.opt = 13.0, label = "Chickpea (Cicer arietinum)",
    stages = c(EM = 150, VEG = 600, FLI = 900, POD = 1200,
               SF = 1550, PM = 1900),
    source = "Soltani & Sinclair (2011) Field Crops Res. 124, 252-260",
    note = "Tbase 0 C reflects cool-season adaptation; desi and kabuli differ."
  ),
  groundnut = list(
    Tbase = 10, Tupper = 32, scale = "Boote R",
    day.type = "neutral", p.opt = 12.5, label = "Groundnut (Arachis hypogaea)",
    stages = c(VE = 130, R1 = 500, R3 = 800, R5 = 1150,
               R6 = 1500, R7 = 1900, R8 = 2300),
    source = "Boote (1982) Peanut Science 9, 35-40",
    note = "Spanish types mature several hundred C d before Virginia types."
  ),
  cowpea = list(
    Tbase = 11, Tupper = 32, scale = "phenophase",
    day.type = "short", p.opt = 12.5, label = "Cowpea (Vigna unguiculata)",
    stages = c(EM = 100, VEG = 350, FLI = 600, POD = 850,
               SF = 1150, PM = 1450),
    source = "Craufurd et al. (1997) Exp. Agric. 33, 267-277",
    note = "Heat-tolerant, short-cycle; many landraces are photoperiodic."
  ),

  ## ---------------- Oilseeds / fibre / roots ----------------
  canola = list(
    Tbase = 5, Tupper = 27, scale = "BBCH-anchored",
    day.type = "long", p.opt = 16.0, label = "Canola / oilseed rape (Brassica napus)",
    stages = c(EM = 120, ROS = 400, STE = 700, BUD = 950,
               FLW = 1250, POD = 1700, PM = 2100),
    source = "Lancashire et al. (1991) Ann. Appl. Biol. 119, 561-601 (BBCH)",
    note = "Spring-type defaults; winter types require vernalisation."
  ),
  sunflower = list(
    Tbase = 6, Tupper = 30, scale = "Schneiter-Miller V/R",
    day.type = "neutral", p.opt = 12.5, label = "Sunflower (Helianthus annuus)",
    stages = c(VE = 90, V8 = 400, R1 = 700, R5 = 1050,
               R6 = 1300, R8 = 1650, R9 = 2000),
    source = "Schneiter & Miller (1981) Crop Sci. 21, 901-903",
    note = "Largely day-neutral, which makes it well behaved across latitudes."
  ),
  cotton = list(
    Tbase = 15.6, Tupper = 32, scale = "nodes/phenophase",
    day.type = "neutral", p.opt = 12.5, label = "Cotton (Gossypium hirsutum)",
    stages = c(EM = 60, SQ = 450, FB = 800, PB = 1150,
               OB = 1600, CO = 2100),
    source = "Oosterhuis (1990) Arkansas Agric. Exp. Stn. Spec. Rep. 145",
    note = paste("Tbase 15.6 C (60 F) is the US industry convention; DD60 heat",
                 "units are this base. Indeterminate growth blurs late stages.")
  ),
  potato = list(
    Tbase = 7, Tupper = 30, scale = "phenophase",
    day.type = "long", p.opt = 14.0, label = "Potato (Solanum tuberosum)",
    stages = c(EM = 200, VEG = 500, TI = 700, TB = 1400,
               MAT = 1900),
    source = "Kooman & Haverkort (1995) in Potato Ecology, Kluwer",
    note = paste("Tuber initiation is promoted by SHORT days in many cultivars,",
                 "opposite to the vegetative long-day response; treat with care.")
  ),
  cassava = list(
    Tbase = 13, Tupper = 30, scale = "phenophase",
    day.type = "neutral", p.opt = 12.5, label = "Cassava (Manihot esculenta)",
    stages = c(EM = 250, BR1 = 1200, BR2 = 2600, BULK = 4500, MAT = 6000),
    source = "Alves (2002) in Cassava: Biology, Production and Utilization, CABI",
    note = "Harvested 8-24 months after planting; thresholds are indicative only."
  ),

  safflower = list(
    Tbase = 5, Tupper = 30, scale = "BBCH-anchored",
    day.type = "long", p.opt = 14.0, label = "Safflower (Carthamus tinctorius)",
    stages = c(EM = 95, ROS = 450, STE = 800, BUD = 1100,
               FL = 1400, SF = 1800, PM = 2200),
    source = paste("Demir (2026) Turk. J. Field Crops 31(2) doi:10.17557/tjfc.1998413;",
                   "Flemmer et al. (2015) Ann. Appl. Biol. 166, 331-339 (BBCH scale)"),
    note = paste("CARDINALS ARE FIELD-VERIFIED: Tbase 5 C / Tupper 30 C from a",
                 "12-year x 5-sowing-date trial; emergence 92.7-96.3 C d was stable",
                 "across sowing dates. Stages after EM are BBCH-anchored ESTIMATES",
                 "and are NOT independently verified -- calibrate before relying",
                 "on them.")
  ),

  ## ---------------- Forage ----------------
  alfalfa = list(
    Tbase = 5, Tupper = 30, scale = "Kalu-Fick MSC",
    day.type = "long", p.opt = 14.0, label = "Alfalfa / lucerne (Medicago sativa)",
    stages = c(EV = 150, LV = 350, EB = 550, LB = 700,
               EF = 850, LF = 1000),
    source = "Kalu & Fick (1981) Crop Sci. 21, 267-271",
    note = "Per REGROWTH CYCLE after each cut, not from sowing."
  )
)

#' Daily thermal time by one of three methods
#'
#' @param tmin,tmax numeric vectors of daily min/max air temperature (C).
#' @param tbase numeric. Base temperature below which development stops.
#' @param tupper numeric. Upper threshold. Interpretation depends on `method`.
#' @param method character. "gdd" caps Tmax at tupper (the standard McMaster &
#'   Wilhelm "Method 1" with an upper cutoff); "baskerville" uses the
#'   Baskerville-Emin sine-wave integration, which is better when the daily
#'   range straddles tbase; "trapezoid" applies a linear decline above tupper
#'   toward an absolute maximum tmaxdev.
#' @param tmaxdev numeric. Absolute temperature at which development ceases,
#'   used only by method = "trapezoid".
#' @return numeric vector of daily thermal time (C d), >= 0.
#' @keywords internal
#' @noRd
.phen_dtt <- function(tmin, tmax, tbase, tupper,
                      method = c("gdd", "baskerville", "trapezoid"),
                      tmaxdev = 45) {
  method <- match.arg(method)

  if (method == "gdd") {
    tmx <- pmin(tmax, tupper)
    tmn <- pmin(pmax(tmin, tbase), tupper)
    tmx <- pmax(tmx, tbase)
    return(pmax((tmx + tmn) / 2 - tbase, 0))
  }

  if (method == "trapezoid") {
    tmean <- (tmin + tmax) / 2
    out <- rep(0, length(tmean))
    ok <- !is.na(tmean)
    rise <- ok & tmean > tbase & tmean <= tupper
    fall <- ok & tmean > tupper & tmean < tmaxdev
    out[rise] <- tmean[rise] - tbase
    out[fall] <- (tupper - tbase) *
      (tmaxdev - tmean[fall]) / (tmaxdev - tupper)
    out[!ok] <- NA_real_
    return(pmax(out, 0))
  }

  # Baskerville-Emin (1969) single-sine with horizontal upper cutoff.
  # Integrates a sine wave fitted through tmin/tmax above tbase, which avoids
  # the bias of the simple mean when tmin < tbase < tmax.
  n <- length(tmin)
  out <- rep(NA_real_, n)
  ok <- !is.na(tmin) & !is.na(tmax)
  tmx <- pmin(tmax[ok], tupper)
  tmn <- pmin(tmin[ok], tupper)
  amp <- (tmx - tmn) / 2
  mid <- (tmx + tmn) / 2
  res <- numeric(length(mid))

  below <- tmx <= tbase
  above <- tmn >= tbase
  strad <- !below & !above

  res[below] <- 0
  res[above] <- mid[above] - tbase
  if (any(strad)) {
    a <- amp[strad]; m <- mid[strad]
    a[a <= 0] <- 1e-8
    theta <- asin(pmin(pmax((tbase - m) / a, -1), 1))
    res[strad] <- ((m - tbase) * (pi / 2 - theta) + a * cos(theta)) / pi
  }
  out[ok] <- pmax(res, 0)
  out
}

#' Photoperiod factor for photothermal time
#'
#' Linear response between a critical and an optimum daylength. For
#' short-day crops (maize, rice, sorghum, soybean) development slows as
#' daylength exceeds `popt`; for long-day crops (wheat) it slows below `popt`.
#'
#' @param daylength numeric. Hours of daylight.
#' @param popt numeric. Daylength at which development is unimpeded (h).
#' @param psens numeric. Sensitivity, C d per hour of deviation, expressed as
#'   a fraction reduction per hour. 0 disables the response.
#' @param type character. "short" or "long" day.
#' @return numeric multiplier in [0, 1].
#' @keywords internal
#' @noRd
.phen_photo <- function(daylength, popt = 12.5, psens = 0.02,
                        type = c("short", "long")) {
  type <- match.arg(type)
  if (is.null(daylength) || psens <= 0) return(rep(1, length(daylength)))
  dev <- if (type == "short") daylength - popt else popt - daylength
  f <- 1 - psens * pmax(dev, 0)
  pmin(pmax(f, 0), 1)
}

#' Daylength (h) from latitude and day of year - CBM model
#'
#' Used only when the weather data carry no daylength column. The package's
#' own .ra_daylength() returns extraterrestrial radiation alongside N; this
#' lighter version avoids a circular dependency on param_radiation().
#'
#' @references Forsythe, W.C. et al. (1995) A model comparison for daylength
#'   as a function of latitude and day of year. Ecol. Modell. 80, 87-95.
#' @keywords internal
#' @noRd
.phen_daylength <- function(doy, lat) {
  p <- 0.8333
  theta <- 0.2163108 + 2 * atan(0.9671396 * tan(0.00860 * (doy - 186)))
  phi <- asin(0.39795 * cos(theta))
  x <- (sin(p * pi / 180) + sin(lat * pi / 180) * sin(phi)) /
    (cos(lat * pi / 180) * cos(phi))
  24 - (24 / pi) * acos(pmin(pmax(x, -1), 1))
}

#' Resolve a stage specification into a named, increasing numeric vector
#' @keywords internal
#' @noRd
.phen_resolve_stages <- function(stages, crop) {
  if (is.null(stages)) {
    if (is.null(crop))
      stop("Supply either 'stages' or 'crop'. Available crop templates: ",
           paste(names(.PHEN_TEMPLATES), collapse = ", "), ".", call. = FALSE)
    crop <- match.arg(tolower(crop), names(.PHEN_TEMPLATES))
    return(.PHEN_TEMPLATES[[crop]]$stages)
  }

  ## ---- cardinals-only mode -------------------------------------------
  ## stages = <single integer N> cuts the season into N equal-thermal-time
  ## intervals labelled TT1..TTN. Use when published cardinal temperatures
  ## exist for a crop but no calibrated stage thresholds do (the common case
  ## outside the major cereals -- e.g. the 117 crops in Paredes et al. 2025).
  ##
  ## This still fixes the calendar-day misalignment that motivates this
  ## module, because TT bins are thermal not chronological. It simply makes
  ## NO CLAIM about where any named stage falls. Inventing thresholds would
  ## produce a 'stage' column that looks authoritative while silently
  ## misaligning every downstream W-matrix column.
  if (is.numeric(stages) && length(stages) == 1L && is.null(names(stages))) {
    n <- stages
    if (is.na(n) || n < 2 || n != as.integer(n))
      stop("When 'stages' is a single number it is the NUMBER of equal ",
           "thermal-time intervals and must be a whole number >= 2. ",
           "Got: ", n, ".", call. = FALSE)
    n <- as.integer(n)
    out <- seq_len(n) / n
    names(out) <- paste0("TT", seq_len(n))
    attr(out, "relative") <- TRUE
    return(out)
  }

  if (!is.numeric(stages) || is.null(names(stages)) || anyNA(names(stages)))
    stop("'stages' must be a NAMED numeric vector of cumulative thermal time ",
         "at the end of each stage, e.g. c(VE = 125, V6 = 475, VT = 1135).",
         call. = FALSE)
  if (is.unsorted(stages, strictly = TRUE))
    stop("'stages' must be strictly increasing cumulative thermal time.",
         call. = FALSE)
  if (any(stages <= 0))
    stop("'stages' must be positive.", call. = FALSE)
  stages
}


#' @title List Built-in Phenology Templates
#'
#' @description
#' \code{phenology_templates()} returns the built-in crop templates as a
#' data.frame; \code{show_phenology()} prints them, optionally for one crop
#' with its stage table, provenance and caveats.
#'
#' Every template is a published tabulation for a REFERENCE cultivar. They
#' give a consistent relative scale across sites -- which is what fixes the
#' calendar-day misalignment in \code{W_matrix()} -- but they are not
#' calibrated predictions for your germplasm. Maturity group alone shifts
#' maize R6 by 300-500 C d. Validate against observed flowering dates before
#' treating stage labels as phenology.
#'
#' @param crop character. Optional single crop; if omitted, all are listed.
#'
#' @return \code{phenology_templates()} a data.frame, one row per crop.
#'   \code{show_phenology()} returns its input invisibly and is called for
#'   the side effect of printing.
#'
#' @examples
#' \dontrun{
#' phenology_templates()
#' show_phenology("safflower")
#' }
#'
#' @export
phenology_templates <- function() {
  do.call(rbind, lapply(names(.PHEN_TEMPLATES), function(k) {
    t <- .PHEN_TEMPLATES[[k]]
    data.frame(crop = k, label = t$label, Tbase = t$Tbase, Tupper = t$Tupper,
               day.type = t$day.type, p.opt = t$p.opt, n_stages = length(t$stages),
               final_stage = names(t$stages)[length(t$stages)],
               final_TT = unname(t$stages[length(t$stages)]),
               scale = t$scale, stringsAsFactors = FALSE)
  }))
}

#' @rdname phenology_templates
#' @export
show_phenology <- function(crop = NULL) {
  if (is.null(crop)) {
    d <- phenology_templates()
    .et_banner("phenology templates",
               paste0(nrow(d), " built-in crops"), TRUE)
    print(d[, c("crop", "Tbase", "Tupper", "day.type", "n_stages",
                "final_stage", "final_TT")], row.names = FALSE)
    cat("\nshow_phenology(\"<crop>\") for stage table, source and caveats.\n")
    cat("stages = <N> gives N equal thermal-time bins with no named stages.\n")
    return(invisible(d))
  }
  crop <- match.arg(tolower(crop), names(.PHEN_TEMPLATES))
  t <- .PHEN_TEMPLATES[[crop]]
  bar <- strrep("-", 63)
  cat(bar, "\n", t$label, "\n", bar, "\n", sep = "")
  cat("  Tbase ", t$Tbase, " C   Tupper ", t$Tupper, " C   ",
      t$day.type, "-day (p.opt ", t$p.opt, " h)\n", sep = "")
  cat("  scale: ", t$scale, "\n\n", sep = "")
  prev <- 0
  for (i in seq_along(t$stages)) {
    cat(sprintf("    %-6s <= %7.0f C d   (+%5.0f)\n",
                names(t$stages)[i], t$stages[i], t$stages[i] - prev))
    prev <- t$stages[i]
  }
  cat("\n  source: ", t$source, "\n", sep = "")
  cat("  NOTE:   ", t$note, "\n", sep = "")
  cat(bar, "\n", sep = "")
  invisible(t)
}


# =========================================================================
# SECTION - env_phenology()
# =========================================================================

#' @title Thermal-Time Phenology: Map Daily Weather onto Developmental Stages
#'
#' @description
#' Accumulates daily thermal time from planting and cuts each environment's
#' season at crop-specific developmental thresholds, so that every downstream
#' interval-based routine compares LIKE DEVELOPMENTAL STATES across sites
#' rather than like calendar dates.
#'
#' The returned \code{stage} column is a drop-in replacement for the
#' \code{interval} produced by \code{time.window} in \code{W_matrix()},
#' \code{summaryWTH()} and \code{env_typing()}.
#'
#' @author Germano Costa Neto
#'
#' @param env.data data.frame. Daily weather, typically a \code{get_weather()}
#'   or \code{processWTH()} output.
#' @param env.id character. Environment-id column. Default \code{"env"}.
#' @param day.id character. Day-counter column. Default \code{"daysFromStart"}.
#' @param Tmax,Tmin character. Daily maximum/minimum temperature columns.
#'   Defaults \code{"T2M_MAX"} / \code{"T2M_MIN"}.
#' @param crop character. Name of a built-in template, e.g. \code{"maize"},
#'   \code{"wheat"}, \code{"soybean"}, \code{"rice"}, \code{"cotton"},
#'   \code{"canola"}. Supplies default \code{stages}, \code{Tbase},
#'   \code{Tupper}, \code{p.opt} and \code{day.type}. Any explicitly supplied
#'   argument wins. Call \code{show_phenology()} for the full list of the 20
#'   templates and their provenance.
#' @param stages named numeric. Cumulative thermal time (C d) at the END of
#'   each stage, strictly increasing, e.g.
#'   \code{c(VE = 125, V6 = 475, VT = 1135, R6 = 2700)}. Overrides \code{crop}.
#' @param Tbase,Tupper numeric. Cardinal temperatures for thermal time.
#'   Override the crop template when supplied.
#' @param method character. Thermal-time method: \code{"gdd"} (default,
#'   McMaster & Wilhelm with upper cutoff), \code{"baskerville"}
#'   (Baskerville-Emin sine integration), or \code{"trapezoid"} (linear
#'   supra-optimal decline to \code{Tmaxdev}).
#' @param Tmaxdev numeric. Absolute developmental maximum, used by
#'   \code{method = "trapezoid"}. Default 45.
#' @param photoperiod boolean. If \code{TRUE}, multiply daily thermal time by a
#'   photoperiod factor, producing photothermal time (PTT).
#' @param lat character or numeric. Latitude column name or a numeric vector,
#'   used to derive daylength when \code{photoperiod = TRUE} and no daylength
#'   column exists. Default \code{"LAT"}.
#' @param daylength.id character. Existing daylength column (h), e.g. \code{"N"}
#'   from \code{param_radiation()}. Used in preference to deriving it.
#' @param p.opt,p.sens numeric. Optimum daylength (h) and sensitivity
#'   (fractional slowdown per hour of deviation). When \code{p.opt} is
#'   \code{NULL} (default) the crop template's value is used -- 12.5 h for
#'   short-day maize, 16 h for long-day wheat, and so on.
#' @param day.type character. \code{"short"}, \code{"long"} or
#'   \code{"neutral"}. \code{NULL} (default) takes the crop template's value.
#'   A \code{"neutral"} crop silently keeps a photoperiod factor of 1 even
#'   when \code{photoperiod = TRUE}.
#' @param planting.id character. Optional column giving planting date, used to
#'   reset accumulation. If \code{NULL}, accumulation starts at the first
#'   record of each environment.
#' @param merge boolean. If \code{TRUE} (default) return \code{env.data} with
#'   the new columns appended; if \code{FALSE} return only the new columns.
#' @param verbose boolean. Print a progress banner.
#'
#' @return
#' A data.frame with these columns added:
#' \describe{
#'   \item{dTT}{daily thermal time (C d), or photothermal time if
#'     \code{photoperiod = TRUE}}
#'   \item{cumTT}{cumulative thermal time from planting}
#'   \item{stage}{ordered factor, the developmental stage of that day}
#'   \item{dstage}{relative position within the stage, in [0, 1)}
#'   \item{stage_day}{day counter within the stage (1-based)}
#' }
#' The attribute \code{"phenology"} carries the stage table, the method and the
#' cardinal temperatures. Environments whose season ends before the last
#' threshold get \code{NA} for the stages never reached, and a warning names
#' them -- silently recycling the final stage would fabricate development that
#' the weather does not support.
#'
#' @details
#' \strong{Which method.} \code{"gdd"} is the convention and what
#' \code{param_temperature()} already computes; use it for comparability.
#' \code{"baskerville"} is materially better where daily range straddles
#' \code{Tbase} -- high-elevation or early-spring sites -- because the simple
#' mean then understates accumulation. \code{"trapezoid"} is the only option
#' that penalises supra-optimal heat, relevant for heat-stress work.
#'
#' \strong{Photothermal time.} Maize, sorghum and rice are quantitative
#' short-day plants: long days delay flowering. With \code{photoperiod = TRUE}
#' the daily increment is scaled by a linear daylength factor, which keeps
#' tropical and temperate plantings on a common developmental scale. The
#' default \code{p.sens = 0.02} is deliberately mild; calibrate it before
#' relying on the absolute stage dates.
#'
#' \strong{Caveat.} These thresholds are cultivar-generic. Maturity group alone
#' moves maize R6 by 300-500 C d. Treat the stage boundaries as a consistent
#' relative scale, not as a validated phenology prediction, unless you have
#' calibrated \code{stages} against observed dates.
#'
#' @examples
#' \dontrun{
#' data("maizeWTH")
#'
#' ## 1. Default maize staging on cumulative GDD
#' ph <- env_phenology(maizeWTH, crop = "maize")
#' table(ph$env, ph$stage)
#'
#' ## 2. Photothermal time, daylength taken from param_radiation()'s N column
#' ph2 <- env_phenology(processWTH(maizeWTH), crop = "maize",
#'                      photoperiod = TRUE, daylength.id = "N")
#'
#' ## 3. Custom stages and the Baskerville-Emin method
#' ph3 <- env_phenology(maizeWTH,
#'                      stages = c(veg = 500, flow = 1200, fill = 2200),
#'                      Tbase = 8, Tupper = 32, method = "baskerville")
#'
#' ## 4. Hand the stage factor to the existing covariable machinery.
#' ##    This is the point of the function: W columns now align on
#' ##    development instead of on calendar date.
#' W <- W_matrix(env.data = ph, var.id = c("T2M", "PRECTOT"),
#'               by.interval = TRUE, time.window = NULL)
#' }
#'
#' @seealso \code{\link{param_temperature}}, \code{\link{W_matrix}},
#'   \code{\link{env_risk_profile}}, \code{\link{summaryWTH}}
#'
#' @references
#' McMaster, G.S. & Wilhelm, W.W. (1997) Growing degree-days: one equation, two
#' interpretations. \emph{Agric. For. Meteorol.} 87(4), 291-300.
#'
#' Baskerville, G.L. & Emin, P. (1969) Rapid estimation of heat accumulation
#' from maximum and minimum temperatures. \emph{Ecology} 50(3), 514-517.
#'
#' Abendroth, L.J. et al. (2011) \emph{Corn Growth and Development}. PMR 1009,
#' Iowa State University Extension.
#'
#' @export
env_phenology <- function(env.data, env.id = "env", day.id = "daysFromStart",
                          Tmax = "T2M_MAX", Tmin = "T2M_MIN",
                          crop = NULL, stages = NULL,
                          Tbase = NULL, Tupper = NULL,
                          method = c("gdd", "baskerville", "trapezoid"),
                          Tmaxdev = 45,
                          photoperiod = FALSE, lat = "LAT",
                          daylength.id = NULL, p.opt = NULL, p.sens = 0.02,
                          day.type = NULL,
                          planting.id = NULL,
                          merge = TRUE, verbose = TRUE) {

  .et_banner("env_phenology", "maps weather onto thermal-time stages", verbose)
  method   <- match.arg(method)
  env.data <- as.data.frame(env.data)

  ## ---- resolve stage table and cardinals --------------------------------
  tmpl <- if (!is.null(crop)) {
    .PHEN_TEMPLATES[[match.arg(tolower(crop), names(.PHEN_TEMPLATES))]]
  } else NULL

  stages <- .phen_resolve_stages(stages, crop)
  if (is.null(Tbase))  Tbase  <- if (!is.null(tmpl)) tmpl$Tbase  else 10
  if (is.null(Tupper)) Tupper <- if (!is.null(tmpl)) tmpl$Tupper else 30

  # Photoperiod defaults follow the crop template, not a global constant: a
  # long-day canola and a short-day soybean must not share p.opt = 12.5.
  if (is.null(p.opt))
    p.opt <- if (!is.null(tmpl)) tmpl$p.opt else 12.5
  if (is.null(day.type)) {
    day.type <- if (!is.null(tmpl)) tmpl$day.type else "short"
  }
  day.type <- match.arg(day.type, c("short", "long", "neutral"))
  if (isTRUE(photoperiod) && identical(day.type, "neutral")) {
    if (verbose)
      message("'", if (!is.null(crop)) crop else "this crop",
              "' is day-neutral in the template; photoperiod factor left at 1.")
    p.sens <- 0
    day.type <- "short"
  }

  if (!(Tbase < Tupper))
    stop("'Tbase' must be lower than 'Tupper'.", call. = FALSE)
  if (method == "trapezoid" && !(Tupper < Tmaxdev))
    stop("method = 'trapezoid' requires Tupper < Tmaxdev.", call. = FALSE)
  ## ---- column checks -----------------------------------------------------
  miss <- setdiff(c(env.id, Tmax, Tmin), names(env.data))
  if (length(miss))
    stop("env_phenology() needs column(s): ", paste(miss, collapse = ", "),
         call. = FALSE)
  if (!day.id %in% names(env.data)) {
    env.data[[day.id]] <- stats::ave(seq_len(nrow(env.data)),
                                     env.data[[env.id]],
                                     FUN = seq_along)
    if (verbose)
      message("No '", day.id, "' column; generated from row order within env.")
  }

  tmn <- .strip_power_na(env.data[[Tmin]])
  tmx <- .strip_power_na(env.data[[Tmax]])
  n_inv <- sum(tmx < tmn, na.rm = TRUE)
  if (n_inv)
    warning(n_inv, " day(s) with Tmax < Tmin; check the input columns.",
            call. = FALSE)

  ## ---- daily thermal time ------------------------------------------------
  .et_step(paste0("computing daily thermal time (method = ", method, ")"), verbose)
  dTT <- .phen_dtt(tmn, tmx, tbase = Tbase, tupper = Tupper,
                   method = method, tmaxdev = Tmaxdev)

  ## ---- optional photoperiod modulation ----------------------------------
  if (isTRUE(photoperiod)) {
    dl <- NULL
    if (!is.null(daylength.id) && daylength.id %in% names(env.data)) {
      dl <- as.numeric(env.data[[daylength.id]])
    } else {
      doy <- if ("DOY" %in% names(env.data)) {
        as.numeric(env.data[["DOY"]])
      } else if ("YYYYMMDD" %in% names(env.data)) {
        as.numeric(format(as.Date(env.data[["YYYYMMDD"]]), "%j"))
      } else NULL
      latv <- if (is.numeric(lat)) {
        rep(lat, length.out = nrow(env.data))
      } else if (lat %in% names(env.data)) {
        as.numeric(env.data[[lat]])
      } else NULL
      if (is.null(doy) || is.null(latv))
        stop("photoperiod = TRUE needs a daylength column ('daylength.id'), ",
             "or both a day-of-year ('DOY'/'YYYYMMDD') and a latitude column ",
             "('lat'). Run param_radiation() first to get 'N'.", call. = FALSE)
      dl <- .phen_daylength(doy, latv)
    }
    pf <- .phen_photo(dl, popt = p.opt, psens = p.sens, type = day.type)
    dTT <- dTT * pf
    .et_step("thermal time modulated by photoperiod (PTT)", verbose)
  }

  ## ---- cumulative accumulation, per environment --------------------------
  env_f <- as.character(env.data[[env.id]])
  ord   <- order(env_f, env.data[[day.id]])
  inv   <- order(ord)

  dTT_o   <- dTT[ord]
  env_o   <- env_f[ord]
  day_o   <- env.data[[day.id]][ord]

  # reset accumulation at planting if a planting column is supplied
  if (!is.null(planting.id)) {
    if (!planting.id %in% names(env.data))
      stop("'planting.id' column '", planting.id, "' not found.", call. = FALSE)
    pl_o <- as.Date(env.data[[planting.id]])[ord]
    dt_o <- if ("YYYYMMDD" %in% names(env.data)) {
      as.Date(env.data[["YYYYMMDD"]])[ord]
    } else NULL
    if (!is.null(dt_o)) dTT_o[!is.na(dt_o) & !is.na(pl_o) & dt_o < pl_o] <- 0
  }

  # NA days contribute nothing but must not poison the cumulative sum
  na_day <- is.na(dTT_o)
  if (any(na_day)) {
    warning(sum(na_day), " day(s) with missing temperature were treated as ",
            "zero thermal time; cumTT is therefore a lower bound for the ",
            "affected environments.", call. = FALSE)
    dTT_o[na_day] <- 0
  }
  cum_o <- stats::ave(dTT_o, env_o, FUN = cumsum)

  ## ---- cardinals-only mode: rescale relative cuts to observed span -------
  ## Relative stages arrive as fractions of the season. Anchor them to the
  ## MAXIMUM cumTT observed across environments so every site is cut on one
  ## common thermal ruler (per-site rescaling would re-introduce exactly the
  ## incomparability this module exists to remove).
  if (isTRUE(attr(stages, "relative"))) {
    span <- suppressWarnings(max(cum_o, na.rm = TRUE))
    if (!is.finite(span) || span <= 0)
      stop("Cannot apply cardinals-only staging: no finite thermal time was ",
           "accumulated. Check the temperature columns and Tbase.",
           call. = FALSE)
    rel   <- stages
    stages <- as.numeric(stages) * span
    names(stages) <- names(rel)
    .et_step(paste0("cardinals-only staging: ", length(stages),
                    " equal TT bins over ", round(span), " C d"), verbose)
  }

  ## ---- cut into stages ---------------------------------------------------
  .et_step("assigning developmental stages", verbose)
  brk  <- c(0, as.numeric(stages))
  labs <- names(stages)
  stage_o <- cut(cum_o, breaks = brk, labels = labs,
                 right = TRUE, include.lowest = TRUE)
  stage_o <- factor(stage_o, levels = labs, ordered = TRUE)

  # relative position within stage, in [0,1)
  lower <- c(0, as.numeric(stages)[-length(stages)])[as.integer(stage_o)]
  upper <- as.numeric(stages)[as.integer(stage_o)]
  dstage <- (cum_o - lower) / (upper - lower)
  dstage[!is.finite(dstage)] <- NA_real_

  key <- paste(env_o, as.integer(stage_o), sep = "\r")
  stage_day <- stats::ave(seq_along(key), key, FUN = seq_along)
  stage_day[is.na(stage_o)] <- NA_integer_

  ## ---- report environments that never reached maturity -------------------
  reached <- tapply(cum_o, env_o, max, na.rm = TRUE)
  short <- names(reached)[reached < max(as.numeric(stages))]
  if (length(short)) {
    warning(length(short), " environment(s) did not accumulate enough thermal ",
            "time to reach '", labs[length(labs)], "': ",
            paste(utils::head(short, 8), collapse = ", "),
            if (length(short) > 8) ", ..." else "",
            ". Days beyond the last reached stage are NA. Extend the series ",
            "or lower the stage thresholds.", call. = FALSE)
  }

  out <- data.frame(dTT = dTT_o[inv], cumTT = cum_o[inv],
                    stage = stage_o[inv], dstage = dstage[inv],
                    stage_day = stage_day[inv],
                    stringsAsFactors = FALSE)

  attr(out, "phenology") <- list(stages = stages, method = method,
                                 Tbase = Tbase, Tupper = Tupper,
                                 Tmaxdev = if (method == "trapezoid") Tmaxdev else NA_real_,
                                 photoperiod = photoperiod,
                                 p.opt = p.opt, p.sens = p.sens,
                                 day.type = day.type,
                                 crop = crop,
                                 scale = if (!is.null(tmpl)) tmpl$scale else "custom",
                                 source = if (!is.null(tmpl)) tmpl$source else "user-supplied")

  if (verbose) {
    message(strrep("-", 63))
    message("Thermal-time stages (", if (photoperiod) "PTT" else "GDD",
            ", Tbase = ", Tbase, " C, Tupper = ", Tupper, " C)")
    for (i in seq_along(stages))
      message(sprintf("  %-6s  <= %7.0f C d   %6d day(s)",
                      labs[i], stages[i], sum(out$stage == labs[i], na.rm = TRUE)))
    message(strrep("-", 63))
  }

  if (isTRUE(merge)) {
    res <- cbind(env.data, out)
    attr(res, "phenology") <- attr(out, "phenology")
    return(res)
  }
  out
}


# =========================================================================
# SECTION - env_risk_profile()
# =========================================================================

#' Classify a numeric vector into stress levels by cardinal breaks
#' @keywords internal
#' @noRd
.risk_classify <- function(x, breaks, labels) {
  cut(x, breaks = c(-Inf, breaks, Inf), labels = labels,
      right = TRUE, include.lowest = TRUE)
}

#' Default stress cardinals for common envirotyping variables
#'
#' Two-sided where the variable has both a cold/dry and a hot/wet failure mode.
#' These are maize-oriented defaults; override via `cardinals`.
#' @keywords internal
#' @noRd
.RISK_CARDINALS <- list(
  T2M_MAX  = list(breaks = c(30, 35),      labels = c("safe", "heat_mild", "heat_severe")),
  T2M_MIN  = list(breaks = c(8, 14),       labels = c("cold_severe", "cold_mild", "safe")),
  T2M      = list(breaks = c(18, 28),      labels = c("cool", "optimal", "hot")),
  PETP     = list(breaks = c(-5, 0),       labels = c("drought_severe", "drought_mild", "wet")),
  PRECTOT  = list(breaks = c(1, 10),       labels = c("dry", "moderate", "wet")),
  VPD      = list(breaks = c(1.5, 2.5),    labels = c("low", "moderate", "high")),
  FRUE     = list(breaks = c(0.5, 0.8),    labels = c("severe", "limiting", "non_limiting")),
  Ks       = list(breaks = c(0.4, 0.8),    labels = c("severe", "moderate", "none"))
)

#' @title Historical Envirotype Risk Profile Across Years
#'
#' @description
#' Quantifies how OFTEN each stress envirotype occurs at a site, by
#' developmental stage, across many historical seasons. Where
#' \code{env_typing()} describes the environments you actually observed,
#' \code{env_risk_profile()} estimates the underlying climatology those
#' environments were drawn from -- the frequency, not the realisation.
#'
#' This is the quantity a breeding programme needs to weight its testing
#' network: a mega-environment that appears in 1 season out of 20 should not
#' carry the same weight as one that appears in 12.
#'
#' @author Germano Costa Neto
#'
#' @param env.data data.frame. Multi-year daily weather for one or more sites,
#'   typically \code{get_weather()} run over a long window, then
#'   \code{processWTH()}. Must contain a site id, a date, and the variables in
#'   \code{var.id}.
#' @param site.id character. Site-id column. Default \code{"env"}. Note this is
#'   a SITE, not a site-year; the site-year is formed internally.
#' @param date.id character. Date column (\code{Date} or \code{YYYYMMDD}).
#'   Default \code{"YYYYMMDD"}.
#' @param var.id vector (character). Variables to profile. Defaults to those
#'   present among the names of \code{cardinals}.
#' @param years vector (numeric). Seasons to include. Defaults to every year
#'   present with a complete season.
#' @param planting character, \code{Date} vector or \code{"DOY"}. Planting date
#'   per season. A single \code{"MM-DD"} string applies the same sowing date
#'   every year; a named vector keyed by site applies a per-site date; a
#'   \code{Date} vector is used as given.
#' @param season.length numeric. Days after planting retained per season when
#'   phenology is not used. Default 150.
#' @param cardinals named list. Per-variable \code{list(breaks =, labels =)}
#'   defining stress classes. Defaults in \code{.RISK_CARDINALS} cover
#'   \code{T2M_MAX}, \code{T2M_MIN}, \code{T2M}, \code{PETP}, \code{PRECTOT},
#'   \code{VPD}, \code{FRUE} and \code{Ks}.
#' @param by.stage boolean. If \code{TRUE} (default) compute risk within
#'   thermal-time stages by calling \code{env_phenology()} per season.
#' @param crop,stages,Tbase,Tupper,method Passed to \code{env_phenology()} when
#'   \code{by.stage = TRUE}.
#' @param min.days integer. Minimum days a site-year must contribute to be kept.
#'   Default 30.
#' @param return.seasons boolean. If \code{TRUE}, also return the per-site-year
#'   table behind the frequencies.
#' @param verbose boolean. Print progress.
#'
#' @return
#' An object of class \code{"env_risk"}: a list with
#' \describe{
#'   \item{frequency}{data.frame of site x stage x variable x class with
#'     \code{n_years}, \code{n_occur}, \code{freq} (probability a season shows
#'     that class at all), and \code{mean_days} (expected days per season)}
#'   \item{severity}{data.frame of site x stage x variable with the mean,
#'     sd and 10th/90th percentile of the variable across seasons}
#'   \item{return_period}{data.frame with \code{1 / freq}, the expected number
#'     of seasons between occurrences of each stress class}
#'   \item{seasons}{per-site-year detail if \code{return.seasons = TRUE}}
#'   \item{call, years, cardinals}{provenance}
#' }
#'
#' @details
#' \strong{Two frequencies are reported and they answer different questions.}
#' \code{freq} is the proportion of SEASONS in which the class occurred at
#' least once -- the probability of encountering that stress. \code{mean_days}
#' is the expected NUMBER OF DAYS per season -- the intensity. A site can have
#' \code{freq = 1.0} for mild heat while averaging only 3 days of it; that is a
#' very different target from one averaging 40 days.
#'
#' \strong{Return period} is simply \code{1 / freq} and inherits all the usual
#' caveats: it is an average recurrence interval estimated from a finite
#' record, not a schedule. With 20 years of data a return period above ~10
#' seasons is barely distinguishable from noise. The function will not print a
#' return period longer than the record.
#'
#' \strong{Stationarity.} Frequencies computed over 1990-2024 describe that
#' period. If the site is warming, they are already an underestimate of future
#' heat risk. Use \code{project_risk()} with a scenario to see the shift.
#'
#' @examples
#' \dontrun{
#' ## 30 seasons of weather for three sites
#' sites <- data.frame(env = c("PIRA", "SETE", "NMAC"),
#'                     lat = c(-22.7, -19.4, -21.2),
#'                     lon = c(-47.6, -44.2, -45.2))
#'
#' wth <- get_weather(env.id = sites$env, lat = sites$lat, lon = sites$lon,
#'                    start.day = "1994-01-01", end.day = "2024-12-31")
#' wth <- processWTH(wth)
#'
#' ## Risk of heat and drought at flowering, sown 15 October each year
#' rp <- env_risk_profile(wth, site.id = "env", planting = "10-15",
#'                        var.id = c("T2M_MAX", "PETP"),
#'                        crop = "maize", by.stage = TRUE)
#' rp
#'
#' ## How often is flowering heat-stressed at each site?
#' subset(rp$frequency, stage == "R1" & class == "heat_severe")
#'
#' ## Weights for a target population of environments
#' w <- tpe_weights(rp, stage = "R1", var = "PETP")
#' }
#'
#' @seealso \code{\link{env_phenology}}, \code{\link{env_typing}},
#'   \code{\link{tpe_weights}}, \code{\link{project_risk}}
#'
#' @references
#' Chapman, S.C. et al. (2000) Using crop simulation to generate genotype by
#' environment interaction effects for sorghum in water-limited environments.
#' \emph{Aust. J. Agric. Res.} 51(2), 209-221.
#'
#' Heinemann, A.B. et al. (2015) Drought impact on rainfed common bean
#' production areas in Brazil. \emph{Agric. For. Meteorol.} 225, 1-12.
#'
#' @importFrom stats sd quantile aggregate
#' @export
env_risk_profile <- function(env.data, site.id = "env", date.id = "YYYYMMDD",
                             var.id = NULL, years = NULL,
                             planting = NULL, season.length = 150,
                             cardinals = NULL, by.stage = TRUE,
                             crop = "maize", stages = NULL,
                             Tbase = NULL, Tupper = NULL, method = "gdd",
                             min.days = 30L, return.seasons = FALSE,
                             verbose = TRUE) {

  .et_banner("env_risk_profile", "multi-year envirotype frequency and risk", verbose)
  env.data <- as.data.frame(env.data)

  if (!site.id %in% names(env.data))
    stop("Site column '", site.id, "' not found.", call. = FALSE)
  if (!date.id %in% names(env.data))
    stop("Date column '", date.id, "' not found. get_weather() returns ",
         "'YYYYMMDD'.", call. = FALSE)

  dates <- as.Date(env.data[[date.id]])
  if (all(is.na(dates)))
    stop("Could not parse '", date.id, "' as dates.", call. = FALSE)
  env.data$.date <- dates
  env.data$.year <- as.numeric(format(dates, "%Y"))
  env.data$.site <- as.character(env.data[[site.id]])

  ## ---- cardinals ---------------------------------------------------------
  if (is.null(cardinals)) cardinals <- .RISK_CARDINALS
  if (is.null(var.id)) {
    var.id <- intersect(names(cardinals), names(env.data))
    if (!length(var.id))
      stop("None of the default risk variables (",
           paste(names(.RISK_CARDINALS), collapse = ", "),
           ") are present. Supply 'var.id' and 'cardinals'.", call. = FALSE)
  }
  miss <- setdiff(var.id, names(env.data))
  if (length(miss))
    stop("Variable(s) not found: ", paste(miss, collapse = ", "),
         ". Did you run processWTH()?", call. = FALSE)
  no_card <- setdiff(var.id, names(cardinals))
  if (length(no_card))
    stop("No cardinals defined for: ", paste(no_card, collapse = ", "),
         ". Add them to 'cardinals' as list(breaks = , labels = ).",
         call. = FALSE)

  ## ---- build site-years --------------------------------------------------
  all_years <- sort(unique(env.data$.year))
  if (is.null(years)) years <- all_years
  years <- intersect(years, all_years)
  if (!length(years)) stop("No requested year is present in the data.", call. = FALSE)

  sites <- sort(unique(env.data$.site))
  .et_step(paste0("assembling ", length(sites), " site(s) x ",
                  length(years), " season(s)"), verbose)

  # resolve planting date for each site-year
  .plant_for <- function(site, yr) {
    if (is.null(planting)) return(as.Date(paste0(yr, "-01-01")))
    if (inherits(planting, "Date")) {
      cand <- planting[format(planting, "%Y") == as.character(yr)]
      return(if (length(cand)) cand[1] else as.Date(NA))
    }
    md <- if (!is.null(names(planting)) && site %in% names(planting)) {
      planting[[site]]
    } else planting[1]
    as.Date(paste0(yr, "-", md))
  }

  season_list <- list()
  for (s in sites) {
    ds <- env.data[env.data$.site == s, , drop = FALSE]
    for (y in years) {
      p0 <- .plant_for(s, y)
      if (is.na(p0)) next
      p1 <- p0 + season.length - 1
      chunk <- ds[!is.na(ds$.date) & ds$.date >= p0 & ds$.date <= p1, , drop = FALSE]
      if (nrow(chunk) < min.days) next
      chunk <- chunk[order(chunk$.date), , drop = FALSE]
      chunk$.season <- y
      chunk$.siteyear <- paste(s, y, sep = "_")
      chunk$daysFromStart <- as.numeric(chunk$.date - p0) + 1
      season_list[[length(season_list) + 1L]] <- chunk
    }
  }
  if (!length(season_list))
    stop("No site-year met 'min.days' = ", min.days,
         ". Check 'planting', 'season.length' and the date range.", call. = FALSE)

  SY <- .rbind_fill(season_list)
  n_sy <- length(unique(SY$.siteyear))
  .et_step(paste0(n_sy, " complete site-year(s) retained"), verbose)

  ## ---- stage assignment --------------------------------------------------
  if (isTRUE(by.stage)) {
    .et_step("assigning thermal-time stages per season", verbose)
    ph <- env_phenology(SY, env.id = ".siteyear", day.id = "daysFromStart",
                        crop = crop, stages = stages,
                        Tbase = Tbase, Tupper = Tupper, method = method,
                        merge = FALSE, verbose = FALSE)
    SY$stage <- ph$stage
    SY <- SY[!is.na(SY$stage), , drop = FALSE]
  } else {
    SY$stage <- factor("whole_season")
  }

  ## ---- classify and tabulate --------------------------------------------
  .et_step("classifying stress envirotypes", verbose)
  freq_rows <- list(); sev_rows <- list(); season_rows <- list()

  for (v in var.id) {
    cd  <- cardinals[[v]]
    if (is.null(cd$breaks) || is.null(cd$labels))
      stop("cardinals[['", v, "']] must be list(breaks = , labels = ).",
           call. = FALSE)
    if (length(cd$labels) != length(cd$breaks) + 1L)
      stop("cardinals[['", v, "']]: 'labels' must be one longer than 'breaks'.",
           call. = FALSE)

    vals <- .strip_power_na(SY[[v]])
    cls  <- .risk_classify(vals, cd$breaks, cd$labels)

    key <- paste(SY$.site, SY$stage, SY$.season, sep = "\r")
    # days in each class per site-stage-season
    tab <- as.data.frame(table(key = key, class = cls),
                         stringsAsFactors = FALSE)
    tab <- tab[tab$Freq >= 0, , drop = FALSE]
    parts <- do.call(rbind, strsplit(tab$key, "\r", fixed = TRUE))
    tab$site <- parts[, 1]; tab$stage <- parts[, 2]; tab$season <- parts[, 3]
    tab$variable <- v
    names(tab)[names(tab) == "Freq"] <- "days"
    season_rows[[v]] <- tab[, c("site", "stage", "season", "variable",
                                "class", "days")]

    # frequency across seasons
    g <- paste(tab$site, tab$stage, tab$class, sep = "\r")
    agg <- data.frame(
      key      = names(tapply(tab$days, g, length)),
      n_years  = as.numeric(tapply(tab$days, g, length)),
      n_occur  = as.numeric(tapply(tab$days, g, function(z) sum(z > 0))),
      mean_days = as.numeric(tapply(tab$days, g, mean)),
      max_days  = as.numeric(tapply(tab$days, g, max)),
      stringsAsFactors = FALSE)
    p <- do.call(rbind, strsplit(agg$key, "\r", fixed = TRUE))
    agg$site <- p[, 1]; agg$stage <- p[, 2]; agg$class <- p[, 3]
    agg$variable <- v
    agg$freq <- agg$n_occur / agg$n_years
    freq_rows[[v]] <- agg[, c("site", "stage", "variable", "class",
                              "n_years", "n_occur", "freq",
                              "mean_days", "max_days")]

    # severity: distribution of the raw variable across seasons
    gs <- paste(SY$.site, SY$stage, sep = "\r")
    sev <- data.frame(
      key  = names(tapply(vals, gs, function(z) mean(z, na.rm = TRUE))),
      mean = as.numeric(tapply(vals, gs, function(z) mean(z, na.rm = TRUE))),
      sd   = as.numeric(tapply(vals, gs, function(z) stats::sd(z, na.rm = TRUE))),
      p10  = as.numeric(tapply(vals, gs, function(z)
        stats::quantile(z, .10, na.rm = TRUE))),
      p90  = as.numeric(tapply(vals, gs, function(z)
        stats::quantile(z, .90, na.rm = TRUE))),
      stringsAsFactors = FALSE)
    ps <- do.call(rbind, strsplit(sev$key, "\r", fixed = TRUE))
    sev$site <- ps[, 1]; sev$stage <- ps[, 2]; sev$variable <- v
    sev_rows[[v]] <- sev[, c("site", "stage", "variable", "mean", "sd", "p10", "p90")]
  }

  frequency <- do.call(rbind, freq_rows)
  severity  <- do.call(rbind, sev_rows)
  rownames(frequency) <- rownames(severity) <- NULL

  # order stages sensibly
  if (isTRUE(by.stage)) {
    lv <- levels(SY$stage)
    frequency$stage <- factor(frequency$stage, levels = lv, ordered = TRUE)
    severity$stage  <- factor(severity$stage,  levels = lv, ordered = TRUE)
    frequency <- frequency[order(frequency$site, frequency$stage,
                                 frequency$variable, frequency$class), ]
  }

  ## ---- return periods ----------------------------------------------------
  rp <- frequency[frequency$freq > 0, , drop = FALSE]
  rp$return_period <- 1 / rp$freq
  rp$reliable <- rp$return_period <= (rp$n_years / 3)
  rp <- rp[, c("site", "stage", "variable", "class", "freq",
               "return_period", "n_years", "reliable")]
  rownames(rp) <- NULL

  out <- list(frequency = frequency,
              severity  = severity,
              return_period = rp,
              seasons = if (isTRUE(return.seasons)) do.call(rbind, season_rows) else NULL,
              years = years, n_siteyears = n_sy,
              cardinals = cardinals[var.id],
              by.stage = by.stage,
              call = match.call())
  class(out) <- "env_risk"
  out
}

#' @export
print.env_risk <- function(x, ...) {
  cat(strrep("-", 63), "\n")
  cat("env_risk_profile() - historical envirotype risk\n")
  cat(strrep("-", 63), "\n")
  cat("Sites .............. ", length(unique(x$frequency$site)), "\n", sep = "")
  cat("Seasons ............ ", length(x$years), " (",
      min(x$years), "-", max(x$years), ")\n", sep = "")
  cat("Site-years ......... ", x$n_siteyears, "\n", sep = "")
  cat("Variables .......... ", paste(unique(x$frequency$variable),
                                     collapse = ", "), "\n", sep = "")
  cat("Stage-resolved ..... [ ", if (isTRUE(x$by.stage)) "x" else " ", " ]\n", sep = "")

  top <- x$frequency[x$frequency$freq > 0, , drop = FALSE]
  # show the most frequent non-benign classes
  benign <- c("safe", "optimal", "none", "non_limiting", "moderate")
  top <- top[!top$class %in% benign, , drop = FALSE]
  if (nrow(top)) {
    top <- top[order(-top$freq, -top$mean_days), , drop = FALSE]
    cat("\nMost frequent stress classes:\n")
    n <- min(10L, nrow(top))
    for (i in seq_len(n))
      cat(sprintf("  %-8s %-6s %-9s %-14s freq %.2f  %5.1f d/season\n",
                  top$site[i], as.character(top$stage[i]), top$variable[i],
                  top$class[i], top$freq[i], top$mean_days[i]))
  }
  cat(strrep("-", 63), "\n")
  cat("$frequency  $severity  $return_period\n")
  invisible(x)
}


# =========================================================================
# SECTION - tpe_weights()
# =========================================================================

#' @title Target Population of Environments (TPE) Weights
#'
#' @description
#' Converts an \code{env_risk_profile()} into per-environment-type weights
#' summing to one, so that selection indices, weighted GxE models and testing
#' network decisions can be made proportional to how often each environment
#' type actually occurs -- rather than how often it happened to be trialled.
#'
#' @param risk An \code{"env_risk"} object from \code{env_risk_profile()}.
#' @param stage character. Restrict to one stage (e.g. \code{"R1"}). If
#'   \code{NULL}, weights are computed across all stages jointly.
#' @param var character. Restrict to one variable. If \code{NULL}, classes from
#'   all variables are pooled, which is only sensible if they are disjoint.
#' @param by character. \code{"class"} (default) weights environment TYPES;
#'   \code{"site"} weights SITES by how representative they are of the TPE.
#' @param measure character. \code{"freq"} (probability of occurrence, default)
#'   or \code{"days"} (expected days, which weights by intensity).
#' @param drop.benign boolean. Exclude classes named \code{safe},
#'   \code{optimal}, \code{none}, \code{non_limiting}. Default \code{FALSE}.
#' @param verbose boolean. Print the weight table.
#'
#' @return A data.frame with the grouping column, the raw measure, and
#'   \code{weight} summing to 1. Attribute \code{"measure"} records which
#'   measure was used.
#'
#' @details
#' \strong{What "site weight" means.} With \code{by = "site"} the weight is the
#' share of the TPE's total stress-occurrence mass that a site contributes. A
#' high weight means the site frequently expresses the stresses of interest,
#' NOT that it is a good discriminator of genotypes -- those are different
#' criteria, and a site can be highly representative yet have poor
#' heritability. Combine these weights with per-site repeatability before
#' cutting a testing location.
#'
#' \strong{Weights are as good as the record.} They inherit every sampling
#' limitation of the underlying profile, and they assume the historical
#' climatology still holds.
#'
#' @examples
#' \dontrun{
#' rp <- env_risk_profile(wth, planting = "10-15", crop = "maize")
#'
#' ## How much of the TPE is each drought class at flowering?
#' tpe_weights(rp, stage = "R1", var = "PETP")
#'
#' ## Which sites carry the most of the TPE's stress mass?
#' tpe_weights(rp, var = "PETP", by = "site", drop.benign = TRUE)
#' }
#'
#' @seealso \code{\link{env_risk_profile}}
#' @export
tpe_weights <- function(risk, stage = NULL, var = NULL,
                        by = c("class", "site"),
                        measure = c("freq", "days"),
                        drop.benign = FALSE, verbose = TRUE) {

  by <- match.arg(by); measure <- match.arg(measure)
  if (!inherits(risk, "env_risk"))
    stop("'risk' must be an object from env_risk_profile().", call. = FALSE)

  f <- risk$frequency
  if (!is.null(stage)) {
    f <- f[as.character(f$stage) %in% stage, , drop = FALSE]
    if (!nrow(f)) stop("No rows for stage '", paste(stage, collapse = "/"),
                       "'. Available: ",
                       paste(unique(as.character(risk$frequency$stage)),
                             collapse = ", "), call. = FALSE)
  }
  if (!is.null(var)) {
    f <- f[f$variable %in% var, , drop = FALSE]
    if (!nrow(f)) stop("No rows for variable '", paste(var, collapse = "/"),
                       "'.", call. = FALSE)
  }
  if (isTRUE(drop.benign))
    f <- f[!f$class %in% c("safe", "optimal", "none", "non_limiting"), ,
           drop = FALSE]
  if (!nrow(f)) stop("Nothing left after filtering.", call. = FALSE)

  val <- if (measure == "freq") f$freq else f$mean_days
  grp <- if (by == "class") f$class else f$site

  tot <- tapply(val, grp, sum, na.rm = TRUE)
  out <- data.frame(group = names(tot), value = as.numeric(tot),
                    stringsAsFactors = FALSE)
  names(out)[1] <- by
  s <- sum(out$value, na.rm = TRUE)
  if (!is.finite(s) || s <= 0)
    stop("All weights are zero; nothing to normalise.", call. = FALSE)
  out$weight <- out$value / s
  out <- out[order(-out$weight), , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "measure") <- measure

  if (isTRUE(verbose)) {
    cat(strrep("-", 48), "\n")
    cat("TPE weights by ", by, " (measure = ", measure, ")\n", sep = "")
    cat(strrep("-", 48), "\n")
    for (i in seq_len(nrow(out)))
      cat(sprintf("  %-18s %6.3f  %s\n", out[[by]][i], out$weight[i],
                  strrep("=", round(out$weight[i] * 30))))
    cat(strrep("-", 48), "\n")
  }
  out
}


# =========================================================================
# SECTION - project_risk()
# =========================================================================

#' @title Compare Envirotype Risk Between Baseline and a Climate Scenario
#'
#' @description
#' Runs \code{env_risk_profile()} on an observed series and on a
#' delta-shifted scenario series, and returns the change in stress frequency
#' per site, stage and class.
#'
#' @param baseline data.frame. Observed multi-year daily weather.
#' @param scenario data.frame. Output of \code{get_climate_scenario()} built
#'   from the same baseline, or any comparably structured series.
#' @param ... Arguments passed to \code{env_risk_profile()} for BOTH runs, so
#'   the two profiles are strictly comparable.
#' @param verbose boolean. Print progress.
#'
#' @return A list of class \code{"env_risk_delta"} with \code{baseline},
#'   \code{scenario} (both \code{"env_risk"} objects) and \code{change}, a
#'   data.frame with \code{freq_base}, \code{freq_scen}, \code{freq_diff} and
#'   \code{days_diff}.
#'
#' @details
#' Because \code{get_climate_scenario(method = "delta")} preserves the observed
#' day-to-day sequence, differences here reflect the MEAN shift only. Rising
#' frequencies of heat classes are real but conservative; unchanged
#' precipitation-class frequencies do not mean rainfall risk is unchanged, only
#' that the mean-preserving perturbation left wet-day counts alone.
#'
#' @examples
#' \dontrun{
#' fut <- get_climate_scenario(env.id = sites$env, lat = sites$lat,
#'                             lon = sites$lon, baseline = wth,
#'                             scenario = "ssp585", period = "2061-2080")
#'
#' chg <- project_risk(baseline = processWTH(wth),
#'                     scenario = processWTH(fut),
#'                     planting = "10-15", crop = "maize",
#'                     var.id = c("T2M_MAX", "PETP"))
#' subset(chg$change, class == "heat_severe")
#' }
#'
#' @seealso \code{\link{env_risk_profile}}, \code{\link{get_climate_scenario}}
#' @export
project_risk <- function(baseline, scenario, ..., verbose = TRUE) {

  .et_banner("project_risk", "baseline vs scenario envirotype risk", verbose)

  .et_step("profiling baseline", verbose)
  rb <- env_risk_profile(baseline, ..., verbose = FALSE)
  .et_step("profiling scenario", verbose)
  rs <- env_risk_profile(scenario, ..., verbose = FALSE)

  kb <- paste(rb$frequency$site, rb$frequency$stage,
              rb$frequency$variable, rb$frequency$class, sep = "\r")
  ks <- paste(rs$frequency$site, rs$frequency$stage,
              rs$frequency$variable, rs$frequency$class, sep = "\r")

  allk <- union(kb, ks)
  p <- do.call(rbind, strsplit(allk, "\r", fixed = TRUE))
  chg <- data.frame(site = p[, 1], stage = p[, 2], variable = p[, 3],
                    class = p[, 4], stringsAsFactors = FALSE)
  chg$freq_base <- rb$frequency$freq[match(allk, kb)]
  chg$freq_scen <- rs$frequency$freq[match(allk, ks)]
  chg$days_base <- rb$frequency$mean_days[match(allk, kb)]
  chg$days_scen <- rs$frequency$mean_days[match(allk, ks)]
  chg[is.na(chg)] <- 0
  chg$freq_diff <- chg$freq_scen - chg$freq_base
  chg$days_diff <- chg$days_scen - chg$days_base
  chg <- chg[order(-abs(chg$freq_diff)), ]
  rownames(chg) <- NULL

  out <- list(baseline = rb, scenario = rs, change = chg)
  class(out) <- "env_risk_delta"

  if (verbose) {
    message(strrep("-", 63))
    message("Largest shifts in stress frequency")
    n <- min(10L, nrow(chg))
    for (i in seq_len(n))
      message(sprintf("  %-8s %-6s %-9s %-14s %+.2f  (%+.1f d)",
                      chg$site[i], chg$stage[i], chg$variable[i],
                      chg$class[i], chg$freq_diff[i], chg$days_diff[i]))
    message(strrep("-", 63))
  }
  out
}

#' @export
print.env_risk_delta <- function(x, ...) {
  cat(strrep("-", 63), "\n")
  cat("project_risk() - baseline vs scenario\n")
  cat(strrep("-", 63), "\n")
  sc <- attr(x$scenario, "scenario")
  cat("Site-years (baseline / scenario) ... ",
      x$baseline$n_siteyears, " / ", x$scenario$n_siteyears, "\n", sep = "")
  up <- x$change[x$change$freq_diff > 0, , drop = FALSE]
  cat("Classes increasing in frequency .... ", nrow(up), "\n", sep = "")
  cat("\n$change  $baseline  $scenario\n")
  invisible(x)
}


# =========================================================================
# T_matrix() - the envirotype counterpart of W_matrix()
# =========================================================================
#
# env_typing() returns FOUR different shapes depending on its arguments (a long
# data.frame, a wide matrix, or 2- and 4-element lists). So downstream code
# cannot write env_kernel(env.data = env_typing(...)) without knowing which
# branch ran. T_matrix() gives envirotypes the same contract W_matrix() gives
# covariables: ALWAYS a q x k numeric matrix, rownames = environments, ready for
# env_kernel(), env_cluster() and the rest of the pipeline. Everything else
# env_typing() produces (descriptions, mining summaries, plots) is preserved as
# ATTRIBUTES, not as list elements that change the return type.
#
# Reference: Costa-Neto, Crossa & Fritsche-Neto (2021), Front. Plant Sci. 12,
# 717552 -- enviromic assembly for genomic prediction of yield plasticity.
# -------------------------------------------------------------------------

#' Pull the typology object out of whatever env_typing() returned.
#'
#' Handles all four return shapes, normalising the variation in ONE place so the
#' rest of T_matrix() sees a single form.
#' @noRd
.tm_unwrap <- function(x) {
  side <- list(envirotype_description = NULL, mining_summary = NULL, plot = NULL)

  if (is.list(x) && !is.data.frame(x) && "typologies" %in% names(x)) {
    for (nm in names(side)) if (!is.null(x[[nm]])) side[[nm]] <- x[[nm]]
    ty <- x$typologies
  } else {
    ty <- x
  }
  list(typologies = ty, side = side)
}

#' Coerce a long typology data.frame to wide env x envirotype.
#'
#' Mirrors the acast() call inside env_typing()'s own format = "wide" branch so
#' the two routes agree exactly: median over duplicates, NA -> 0.
#' @noRd
.tm_long_to_wide <- function(d, env.col = "env", type.col = "env.variable",
                             value.col = "Freq") {
  miss <- setdiff(c(env.col, type.col, value.col), names(d))
  if (length(miss))
    stop("Long typology frame is missing column(s): ",
         paste(miss, collapse = ", "),
         ".\n  Expected the output of env_typing(format = 'long').",
         call. = FALSE)

  d <- d[!is.na(d[[type.col]]), , drop = FALSE]
  if (!nrow(d)) stop("Typology frame has no non-missing envirotypes.", call. = FALSE)

  envs  <- sort(unique(as.character(d[[env.col]])))
  types <- unique(as.character(d[[type.col]]))

  M <- matrix(NA_real_, length(envs), length(types), dimnames = list(envs, types))
  key <- paste(as.character(d[[env.col]]), as.character(d[[type.col]]), sep = "\r")
  agg <- tapply(as.numeric(d[[value.col]]), key, stats::median, na.rm = TRUE)
  kk  <- strsplit(names(agg), "\r", fixed = TRUE)
  ri  <- match(vapply(kk, `[`, "", 1L), envs)
  ci  <- match(vapply(kk, `[`, "", 2L), types)
  M[cbind(ri, ci)] <- as.numeric(agg)
  M[is.na(M)] <- 0
  M
}

#' Row-normalise so each environment's frequencies sum to 1.
#' @noRd
.tm_row_ratio <- function(M) {
  rs <- rowSums(M, na.rm = TRUE)
  bad <- !is.finite(rs) | rs <= 0
  if (any(bad))
    warning(sum(bad), " environment(s) have zero total frequency; ",
            "left as zero rows: ",
            paste(rownames(M)[bad], collapse = ", "), call. = FALSE)
  rs[bad] <- 1
  M / rs
}

#' @title Envirotype Frequency Matrix, Ready for the Kernel Pipeline
#'
#' @description
#' The envirotype counterpart of \code{\link{W_matrix}}. Runs
#' \code{\link{env_typing}} and returns, unconditionally, an
#' environment x envirotype numeric matrix suitable for \code{\link{env_kernel}}
#' and \code{\link{env_cluster}} -- regardless of which of the four
#' \code{env_typing()} branches ran. Where \code{W_matrix} summarises continuous
#' covariables, \code{T_matrix} summarises the \emph{frequencies} of discrete
#' environmental types (envirotypes), the "enviromic assembly" building block of
#' Costa-Neto, Crossa & Fritsche-Neto (2021).
#'
#' @param env.data data.frame of processed weather data, or an object already
#'   returned by \code{\link{env_typing}} (any of its four shapes), in which case
#'   it is reshaped rather than recomputed.
#' @param var.id,env.id,cardinals,days.id,time.window,names.window,quantiles
#'   Passed to \code{\link{env_typing}}.
#' @param id.names,by.interval,joint,joint.by.interval Passed to
#'   \code{\link{env_typing}}.
#' @param envirotype_mining,k.range,nstart,iter.max,seed Passed to
#'   \code{\link{env_typing}}.
#' @param combine.features,combine.max.order Passed to \code{\link{env_typing}}.
#' @param ratio logical. Convert counts to within-environment relative
#'   frequencies. Default \code{TRUE} -- see Details.
#' @param center,scale logical. Passed to the shared scaler. \strong{Both default
#'   to \code{FALSE}}, unlike \code{\link{W_matrix}} -- see Details.
#' @param sd.tol,tol,QC,impute Passed to the shared scaler.
#' @param drop.constant logical. Drop envirotype columns that are constant across
#'   environments (they contribute nothing to a kernel). Default \code{TRUE}.
#' @param verbose logical.
#'
#' @details
#' \strong{Why \code{ratio = TRUE} by default.} Raw envirotype counts are days,
#' so an environment with a longer season has larger counts in every column. A
#' kernel built on raw counts is then partly a kernel on season length.
#' Row-normalising makes each environment a composition -- the share of its
#' season spent in each envirotype -- which is what "envirotype frequency" is
#' meant to convey.
#'
#' \strong{Why \code{center = scale = FALSE} by default.} This deliberately
#' differs from \code{\link{W_matrix}}, where \code{center = scale = TRUE}. After
#' \code{ratio = TRUE} the columns are already on a common [0, 1] scale and each
#' row sums to 1. Scaling a composition to unit variance per column inflates rare
#' envirotypes -- a type present in 2\% of days in one environment and 0\%
#' elsewhere becomes a high-leverage column. If you want the \code{W_matrix}
#' convention, set them explicitly; the choice is yours but it should be a choice.
#'
#' \strong{Compositional caveat.} With \code{ratio = TRUE} the rows sum to 1, so
#' the columns are linearly dependent and the matrix is rank-deficient by exactly
#' one. This is harmless for \code{\link{env_kernel}} (a Gaussian/linear kernel is
#' still positive semi-definite) but will produce one zero eigenvalue. Do not
#' "fix" it by dropping a column -- that makes the result depend on which column
#' you dropped.
#'
#' @return A numeric matrix of class \code{c("T_matrix", "matrix", "array")},
#'   with \code{rownames} = environments and \code{colnames} = envirotypes.
#'   Attributes: \code{"typologies"} (the long frame),
#'   \code{"envirotype_description"}, \code{"mining_summary"}, \code{"plot"},
#'   \code{"ratio"}, \code{"dropped"}, and the \code{.env_w_scale()} attributes
#'   when scaling is requested.
#'
#' @examples
#' \donttest{
#' data(maizeWTH)
#' Tm <- T_matrix(maizeWTH, var.id = c("T2M", "PRECTOT"), env.id = "env")
#' dim(Tm)
#' K <- env_kernel(env.data = Tm, is.scaled = TRUE)
#' }
#'
#' @seealso \code{\link{W_matrix}} for the covariable counterpart,
#'   \code{\link{env_typing}} for the underlying typology,
#'   \code{\link{env_kernel}} and \code{\link{env_cluster}} for the consumers.
#'
#' @references
#' Costa-Neto, G., Crossa, J., & Fritsche-Neto, R. (2021). Enviromic assembly
#' increases accuracy and reduces costs of the genomic prediction for yield
#' plasticity in maize. \emph{Frontiers in Plant Science} 12, 717552.
#'
#' @importFrom stats sd median
#' @export
T_matrix <- function(env.data, var.id = NULL, env.id = NULL,
                     cardinals = NULL, days.id = NULL,
                     time.window = NULL, names.window = NULL,
                     quantiles = NULL, id.names = NULL,
                     by.interval = FALSE,
                     joint = FALSE, joint.by.interval = TRUE,
                     envirotype_mining = FALSE, k.range = 2:10,
                     nstart = 25, iter.max = 100, seed = 1234,
                     combine.features = FALSE, combine.max.order = 2L,
                     ratio = TRUE,
                     center = FALSE, scale = FALSE, sd.tol = 10,
                     tol = 1e-3, QC = FALSE,
                     impute = c("none", "mean", "drop"),
                     drop.constant = TRUE,
                     verbose = TRUE) {

  impute <- match.arg(impute)
  if (exists(".et_banner"))
    .et_banner("T_matrix", "builds the envirotype (T) matrix", verbose)

  # ---- 1. obtain typologies ------------------------------------------------
  already <- is.matrix(env.data) ||
    (is.list(env.data) && !is.data.frame(env.data) &&
       "typologies" %in% names(env.data)) ||
    (is.data.frame(env.data) && all(c("env", "env.variable", "Freq") %in% names(env.data)))

  if (already) {
    if (exists(".et_step")) .et_step("reshaping an existing typology object", verbose)
    raw <- env.data
  } else {
    if (is.null(var.id) || is.null(env.id))
      stop("'var.id' and 'env.id' are required when 'env.data' is raw weather ",
           "data.", call. = FALSE)
    if (exists(".et_step")) .et_step("mining envirotypes with env_typing()", verbose)
    raw <- env_typing(env.data = env.data, var.id = var.id, env.id = env.id,
                      cardinals = cardinals, days.id = days.id,
                      time.window = time.window, names.window = names.window,
                      quantiles = quantiles, id.names = id.names,
                      by.interval = by.interval, format = "long",
                      ratio = FALSE,                 # normalise here, not there
                      joint = joint, joint.by.interval = joint.by.interval,
                      envirotype_mining = envirotype_mining,
                      k.range = k.range, nstart = nstart, iter.max = iter.max,
                      seed = seed, combine.features = combine.features,
                      combine.max.order = combine.max.order,
                      plot_envirotypes = FALSE, verbose = FALSE)
  }

  # ---- 2. normalise the FOUR possible shapes to one ------------------------
  uw   <- .tm_unwrap(raw)
  ty   <- uw$typologies
  side <- uw$side

  if (is.matrix(ty)) {
    M    <- ty
    long <- NULL
  } else if (is.data.frame(ty)) {
    long <- ty
    M    <- .tm_long_to_wide(ty)
  } else {
    stop("Unrecognised typology object of class <",
         paste(class(ty), collapse = "/"),
         ">. Expected a data.frame, a matrix, or a list with $typologies.",
         call. = FALSE)
  }

  storage.mode(M) <- "double"
  if (!nrow(M) || !ncol(M))
    stop("Typology matrix is empty (", nrow(M), " x ", ncol(M), ").",
         call. = FALSE)

  # ---- 3. compositional normalisation --------------------------------------
  if (isTRUE(ratio)) {
    if (exists(".et_step"))
      .et_step("normalising envirotype counts to within-environment shares", verbose)
    M <- .tm_row_ratio(M)
  }

  # ---- 4. drop uninformative columns ---------------------------------------
  dropped <- character(0)
  if (isTRUE(drop.constant) && nrow(M) > 1L) {
    cs <- apply(M, 2, function(v) stats::sd(v, na.rm = TRUE))
    bad <- !is.finite(cs) | cs < tol
    if (any(bad)) {
      dropped <- colnames(M)[bad]
      if (all(bad))
        stop("All ", ncol(M), " envirotype column(s) are constant across ",
             "environments, so the matrix carries no information for a ",
             "kernel.\n  This usually means every environment was assigned ",
             "the same typology. Check 'cardinals'/'quantiles', or widen ",
             "'var.id'.", call. = FALSE)
      M <- M[, !bad, drop = FALSE]
      if (verbose)
        message("  dropped ", length(dropped), " constant envirotype column(s).")
    }
  }

  # ---- 5. optional scaling, via the SHARED scaler --------------------------
  # Use the package's own .env_w_scale() so T and W carry identical attributes.
  if (isTRUE(center) || isTRUE(scale) || isTRUE(QC) || impute != "none") {
    if (exists(".env_w_scale")) {
      if (exists(".et_step")) .et_step("centring/scaling via .env_w_scale()", verbose)
      M <- .env_w_scale(env.data = M, center = center, scale = scale,
                        sd.tol = sd.tol, tol = tol, QC = QC,
                        impute = impute, verbose = verbose)
    } else {
      M <- scale(M, center = center, scale = scale)
    }
  }

  # ---- 6. one contract, always ---------------------------------------------
  attr(M, "typologies")             <- long
  attr(M, "envirotype_description") <- side$envirotype_description
  attr(M, "mining_summary")         <- side$mining_summary
  attr(M, "plot")                   <- side$plot
  attr(M, "ratio")                  <- ratio
  attr(M, "dropped")                <- dropped
  attr(M, "source")                 <- "T_matrix"
  class(M) <- c("T_matrix", "matrix", "array")

  if (verbose) print(M)
  M
}

#' @title Print an Envirotype (T) Matrix
#' @param x a \code{T_matrix} object.
#' @param ... ignored.
#' @return \code{x}, invisibly.
#' @seealso \code{\link{T_matrix}}
#' @export
print.T_matrix <- function(x, ...) {
  cat("<T_matrix>\n")
  cat(sprintf("  environments x envirotypes . %d x %d\n", nrow(x), ncol(x)))
  cat(sprintf("  row-normalised (ratio) ..... %s\n", attr(x, "ratio")))
  if (length(attr(x, "dropped")))
    cat(sprintf("  dropped constant columns ... %d\n", length(attr(x, "dropped"))))
  if (isTRUE(attr(x, "ratio"))) {
    rs <- range(rowSums(unclass(x)))
    cat(sprintf("  row sums ................... %.4f - %.4f\n", rs[1], rs[2]))
  }
  if (!is.null(attr(x, "mining_summary")))
    cat("  mining summary ............. present\n")
  if (!is.null(attr(x, "envirotype_description")))
    cat("  envirotype description ..... present\n")
  invisible(x)
}