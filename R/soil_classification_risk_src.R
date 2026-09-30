#==================================================================================================
# soil_classification.R
#
# Gaussian mixture model clustering of soil profiles from get_soil(), for
# delineating soil zones (management zones, risk zones, texture groupings).
#
# Method follows the general GMM zoning workflow described in:
#
#   Zhang, F., Jia, Z., Wu, S., Chen, C., Chen, X., Zheng, C., Xu, M. (2025)
#   "Machine Learning and Gaussian Mixture Model for Delineating Soil
#   Cadmium Risk Zones." Ecosystem Health and Sustainability.
#   doi:10.34133/ehs.0402
#
# IMPORTANT -- READ BEFORE USE. See "Relationship to the source paper" in the
# soil_classification() documentation. The publisher blocks automated access,
# so the paper's exact preprocessing, covariance families and model-selection
# settings could NOT be verified. This implements the standard GMM zoning
# pipeline the title describes. It is not a byte-for-byte reproduction, and it
# cannot reproduce the paper's cadmium work because SoilGrids carries no heavy
# metals (see "What this cannot do").
#
# Companion to env_data_collection.R and water_balance.R.
#
# --------------------------------------------------------------------------
# METHOD REFERENCES
# --------------------------------------------------------------------------
# Gaussian mixture models and BIC model selection
#   Fraley, C. & Raftery, A.E. (2002) Model-based clustering, discriminant
#     analysis and density estimation. J. Am. Stat. Assoc. 97(458), 611-631.
#     doi:10.1198/016214502760047131
#   Scrucca, L., Fop, M., Murphy, T.B. & Raftery, A.E. (2016) mclust 5:
#     clustering, classification and density estimation using Gaussian finite
#     mixture models. The R Journal 8(1), 205-233.
#     doi:10.32614/RJ-2016-021
#   Schwarz, G. (1978) Estimating the dimension of a model. Ann. Stat. 6(2),
#     461-464. doi:10.1214/aos/1176344136
#
# Compositional data (why texture needs a log-ratio transform)
#   Aitchison, J. (1982) The statistical analysis of compositional data.
#     J. R. Stat. Soc. B 44(2), 139-177.
#   Egozcue, J.J., Pawlowsky-Glahn, V., Mateu-Figueras, G. & Barcelo-Vidal, C.
#     (2003) Isometric logratio transformations for compositional data
#     analysis. Mathematical Geology 35(3), 279-300.
#     doi:10.1023/A:1023818214614
#   Filzmoser, P., Hron, K. & Reimann, C. (2009) Univariate statistical
#     analysis of environmental (compositional) data. Science of the Total
#     Environment 407(23), 6100-6108. doi:10.1016/j.scitotenv.2009.08.008
#
# Management-zone delineation in precision agriculture
#   Fridgen, J.J. et al. (2004) Management Zone Analyst (MZA): software for
#     subfield management zone delineation. Agronomy Journal 96(1), 100-108.
#     doi:10.2134/agronj2004.1000
#   Corwin, D.L. & Lesch, S.M. (2003) Application of soil electrical
#     conductivity to precision agriculture. Agronomy Journal 95(3), 455-471.
#     doi:10.2134/agronj2003.4550
#
# Soil data source
#   Poggio, L. et al. (2021) SoilGrids 2.0: producing soil information for the
#     globe with quantified spatial uncertainty. SOIL 7, 217-240.
#     doi:10.5194/soil-7-217-2021
#==================================================================================================


# =========================================================================
# SECTION 0 - internal helpers (not exported)
# =========================================================================

#' Centred log-ratio transform for compositional variables
#'
#' Soil texture fractions (clay, sand, silt) are compositional: they sum to a
#' constant, so they live on a simplex, not in Euclidean space. Fitting a
#' Gaussian directly to them is a modelling error -- the closure constraint
#' induces spurious negative correlation between the parts.
#'
#' The clr transform maps the simplex to real space:
#' \deqn{clr(x)_i = \log(x_i / g(x))}
#' where \eqn{g(x)} is the geometric mean of the composition.
#'
#' CAUTION. clr coordinates sum to zero at every site by construction, so a
#' clr'd composition of \eqn{D} parts spans only \eqn{D-1} dimensions and its
#' covariance matrix is SINGULAR. That is fatal for a Gaussian mixture. This
#' function is retained because clr coordinates are directly interpretable
#' (each is the log-ratio of one part to the geometric mean), but the fitting
#' pipeline uses \code{.ilr()} below. See Aitchison (1982), Egozcue et al.
#' (2003).
#'
#' @param x numeric matrix of compositional parts. Non-positive values are
#'   replaced before the log.
#' @return matrix of clr coordinates, same dimensions as \code{x}, rows summing
#'   to zero.
#' @keywords internal
#' @noRd
.clr <- function(x) {
  x <- as.matrix(x)
  if (any(x <= 0, na.rm = TRUE)) {
    # Multiplicative zero replacement: shift zeros off the simplex boundary
    # without destroying the ratios between the remaining parts. The 0.65
    # multiplier follows the usual recommendation for rounded zeros.
    small <- min(x[x > 0], na.rm = TRUE) * 0.65
    x[x <= 0] <- small
  }
  gm <- exp(rowMeans(log(x), na.rm = TRUE))
  log(x / gm)
}

#' Isometric log-ratio transform for compositional variables
#'
#' Maps \eqn{D} compositional parts to \eqn{D-1} orthonormal real coordinates
#' by applying a Helmert basis to the clr coordinates. Unlike clr the result is
#' full rank, so a Gaussian mixture can estimate a non-singular covariance.
#'
#' This is not a cosmetic choice. Applying clr to clay/sand/silt returns three
#' columns that sum to zero at every site; the covariance then has a zero
#' eigenvalue, EM either fails outright or silently inflates the likelihood of
#' over-complex models, and BIC selects the wrong number of clusters. This was
#' an actual defect in an earlier version of this file: ground-truth recovery
#' returned 2 clusters where 3 were planted, and the feature matrix had rank 8
#' of 10. Switching to ilr restored full rank and exact recovery.
#'
#' ilr preserves Aitchison distances, so distances between sites are unchanged
#' relative to clr -- only the redundant dimension is removed.
#'
#' @param x numeric matrix of compositional parts (at least 2 columns).
#' @return matrix with \code{ncol(x) - 1} orthonormal ilr coordinates.
#' @references Egozcue, J.J. et al. (2003) Mathematical Geology 35(3), 279-300.
#' @keywords internal
#' @noRd
.ilr <- function(x) {
  x <- as.matrix(x)
  D <- ncol(x)
  if (D < 2)
    stop("ilr needs at least two compositional parts.", call. = FALSE)
  # Helmert sub-matrix: D-1 orthonormal contrasts spanning the clr hyperplane.
  V <- matrix(0, nrow = D, ncol = D - 1)
  for (i in seq_len(D - 1)) {
    V[seq_len(i), i] <- 1 / i
    V[i + 1, i]      <- -1
    V[, i]           <- V[, i] * sqrt(i / (i + 1))
  }
  .clr(x) %*% V
}

#' Scale columns to zero mean and unit variance, dropping degenerate ones
#'
#' GMM covariance estimation fails on constant columns. Rather than let
#' \code{Mclust} error out with an opaque message, drop them and warn.
#'
#' A warning rather than a message: dropping a variable changes the clustering,
#' so it must be visible even when \code{verbose = FALSE}.
#'
#' @param m numeric matrix.
#' @return scaled matrix with degenerate columns removed.
#' @keywords internal
#' @noRd
.safe_scale <- function(m) {
  m <- as.matrix(m)
  sds <- apply(m, 2, stats::sd, na.rm = TRUE)
  bad <- is.na(sds) | sds < .Machine$double.eps^0.5
  if (any(bad)) {
    warning("Dropping ", sum(bad), " constant/degenerate variable(s): ",
            paste(colnames(m)[bad], collapse = ", "), call. = FALSE)
    m <- m[, !bad, drop = FALSE]
  }
  if (!ncol(m))
    stop("No variable has any variation; cannot cluster.", call. = FALSE)
  scale(m)
}

#' Order clusters by a composite risk/severity score
#'
#' Mixture component labels are arbitrary -- re-running, or changing the seed,
#' can permute them. For output to be interpretable and reproducible, clusters
#' are renumbered by an ordering score so that cluster 1 always means the same
#' thing (for example, always the lowest-SOC zone).
#'
#' @param cl integer vector of raw cluster labels.
#' @param scores numeric vector, one score per cluster.
#' @return integer vector of renumbered labels.
#' @keywords internal
#' @noRd
.order_clusters <- function(cl, scores) {
  ord <- order(scores)
  remap <- integer(length(scores))
  remap[ord] <- seq_along(ord)
  remap[cl]
}

#' Build the feature matrix from wide get_soil() output
#'
#' Parses \code{property|depth} column names, filters by requested variables
#' and depths, converts SoilGrids integers to conventional units, and replaces
#' texture parts with ilr coordinates.
#'
#' @keywords internal
#' @noRd
.soil_features <- function(soil.data, env.id, variables, depths,
                           convert.units, clr.texture, verbose) {

  nm  <- names(soil.data)
  key <- grep("\\|", nm, value = TRUE)
  if (!length(key))
    stop("No 'property|depth' columns found. Pass the output of ",
         "get_soil(wide = TRUE).", call. = FALSE)

  parts <- do.call(rbind, strsplit(key, "|", fixed = TRUE))
  meta  <- data.frame(col = key, property = parts[, 1], layer = parts[, 2],
                      stringsAsFactors = FALSE)

  if (!is.null(variables)) {
    unknown <- setdiff(variables, unique(meta$property))
    if (length(unknown))
      stop("Requested variable(s) not present in soil.data: ",
           paste(unknown, collapse = ", "), call. = FALSE)
    meta <- meta[meta$property %in% variables, , drop = FALSE]
  }
  if (!is.null(depths)) {
    unknown <- setdiff(depths, unique(meta$layer))
    if (length(unknown))
      stop("Requested depth(s) not present in soil.data: ",
           paste(unknown, collapse = ", "), call. = FALSE)
    meta <- meta[meta$layer %in% depths, , drop = FALSE]
  }
  if (!nrow(meta))
    stop("No soil columns left after filtering by 'variables'/'depths'.",
         call. = FALSE)

  X <- as.matrix(soil.data[, meta$col, drop = FALSE])
  storage.mode(X) <- "double"

  # SoilGrids integer -> conventional units. Without this, properties enter on
  # wildly different scales; standardising fixes the per-column scale but would
  # still distort the ratios between texture parts before the log-ratio step.
  if (isTRUE(convert.units)) {
    if (!exists(".SG_FACTORS"))
      stop("Conversion requires .SG_FACTORS from water_balance.R. ",
           "Source that file, or set convert.units = FALSE.", call. = FALSE)
    f <- .SG_FACTORS[meta$property]
    unknown <- unique(meta$property[is.na(f)])
    if (length(unknown) && verbose)
      message("No conversion factor for: ", paste(unknown, collapse = ", "),
              " (left unconverted).")
    f[is.na(f)] <- 1
    X <- sweep(X, 2, as.numeric(f), "/")
  }

  # Texture fractions are compositional. Transform per depth before they meet a
  # Gaussian, using ilr not clr: clr on D parts returns D columns summing to
  # zero, which is singular by construction. ilr returns D-1 full-rank coords.
  tex <- c("clay", "sand", "silt")
  colnames(X) <- meta$col
  if (isTRUE(clr.texture)) {
    keep_cols <- rep(TRUE, ncol(X))
    extra <- list()
    for (d in unique(meta$layer)) {
      sel <- which(meta$layer == d & meta$property %in% tex)
      if (length(sel) >= 2) {
        Z <- .ilr(X[, sel, drop = FALSE])
        colnames(Z) <- paste0("ilr", seq_len(ncol(Z)), "_texture|", d)
        extra[[length(extra) + 1L]] <- Z
        keep_cols[sel] <- FALSE
        if (verbose)
          message("ilr transform: ", length(sel), " texture parts at depth ",
                  d, " -> ", ncol(Z), " coordinate(s).")
      }
    }
    if (length(extra))
      X <- cbind(X[, keep_cols, drop = FALSE], do.call(cbind, extra))
  }

  rownames(X) <- as.character(soil.data[[env.id]])
  attr(X, "meta") <- meta
  X
}


# =========================================================================
# SECTION 1 - soil_classification()
# =========================================================================

#' @title Delineate Soil Zones with a Gaussian Mixture Model
#'
#' @description
#' Clusters soil profiles returned by \code{get_soil} into a small number
#' of zones using a Gaussian mixture model, selecting both the number of
#' clusters and the covariance structure by BIC. Returns hard cluster labels,
#' the full matrix of posterior membership probabilities, and a per-site
#' uncertainty measure.
#'
#' Unlike k-means, a mixture model gives every site a probability of belonging
#' to every zone. A site on a boundary is genuinely ambiguous, and that
#' ambiguity is reported rather than hidden behind a hard label.
#'
#' @param soil.data data.frame. Wide-format output of
#'   \code{get_soil(wide = TRUE)}, one row per site, with columns named
#'   \code{property|depth} (for example \code{"clay|0_5cm"}).
#' @param env.id character. Name of the site-id column. Default \code{"env"}.
#' @param variables character. Soil properties to cluster on, for example
#'   \code{c("clay", "sand", "soc")}. \code{NULL} (default) uses every property
#'   present. Filtering is strongly recommended -- see \strong{Choosing
#'   variables} below.
#' @param depths character. Depth layers to use, for example
#'   \code{c("0_5cm", "30_60cm")}. \code{NULL} (default) uses all layers
#'   present. Adjacent SoilGrids layers are highly correlated, so using all six
#'   rarely helps.
#' @param G integer vector. Candidate numbers of clusters. Default \code{1:9}.
#'   \code{G = 1} is included deliberately: if BIC selects it, the data provide
#'   no evidence for more than one zone, which is a legitimate and useful
#'   finding.
#' @param model.names character. \code{mclust} covariance families to try, for
#'   example \code{"EII"} (spherical, equal volume) or \code{"VVV"}
#'   (ellipsoidal, all free). \code{NULL} (default) lets \code{mclust} try every
#'   family it can fit. Restrict this when sites are few relative to variables.
#' @param scale.data logical. Standardise variables to zero mean and unit
#'   variance before fitting. Default \code{TRUE}, and strongly recommended:
#'   soil properties differ by orders of magnitude, and an unscaled fit is
#'   dominated by whichever variable happens to have the largest units.
#' @param clr.texture logical. Apply an isometric log-ratio (ilr) transform to
#'   texture fractions (clay/sand/silt) at each depth. Default \code{TRUE}.
#' @param convert.units logical. Apply SoilGrids conversion factors, turning
#'   mapped integers into conventional units. Requires \code{.SG_FACTORS} from
#'   \code{water_balance.R}. Default \code{TRUE}.
#' @param use.pca logical. Reduce the feature matrix to principal components
#'   before fitting. Default \code{FALSE}. Useful when variables are many and
#'   collinear.
#' @param pca.var numeric in (0, 1]. Proportion of variance the retained
#'   components must explain when \code{use.pca = TRUE}. Default \code{0.95}.
#' @param risk.vars character. Variables used to order the cluster labels, for
#'   example \code{"soc"}. \code{NULL} (default) orders by the mean of all
#'   scaled features. See \strong{Label ordering} below.
#' @param risk.direction \code{"high"} or \code{"low"}. Whether high values of
#'   \code{risk.vars} should sort last (\code{"high"}, the default, so cluster 1
#'   is the lowest) or first.
#' @param uncertainty.threshold numeric in (0, 1). A site whose maximum
#'   posterior probability falls below this is flagged \code{uncertain}.
#'   Default \code{0.7}.
#' @param seed integer or \code{NULL}. Seed for reproducible EM starts. Default
#'   \code{1}. Cluster \emph{labels} are seed-invariant because of the ordering
#'   step, though the underlying fit may vary.
#' @param verbose logical. Print progress. Default \code{TRUE}.
#'
#' @details
#' \strong{What the model does.} A Gaussian mixture assumes the sites are drawn
#' from \eqn{G} multivariate normal components with unknown means, covariances
#' and mixing weights, fitted by EM. \code{mclust} additionally searches over
#' constrained covariance families (spherical, diagonal, ellipsoidal; equal or
#' varying across components), and BIC selects both \eqn{G} and the family in
#' one criterion. See Fraley & Raftery (2002) and Scrucca et al. (2016).
#'
#' \strong{Compositional data.} Texture fractions sum to a constant. Feeding
#' them to a Gaussian untransformed is a category error: the closure constraint
#' forces spurious negative correlations. The default \code{clr.texture = TRUE}
#' applies an \emph{isometric} log-ratio transform per depth, mapping \eqn{D}
#' parts to \eqn{D-1} orthonormal coordinates. A plain clr transform is
#' deliberately not used for fitting: its coordinates sum to zero, leaving the
#' covariance singular and causing BIC to select the wrong number of clusters.
#' Turn this off only if your variables are not compositional.
#'
#' \strong{Label ordering.} Mixture component labels are arbitrary and permute
#' between runs. Clusters are therefore renumbered by \code{risk.vars} so that
#' cluster 1 is always the low end and cluster \eqn{G} the high end. Without
#' this, "cluster 3" would mean something different on every run and results
#' could not be compared or cached.
#'
#' \strong{Choosing variables.} A mixture estimates a mean vector and a
#' covariance per component, so the parameter count grows quickly with the
#' number of features. Passing every property at every depth gives 20+ highly
#' correlated columns and will overwhelm a modest number of sites. The function
#' warns when \code{n < 10 * p}. Prefer a handful of properties at two
#' contrasting depths, or set \code{use.pca = TRUE}.
#'
#' \strong{Sample size.} This is the single most common way to get a
#' meaningless result. Clustering a handful of sites produces a model in which
#' every site owns its own component: the output reports
#' \code{max_posterior = 1} and \code{uncertainty = 0} for every row, which
#' looks like total confidence and means nothing. Zone delineation samples a
#' grid of points across the area of interest -- typically hundreds -- not one
#' point per farm. Check \eqn{n} against \eqn{p} before believing any output.
#'
#' \strong{Reading the BIC gap.} \code{print()} reports the BIC difference
#' between the selected model and the runner-up, labelled weak (< 2), positive
#' (2-6), strong (6-10) or very strong (> 10) following the usual convention for
#' BIC differences. A weak gap means several partitions fit about equally well
#' and the zone count should not be presented as settled.
#'
#' \strong{Relationship to the source paper.} The paper behind this function
#' (Zhang et al., doi:10.34133/ehs.0402) could not be read: the publisher blocks
#' automated access. This implements the standard GMM zoning workflow its title
#' describes -- log-ratio transform, standardise, BIC-selected mixture,
#' posterior-based uncertainty -- but the paper's exact preprocessing,
#' covariance families and thresholds are unverified. Treat this as a
#' methodologically conventional starting point, not a reproduction.
#'
#' \strong{What this cannot do.} SoilGrids contains no cadmium, and no heavy
#' metal of any kind. The source paper's subject is cadmium risk; this function
#' can only cluster the properties SoilGrids actually provides (texture, organic
#' carbon, pH, CEC, water retention). Any "risk" interpretation rests entirely
#' on whichever proxies you pass to \code{risk.vars}, and the paper's machine
#' learning step -- predicting the contaminant before zoning -- has no
#' counterpart here. Do not present output from this function as a contaminant
#' risk map.
#'
#' @return An object of class \code{soil_gmm}, a list with:
#'   \describe{
#'     \item{\code{classification}}{data.frame, one row per input site in input
#'       order: \code{env}, \code{cluster}, \code{max_posterior},
#'       \code{uncertainty} (= 1 - max_posterior), and logical
#'       \code{uncertain}.}
#'     \item{\code{posterior}}{matrix of posterior membership probabilities,
#'       sites x clusters, each row summing to 1.}
#'     \item{\code{profiles}}{data.frame of per-cluster means in original units,
#'       plus \code{n_sites}.}
#'     \item{\code{model}}{list with \code{G}, \code{modelName},
#'       \code{bic}, \code{loglik}.}
#'     \item{\code{BIC}}{the full BIC table over \eqn{G} and covariance family.}
#'     \item{\code{features}}{the feature matrix actually clustered.}
#'     \item{\code{pca}}{the \code{prcomp} object, or \code{NULL}.}
#'     \item{\code{settings}}{the call arguments, for reproducibility.}
#'   }
#'
#' @section Output granularity:
#' The function classifies \emph{rows}, one output row per input row, in input
#' order. If each row is a grid point, you get a cluster per grid point; a farm
#' spanning several rows may well split across zones, and that split is usually
#' the informative result. Join back to your own table on \code{env}.
#'
#' @references
#' Zhang, F., Jia, Z., Wu, S., Chen, C., Chen, X., Zheng, C., Xu, M. (2025)
#' Machine Learning and Gaussian Mixture Model for Delineating Soil Cadmium
#' Risk Zones. \emph{Ecosystem Health and Sustainability}.
#' \doi{10.34133/ehs.0402}
#'
#' Fraley, C. & Raftery, A.E. (2002) Model-based clustering, discriminant
#' analysis and density estimation. \emph{Journal of the American Statistical
#' Association} 97(458), 611-631. \doi{10.1198/016214502760047131}
#'
#' Scrucca, L., Fop, M., Murphy, T.B. & Raftery, A.E. (2016) mclust 5:
#' clustering, classification and density estimation using Gaussian finite
#' mixture models. \emph{The R Journal} 8(1), 205-233.
#' \doi{10.32614/RJ-2016-021}
#'
#' Aitchison, J. (1982) The statistical analysis of compositional data.
#' \emph{Journal of the Royal Statistical Society B} 44(2), 139-177.
#'
#' Egozcue, J.J., Pawlowsky-Glahn, V., Mateu-Figueras, G. & Barcelo-Vidal, C.
#' (2003) Isometric logratio transformations for compositional data analysis.
#' \emph{Mathematical Geology} 35(3), 279-300. \doi{10.1023/A:1023818214614}
#'
#' Schwarz, G. (1978) Estimating the dimension of a model. \emph{Annals of
#' Statistics} 6(2), 461-464. \doi{10.1214/aos/1176344136}
#'
#' Poggio, L. et al. (2021) SoilGrids 2.0. \emph{SOIL} 7, 217-240.
#' \doi{10.5194/soil-7-217-2021}
#'
#' @seealso \code{get_soil} for the input, \code{water_balance}
#'   for per-site indices that need no cross-site sample,
#'   \code{Mclust} for the underlying model.
#'
#' @examples
#' \dontrun{
#' library(mclust)
#'
#' ## ---------------------------------------------------------------
#' ## 1. Minimal use
#' ## ---------------------------------------------------------------
#' soil <- get_soil(env.id = c("NM", "SO", "IP"),
#'                  lat = c(-13.05, -12.32, -21.98),
#'                  lon = c(-56.08, -55.71, -47.88),
#'                  variables.names = c("clay", "sand", "silt", "soc", "phh2o"))
#'
#' fit <- soil_classification(soil)
#' fit                                  # summary, BIC gap, cluster profiles
#' head(fit$classification)             # per-site cluster + uncertainty
#'
#'
#' ## ---------------------------------------------------------------
#' ## 2. Realistic zoning: sample a GRID, not one point per farm
#' ## ---------------------------------------------------------------
#' ## A mixture cannot be identified from a handful of rows. Sample the
#' ## area of interest, here 40 points around each of 5 locations.
#' sites <- data.frame(
#'   env = c("SOR", "LON", "SLG", "PAL", "BAR"),
#'   lat = c(-12.5453, -23.3045, -28.4083, -10.1689, -12.1530),
#'   lon = c(-55.7113, -51.1696, -54.9608, -48.3317, -44.9900))
#'
#' g <- expand.grid(k = 1:40, i = seq_len(nrow(sites)))
#' pts <- data.frame(
#'   env  = paste0(sites$env[g$i], "_", sprintf("%03d", g$k)),
#'   farm = sites$env[g$i],
#'   lat  = sites$lat[g$i] + runif(nrow(g), -0.05, 0.05),
#'   lon  = sites$lon[g$i] + runif(nrow(g), -0.05, 0.05))
#'
#' soil_grid <- get_soil(env.id = pts$env, lat = pts$lat, lon = pts$lon,
#'                       variables.names = c("clay", "sand", "silt",
#'                                           "soc", "phh2o"))
#'
#' fit <- soil_classification(soil_grid,
#'                            variables = c("clay", "sand", "silt",
#'                                          "soc", "phh2o"),
#'                            depths    = c("0_5cm", "30_60cm"),
#'                            risk.vars = "soc")
#' fit
#'
#' ## Do farms split across zones? Usually the interesting question.
#' table(farm = pts$farm, cluster = fit$classification$cluster)
#'
#'
#' ## ---------------------------------------------------------------
#' ## 3. Posterior probabilities and ambiguous sites
#' ## ---------------------------------------------------------------
#' round(head(fit$posterior), 3)        # membership in EVERY zone
#'
#' ## Sites the model is least sure about -- these sit on a boundary.
#' amb <- fit$classification[fit$classification$uncertain, ]
#' amb[order(amb$max_posterior), ]
#'
#' ## Tighten or relax the flag without refitting the model:
#' fit2 <- soil_classification(soil_grid, uncertainty.threshold = 0.9)
#' sum(fit2$classification$uncertain)
#'
#'
#' ## ---------------------------------------------------------------
#' ## 4. Ordering zones so cluster 1 always means the same thing
#' ## ---------------------------------------------------------------
#' lo  <- soil_classification(soil_grid, risk.vars = "soc",
#'                            risk.direction = "high")   # 1 = lowest SOC
#' hi  <- soil_classification(soil_grid, risk.vars = "soc",
#'                            risk.direction = "low")    # 1 = highest SOC
#' lo$profiles[, c("cluster", "n_sites")]
#'
#'
#' ## ---------------------------------------------------------------
#' ## 5. Many correlated variables: reduce first
#' ## ---------------------------------------------------------------
#' fit_pca <- soil_classification(soil_grid, use.pca = TRUE, pca.var = 0.95)
#' fit_pca$pca                          # how many components were kept
#'
#' ## Or restrict the covariance family when sites are few:
#' fit_par <- soil_classification(soil_grid, model.names = c("EII", "VII",
#'                                                           "EEI", "VVI"))
#'
#'
#' ## ---------------------------------------------------------------
#' ## 6. Joining results back to your own table
#' ## ---------------------------------------------------------------
#' out <- merge(pts, fit$classification, by = "env")
#' head(out[, c("env", "farm", "lat", "lon", "cluster",
#'              "max_posterior", "uncertain")])
#'
#' ## Dominant zone per farm, and whether the farm is homogeneous
#' with(out, table(farm, cluster))
#'
#'
#' ## ---------------------------------------------------------------
#' ## 7. Is the cluster count actually supported?
#' ## ---------------------------------------------------------------
#' fit$BIC                              # full table over G and family
#' fit$model$G                          # selected number of zones
#' ## print(fit) reports the gap to the runner-up; a weak gap (< 2) means
#' ## the zone count is not settled, whatever the point estimate says.
#' }
#'
#' @importFrom stats sd prcomp complete.cases
#' @importFrom utils head
#' @export
soil_classification <- function(soil.data, env.id = "env",
                                variables = NULL, depths = NULL,
                                G = 1:9, model.names = NULL,
                                scale.data = TRUE, clr.texture = TRUE,
                                convert.units = TRUE,
                                use.pca = FALSE, pca.var = 0.95,
                                risk.vars = NULL,
                                risk.direction = c("high", "low"),
                                uncertainty.threshold = 0.7,
                                seed = 1, verbose = TRUE) {

  .et_banner("soil_classification",
             "delineates soil zones with Gaussian mixtures", verbose)
  if (!requireNamespace("mclust", quietly = TRUE))
    stop("soil_classification() requires the 'mclust' package.\n",
         "  install.packages(\"mclust\")", call. = FALSE)

  # mclust::Mclust() calls mclustBIC() and other helpers WITHOUT namespace
  # qualification, so loading the namespace is not enough -- the package must be
  # attached or those calls fail with "could not find function mclustBIC".
  # Attach for the duration of this call, then restore the search path.
  if (!"package:mclust" %in% search()) {
    suppressPackageStartupMessages(attachNamespace(loadNamespace("mclust")))
    on.exit(try(detach("package:mclust", unload = FALSE, character.only = TRUE),
                silent = TRUE), add = TRUE)
  }

  risk.direction <- match.arg(risk.direction)

  if (!is.data.frame(soil.data))
    stop("'soil.data' must be a data.frame.", call. = FALSE)
  if (!env.id %in% names(soil.data))
    stop("Site column '", env.id, "' not found in soil.data.", call. = FALSE)
  if (anyDuplicated(soil.data[[env.id]]))
    stop("Duplicated values in '", env.id, "'. soil_classification() expects ",
         "one row per site.", call. = FALSE)
  if (!is.numeric(uncertainty.threshold) || length(uncertainty.threshold) != 1 ||
      uncertainty.threshold <= 0 || uncertainty.threshold >= 1)
    stop("'uncertainty.threshold' must lie in (0, 1).", call. = FALSE)
  if (!is.numeric(pca.var) || length(pca.var) != 1 ||
      pca.var <= 0 || pca.var > 1)
    stop("'pca.var' must lie in (0, 1].", call. = FALSE)

  if (!is.null(seed)) set.seed(seed)

  X <- .soil_features(soil.data, env.id, variables, depths,
                      convert.units, clr.texture, verbose)
  X_raw <- X

  if (anyNA(X)) {
    keep <- stats::complete.cases(X)
    if (!any(keep))
      stop("Every site has at least one missing soil value.", call. = FALSE)
    if (sum(keep) < 2)
      stop("Fewer than two sites have complete soil data.", call. = FALSE)
    # Warning, not message: silently clustering a subset of the sites the user
    # asked for is exactly the kind of thing that must not be quiet.
    warning("Dropping ", sum(!keep), " site(s) with missing soil values: ",
            paste(utils::head(as.character(soil.data[[env.id]])[!keep], 5),
                  collapse = ", "),
            if (sum(!keep) > 5) ", ..." else "", call. = FALSE)
    X <- X[keep, , drop = FALSE]
    soil.data <- soil.data[keep, , drop = FALSE]
    X_raw <- X_raw[keep, , drop = FALSE]
  }

  n <- nrow(X)
  if (n < 3)
    stop("Need at least 3 sites to fit a mixture model; got ", n, ".",
         call. = FALSE)

  if (isTRUE(scale.data)) X <- .safe_scale(X)

  pca <- NULL
  if (isTRUE(use.pca)) {
    pca <- stats::prcomp(X, center = TRUE, scale. = FALSE)
    cum <- cumsum(pca$sdev^2) / sum(pca$sdev^2)
    k   <- max(1L, which(cum >= pca.var)[1])
    X   <- pca$x[, seq_len(k), drop = FALSE]
    if (verbose)
      message("PCA: ", k, " component(s) retained (",
              round(100 * cum[k], 1), "% of variance).")
  }

  p <- ncol(X)
  if (n < 10 * p && verbose) {
    message("Note: ", n, " sites for ", p, " variables. A GMM estimates a ",
            "covariance per component;")
    message("      consider variables/depths filtering, use.pca = TRUE, or a ",
            "parsimonious model.names.")
  }

  G <- G[G >= 1 & G <= n - 1]
  if (!length(G))
    stop("No valid cluster count: with ", n, " sites, G must lie in 1..",
         n - 1, ".", call. = FALSE)

  if (verbose) {
    message(strrep("-", 63))
    message("soil_classification() - Gaussian mixture model")
    message(strrep("-", 63))
    message("Sites .............................. ", n)
    message("Features ........................... ", p)
    message("Candidate clusters ................. ", min(G), "-", max(G))
  }

  fit <- try(mclust::Mclust(X, G = G, modelNames = model.names,
                            verbose = FALSE), silent = TRUE)
  if (inherits(fit, "try-error") || is.null(fit)) {
    msg  <- if (inherits(fit, "try-error")) trimws(as.character(fit)) else "NULL model"
    hint <- if (grepl("could not find function", msg))
      paste0("\n  This looks like an mclust namespace problem. Try running ",
             "library(mclust) before soil_classification().")
    else ""
    stop("Mclust failed to fit any model. With few sites relative to ",
         "variables, try use.pca = TRUE or model.names = \"EII\".\n",
         "  Original error: ", msg, hint, call. = FALSE)
  }

  z  <- fit$z
  if (is.null(dim(z))) z <- matrix(z, ncol = 1)
  cl <- as.integer(fit$classification)
  Gf <- fit$G

  # ---- order clusters so labels mean something --------------------------
  score_mat <- if (!is.null(risk.vars)) {
    # Match against raw names, the property prefix, and the ilr texture label,
    # so risk.vars = "clay" still selects the texture coordinates at each depth.
    base_nm <- sub("\\|.*$", "", colnames(X_raw))
    is_tex  <- grepl("^ilr[0-9]+_texture", colnames(X_raw))
    sel <- which(colnames(X_raw) %in% risk.vars |
                 base_nm %in% risk.vars |
                 (is_tex & any(c("clay", "sand", "silt") %in% risk.vars)))
    if (!length(sel))
      stop("None of 'risk.vars' matched the clustered variables: ",
           paste(risk.vars, collapse = ", "), call. = FALSE)
    scale(X_raw[, sel, drop = FALSE])
  } else {
    scale(X_raw)
  }
  cl_score <- tapply(rowMeans(score_mat, na.rm = TRUE), cl, mean)
  if (risk.direction == "low") cl_score <- -cl_score

  cl_new <- .order_clusters(cl, as.numeric(cl_score))
  perm   <- order(as.numeric(cl_score))
  z      <- z[, perm, drop = FALSE]
  colnames(z) <- paste0("p_cluster", seq_len(Gf))
  rownames(z) <- as.character(soil.data[[env.id]])

  maxp <- apply(z, 1, max)

  classification <- data.frame(
    env           = soil.data[[env.id]],
    cluster       = cl_new,
    max_posterior = round(maxp, 4),
    uncertainty   = round(1 - maxp, 4),
    uncertain     = maxp < uncertainty.threshold,
    stringsAsFactors = FALSE
  )
  names(classification)[1] <- env.id

  # ---- per-cluster profiles in original units ---------------------------
  profiles <- data.frame(cluster = seq_len(Gf),
                         n_sites = as.integer(table(
                           factor(cl_new, levels = seq_len(Gf)))))
  means <- t(sapply(seq_len(Gf), function(k) {
    idx <- cl_new == k
    if (!any(idx)) rep(NA_real_, ncol(X_raw))
    else colMeans(X_raw[idx, , drop = FALSE], na.rm = TRUE)
  }))
  colnames(means) <- colnames(X_raw)
  profiles <- cbind(profiles, as.data.frame(means))

  if (verbose) {
    message("Selected model ..................... ", fit$modelName,
            " with ", Gf, " cluster(s)")
    message("BIC ................................ ", round(fit$bic, 2))
    message("Cluster sizes ...................... ",
            paste(paste0(seq_len(Gf), ":", profiles$n_sites), collapse = "  "))
    nu <- sum(classification$uncertain)
    if (nu)
      message("Uncertain sites (p < ",
              uncertainty.threshold, ") ......... ", nu, " of ", n)
    message("Done.")
  }

  structure(list(
    classification = classification,
    posterior      = z,
    profiles       = profiles,
    model          = list(G = Gf, modelName = fit$modelName,
                          bic = fit$bic, loglik = fit$loglik),
    BIC            = fit$BIC,
    features       = X_raw,
    pca            = pca,
    settings       = list(env.id = env.id, variables = variables,
                          depths = depths, G = G, model.names = model.names,
                          scale.data = scale.data, clr.texture = clr.texture,
                          convert.units = convert.units, use.pca = use.pca,
                          risk.vars = risk.vars,
                          risk.direction = risk.direction,
                          uncertainty.threshold = uncertainty.threshold,
                          seed = seed)
  ), class = "soil_gmm")
}


# =========================================================================
# SECTION 2 - print method
# =========================================================================

#' @title Print a soil_gmm object
#'
#' @description
#' Compact summary of a fitted soil mixture model: the selected model and
#' cluster count, the BIC gap to the runner-up, cluster sizes, the number of
#' uncertain sites, and the per-cluster profile means in original units.
#'
#' The BIC gap is reported because a point estimate of the zone count is
#' misleading on its own. Differences are graded weak (< 2), positive (2-6),
#' strong (6-10) and very strong (> 10). A weak gap means several partitions fit
#' about equally well.
#'
#' @param x an object of class \code{soil_gmm}.
#' @param ... ignored.
#' @return \code{x}, invisibly.
#' @examples
#' \dontrun{
#' fit <- soil_classification(soil_grid, risk.vars = "soc")
#' print(fit)
#' }
#' @export
print.soil_gmm <- function(x, ...) {
  cat("Gaussian mixture soil classification\n")
  cat("------------------------------------\n")
  cat("Sites           :", nrow(x$classification), "\n")
  cat("Features        :", ncol(x$features), "\n")
  cat("Model           : ", x$model$modelName, " with ", x$model$G,
      " cluster(s)\n", sep = "")
  cat("BIC             :", round(x$model$bic, 2), "\n")

  b <- x$BIC
  if (!is.null(b)) {
    vals <- sort(as.numeric(b[is.finite(b)]), decreasing = TRUE)
    if (length(vals) >= 2) {
      gap <- vals[1] - vals[2]
      lab <- if (gap < 2) "weak" else if (gap < 6) "positive" else
             if (gap < 10) "strong" else "very strong"
      cat("BIC gap to 2nd  :", round(gap, 2), " (", lab, ")\n", sep = "")
    }
  }

  cat("Cluster sizes   :",
      paste(paste0(x$profiles$cluster, "=", x$profiles$n_sites),
            collapse = "  "), "\n")

  nu <- sum(x$classification$uncertain)
  cat("Uncertain sites :", nu, "of", nrow(x$classification),
      "(max posterior <", x$settings$uncertainty.threshold, ")\n")

  if (nu) {
    cat("\nSites needing review:\n")
    amb <- x$classification[x$classification$uncertain, , drop = FALSE]
    print(utils::head(amb[order(amb$max_posterior), ], 10), row.names = FALSE)
  }

  cat("\nCluster profiles (original units):\n")
  print(x$profiles, row.names = FALSE)
  invisible(x)
}
