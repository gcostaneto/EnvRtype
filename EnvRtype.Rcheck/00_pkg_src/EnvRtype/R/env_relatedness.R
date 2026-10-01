#==================================================================================================
# env_relatedness.R
#
# Environmental relatedness analysis for reaction-norm and enviromic-assisted models.
#
# Module contents
#   env_kernel()  - environmental relatedness kernels (linear / Gaussian / arc-cosine deep)
#   env_expand()  - expand an environment-level W matrix onto the phenotypic records
#
# Internal helpers (not exported) are prefixed with '.env_'. The shared helper
# .env_as_matrix() lives in env_characterizaton_src.R and is available package-wide,
# since R sources every file in R/ into the same namespace.
#
# Author: G Costa-Neto and contributors
#==================================================================================================


# =========================================================================
# SECTION 1 - env_kernel()
# =========================================================================

#' @title Easily Building of Environmental Relatedness Kernels
#'
#' @description
#' Returns environmental kinships for reaction-norm models. The output is a list containing
#' \code{varCov}, the relatedness among environmental covariables, and \code{envCov}, the
#' relatedness among environments. Linear (GBLUP-like), nonlinear Gaussian and arc-cosine
#' deep kernels are supported.
#'
#' @author Germano Costa Neto
#'
#' @param env.data matrix. Environmental variables (or markers) per environment (or per
#'   genotype-environment combination), typically a \code{W_matrix} output. An envirotype
#'   frequency matrix from \code{\link{T_matrix}} is accepted identically (rows =
#'   environments); pass \code{is.scaled = TRUE} since \code{T_matrix} rows are already
#'   normalised compositions.
#' @param Y data.frame. Phenotypic data set containing environment id, genotype id and trait
#'   value. Only used when \code{merge = TRUE}.
#' @param is.scaled boolean. If the environmental data are already mean-centred and scaled
#'   (default \code{TRUE}), assuming \eqn{x \sim N(0,1)}.
#' @param sd.tol numeric. Maximum standard deviation tolerated in quality control. Columns
#'   above this value are eliminated. Default 10.
#' @param digits numeric. Number of digits used for rounding. Default 5.
#' @param tol numeric. Numerical tolerance. Default 1E-3.
#' @param merge boolean. If \code{TRUE}, the environmental covariables are merged with \code{Y}
#'   to build an \eqn{n \times n} kernel.
#' @param Z_E matrix. Model matrix for environments, used when \code{merge = TRUE}. If
#'   \code{NULL} it is built from \code{Y}.
#' @param stages vector (character). Names of each stage or time interval. If not \code{NULL},
#'   one kernel is produced per development stage.
#' @param env.id character. Identification of the environment column. Default \code{'env'}.
#' @param gaussian boolean. If \code{TRUE}, uses the Gaussian kernel parametrisation
#'   \eqn{envCov = exp(-h \cdot d / q)}.
#' @param h.gaussian numeric. Bandwidth \eqn{h} used when \code{gaussian = TRUE}.
#' @param deep.kernel boolean. If \code{TRUE}, an arc-cosine deep kernel is built from the
#'   environmental data instead of the linear/Gaussian kernel. Default \code{FALSE}.
#' @param deep.layers integer. Number of hidden layers when \code{deep.kernel = TRUE} (>= 1).
#' @param QC boolean. If \code{TRUE}, applies quality control based on \code{sd.tol}
#'   regardless of \code{is.scaled}. Default \code{FALSE}.
#' @param verbose boolean. If \code{TRUE} (default) prints quality-control messages.
#'
#' @return
#' A list with two elements: \code{varCov}, the covariable x covariable relatedness, and
#' \code{envCov}, the environment x environment relatedness. When \code{stages} is supplied,
#' each element is itself a named list of kernels, one per stage.
#'
#' @details
#' Three kernel methods are available.
#'
#' \strong{Linear (GB).} \eqn{K = WW' / \mathrm{trace}(WW')/n}, the environmental analogue of
#' the GBLUP kernel of VanRaden (2008).
#'
#' \strong{Gaussian (GK).} \eqn{K = \exp(-h\,d/q)}, with \eqn{d} the squared Euclidean distance
#' and \eqn{q} its median, so the bandwidth is scale-free.
#'
#' \strong{Arc-cosine deep (DK).} Emulates a deep neural network with \code{deep.layers} hidden
#' layers (Cuevas et al. 2019); no bandwidth parameter is required.
#'
#' Note that quality control is now controlled by the explicit \code{QC} argument. In earlier
#' versions it ran only when \code{is.scaled = FALSE}, so with the default
#' \code{is.scaled = TRUE} no filtering happened despite \code{sd.tol} being documented as
#' active.
#'
#' @examples
#' \donttest{
#' data('maizeYield'); data("maizeWTH")
#'
#' ## Environmental covariable matrix (environments x covariables)
#' W.cov <- W_matrix(env.data = maizeWTH[maizeWTH$daysFromStart < 100, ],
#'                   var.id = c("T2M", "T2M_MAX", "PRECTOT"), statistic = "mean")
#'
#' ## 1. Linear (co)variance kernel
#' K.lin <- env_kernel(env.data = W.cov, gaussian = FALSE)
#' round(K.lin$envCov, 2)
#'
#' ## 2. Nonlinear Gaussian kernel
#' K.gau <- env_kernel(env.data = W.cov, gaussian = TRUE)
#' round(K.gau$envCov, 2)
#'
#' ## 3. Arc-cosine deep kernel with 2 hidden layers
#' K.deep <- env_kernel(env.data = W.cov, deep.kernel = TRUE, deep.layers = 2)
#' round(K.deep$envCov, 2)
#'
#' ## 4. One kernel per development stage
#' stages   <- c('VE', 'V1_V6', 'V6_VT', 'VT_R1')
#' interval <- c(0, 7, 30, 65, 90)
#' W.stage  <- W_matrix(env.data = maizeWTH, var.id = c('FRUE', 'PETP'),
#'                      by.interval = TRUE, time.window = interval,
#'                      names.window = stages)
#' K.stage <- env_kernel(env.data = W.stage, stages = stages, gaussian = TRUE)
#' names(K.stage$envCov)
#'
#' ## 5. Feeding a reaction-norm model
#' KE <- list(W = env_kernel(env.data = W.cov)$envCov)
#' }
#'
#' @seealso \code{W_matrix}, \code{get_kernel}, \code{env_cluster},
#'   \code{\link{T_matrix}}
#'
#' @references
#' Cuevas J. et al. (2019). Deep kernel for genomic and near infrared predictions in
#' multi-environment breeding trials. \emph{G3} 9(9), 2913-2924.
#'
#' Costa-Neto G., Fritsche-Neto R., Crossa J. (2021). Nonlinear kernels, dominance, and
#' envirotyping data increase the accuracy of genome-based prediction in multi-environment
#' trials. \emph{Heredity} 126, 92-106.
#'
#' @importFrom stats sd dist model.matrix median cor
#' @export
env_kernel <- function(env.data, Y = NULL, is.scaled = TRUE, sd.tol = 10, digits = 5,
                       tol = 1E-3, merge = FALSE, Z_E = NULL, stages = NULL,
                       env.id = "env", gaussian = FALSE, h.gaussian = NULL,
                       deep.kernel = FALSE, deep.layers = 1,
                       QC = FALSE, verbose = TRUE) {

  .et_banner("env_kernel", "builds environmental relatedness kernels", verbose)
  if (isTRUE(gaussian) && isTRUE(deep.kernel))
    stop("Choose either 'gaussian = TRUE' or 'deep.kernel = TRUE', not both.",
         call. = FALSE)

  args <- list(Y = Y, is.scaled = is.scaled, sd.tol = sd.tol, tol = tol,
               merge = merge, Z_E = Z_E, env.id = env.id, gaussian = gaussian,
               h.gaussian = h.gaussian, deep.kernel = deep.kernel,
               deep.layers = deep.layers, QC = QC, verbose = verbose)

  if (is.null(stages)) {
    .et_step("computing a single environmental kernel", verbose)
    K_E <- do.call(.env_kernel_0, c(list(env.data = env.data), args))
    return(list(varCov = round(K_E$varCov, digits),
                envCov = round(K_E$envCov, digits)))
  }

  .et_step("computing one kernel per developmental stage", verbose)
  K_E <- K_W <- vector("list", length(stages))
  for (i in seq_along(stages)) {
    id <- grep(colnames(env.data), pattern = stages[i])
    if (!length(id))
      stop("No covariable matched stage '", stages[i], "' in the column names of env.data.",
           call. = FALSE)
    .k <- do.call(.env_kernel_0,
                  c(list(env.data = env.data[, id, drop = FALSE]), args))
    K_E[[i]] <- round(.k$envCov, digits)
    K_W[[i]] <- round(.k$varCov, digits)
  }
  names(K_E) <- names(K_W) <- stages
  list(varCov = K_W, envCov = K_E)
}

#' Single-kernel engine behind env_kernel()
#'
#' @inheritParams env_kernel
#' @return A list with \code{varCov} and \code{envCov}.
#' @keywords internal
#' @noRd
.env_kernel_0 <- function(env.data, Y = NULL, is.scaled = TRUE, sd.tol = 10,
                          tol = 1E-3, merge = FALSE, Z_E = NULL, env.id = "env",
                          gaussian = FALSE, h.gaussian = NULL,
                          deep.kernel = FALSE, deep.layers = 1,
                          QC = FALSE, verbose = TRUE) {

  env.data <- .env_as_matrix(env.data, "env.data")

  # ---- optional centring/scaling and quality control ----
  if (isFALSE(is.scaled)) {
    Amean <- sweep(env.data, 2, colMeans(env.data, na.rm = TRUE), "-")
    sdA   <- apply(Amean, 2, stats::sd, na.rm = TRUE)
    sdA[!is.finite(sdA) | sdA < tol] <- tol
    env.data <- sweep(Amean, 2, sdA, "/")
    QC <- TRUE
  }

  if (isTRUE(QC)) {
    sdA <- apply(env.data, 2, stats::sd, na.rm = TRUE)
    t   <- ncol(env.data)
    removed <- names(sdA[sdA > sd.tol | !is.finite(sdA)])
    if (length(removed)) {
      keep <- !colnames(env.data) %in% removed
      if (!any(keep))
        stop("Quality control removed every covariable; relax 'sd.tol'.", call. = FALSE)
      env.data <- env.data[, keep, drop = FALSE]
    }
    if (isTRUE(verbose)) {
      message(strrep("-", 48))
      message("Removed envirotype markers: ", length(removed), " from ", t)
      if (length(removed)) message(paste(removed, collapse = "\n"))
      message(strrep("-", 48))
    }
  }

  # ---- expand to n x n using the phenotypic table, if requested ----
  if (isTRUE(merge)) {
    if (is.null(Z_E)) {
      if (is.null(Y))
        stop("'merge = TRUE' requires either 'Y' or a model matrix 'Z_E'.", call. = FALSE)
      .DF <- data.frame(env = as.factor(Y[, env.id]))
      Z_E <- stats::model.matrix(~ 0 + env, .DF)
    }
    env.data <- Z_E %*% env.data
  }

  # ---- kernel construction ----
  if (isTRUE(gaussian)) {
    O <- .env_gaussian_kernel(env.data, h = h.gaussian)
    H <- .env_gaussian_kernel(t(env.data), h = h.gaussian)
  } else if (isTRUE(deep.kernel)) {
    nl <- max(1L, as.integer(deep.layers))
    O <- .env_arccos_kernel(env.data)
    H <- .env_arccos_kernel(t(env.data))
    if (nl > 1L) {
      O <- .env_arccos_deep(O, nl = nl - 1L)
      H <- .env_arccos_deep(H, nl = nl - 1L)
    }
  } else {
    O <- .env_GB_kernel(env.data)
    H <- .env_GB_kernel(t(env.data))
  }

  list(varCov = H, envCov = O)
}

#' Linear (GBLUP-like) environmental kernel
#'
#' \eqn{K = WW' / (\mathrm{trace}(WW') / n)}, the environmental analogue of the
#' VanRaden (2008) genomic relationship matrix used throughout \pkg{EnvRtype}.
#'
#' @param x numeric matrix (entities to relate in rows).
#' @return A symmetric relationship matrix.
#' @keywords internal
#' @noRd
.env_GB_kernel <- function(x) {
  x <- as.matrix(x)
  XX <- tcrossprod(x)
  denom <- sum(diag(XX)) / nrow(XX)
  if (!is.finite(denom) || denom == 0) denom <- 1
  XX / denom
}

#' First-layer arc-cosine kernel (Cho & Saul 2009; Cuevas et al. 2019)
#'
#' @param x numeric matrix (entities to relate in rows).
#' @return A symmetric arc-cosine (n = 1) kernel.
#' @keywords internal
#' @noRd
.env_arccos_kernel <- function(x) {
  x <- as.matrix(x)
  XX <- tcrossprod(x)
  nrm <- sqrt(diag(XX))
  denom <- outer(nrm, nrm)
  denom[denom == 0] <- .Machine$double.eps
  cost <- pmin(pmax(XX / denom, -1), 1)
  theta <- acos(cost)
  (1 / pi) * denom * (sin(theta) + (pi - theta) * cost)
}

#' Additional hidden layers of the arc-cosine deep kernel
#'
#' Applies the arc-cosine recursion \code{nl} further times to an existing
#' first-layer kernel, emulating a deep network of \code{nl + 1} layers.
#'
#' @param GC a first-layer arc-cosine kernel.
#' @param nl integer. Number of extra layers (>= 0).
#' @return The deep arc-cosine kernel.
#' @keywords internal
#' @noRd
.env_arccos_deep <- function(GC, nl = 1L) {
  for (i in seq_len(nl)) {
    nrm <- sqrt(diag(GC))
    denom <- outer(nrm, nrm)
    denom[denom == 0] <- .Machine$double.eps
    cost <- pmin(pmax(GC / denom, -1), 1)
    theta <- acos(cost)
    GC <- (1 / pi) * denom * (sin(theta) + (pi - theta) * cost)
  }
  GC
}

#' Gaussian kernel from a numeric matrix
#'
#' @param x numeric matrix (rows are the entities to relate).
#' @param h numeric bandwidth; 1 when \code{NULL}.
#' @return A Gaussian relationship matrix.
#' @keywords internal
#' @noRd
.env_gaussian_kernel <- function(x, h = NULL) {
  d <- as.matrix(stats::dist(x, upper = TRUE, diag = TRUE))^2
  q <- stats::median(d)
  if (!is.finite(q) || q == 0) q <- stats::quantile(d, .05, na.rm = TRUE)
  if (!is.finite(q) || q == 0) q <- 1
  if (is.null(h)) h <- 1
  exp(-h * d / q)
}


# =========================================================================
# SECTION 2 - env_expand()
# =========================================================================

#' Merge an environmental covariable matrix onto a phenotypic table
#'
#' Expands an environment-level \strong{W} matrix to the \eqn{n} observations of a
#' phenotypic data.frame, so that each row of the output corresponds to one
#' genotype-in-environment record.
#'
#' @param env.data matrix. Environment x covariable matrix with environments as rownames.
#' @param df.pheno data.frame. Phenotypic records containing the environment column.
#' @param env.id character. Name of the environment column in \code{df.pheno}.
#' @param skip integer. Number of leading columns of the merged table to drop
#'   (the phenotypic id columns). Default 3.
#' @param verbose boolean. If \code{TRUE} (default) prints a progress banner.
#'
#' @return A numeric matrix of \eqn{n} rows (observations) by \eqn{k} covariables.
#'
#' @examples
#' \donttest{
#' data("maizeYield"); data("maizeWTH")
#' W <- W_matrix(env.data = maizeWTH[maizeWTH$daysFromStart < 100, ],
#'               var.id = c("T2M", "PRECTOT"), statistic = "mean")
#' Wn <- env_expand(env.data = W, df.pheno = maizeYield, env.id = "env")
#' dim(Wn)
#' }
#'
#' @seealso \code{W_matrix}, \code{env_kernel}
#' @export
env_expand <- function(env.data, df.pheno, env.id = "env", skip = 3, verbose = TRUE) {
  .et_banner("env_expand", "expands W to phenotypic observations", verbose)
  df.pheno <- data.frame(df.pheno)
  env.data <- data.frame(env.data)
  env.data[[env.id]] <- rownames(env.data)
  .et_step("merging covariables onto phenotypic records", verbose)
  out <- merge(df.pheno, env.data, by = env.id)
  as.matrix(out[, -seq_len(skip), drop = FALSE])
}
