## ############################################################################
## env_prediction.R
##
## Envirotype-informed genomic prediction: kernel construction, model fitting,
## and prediction at environments that were never planted -- new LOCATIONS,
## new PLANTING DATES, or both.
##
## This file is a CONCATENATION of five tested source files. The code is
## byte-identical to the versions that passed testing; only section banners
## have been inserted. Nothing was re-typed.
## # Organized by F. Pontes
## ----------------------------------------------------------------------------
## TYPICAL WORKFLOW
## ----------------------------------------------------------------------------
##   K   <- get_kernel(K_G = list(G = K_A), K_E = list(W = K_En),
##                     data = df, model = 'RNMM')
##   fit <- kernel_model(y = 'value', data = df, random = K,
##                       env = 'env', gid = 'gid', keep_effects = TRUE)
##          #                                      ^^^^^^^^^^^^^^^^^^^
##          # REQUIRED for any scan: the scan needs per-draw effect chains,
##          # not posterior means, or its intervals would be fabricated.
##
##   ## new locations ...
##   sc  <- scan_untested_envs(fit, W_train = Wtr, W_new = Wnew, K_G = K_A)
##   map_scan(sc, coords, type = 'coverage')      # where can I predict?
##   map_scan(sc, coords, genotypes = c('G01'))   # what will it yield?
##
##   ## ... or new planting dates
##   plot_planting_window(sc, dates = sowing_dates)
##   best_planting_date(sc, dates = sowing_dates)
##
##   ## ... or both (one panel per date)
##   map_scan(sc, coords, date_col = 'date', facet_by = 'date',
##            genotypes = 'G01')
##
##   ## large grids: chunked, memory-bounded
##   gs <- grid_scan(fit, W_train = Wtr, W_new = W_grid, K_G = K_A)
##
## ----------------------------------------------------------------------------
## DEPENDENCIES
## ----------------------------------------------------------------------------
##   required : BGGE (fitting)
##   optional : ggplot2  -- all plots except engine = 'base'
##              sf       -- shapefile overlays in map_scan(shp = )
##   base R only otherwise. No rlang dependency.
##
## ----------------------------------------------------------------------------
## WHAT THE ENVELOPE DIAGNOSTICS DO AND DO NOT TELL YOU
## ----------------------------------------------------------------------------
## The scan reports WHERE a site sits geometrically relative to the training
## covariate cloud. That is NOT a statement that the prediction is wrong.
## Calibration against 60 held-out environments gave correlations with accuracy
## of only r ~ -0.37 to -0.66; an earlier predictable/unreliable verdict was
## removed because the thresholds could not support it. Treat 'interpolation'
## as geometry, never as a correctness guarantee.
##
## Position thresholds have NOT been retuned against real covariate data --
## see the caveat in PART 3. weight_negativity is the cleanest monotone signal
## observed in testing but was NOT part of the 60-environment calibration, so
## its relationship to accuracy is unknown.
##
## ----------------------------------------------------------------------------
## VERIFICATION STATUS AT TIME OF ASSEMBLY
## ----------------------------------------------------------------------------
##   verified : all numerical paths, envelope modes, chunking (chunked ==
##              unchunked), guards, base-graphics engine, and ggplot2 CALL
##              STRUCTURE via a recording mock.
##   NOT verified : actual ggplot2 rendering, and the sf overlay path, which
##              has never been executed. The development sandbox has no
##              ggplot2/sf and cannot install packages. Run
##              test_scan_maps_LOCAL.R on a machine with both to close this.
## ############################################################################



## ############################################################################
## PART 1 -- KERNEL CONSTRUCTION      get_kernel()
## ############################################################################

#' @title Envirotype-informed Kernels for Multi-Environment Genomic Prediction
#'
#' @description
#' Builds the genomic, envirotypic and genotype-by-environment kernels required
#' by \code{\link{kernel_model}}. Given a genomic relationship matrix
#' (\code{p x p} genotypes) and optionally an environmental relationship matrix
#' (\code{q x q} environments), it returns one kernel per model term, already
#' spectrally decomposed and ready to sample from.
#'
#' Two computational engines are available and produce numerically identical
#' results:
#' \describe{
#'   \item{\code{engine = "factored"} (default)}{Never allocates the
#'     \code{n x n} kernel. Each kernel is represented by its eigenvectors
#'     \code{U} (\code{n x r}) and eigenvalues \code{s}, obtained from the small
#'     \code{p x p} / \code{q x q} inputs. Memory is O(n*r) instead of O(n^2),
#'     which is what makes large multi-environment trials feasible.}
#'   \item{\code{engine = "dense"}}{Builds the explicit \code{n x n} matrix and
#'     attaches it as \code{$Kernel}. Slower and quadratic in memory, but needed
#'     if you want to inspect or post-process the kernel itself.}
#' }
#'
#' Set \code{verbose = TRUE} to trace each step with timings;
#' \code{verbose = FALSE} (default) is entirely silent.
#'
#' @section Model structures:
#' \tabular{lll}{
#'   \strong{model} \tab \strong{terms}                \tab \strong{kernels} \cr
#'   \code{MM}      \tab y = mu + G                    \tab KG_G \cr
#'   \code{MDs}     \tab y = mu + G + GE               \tab KG_G, KGE_GE \cr
#'   \code{EMM}     \tab y = mu + E + G                \tab KE_W, KG_G \cr
#'   \code{EMDs}    \tab y = mu + E + G + GE           \tab KE_W, KG_G, KGE_GE \cr
#'   \code{RNMM}    \tab y = mu + E + G + GxW          \tab KE_W, KG_G, KGE_GW \cr
#'   \code{RNMDs}   \tab y = mu + E + G + GE + GxW     \tab KE_W, KG_G, KGE_GE, KGE_GW
#' }
#' \code{GE} is a block-diagonal (within-environment) interaction; \code{GxW} is
#' the reaction-norm interaction built from the environmental kernel.
#'
#' @section Why the factored engine is exact:
#' Every kernel here is a congruence transform of a small matrix. With
#' \code{ig} mapping observations to genotypes and \code{ie} to environments,
#' \deqn{K_G = Z_g K_g Z_g' = K_g[ig, ig], \quad K_E = K_e[ie, ie]}
#' so the eigendecomposition of the \code{n x n} kernel follows from the small
#' one in O(p^3) rather than O(n^3). Writing \eqn{K = A A'} with
#' \eqn{A = V[ig,]\sqrt{D}}, the nonzero eigenpairs come from the small Gram
#' matrix \eqn{A'A}, giving an orthonormal \code{U} and exact eigenvalues.
#'
#' The interaction kernel is an elementwise product of two expansions, which is
#' the expansion of a Kronecker product: with \eqn{j = (ie-1) n_g + ig},
#' \eqn{K_{GE} = (K_e \otimes K_g)[j, j]}. On a complete balanced grid its
#' eigenpairs are closed-form (eigenvalues \eqn{d_e \otimes d_g}, eigenvectors a
#' Khatri-Rao product) so no large eigendecomposition is ever required.
#'
#' @section Original Version:
#' Costa-Neto et al (2021)
#'
#' EnvRtype v.1.2.3, Sep 2026
#'
#' @param K_G list of genomic relationship matrices (\code{p x p}), named. If
#'   \code{NULL}, an identity genotype kernel is used.
#' @param K_E list of environmental relationship matrices (\code{q x q}), named.
#'   If \code{NULL}, benchmark genomic-only models are built and \code{model} is
#'   restricted to \code{MM} / \code{MDs}.
#' @param data data.frame with the environment, genotype and phenotype columns.
#' @param model character. One of \code{"MM"}, \code{"MDs"}, \code{"EMM"},
#'   \code{"EMDs"}, \code{"RNMM"}, \code{"RNMDs"}. See \emph{Model structures}.
#' @param intercept.random logical. Append an identity genotype kernel
#'   \code{KG_Gi}, separating marker-predictable genetic variance from
#'   genotype-specific deviation the markers do not capture. Requires
#'   replication (a genotype seen in more than one environment) to be
#'   identifiable; with one observation per genotype it is confounded with the
#'   residual and a warning is issued. Default \code{FALSE}.
#' @param engine character. \code{"factored"} (default) or \code{"dense"}.
#' @param env,gid,y character. Column names in \code{data}.
#' @param tol numeric. Eigenvalue cutoff, relative to the largest eigenvalue.
#'   Must match \code{kernel_model(tol =)}. Default \code{1e-10}.
#' @param keep_var numeric in (0, 1]. Retain only the leading eigenvalues
#'   carrying this proportion of each kernel's trace. Sampler cost scales with
#'   total rank, so e.g. \code{0.99} is often a large speedup for negligible
#'   change in predictions. Default \code{1} (keep everything).
#' @param scale_kernels logical. Normalise each kernel to
#'   \code{mean(diag(K)) = 1} before decomposing. Must match
#'   \code{kernel_model(scale_kernels =)}. Default \code{TRUE}.
#' @param verbose logical. Narrate each step with timings. Default \code{FALSE}.
#'
#' @return A named list of kernels, class \code{"kernel_set"}. Each element is a
#'   list with:
#'   \describe{
#'     \item{\code{Type}}{\code{"D"} (dense) or \code{"BD"} (block-diagonal).}
#'     \item{\code{decomp}}{list of blocks, each
#'       \code{list(s, U, nr, deltav, type, pos, scaled, tol)}.}
#'     \item{\code{diag_mean}}{\code{mean(diag(K))}, computed in O(n*r).}
#'     \item{\code{Kernel}}{the \code{n x n} matrix -- only when
#'       \code{engine = "dense"}.}
#'   }
#'   Attributes: \code{"ne"} (observations per environment), \code{"kd_meta"}
#'   (decomposition settings, so \code{kernel_model()} can validate the cache),
#'   and \code{"obs"} (index vectors and level names).
#'
#' @section Reusing kernels:
#' Kernels depend only on genotypes and environments, never on the phenotype.
#' Masking \code{y} for a cross-validation fold does \strong{not} invalidate
#' them, so build once and reuse across folds, chains and models.
#'
#' @author Germano Costa Neto. Refactored version.
#'
#' @references
#' Costa-Neto G, Galli G, Carvalho HF, Crossa J, Fritsche-Neto R (2021).
#' EnvRtype: a software to interplay enviromics and quantitative genomics in
#' agriculture. \emph{G3}, 11(4), jkab040.
#'
#' @seealso \code{\link{kernel_model}}, \code{\link{decompose_kernels}}
#'
#' @examples
#' \dontrun{
#' data("maizeYield"); data("maizeG"); data("maizeWTH")
#'
#' ECs <- W_matrix(env.data = maizeWTH[maizeWTH$daysFromStart < 100, ],
#'                 var.id = c("FRUE", "PETP", "SRAD", "T2M_MAX"),
#'                 statistic = "mean")
#' KE <- list(W = env_kernel(env.data = ECs)[[2]])
#' KG <- list(G = maizeG)
#'
#' ## factored engine (default): no n x n matrix is ever allocated
#' K <- get_kernel(K_G = KG, K_E = KE, data = maizeYield, model = "RNMM",
#'                 env = "env", gid = "gid", y = "value", verbose = TRUE)
#'
#' fit <- kernel_model(y = "value", data = maizeYield, random = K,
#'                     env = "env", gid = "gid",
#'                     iterations = 5000, burnin = 1000, verbose = TRUE)
#'
#' ## trim the eigenvalue tail: usually a large speedup, tiny accuracy cost
#' K99 <- get_kernel(K_G = KG, K_E = KE, data = maizeYield, model = "RNMM",
#'                   keep_var = 0.99)
#'
#' ## separate marker-predictable from residual genetic variance
#' Ki <- get_kernel(K_G = KG, K_E = KE, data = maizeYield, model = "RNMM",
#'                  intercept.random = TRUE)
#' }
#'
#' @importFrom stats var
#' @export
get_kernel <- function(K_G = NULL,
                       K_E = NULL,
                       data = NULL,
                       model = "RNMM",
                       intercept.random = FALSE,
                       engine = c("factored", "dense"),
                       env = "env",
                       gid = "gid",
                       y   = "value",
                       tol = 1e-10,
                       keep_var = 1,
                       scale_kernels = TRUE,
                       verbose = FALSE) {

  engine <- match.arg(engine)
  model  <- match.arg(model, c("MM", "MDs", "EMM", "EMDs", "RNMM", "RNMDs"))

  ## ---- verbose step tracker -------------------------------------------------
  ## Mirrors kernel_model(): every call is a no-op when verbose is FALSE/0, so
  ## the function emits nothing at all.
  .v      <- as.numeric(verbose) != 0
  .t_all  <- proc.time()[3]
  .t_step <- .t_all
  .step <- function(txt) {
    if (!.v) return(invisible(NULL))
    .t_step <<- proc.time()[3]
    cat(sprintf("  %-52s", paste0(txt, " ...")))
    utils::flush.console()
    invisible(NULL)
  }
  .done <- function(extra = NULL) {
    if (!.v) return(invisible(NULL))
    cat(sprintf(" %7.2fs%s\n", proc.time()[3] - .t_step,
                if (is.null(extra)) "" else paste0("  ", extra)))
    invisible(NULL)
  }

  if (.v) {
    cat("\n=========================================================\n")
    cat(" get_kernel: envirotype-informed kernels for multi-\n")
    cat("   environment genomic prediction. Builds the G, E and\n")
    cat("   GxE relationship kernels for the chosen model and\n")
    cat("   returns them spectrally decomposed, ready for\n")
    cat("   kernel_model().\n")
    cat(" Original Version: Costa-Neto et al (2021)\n")
    cat(" EnvRtype v.1.2.3, Sep 2026\n")
    cat("=========================================================\n")
    cat(sprintf(" model = %-7s engine = %s\n", model, engine))
    cat("---------------------------------------------------------\n")
  }

  ## ---- validation -----------------------------------------------------------
  .step("[1/4] Validating inputs")
  if (is.null(data)) stop("`data` is required.")
  if (!is.data.frame(data)) stop("`data` must be a data.frame.")
  miss <- setdiff(c(env, gid, y), names(data))
  if (length(miss))
    stop("Columns not found in `data`: ", paste(miss, collapse = ", "))
  if (!is.null(keep_var) && (keep_var <= 0 || keep_var > 1))
    stop("`keep_var` must be in (0, 1].")

  Y <- data.frame(env = as.character(data[[env]]),
                  gid = as.character(data[[gid]]),
                  y   = data[[y]], stringsAsFactors = FALSE)

  gl <- sort(unique(Y$gid))
  el <- unique(Y$env)                  # first-appearance order, so that the
  ig <- match(Y$gid, gl)               # block-diagonal layout matches the data
  ie <- match(Y$env, el)
  n  <- nrow(Y); ng <- length(gl); q <- length(el)
  ne_blocks <- as.vector(table(factor(Y$env, levels = el)))

  ## Block-diagonal kernels are only genuinely block-diagonal when rows are
  ## grouped by environment. Silently producing the wrong model here was a
  ## defect in the original; we refuse instead.
  if (model %in% c("MDs", "EMDs", "RNMDs")) {
    if (any(duplicated(rle(ie)$values)))
      stop("Rows of `data` are not grouped by environment, but model '", model,
           "' builds a block-diagonal GxE kernel which requires contiguous\n",
           "  environment blocks. Fix: data <- data[order(data[['", env,
           "']], data[['", gid, "']]), ]")
  }
  if (model %in% c("EMM", "EMDs", "RNMM", "RNMDs") && is.null(K_E))
    stop("model '", model, "' requires `K_E`, but none was supplied.")
  .done(sprintf("%d obs | %d gid | %d env", n, ng, q))

  ## ---- align the small input matrices --------------------------------------
  .step("[2/4] Aligning relationship matrices")
  .align <- function(M, levels_, what) {
    M <- as.matrix(M)
    rn <- rownames(M)
    if (is.null(rn)) {
      if (nrow(M) != length(levels_))
        stop(what, " matrix has no dimnames and its dimension (", nrow(M),
             ") does not match the number of ", what, " (", length(levels_), ").")
      dimnames(M) <- list(levels_, levels_)
      return(M)
    }
    j <- match(levels_, rn)
    if (anyNA(j))
      stop("These ", what, " are missing from the supplied matrix: ",
           paste(utils::head(levels_[is.na(j)], 5), collapse = ", "),
           if (sum(is.na(j)) > 5) ", ..." else "")
    M[j, j, drop = FALSE]
  }

  if (!is.null(K_G) && !is.list(K_G)) K_G <- list(G = K_G)
  if (!is.null(K_E) && !is.list(K_E)) K_E <- list(W = K_E)
  KG <- if (!is.null(K_G)) .align(K_G[[1]], gl, "genotypes") else diag(ng)
  KE <- if (!is.null(K_E)) .align(K_E[[1]], el, "environments") else NULL
  if (model %in% c("MM", "MDs")) KE <- NULL
  nmG <- if (!is.null(K_G) && !is.null(names(K_G))) names(K_G)[1] else "G"
  nmE <- if (!is.null(K_E) && !is.null(names(K_E))) names(K_E)[1] else "W"
  .done(sprintf("K_G %d x %d%s", nrow(KG), ncol(KG),
                if (is.null(KE)) "" else sprintf(", K_E %d x %d", nrow(KE), ncol(KE))))

  ## ---- factorise ------------------------------------------------------------
  .step("[3/4] Factorising kernels")
  fk <- list()
  if (!is.null(KE))
    fk[[paste0("KE_", nmE)]] <- .fk_expand(KE, ie, tol, keep_var, scale_kernels)
  fk[[paste0("KG_", nmG)]]   <- .fk_expand(KG, ig, tol, keep_var, scale_kernels)
  if (model %in% c("MDs", "EMDs", "RNMDs"))
    fk[["KGE_GE"]] <- .fk_kron(KG, ig, diag(q), ie, tol, keep_var, scale_kernels)
  if (!is.null(KE) && model %in% c("RNMM", "RNMDs"))
    fk[[paste0("KGE_", nmG, nmE)]] <-
      .fk_kron(KG, ig, KE, ie, tol, keep_var, scale_kernels)

  if (isTRUE(intercept.random)) {
    fk[["KG_Gi"]] <- .fk_expand(diag(ng), ig, tol, keep_var, scale_kernels)
    if (n == ng)
      warning("intercept.random = TRUE with one observation per genotype: ",
              "KG_Gi is confounded with the residual and is not identifiable.")
  }
  tot_r <- sum(vapply(fk, function(k) k$r, numeric(1)))
  .done(sprintf("%d kernels, total rank %d", length(fk), tot_r))

  ## ---- assemble in kernel_model()'s expected format ------------------------
  ## The field list and the LENGTH of $pos are load-bearing: kernel_model()
  ## tests `!is.na(d$pos[1])` to decide between the dense and block-diagonal
  ## branches, so $pos must be a length-1 NA for Type "D".
  .step("[4/4] Assembling output")
  K <- lapply(names(fk), function(nm) {
    f <- fk[[nm]]; U <- f$U(); s <- f$s
    dm <- mean(rowSums((U^2) * rep(s, each = nrow(U))))
    if (!is.finite(dm) || dm <= 0) dm <- 1
    out <- list(Type = "D", diag_mean = dm,
                decomp = list(list(s = s, U = U, nr = length(s),
                                   deltav = 1 / s, type = "D", pos = NA,
                                   scaled = isTRUE(scale_kernels), tol = tol)))
    if (identical(engine, "dense")) out$Kernel <- tcrossprod(U * rep(sqrt(s), each = nrow(U)))
    out
  })
  names(K) <- names(fk)

  attr(K, "ne") <- ne_blocks
  attr(K, "kd_meta") <- list(tol = tol, scale_kernels = scale_kernels,
                             ne = ne_blocks, n = n,
                             keep_var = stats::setNames(rep(keep_var, length(K)),
                                                        names(K)),
                             engine = engine,
                             elapsed = proc.time()[3] - .t_all, version = 3L)
  attr(K, "obs") <- list(ig = ig, ie = ie, gl = gl, el = el, n = n)
  class(K) <- c("kernel_set", "list")
  .done(if (identical(engine, "dense"))
          sprintf("%d dense matrices built", length(K)) else "no n x n matrix built")

  if (.v) {
    cat("---------------------------------------------------------\n")
    cat(sprintf(" Kernels: %s\n", paste(names(K), collapse = ", ")))
    cat(sprintf(" Total rank %d / %d observations (%.1f%%)\n",
                tot_r, n, 100 * tot_r / n))
    if (identical(engine, "factored"))
      cat(sprintf(" Memory: ~%.1f MB factored vs %.1f MB dense (%.0fx less)\n",
                  8 * n * tot_r / 1e6, 8 * n^2 * length(K) / 1e6,
                  (n * length(K)) / max(tot_r, 1)))
    cat(sprintf(" Total running time: %.2fs\n", proc.time()[3] - .t_all))
    cat("=========================================================\n\n")
  }
  K
}


## ===========================================================================
## Internal factorisation engine
## ===========================================================================

#' Exact orthonormal decomposition of K = Kb[idx, idx]
#'
#' Never forms the \code{n x n} matrix. With \code{Kb = V D V'},
#' \code{K = A A'} for \code{A = V[idx,] sqrt(D)}; the nonzero eigenpairs of
#' \code{A A'} follow from the small Gram matrix \code{A'A = W L W'} as
#' eigenvalues \code{L} and eigenvectors \code{A W L^-1/2}.
#'
#' This routes through the Gram matrix rather than returning \code{V[idx,]}
#' directly. \code{V[idx,]} is \emph{not} orthonormal when \code{idx} repeats --
#' its columns have norm \code{sqrt(multiplicity)} -- and although
#' \code{K = U S U'} still reconstructs exactly, \code{kernel_model()} uses
#' \code{U} and \code{deltav = 1/s} separately, so the error does not cancel.
#'
#' @param Kb small symmetric matrix (\code{p x p}).
#' @param idx integer vector mapping observations to rows of \code{Kb}.
#' @param tol relative eigenvalue cutoff.
#' @param keep_var proportion of the trace to retain.
#' @param scale normalise to \code{mean(diag(K)) = 1}.
#' @return list with \code{s}, \code{r}, \code{n} and \code{U()}.
#' @keywords internal
.fk_expand <- function(Kb, idx, tol = 1e-10, keep_var = 1, scale = TRUE) {
  n <- length(idx)
  if (isTRUE(scale)) {
    f <- mean(diag(Kb)[idx])
    if (!is.finite(f) || f <= 0) f <- 1
    Kb <- Kb / f
  }
  e  <- eigen(Kb, symmetric = TRUE)
  rM <- sum(e$values > tol * max(e$values, 1))
  if (rM == 0L) stop("Base matrix has no eigenvalue above tol.")

  A <- e$vectors[idx, seq_len(rM), drop = FALSE] *
       rep(sqrt(e$values[seq_len(rM)]), each = n)

  eg <- eigen(crossprod(A), symmetric = TRUE)
  r  <- sum(eg$values > tol)
  if (r == 0L) stop("Expanded kernel is numerically zero.")
  s <- eg$values[seq_len(r)]

  if (!is.null(keep_var) && keep_var < 1) {
    rr <- which(cumsum(s) / sum(s) >= keep_var)[1]
    if (!is.na(rr)) { r <- max(2L, min(as.integer(rr), r)); s <- s[seq_len(r)] }
  }
  U <- A %*% (eg$vectors[, seq_len(r), drop = FALSE] *
              rep(1 / sqrt(s), each = rM))
  list(s = s, r = r, n = n, U = function() U)
}

#' Exact decomposition of a Hadamard (GxE) kernel via its Kronecker form
#'
#' The elementwise product of two expansions is the expansion of a Kronecker
#' product: with \code{j = (ie-1)*ng + ig}, \code{K = (Ke \%x\% Kg)[j, j]}.
#'
#' On a complete balanced grid the eigenpairs are closed-form -- eigenvalues
#' \code{outer(dg, de)} and eigenvectors the Khatri-Rao product of rows -- so no
#' large eigendecomposition is needed and truncation is applied \emph{before}
#' the eigenvectors are materialised. On unbalanced designs that construction is
#' not orthonormal, so the Gram route is used instead.
#'
#' @param Kg,Ke small genomic and environmental matrices.
#' @param ig,ie observation index vectors.
#' @param tol,keep_var,scale as in \code{.fk_expand}.
#' @return list with \code{s}, \code{r}, \code{n} and \code{U()}.
#' @keywords internal
.fk_kron <- function(Kg, ig, Ke, ie, tol = 1e-10, keep_var = 1, scale = TRUE) {
  eg <- eigen(Kg, symmetric = TRUE); ee <- eigen(Ke, symmetric = TRUE)
  rg <- sum(eg$values > tol * max(eg$values, 1))
  re <- sum(ee$values > tol * max(ee$values, 1))
  if (rg == 0L || re == 0L) stop("A factor of the GxE kernel is numerically zero.")
  ng <- nrow(Kg); nq <- nrow(Ke); n <- length(ig)

  sv <- as.vector(outer(eg$values[seq_len(rg)], ee$values[seq_len(re)]))
  gi <- rep(seq_len(rg), times = re); ei <- rep(seq_len(re), each = rg)
  o  <- order(sv, decreasing = TRUE); sv <- sv[o]; gi <- gi[o]; ei <- ei[o]
  keep <- sv > tol; sv <- sv[keep]; gi <- gi[keep]; ei <- ei[keep]
  if (!length(sv)) stop("GxE kernel is numerically zero.")

  if ((n == ng * nq) && !anyDuplicated((ie - 1) * ng + ig)) {
    ## balanced: closed form, truncate before building eigenvectors
    if (!is.null(keep_var) && keep_var < 1) {
      rr <- which(cumsum(sv) / sum(sv) >= keep_var)[1]
      if (!is.na(rr)) {
        rr <- max(2L, min(as.integer(rr), length(sv)))
        sv <- sv[seq_len(rr)]; gi <- gi[seq_len(rr)]; ei <- ei[seq_len(rr)]
      }
    }
    U <- eg$vectors[ig, gi, drop = FALSE] * ee$vectors[ie, ei, drop = FALSE]
    if (isTRUE(scale)) {
      f <- mean(rowSums((U^2) * rep(sv, each = n)))
      if (!is.finite(f) || f <= 0) f <- 1
      sv <- sv / f
    }
    return(list(s = sv, r = length(sv), n = n, U = function() U))
  }

  ## unbalanced: Khatri-Rao factor then Gram, to restore orthonormality
  Ag <- eg$vectors[ig, seq_len(rg), drop = FALSE] *
        rep(sqrt(eg$values[seq_len(rg)]), each = n)
  Ae <- ee$vectors[ie, seq_len(re), drop = FALSE] *
        rep(sqrt(ee$values[seq_len(re)]), each = n)
  A <- Ag[, rep(seq_len(rg), times = re), drop = FALSE] *
       Ae[, rep(seq_len(re), each  = rg), drop = FALSE]
  if (isTRUE(scale)) {
    f <- mean(rowSums(A^2))
    if (!is.finite(f) || f <= 0) f <- 1
    A <- A / sqrt(f)
  }
  eG <- eigen(crossprod(A), symmetric = TRUE)
  r  <- sum(eG$values > tol)
  if (r == 0L) stop("GxE kernel is numerically zero.")
  s <- eG$values[seq_len(r)]
  if (!is.null(keep_var) && keep_var < 1) {
    rr <- which(cumsum(s) / sum(s) >= keep_var)[1]
    if (!is.na(rr)) { r <- max(2L, min(as.integer(rr), r)); s <- s[seq_len(r)] }
  }
  U <- A %*% (eG$vectors[, seq_len(r), drop = FALSE] *
              rep(1 / sqrt(s), each = ncol(A)))
  list(s = s, r = r, n = n, U = function() U)
}


## ===========================================================================
## S3 methods and utilities
## ===========================================================================

#' Print a kernel set
#'
#' @param x object of class \code{"kernel_set"} from \code{\link{get_kernel}}.
#' @param ... ignored.
#' @export
print.kernel_set <- function(x, ...) {
  m <- attr(x, "kd_meta"); o <- attr(x, "obs")
  cat("<kernel_set>", length(x), "kernels |", m$n, "observations\n")
  if (!is.null(o))
    cat("  ", length(o$gl), " genotypes x ", length(o$el), " environments\n", sep = "")
  cat("  engine:", m$engine %||% "dense",
      "| scaled:", m$scale_kernels, "| tol:", format(m$tol, scientific = TRUE), "\n")
  tb <- data.frame(kernel = names(x),
                   rank = vapply(x, function(k)
                     sum(vapply(k$decomp, function(d) d$nr, numeric(1))), numeric(1)),
                   diag_mean = round(vapply(x, function(k)
                     k$diag_mean %||% NA_real_, numeric(1)), 4),
                   row.names = NULL)
  tb$pct_of_n <- sprintf("%.1f%%", 100 * tb$rank / m$n)
  print(tb, row.names = FALSE)
  cat("  total rank:", sum(tb$rank), "\n")
  invisible(x)
}


## ===========================================================================
## Spectral-decomposition machinery
##
## Dense kernels supplied by the user (list(Kernel = <n x n>, Type = "D"|"BD"))
## are turned into cached spectral decompositions here. Kernels built as a
## congruence transform of a SMALL matrix (KG = Kg[ig, ig], KE = Ke[ie, ie])
## carry that base matrix, so the eigendecomposition costs O(p^3) instead of
## O(n^3). The Hadamard GxE kernel destroys low-rank structure and falls back
## to a full eigen() unless the (genotype, environment) grid is complete, in
## which case the Kronecker spectrum is available in closed form.
##
## Factored kernels from get_kernel() arrive pre-decomposed (they carry
## $decomp / $diag_mean and no $Kernel) and never pass through here.
## ===========================================================================

#' Eigendecompose a symmetric matrix, dropping negligible eigenvalues.
#' eigen(symmetric = TRUE) returns DECREASING eigenvalues, so the retained set
#' is a leading prefix and we slice rather than which().
#' @keywords internal
.kd_eig <- function(K, tol = 1e-10) {
  ei <- eigen(K, symmetric = TRUE)
  r  <- sum(ei$values > tol)
  if (r == 0L)
    stop("Kernel has no eigenvalue above tol = ", tol, "; it is numerically zero.")
  s <- ei$values[seq_len(r)]
  list(s = s, U = ei$vectors[, seq_len(r), drop = FALSE], nr = r, deltav = 1 / s)
}

#' Exact spectral decomposition of K = M[idx, idx] without forming n x n eigen.
#'
#' M = V D V'  =>  K = (V[idx,]) D (V[idx,])'.
#' Let A = V[idx,] sqrt(D) so K = A A'. The nonzero eigenpairs of A A' follow
#' from the small Gram matrix A'A = W L W': eigenvalues L, eigenvectors A W L^-1/2.
#' Exact (not Nystrom): verified to 1e-14 against full eigen().
#' @keywords internal
.kd_eig_expanded <- function(M, idx, tol = 1e-10) {
  eiM <- eigen(M, symmetric = TRUE)
  rM  <- sum(eiM$values > tol)
  if (rM == 0L) stop("Base matrix has no eigenvalue above tol.")

  A <- eiM$vectors[idx, seq_len(rM), drop = FALSE] *
       rep(sqrt(eiM$values[seq_len(rM)]), each = length(idx))

  eiG <- eigen(crossprod(A), symmetric = TRUE)
  r   <- sum(eiG$values > tol)
  if (r == 0L) stop("Expanded kernel is numerically zero.")

  s <- eiG$values[seq_len(r)]
  U <- A %*% (eiG$vectors[, seq_len(r), drop = FALSE] * rep(1 / sqrt(s), each = rM))
  list(s = s, U = U, nr = r, deltav = 1 / s)
}

#' Exact spectral decomposition of a Hadamard GxE kernel via its Kronecker form.
#'
#' The GxE kernel is built as an ELEMENTWISE product of two expanded kernels:
#'   Kern[a,b] = Kg[ig[a], ig[b]] * Ke[ie[a], ie[b]]
#' The Hadamard product of two expansions IS the expansion of the Kronecker
#' product. With B = Ke (x) Kg and the composite index j = (ie - 1) * ng + ig,
#' Kern = B[j, j] exactly, so .kd_eig_expanded()'s logic applies -- we never
#' need the n x n eigen and never even form B. When the grid is complete the
#' eigenpairs are known in closed form (eigenvalues = d_e (x) d_g, eigenvectors
#' = rows of Ue (x) Ug); otherwise we go through the small Gram matrix.
#' @keywords internal
.kd_eig_hadamard <- function(Kg, ig, Ke, ie, tol = 1e-10, keep_var = 1,
                             verbose = FALSE) {
  eg <- eigen(Kg, symmetric = TRUE)
  ee <- eigen(Ke, symmetric = TRUE)
  rg <- sum(eg$values > tol * max(eg$values, 1))
  re <- sum(ee$values > tol * max(ee$values, 1))
  if (rg == 0L || re == 0L)
    stop("A factor of the GxE kernel is numerically zero.")

  ng <- nrow(Kg); nq <- nrow(Ke); n <- length(ig)

  ## FAST PATH: complete balanced grid -- closed-form Kronecker spectrum.
  key <- (ie - 1) * ng + ig
  balanced <- (n == ng * nq) && !anyDuplicated(key)

  if (balanced) {
    sv <- as.vector(outer(eg$values[seq_len(rg)], ee$values[seq_len(re)]))
    gi <- rep(seq_len(rg), times = re)     # which Ug column
    ei <- rep(seq_len(re), each  = rg)     # which Ue column
    ord <- order(sv, decreasing = TRUE)
    sv <- sv[ord]; gi <- gi[ord]; ei <- ei[ord]

    keep <- sv > tol
    if (!any(keep)) stop("GxE kernel is numerically zero.")
    sv <- sv[keep]; gi <- gi[keep]; ei <- ei[keep]

    ## truncate FIRST -- avoids building eigenvectors we would discard
    if (!is.null(keep_var) && keep_var < 1) {
      cum <- cumsum(sv) / sum(sv)
      rr  <- which(cum >= keep_var)[1]
      if (!is.na(rr)) {
        rr <- max(2L, min(as.integer(rr), length(sv)))
        sv <- sv[seq_len(rr)]; gi <- gi[seq_len(rr)]; ei <- ei[seq_len(rr)]
      }
    }
    r <- length(sv)
    if (isTRUE(verbose))
      cat(sprintf("[closed-form Kronecker: %d x %d -> rank %d, no large eigen] ",
                  rg, re, r))

    ## U[a, k] = Ug[ig[a], gi[k]] * Ue[ie[a], ei[k]]  -- orthonormal already
    U <- eg$vectors[ig, gi, drop = FALSE] * ee$vectors[ie, ei, drop = FALSE]
    return(list(s = sv, U = U, nr = r, deltav = 1 / sv,
                truncated = !is.null(keep_var) && keep_var < 1 && r < rg*re))
  }

  ## GENERAL PATH: unbalanced / repeated cells. Build the Khatri-Rao factor A
  ## and go through the small Gram matrix. Exact, but costs O((rg*re)^3).
  if (isTRUE(verbose)) cat("[unbalanced: Gram route] ")
  Ag <- eg$vectors[ig, seq_len(rg), drop = FALSE] *
        rep(sqrt(eg$values[seq_len(rg)]), each = length(ig))
  Ae <- ee$vectors[ie, seq_len(re), drop = FALSE] *
        rep(sqrt(ee$values[seq_len(re)]), each = length(ie))
  A <- Ag[, rep(seq_len(rg), times = re), drop = FALSE] *
       Ae[, rep(seq_len(re), each  = rg), drop = FALSE]

  eiG <- eigen(crossprod(A), symmetric = TRUE)
  r   <- sum(eiG$values > tol)
  if (r == 0L) stop("GxE kernel is numerically zero.")
  s <- eiG$values[seq_len(r)]
  U <- A %*% (eiG$vectors[, seq_len(r), drop = FALSE] *
              rep(1 / sqrt(s), each = ncol(A)))
  list(s = s, U = U, nr = r, deltav = 1 / s)
}

#' Truncate a decomposition to the leading eigenvalues carrying `keep_var` of
#' the trace. The Gibbs inner loop is O(T * k * n * r), so dropping the flat
#' eigenvalue tail speeds up EVERY iteration of EVERY fit.
#' @keywords internal
.kd_truncate <- function(d, keep_var) {
  if (is.null(keep_var) || keep_var >= 1) return(d)
  tot <- sum(d$s)
  if (!is.finite(tot) || tot <= 0) return(d)
  r <- which(cumsum(d$s) / tot >= keep_var)[1]
  if (is.na(r)) return(d)
  r <- max(2L, min(as.integer(r), d$nr))
  if (r >= d$nr) return(d)
  d$s <- d$s[seq_len(r)]
  d$U <- d$U[, seq_len(r), drop = FALSE]
  d$nr <- r
  d$deltav <- 1 / d$s
  d$truncated <- TRUE
  d
}

#' Decompose one kernel element (dense or block-diagonal).
#' Scaling is applied BEFORE decomposition so it matches kernel_model().
#' @keywords internal
.kd_decompose_one <- function(k, ne = NULL, tol = 1e-10, scale_kernel = TRUE,
                              keep_var = 1, verbose = FALSE) {
  K <- k$Kernel

  f <- 1
  if (isTRUE(scale_kernel)) {
    f <- mean(diag(K))
    if (!is.finite(f) || f <= 0) f <- 1
  }

  if (identical(k$Type, "D")) {
    d <- if (!is.null(k$kron)) {
      ## GxE Hadamard: closed-form Kronecker spectrum when the grid is
      ## complete. keep_var is passed IN so truncation happens before the
      ## eigenvectors are materialised, not after.
      .kd_eig_hadamard(k$kron$Kg / f, k$kron$ig,
                       k$kron$Ke,     k$kron$ie, tol = tol,
                       keep_var = keep_var, verbose = verbose)
    } else if (!is.null(k$base) && !is.null(k$idx)) {
      .kd_eig_expanded(k$base / f, k$idx, tol = tol)     # low-rank route
    } else {
      .kd_eig(K / f, tol = tol)                          # full eigen
    }
    ## .kd_truncate is a no-op when the fast path already truncated
    if (!isTRUE(d$truncated)) d <- .kd_truncate(d, keep_var)
    d$type <- "D"; d$pos <- NA
    d$scaled <- isTRUE(scale_kernel); d$tol <- tol
    return(list(d))
  }

  if (identical(k$Type, "BD")) {
    if (is.null(ne) || length(ne) <= 1)
      stop("Type 'BD' requires `ne`, the size of each diagonal block.")
    posf <- cumsum(ne)
    posi <- cumsum(c(1, ne[-length(ne)]))
    out <- vector("list", length(ne)); cnt <- 0L
    for (j in seq_along(ne)) {
      rng <- posi[j]:posf[j]
      d <- tryCatch(.kd_eig(K[rng, rng, drop = FALSE] / f, tol = tol),
                    error = function(e) NULL)
      if (is.null(d)) next
      d <- .kd_truncate(d, keep_var)
      cnt <- cnt + 1L
      d$type <- "BD"; d$pos <- c(posi[j], posf[j])
      d$scaled <- isTRUE(scale_kernel); d$tol <- tol
      out[[cnt]] <- d
    }
    if (cnt == 0L) stop("All blocks of a BD kernel were numerically zero.")
    return(out[seq_len(cnt)])
  }

  stop("Kernel `Type` must be 'D' or 'BD'; got '", k$Type, "'.")
}

#' Attach spectral decompositions to a list of dense kernels
#'
#' Turns each dense kernel (\code{list(Kernel = <n x n>, Type = "D"|"BD")}) into
#' a cached spectral decomposition stored on \code{$decomp}. Call once and reuse
#' across any number of \code{\link{kernel_model}} fits: kernels depend only on
#' genotypes and environments -- never on the phenotype vector -- so masking
#' \code{y} for a cross-validation fold does NOT invalidate the decomposition.
#'
#' Factored kernels from \code{\link{get_kernel}} arrive pre-decomposed and do
#' not need this function.
#'
#' @param K list of kernels.
#' @param ne per-environment observation counts; defaults to \code{attr(K, "ne")}.
#' @param tol eigenvalue cutoff; must match \code{kernel_model(tol=)}.
#' @param scale_kernels must match \code{kernel_model(scale_kernels=)}.
#' @param keep_var numeric in (0, 1]. Retain only the leading eigenvalues
#'   carrying this proportion of each kernel's trace. 1 (default) keeps
#'   everything. A named vector or list keyed by kernel name / regex applies
#'   per kernel, e.g. \code{c(KGE_ = 0.95)}.
#' @param verbose report rank, timing and route per kernel.
#' @return \code{K} with \code{$decomp} on each element and a \code{"kd_meta"} attribute.
#' @seealso \code{\link{get_kernel}}, \code{\link{kernel_model}}, \code{\link{undecompose_kernels}}
#' @export
decompose_kernels <- function(K, ne = NULL, tol = 1e-10,
                              scale_kernels = TRUE, keep_var = 1,
                              verbose = TRUE) {
  if (!is.list(K) || !length(K)) stop("`K` must be a non-empty list of kernels.")
  if (is.null(ne)) ne <- attr(K, "ne")

  types <- vapply(K, function(k) k$Type %||% NA_character_, character(1))
  if (any(types == "BD", na.rm = TRUE) && (is.null(ne) || length(ne) <= 1))
    stop("Kernels of Type 'BD' are present but `ne` (block sizes) was not ",
         "supplied and is not stored on the kernel list.")

  ## resolve keep_var to one value per kernel
  kv <- rep(1, length(K)); names(kv) <- names(K)
  if (!is.null(keep_var)) {
    if (is.null(names(keep_var)) && length(keep_var) == 1L) {
      kv[] <- keep_var
    } else {
      for (pat in names(keep_var)) {
        hit <- grep(pat, names(K))
        if (length(hit)) kv[hit] <- keep_var[[pat]]
      }
    }
  }
  if (any(kv <= 0 | kv > 1)) stop("`keep_var` must be in (0, 1].")

  keep_attr <- attributes(K)[c("ne")]
  t0 <- proc.time()[3]
  if (isTRUE(verbose)) {
    cat("\n[3/3] Spectral decomposition\n")
    cat("      route: 'kronecker' = via the two small factors (fastest)\n")
    cat("             'low-rank'  = via the p x p base matrix\n")
    cat("             'full eigen'= dense n x n (slowest; unavoidable)\n")
  }
  for (j in seq_along(K)) {
    tj <- proc.time()[3]
    route <- if (!is.null(K[[j]]$kron)) "kronecker"
             else if (!is.null(K[[j]]$base)) "low-rank"
             else "full eigen"
    n_j <- nrow(K[[j]]$Kernel)
    if (isTRUE(verbose))
      cat(sprintf("      %-12s dim %5d  decomposing via %-11s ...",
                  names(K)[j], n_j, route))
    K[[j]]$decomp <- .kd_decompose_one(K[[j]], ne = ne, tol = tol,
                                       scale_kernel = scale_kernels,
                                       keep_var = kv[j], verbose = verbose)
    if (isTRUE(verbose)) {
      r <- sum(vapply(K[[j]]$decomp, function(d) d$nr, numeric(1)))
      trunc <- any(vapply(K[[j]]$decomp,
                          function(d) isTRUE(d$truncated), logical(1)))
      cat(sprintf(" rank %5d / %d  (%.1f%%)  %6.2fs%s\n",
                  r, n_j, 100*r/n_j, proc.time()[3] - tj,
                  if (trunc) sprintf("  [truncated @ %.4f]", kv[j]) else ""))
    }
  }
  attr(K, "ne") <- keep_attr$ne %||% ne
  attr(K, "kd_meta") <- list(tol = tol, scale_kernels = scale_kernels,
                             ne = ne, n = nrow(K[[1]]$Kernel),
                             keep_var = kv,
                             elapsed = proc.time()[3] - t0, version = 2L)
  if (isTRUE(verbose)) {
    tot_r <- sum(vapply(K, function(k)
      sum(vapply(k$decomp, function(d) d$nr, numeric(1))), numeric(1)))
    cat(sprintf("      -> %d kernels decomposed in %.2fs | total rank %d\n",
                length(K), proc.time()[3] - t0, tot_r))
    cat(sprintf("      Gibbs cost per iteration scales with total rank.\n"))
  }
  K
}

#' Are cached decompositions valid for the requested settings?
#'
#' A mismatch in \code{tol} or \code{scale_kernels} silently changes the model,
#' so \code{kernel_model()} re-decomposes rather than trust a stale cache.
#'
#' @param K a kernel list.
#' @param tol,scale_kernels,n settings to validate against.
#' @return logical.
#' @keywords internal
.kd_valid <- function(K, tol, scale_kernels, n) {
  m <- attr(K, "kd_meta")
  if (is.null(m)) return(FALSE)
  if (!isTRUE(all.equal(m$tol, tol))) return(FALSE)
  if (!identical(isTRUE(m$scale_kernels), isTRUE(scale_kernels))) return(FALSE)
  if (!identical(as.integer(m$n), as.integer(n))) return(FALSE)
  all(vapply(K, function(k) !is.null(k$decomp), logical(1)))
}

#' Drop cached decompositions
#'
#' Use after modifying a kernel in place, to force recomputation.
#'
#' @param K a kernel list from \code{\link{get_kernel}}.
#' @return \code{K} with decompositions and metadata removed.
#' @export
undecompose_kernels <- function(K) {
  for (j in seq_along(K)) K[[j]]$decomp <- NULL
  attr(K, "kd_meta") <- NULL
  K
}

`%||%` <- function(a, b) if (is.null(a)) b else a


## ############################################################################
## PART 2 -- MODEL FITTING            kernel_model() and friends
## ############################################################################

#' @title Kernel Models for Predicting Phenotypes across Multi-Environment Conditions
#'
#' @name kernel_model
#' @rdname kernel_model
#'
#' @description
#' Fits Bayesian linear mixed models for multi-environment trials (MET) using
#' the kernels produced by \code{get_kernel()}. The model is
#' \deqn{y = 1\mu + X\beta + \sum_j u_j + e}
#' with \eqn{u_j \sim N(0, K_j \sigma^2_j)} and \eqn{e \sim N(0, I\sigma^2_e)}.
#' Sampling runs in the spectral basis of each kernel, so every kernel is used
#' through its eigendecomposition rather than as an n x n matrix.
#'
#' The function returns genotype predictions, per-kernel variance components,
#' broad- and narrow-sense heritabilities, and -- with
#' \code{keep_effects = TRUE} -- posterior effect chains and breeder-facing
#' report tables (per-environment rankings, stability, selection candidates).
#'
#' Set \code{verbose = TRUE} to trace the six execution steps with per-step and
#' total timings; \code{verbose = FALSE} (default) is entirely silent.
#'
#' @section Original Version:
#' Costa-Neto et al (2021)
#'
#' EnvRtype v.1.2.3, Sep 2026
#'
#' @author Germano Costa Neto. Gibbs sampler adapted from Granato et al. (2018),
#'   BGGE package. Refactored version.
#'
#' @references
#' Costa-Neto G, Galli G, Carvalho HF, Crossa J, Fritsche-Neto R (2021).
#' EnvRtype: a software to interplay enviromics and quantitative genomics in
#' agriculture. \emph{G3}, 11(4), jkab040.
#'
#' Granato I, Cuevas J, Luna-Vazquez F, Crossa J, Montesinos-Lopez O,
#' Burgueno J, Fritsche-Neto R (2018). BGGE: a new package for genomic-enabled
#' prediction incorporating genotype x environment interaction models.
#' \emph{G3}, 8(9), 3039-3047.
#'
#' @section Performance changes vs. the original:
#' \itemize{
#'   \item \strong{Cached eigendecomposition.} If \code{random} carries \code{$decomp}
#'     (attached by \code{get_kernel(decompose = TRUE)}) and the cache matches
#'     \code{tol} / \code{scale_kernels} / \code{n}, \code{eigen()} is skipped
#'     entirely. Cross-validation folds and multiple chains reuse one decomposition.
#'   \item \strong{No stored transposes.} The original kept both \code{U} and
#'     \code{t(U)} per kernel (doubling eigenvector memory) and computed
#'     \code{crossprod(tU, b)}. Now only \code{U} is stored and the back-transform
#'     is \code{U \%*\% b}.
#'   \item \strong{Optional chain storage.} Full per-draw effect matrices
#'     (\code{k x nCum x n} doubles) are allocated only when needed
#'     (\code{keep_effects = TRUE} or \code{yHat.CI}). Otherwise running sums
#'     reduce memory from O(k*nCum*n) to O(k*n).
#'   \item \strong{Vectorised post-processing.} Row-wise \code{apply()} over draw
#'     matrices replaced by \code{rowSums}/\code{rowMeans} identities; \code{sweep()}
#'     replaced by recycled arithmetic.
#'   \item \strong{Scalar chains in a matrix} rather than a nested list, removing
#'     per-iteration copy-on-modify of the chain container.
#' }
#'
#' @param y character. Name of the phenotype column in \code{data}.
#' @param data data.frame with the environment, genotype and phenotype columns.
#'   Rows must be grouped by environment if any kernel is block-diagonal.
#' @param random list of kernels from \code{\link{get_kernel}}. Each element is a
#'   list carrying \code{$decomp} (the spectral decomposition) and \code{$Type}
#'   (\code{"D"} or \code{"BD"}); \code{$Kernel} is optional and is not required
#'   by the sampler. Factored kernels (no \code{$Kernel}) and dense kernels may
#'   be mixed in one list. Names should follow the \code{KG_} / \code{KE_} /
#'   \code{KGE_} convention, which drives the variance-component roles.
#' @param fixed matrix (\code{n x p}) of fixed effects, or \code{NULL}.
#' @param env,gid character. Environment and genotype column names in \code{data}.
#' @param verbose logical or 0/1. \code{TRUE} traces the six execution steps with
#'   per-step and total timings; \code{FALSE} (default) is entirely silent --
#'   no \code{cat()}, no \code{message()}. Results are identical either way.
#' @param iterations integer. Total Gibbs iterations. Default \code{1000}.
#' @param burnin integer. Iterations discarded before collecting draws.
#'   Default \code{200}.
#' @param thining integer. Keep every \code{thining}-th post-burn-in draw.
#'   Default \code{10}.
#' @param tol numeric. Eigenvalue cutoff. Must match \code{get_kernel(tol =)} or
#'   the cached decomposition is rejected and recomputed. Default \code{1e-10}.
#' @param R2 numeric in (0, 1). Prior proportion of phenotypic variance explained
#'   by the random terms. Default \code{0.5}.
#' @param digits integer. Rounding for the variance-component tables.
#' @param seed integer or \code{NULL}. Sets the RNG state for reproducibility.
#' @param probs length-2 numeric. Quantiles for posterior prediction intervals.
#'   Default \code{c(0.025, 0.975)}.
#' @param variance \code{"realized"} (default) or \code{"scale"}. See \emph{Details}.
#' @param scale_kernels logical. Normalise each kernel to \code{mean(diag(K)) = 1}.
#'   Must match \code{get_kernel(scale_kernels =)} for the cache to be reused.
#'   Default \code{TRUE}.
#' @param keep_effects logical. Retain per-draw effect chains. Required for
#'   \code{varcomp_method = "postmean"} and for the breeder-facing report tables
#'   (\code{$predictions}, \code{$genotype_effects}, \code{$env_summary}).
#'   Costs O(k * draws * n) memory. Default \code{FALSE}.
#' @param compute_CI logical. Compute posterior prediction intervals for
#'   \code{yHat}. Requires storing effect chains. Default \code{TRUE}.
#' @param R2_split \code{"per_kernel"} (default) or \code{"shared"}. Controls how
#'   the prior variance \code{R2} is distributed across kernels. See
#'   \emph{Prior parameterisation}.
#' @param prior_scale \code{"eigen"} (default) or \code{"diag"}. Whether the prior
#'   scale is divided by the mean eigenvalue of each kernel, making kernels of
#'   differing rank comparable. See \emph{Prior parameterisation}.
#' @param varcomp_method \code{"postmean"} (default) or \code{"gibbs"}. The
#'   estimator used for \code{$varcomp}. \code{"postmean"} needs
#'   \code{keep_effects = TRUE}; when unavailable it falls back to
#'   \code{"gibbs"} (reported in the trace when \code{verbose = TRUE}).
#'
#' @details
#' Fits \deqn{y = 1\mu + X\beta + \sum_j u_j + e}
#' with \eqn{u_j \sim N(0, K_j \sigma^2_{u_j})} and
#' \eqn{e \sim N(0, I\sigma^2_e)}, by Gibbs sampling in the spectral basis of
#' each kernel. Writing \eqn{K_j = U_j S_j U_j'}, the effects are sampled as
#' \eqn{u_j = U_j b_j} with \eqn{b_j \sim N(0, S_j \sigma^2_{u_j})}, so the cost
#' per iteration scales with the retained \strong{rank} of each kernel, not with
#' \eqn{n^2}. Trimming the eigenvalue tail via \code{get_kernel(keep_var =)} is
#' therefore the most direct way to speed up sampling.
#'
#' \code{variance = "scale"} reports the raw kernel-scale parameter
#' \code{sigb[j]} (original BGGE behaviour, inflated for low-rank incidence
#' kernels such as \code{KE_}). \code{variance = "realized"} reports
#' \code{var(u_j)} averaged over draws -- what most users mean by a variance
#' component. The raw scale is always returned as \code{$sigb}.
#'
#' @section Prior parameterisation:
#' Two defects in the original parameterisation are addressed by arguments:
#'
#' \strong{\code{R2_split}}. The original gave \emph{every} kernel the full
#' \code{R2} share, so with \code{nk} kernels the priors collectively assert
#' \code{nk * R2} of the phenotypic variance -- with \code{RNMDs} (4 kernels,
#' \code{R2 = 0.5}) that is twice \code{var(y)}. \code{"per_kernel"} (default)
#' uses \code{R2 / nk} so total prior mass stays at \code{R2} regardless of model
#' size, keeping model comparisons on equal footing. \code{"shared"} restores the
#' original behaviour.
#'
#' \strong{\code{prior_scale}}. Kernels normalised to \code{mean(diag(K)) = 1}
#' are equal in scale but not in effective dimension: a rank-5 and a rank-750
#' kernel with the same mean diagonal carry very different total variance.
#' \code{"eigen"} (default) divides the prior scale by \code{mean(s)}
#' (= trace / rank), making priors comparable across kernels of differing rank.
#'
#' @return An object of class \code{"kernel_model"}, a list with:
#' \describe{
#'   \item{\code{yHat}}{numeric vector of fitted/predicted values, length n.}
#'   \item{\code{yHat.CI}}{matrix of posterior prediction intervals
#'     (\code{NULL} unless \code{compute_CI = TRUE}).}
#'   \item{\code{varE}}{posterior mean residual variance.}
#'   \item{\code{random}}{per-kernel posterior effects and variances.}
#'   \item{\code{fit}}{the full sampler object: chains, kernel names, MCMC
#'     settings and -- with \code{keep_effects = TRUE} -- \code{Uchain}.
#'     (Named \code{$BGGE} in versions before v1.2.3.)}
#'   \item{\code{meta}}{the aligned env / gid / y data frame.}
#'   \item{\code{VarComp}}{variance-component table with confidence intervals.}
#'   \item{\code{varcomp}}{validated components and heritabilities, class
#'     \code{"km_varcomp"}; see \code{varcomp_summary}.}
#'   \item{\code{varcomp_gibbs}}{the same computed by the Gibbs-chain estimator,
#'     for comparison.}
#'   \item{\code{runtime}}{sampler wall time in seconds.}
#' }
#' With \code{keep_effects = TRUE} it additionally carries
#' \code{predictions}, \code{genetic_cor_env}, \code{genotype_effects},
#' \code{env_summary} and \code{variance_proportions}.
#'
#' @examples
#' \dontrun{
#' data("maizeYield"); data("maizeG"); data("maizeWTH")
#'
#' ECs <- W_matrix(env.data = maizeWTH[maizeWTH$daysFromStart < 100, ],
#'                 var.id = c("FRUE", "PETP", "SRAD", "T2M_MAX"),
#'                 statistic = "mean")
#' K <- get_kernel(K_G = list(G = maizeG),
#'                 K_E = list(W = env_kernel(env.data = ECs)[[2]]),
#'                 data = maizeYield, model = "RNMM",
#'                 env = "env", gid = "gid", y = "value")
#'
#' ## silent fit
#' fit <- kernel_model(y = "value", data = maizeYield, random = K,
#'                     env = "env", gid = "gid",
#'                     iterations = 5000, burnin = 1000)
#' fit
#' fit$varcomp
#'
#' ## traced fit, with effect chains for the report tables
#' fit2 <- kernel_model(y = "value", data = maizeYield, random = K,
#'                      env = "env", gid = "gid",
#'                      iterations = 5000, burnin = 1000,
#'                      keep_effects = TRUE, verbose = TRUE)
#' head(fit2$genotype_effects)
#'
#' ## leave-one-environment-out: kernels do not depend on y, so build once
#' for (e in unique(maizeYield$env)) {
#'   tr <- maizeYield; tr$value[tr$env == e] <- NA
#'   f <- kernel_model(y = "value", data = tr, random = K,
#'                     env = "env", gid = "gid", iterations = 2000, burnin = 500)
#' }
#' }
#'
#' @seealso \code{\link{get_kernel}}, \code{varcomp_summary},
#'   \code{kernel_model_mc}, \code{kernel_cv}
#' @importFrom stats rnorm rgamma var sd qchisq quantile cov cov2cor model.matrix
#' @importFrom utils flush.console head
#' @export


## ===========================================================================
## Variance-component machinery
##
## Established by simulation against known truth:
##   * Gibbs "realized" variances inflate sigma2_GE (+70% to +90% at q=8) and
##     deflate sigma2_e correspondingly, because the GxE Hadamard kernel is
##     ~59% identity-like and competes with the residual.
##   * REJECTED after testing: weaker prior, 3x iterations, kernel bending,
##     double-centering, dropping G, heteroscedastic reweighting.
##   * WORKS: component variances from POSTERIOR-MEAN effects, residual by
##     SUBTRACTION from var(y), varG retained in the h2 denominator.
##
## Validation (60 genotypes x 8 environments, seed 11):
##   quantity   truth   gibbs   postmean
##   varGE      0.380   0.721   0.393
##   residual   0.932   0.532   0.970
##   h2 entry   0.352   0.534   0.365
##   rg         0.595   0.480   0.597
## ===========================================================================

## Local copy: kernel_model.R must be sourceable without get_kernel.R.
if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a

#' Map kernel names to G / E / GE roles. GxE is tested first so "GE" is never
#' caught by a "^G" rule.
#' @keywords internal
.km_role_of <- function(nms) {
  role <- rep(NA_character_, length(nms))
  for (i in seq_along(nms)) {
    n <- nms[i]
    if (grepl("GE", n))                                        role[i] <- "GE"
    else if (grepl("^KE_|^E$|env", n, ignore.case = TRUE))     role[i] <- "E"
    else if (grepl("^KG_|^G$|gid|geno", n, ignore.case = TRUE)) role[i] <- "G"
  }
  role
}

#' Realized variance of each random term.
#' @keywords internal
.km_component_variances <- function(fit, method = c("postmean", "gibbs")) {
  method <- match.arg(method)
  out <- c(varG = NA_real_, varE = NA_real_, varGE = NA_real_)

  if (method == "gibbs") {
    vc <- fit$VarComp
    if (is.null(vc)) return(out)
    nm <- if ("K" %in% names(vc)) as.character(vc$K) else rownames(vc)
    v  <- as.numeric(vc$Var)
    role <- .km_role_of(nm)
    for (r in c("G", "E", "GE")) {
      i <- which(role == r)
      if (length(i)) out[[paste0("var", r)]] <- v[i[1]]
    }
    return(out)
  }

  U <- fit$fit$Uchain %||% fit$BGGE$Uchain   # $BGGE accepted for old objects
  if (is.null(U) || !length(U)) return(out)
  pm <- lapply(U, colMeans)
  role <- .km_role_of(names(pm))
  for (r in c("G", "E", "GE")) {
    i <- which(role == r)
    if (length(i)) out[[paste0("var", r)]] <- stats::var(Reduce(`+`, pm[i]))
  }
  out
}

#' Validated variance components, heritabilities and genetic correlation.
#'
#' var_total = var(y) (PHENOTYPIC; var(yHat) carries no residual and yields a
#' NEGATIVE residual by subtraction).
#' var_resid = var_total - (varG + varE + varGE)
#' h2_entry  = varG / (varG + varE/q + varGE/(p*q) + var_resid)
#'   varG is retained in the denominator; omitting it produced h2 > 1 even with
#'   true inputs (observed 1.025 and 0.983).
#'
#' @param fit a fitted \code{kernel_model} object (or a compatible list carrying
#'   variance components).
#' @param method character. \code{"postmean"} (default) uses posterior-mean effects;
#'   \code{"gibbs"} uses the stored Gibbs variance components.
#' @param y numeric. Observed phenotypes used as the total-variance reference.
#'   Taken from \code{fit$meta$y} when \code{NULL}.
#' @param p integer. Number of genotypes. Inferred from \code{fit$meta$gid} when \code{NULL}.
#' @param q integer. Number of environments. Inferred from \code{fit$meta$env} when \code{NULL}.
#' @param digits integer. Number of decimal places in the summary table. Default 4.
#' @export
varcomp_summary <- function(fit, method = c("postmean", "gibbs"),
                            y = NULL, p = NULL, q = NULL, digits = 4) {
  method <- match.arg(method)
  cmp   <- .km_component_variances(fit, method = method)
  varG  <- unname(cmp[["varG"]]); varE <- unname(cmp[["varE"]])
  varGE <- unname(cmp[["varGE"]])

  yobs <- if (!is.null(y)) as.numeric(y)
          else if (!is.null(fit$meta) && !is.null(fit$meta$y)) as.numeric(fit$meta$y)
          else stop("Cannot find observed phenotype; pass it via `y =`.")
  var_total <- stats::var(yobs, na.rm = TRUE)

  if (is.null(p)) p <- if (!is.null(fit$meta$gid)) length(unique(fit$meta$gid)) else NA_integer_
  if (is.null(q)) q <- if (!is.null(fit$meta$env)) length(unique(fit$meta$env)) else NA_integer_

  var_resid <- var_total - sum(c(varG, varE, varGE), na.rm = TRUE)
  negative_residual <- isTRUE(var_resid < 0)
  if (negative_residual) {
    warning("Residual variance negative (", round(var_resid, 4), "): genetic ",
            "components exceed var(y). Try method = 'postmean' or ",
            "truncate_gxe_kernel(). Floored at a small positive value.")
    var_resid <- max(var_total * 1e-6, 1e-8)
  }

  g_ <- function(x) if (is.na(x)) 0 else x
  denom_entry <- g_(varG) +
    (if (!is.na(varE)  && !is.na(q) && q > 0) varE / q else 0) +
    (if (!is.na(varGE) && !is.na(p) && !is.na(q) && p * q > 0) varGE / (p * q) else 0) +
    var_resid
  h2_entry <- if (denom_entry > 0) g_(varG) / denom_entry else NA_real_
  denom_plot <- g_(varG) + g_(varGE) + var_resid
  h2_plot <- if (denom_plot > 0) g_(varG) / denom_plot else NA_real_
  rg <- if (!is.na(varG) && !is.na(varGE) && (varG + varGE) > 0)
    varG / (varG + varGE) else NA_real_

  tbl <- data.frame(
    Component = c("Genotype (G)", "Environment (E)", "GxE", "Residual",
                  "Total (var(y))"),
    Variance = round(c(varG, varE, varGE, var_resid, var_total), digits),
    Prop = round(c(varG, varE, varGE, var_resid, var_total) / var_total, digits),
    stringsAsFactors = FALSE)

  out <- list(table = tbl, varG = varG, varE = varE, varGE = varGE,
              var_resid = var_resid, var_total = var_total,
              h2_entry_mean = h2_entry, h2_plot = h2_plot, rg = rg,
              p = p, q = q, method = method,
              negative_residual = negative_residual,
              few_environments = isTRUE(!is.na(q) && q < 20))
  class(out) <- "km_varcomp"
  out
}

#' @export
print.km_varcomp <- function(x, ...) {
  cat("------------------------------------------------------------\n")
  cat("Variance components  (method = '", x$method, "')\n", sep = "")
  cat("  genotypes p = ", x$p, "   environments q = ", x$q, "\n", sep = "")
  cat("------------------------------------------------------------\n")
  print(x$table, row.names = FALSE)
  cat("------------------------------------------------------------\n")
  cat(sprintf("  h2 (entry-mean) : %.4f\n", x$h2_entry_mean))
  cat(sprintf("  h2 (plot-level) : %.4f\n", x$h2_plot))
  cat(sprintf("  rg  vG/(vG+vGE) : %.4f\n", x$rg))
  cat("------------------------------------------------------------\n")
  if (isTRUE(x$negative_residual))
    cat("  ! residual was negative before flooring - components inflated\n")
  if (isTRUE(x$few_environments))
    cat("  ! q < 20 environments: sigma2_GE is only approximate.\n",
        "    Validated GxE bias: +0.90 at q=8 vs +0.09 at q=20.\n", sep = "")
  invisible(x)
}

#' Truncate the GxE kernel to its structural rank.
#'
#' The GxE kernel has a flat eigenvalue tail that behaves like an identity matrix
#' and competes with the residual. Dropping it reduces sigma2_GE inflation (mean
#' abs bias 0.333 -> 0.130) at essentially no cost in cross-validated prediction
#' (CV2 0.678 full vs 0.669 truncated).
#'
#' CAUTION: truncation distorts the G:GxE ratio -- rg moved to 0.72 vs a true
#' 0.595 in testing. Use only when sigma2_GE specifically matters.
#'
#' Any cached decomposition is dropped, since the kernel changes.
#'
#' @param random a named list of kernels (the \code{random} argument of \code{kernel_model}).
#' @param keep_prop numeric in (0, 1]. Proportion of the GxE eigenvalue mass to retain.
#'   Default 0.75.
#' @param pattern character. Regular expression selecting the GxE kernels to truncate.
#'   Default \code{"GE"}.
#' @export
truncate_gxe_kernel <- function(random, keep_prop = 0.75, pattern = "GE") {
  if (keep_prop >= 1) return(random)
  if (keep_prop <= 0) stop("keep_prop must be in (0, 1].")
  sel <- grep(pattern, names(random))
  if (!length(sel)) {
    warning("No kernel matched '", pattern, "'; returning input unchanged.")
    return(random)
  }
  for (j in sel) {
    K  <- random[[j]]$Kernel
    ei <- eigen(K, symmetric = TRUE)
    v  <- pmax(ei$values, 0)
    if (sum(v) <= 0) next
    r <- which(cumsum(v) / sum(v) >= keep_prop)[1]
    r <- max(2L, min(r, length(v) - 1L))
    Kt <- ei$vectors[, 1:r, drop = FALSE] %*% diag(v[1:r], r) %*%
          t(ei$vectors[, 1:r, drop = FALSE])
    Kt <- Kt + diag(1e-8, nrow(Kt))
    dimnames(Kt) <- dimnames(K)
    random[[j]]$Kernel <- Kt
    random[[j]]$decomp <- NULL          # invalidate stale cache
    random[[j]]$base <- random[[j]]$idx <- NULL
    message(sprintf("truncate_gxe_kernel: '%s' rank %d -> %d (%.0f%% of trace)",
                    names(random)[j], length(v), r, 100 * keep_prop))
  }
  attr(random, "kd_meta") <- NULL
  random
}


## ===========================================================================
## kernel_model
## ===========================================================================

#' @rdname kernel_model
#' @export
kernel_model <- function(y, data = NULL, random = NULL, fixed = NULL, env, gid,
                         verbose = FALSE, iterations = 1E3, burnin = 2E2,
                         thining = 10, tol = 1e-10, R2 = 0.5, digits = 4,
                         seed = NULL, probs = c(0.025, 0.975),
                         variance = c("realized", "scale"),
                         scale_kernels = TRUE, keep_effects = FALSE,
                         compute_CI = TRUE,
                         R2_split = c("per_kernel", "shared"),
                         prior_scale = c("eigen", "diag"),
                         varcomp_method = c("postmean", "gibbs")) {

  variance <- match.arg(variance)
  varcomp_method <- match.arg(varcomp_method)
  r2_split <- match.arg(R2_split)
  pscale   <- match.arg(prior_scale)

  ## ---- verbose step tracker -------------------------------------------------
  ## One helper drives all progress output. When verbose is FALSE/0 every call
  ## is a no-op, so the function is completely silent -- no cat(), no message().
  ## Each step records its own elapsed time; the total is reported at the end.
  .v      <- as.numeric(verbose) != 0
  .t_all  <- proc.time()[3]
  .t_step <- .t_all
  .steps  <- list()
  .step <- function(txt) {                     # announce the start of a step
    if (!.v) return(invisible(NULL))
    .t_step <<- proc.time()[3]
    cat(sprintf("  %-52s", paste0(txt, " ...")))
    utils::flush.console()
    invisible(NULL)
  }
  .done <- function(extra = NULL) {            # close it and print the timing
    if (!.v) return(invisible(NULL))
    el <- proc.time()[3] - .t_step
    .steps[[length(.steps) + 1L]] <<- el
    cat(sprintf(" %7.2fs%s\n", el, if (is.null(extra)) "" else paste0("  ", extra)))
    invisible(NULL)
  }
  .note <- function(txt) {                     # sub-detail, no timing
    if (.v) cat(sprintf("      %s\n", txt))
    invisible(NULL)
  }

  if (.v) {
    cat("\n=========================================================\n")
    cat(" kernel_model: Bayesian genomic prediction for multi-\n")
    cat("   environment trials. Fits y = 1mu + Xb + sum_j u_j + e\n")
    cat("   with u_j ~ N(0, K_j sigma2_j) by Gibbs sampling in the\n")
    cat("   spectral basis of each kernel, returning predictions,\n")
    cat("   variance components and heritabilities.\n")
    cat(" Original Version: Costa-Neto et al (2021)\n")
    cat(" EnvRtype v.1.2.3, Sep 2026\n")
    cat("=========================================================\n")
    cat(" Gibbs sampler after Granato et al (2018) G3\n")
    cat("---------------------------------------------------------\n")
  }

  .step("[1/6] Validating inputs")

  ## ---- validation ----------------------------------------------------------
  if (is.null(data)) stop("`data` is required.")
  if (is.null(random)) stop("Missing the list of kernels for random effects (`random`).")
  if (!is.data.frame(data)) stop("`data` must be a data.frame.")

  miss_cols <- setdiff(c(env, gid, y), names(data))
  if (length(miss_cols))
    stop("Columns not found in `data`: ", paste(miss_cols, collapse = ", "))

  if (!is.list(random) || length(random) == 0)
    stop("`random` must be a non-empty list of kernels.")
  ## ------------------------------------------------------------------------
  ## Kernels may arrive in either of two forms:
  ##   DENSE    $Kernel (n x n) + $Type          -- as produced by get_kernel()
  ##   FACTORED $decomp + $Type + $diag_mean     -- from get_kernel_test()
  ## The Gibbs sampler below reads ONLY $decomp ($U, $s, $deltav, $nr), so a
  ## dense matrix is never required by the algorithm itself -- it was only
  ## ever used for validation and for the mean(diag(K)) prior scale. Allowing
  ## the factored form removes the n x n allocation entirely: at 2000
  ## genotypes x 30 environments that is 29 GB per kernel avoided.
  ## ------------------------------------------------------------------------
  bad_struct <- vapply(random, function(k)
    !is.list(k) || is.null(k$Type) ||
      (is.null(k$Kernel) && is.null(k$decomp)), logical(1))
  if (any(bad_struct))
    stop("Each element of `random` must be a list with `$Type` and either ",
         "`$Kernel` (dense) or `$decomp` (factored). ",
         "Offending element(s): ", paste(which(bad_struct), collapse = ", "))

  is_factored <- vapply(random, function(k)
    is.null(k$Kernel) && !is.null(k$decomp), logical(1))

  n_obs <- nrow(data)
  bad_dim <- vapply(seq_along(random), function(j) {
    k <- random[[j]]
    if (is_factored[j]) {
      ## factored: every eigenvector block must have n rows
      nr <- vapply(k$decomp, function(e)
        if (is.null(e$U)) NA_integer_ else nrow(e$U), integer(1))
      !all(is.na(nr) | nr == n_obs)
    } else {
      K <- k$Kernel
      !is.matrix(K) || nrow(K) != n_obs || ncol(K) != n_obs
    }
  }, logical(1))
  if (any(bad_dim))
    stop("Every kernel must correspond to n = nrow(data) = ", n_obs,
         " observations. Offending kernel(s): ",
         paste(which(bad_dim), collapse = ", "))

  if (!is.numeric(iterations) || iterations < 1)
    stop("`iterations` must be a positive integer.")
  if (!is.numeric(thining) || thining < 1) stop("`thining` must be >= 1.")
  if (burnin >= iterations)
    stop("`burnin` (", burnin, ") must be smaller than `iterations` (", iterations, ").")
  if (length(probs) != 2 || any(probs < 0) || any(probs > 1) || probs[1] >= probs[2])
    stop("`probs` must be two increasing values in [0, 1].")
  if (!is.null(fixed)) {
    fixed <- as.matrix(fixed)
    if (nrow(fixed) != n_obs)
      stop("`fixed` must have nrow(data) = ", n_obs, " rows.")
  }

  .done(sprintf("%d obs, %d kernel%s", nrow(data), length(random),
                if (length(random) == 1L) "" else "s"))

  if (!is.null(seed)) set.seed(as.integer(seed))

  .step("[2/6] Assembling response and design")
  Y <- data.frame(env = data[[env]], gid = data[[gid]], y = data[[y]])
  ne <- as.vector(table(factor(Y$env, levels = unique(as.character(Y$env)))))
  .done(sprintf("%d env, %d gid", length(unique(Y$env)), length(unique(Y$gid))))
  ## ---- decomposition: reuse cache or build it ------------------------------
  ## Factored kernels arrive pre-decomposed, so the eigen step is skipped for
  ## them. Mixed lists (some dense, some factored) are supported: only the
  ## dense members are sent through decompose_kernels().
  .step("[3/6] Spectral decomposition")
  if (all(is_factored)) {
    .done("all kernels pre-factored; eigen skipped")
  } else {
    cache_ok <- .kd_valid(random, tol = tol, scale_kernels = scale_kernels, n = n_obs)
    if (!cache_ok) {
      stale <- !is.null(attr(random, "kd_meta"))
      dense_idx <- which(!is_factored)
      sub <- random[dense_idx]
      ## decompose_kernels() prints its own multi-line route report; keep it
      ## silent here so the step tracker owns the output format.
      sub <- decompose_kernels(sub, ne = ne, tol = tol,
                               scale_kernels = scale_kernels,
                               verbose = FALSE)
      for (i in seq_along(dense_idx)) random[[dense_idx[i]]] <- sub[[i]]
      tot_r <- sum(vapply(sub, function(k)
        sum(vapply(k$decomp, function(d) d$nr, numeric(1))), numeric(1)))
      .done(sprintf("%d dense kernel%s, total rank %d%s",
                    length(dense_idx), if (length(dense_idx) == 1L) "" else "s",
                    tot_r, if (stale) "  [cache stale, recomputed]" else ""))
    } else {
      .done("cached decomposition reused; eigen skipped")
    }
  }

  ## kernel scaling is folded into $decomp; keep scale factors for priors.
  ## A factored kernel cannot expose diag(K) without forming it, so
  ## get_kernel_test() records mean(diag(K)) as $diag_mean at build time --
  ## computable as sum over eigen-directions at O(n*r), never O(n^2).
  .step("[4/6] Computing prior scales")
  kscale <- vapply(random, function(k) {
    f <- if (!is.null(k$Kernel)) mean(diag(k$Kernel))
         else if (!is.null(k$diag_mean)) k$diag_mean else 1
    if (isTRUE(scale_kernels) && is.finite(f) && f > 0) f else 1
  }, numeric(1))
  .done(sprintf("R2 = %.2f split '%s', prior_scale '%s'", R2, r2_split, pscale))

  ## =========================================================================
  ## Gibbs sampler
  ## =========================================================================
  BGGE <- function(y, K, XF = NULL, ite = 1000, burn = 200, thin = 3,
                   verbose = FALSE, R2 = 0.5, probs = c(0.025, 0.975),
                   variance = "realized", keep_effects = FALSE,
                   compute_CI = TRUE, kscale = NULL) {

    dcondsigb <- function(b, deltav, n, nu, Sc) {
      z <- sum(b * b * deltav)
      1 / rgamma(1, (n + nu) / 2, (z + nu * Sc) / 2)
    }
    dcondsigsq <- function(Aux, n, nu, Sce)
      1 / rgamma(1, (n + nu) / 2, crossprod(Aux) / 2 + Sce / 2)
    rmvnor_chol <- function(n, media, cholR) media + crossprod(cholR, rnorm(n))

    y <- as.numeric(y)
    yNA <- is.na(y); whichNa <- which(yNA); nNa <- length(whichNa)
    mu <- mean(y, na.rm = TRUE)
    n <- length(y)
    if (nNa > 0) y[whichNa] <- mu

    if (is.null(names(K))) names(K) <- paste0("K", seq_along(K))
    nk <- length(K)
    typeM <- vapply(K, function(x) x$Type, character(1))

    ## decompositions (already computed and cached upstream)
    Ei <- lapply(K, function(k) k$decomp)

    ## flatten BD bookkeeping once
    bd_meta <- vector("list", nk)
    for (j in seq_len(nk)) {
      if (length(Ei[[j]]) > 1) {
        blocks <- Ei[[j]]
        neiv <- vapply(blocks, function(e) e$nr, numeric(1))
        bd_meta[[j]] <- list(
          nsk  = length(blocks),
          pos  = t(vapply(blocks, function(e) e$pos, numeric(2))),
          posf = cumsum(neiv),
          posi = cumsum(c(1, neiv[-length(neiv)])),
          s    = unlist(lapply(blocks, function(e) e$s), use.names = FALSE),
          deltav = unlist(lapply(blocks, function(e) e$deltav), use.names = FALSE))
      } else if (!is.na(Ei[[j]][[1]]$pos[1])) {
        e1 <- Ei[[j]][[1]]
        bd_meta[[j]] <- list(single_rng = e1$pos[1]:e1$pos[2])
      }
    }

    if (!is.null(XF)) {
      tXX  <- solve(crossprod(XF))
      Bet  <- tXX %*% crossprod(XF, y)
      nBet <- length(Bet)
      cholTXX <- chol(tXX)
    }

    nCum <- sum(seq_len(ite) %% thin == 0)
    saved_iter <- which(seq_len(ite) %% thin == 0)
    draw <- saved_iter > burn
    if (sum(draw) < 2)
      stop("`burnin` too large relative to `iterations`/`thining`: ",
           "no post-burn-in samples remain.")

    ## scalar chains as a matrix (no nested-list copy-on-modify)
    chain_mu   <- numeric(nCum)
    chain_varE <- numeric(nCum)
    chain_varU <- matrix(0, nrow = nCum, ncol = nk,
                         dimnames = list(NULL, names(K)))

    ## Full per-draw effect chains only when actually needed.
    store_chains <- isTRUE(keep_effects) || isTRUE(compute_CI)
    if (store_chains) {
      cpred <- vector("list", nk); names(cpred) <- names(K)
      for (j in seq_len(nk)) cpred[[j]] <- matrix(NA_real_, nrow = nCum, ncol = n)
    } else {
      cpred <- NULL
      u_sum <- lapply(seq_len(nk), function(j) numeric(n))
      u_sq  <- lapply(seq_len(nk), function(j) numeric(n))
      vu_dr <- matrix(NA_real_, nrow = nCum, ncol = nk)
    }

    nu  <- 3
    vy  <- var(y, na.rm = TRUE)
    Sce <- (nu + 2) * (1 - R2) * vy

    ## ---- prior scale per kernel -------------------------------------------
    ## Two defects in the original parameterisation:
    ##
    ## (1) R2_SPLIT. `Sc <- rep((nu+2)*R2*vy, nk)` gives EVERY kernel the full
    ##     R2 share, so with nk kernels the priors collectively assert nk * R2
    ##     of the variance. With RNMDs (nk = 4) that is 2x var(y). Splitting
    ##     R2 across kernels keeps the total prior mass at R2 regardless of how
    ##     many kernels the model has.
    ##
    ## (2) PRIOR_SCALE. Kernels are normalised to mean(diag(K)) = 1, which
    ##     equalises scale but NOT effective dimension. A rank-5 kernel and a
    ##     rank-749 kernel with the same mean diagonal carry very different
    ##     total variance (sum of eigenvalues), so an identical Sc means
    ##     something different for each. Dividing by mean(s) (= trace / rank)
    ##     makes the prior comparable across kernels of differing rank.
    R2_j <- if (identical(r2_split, "per_kernel")) R2 / nk else R2
    Sc <- numeric(nk)
    for (j in seq_len(nk)) {
      denom <- 1
      if (identical(pscale, "eigen")) {
        sj <- unlist(lapply(Ei[[j]], function(e) e$s), use.names = FALSE)
        mj <- mean(sj)
        if (is.finite(mj) && mj > 0) denom <- mj
      }
      Sc[j] <- (nu + 2) * R2_j * vy / denom
    }

    tau   <- 0.01
    u     <- lapply(seq_len(nk), function(j) rnorm(n, 0, 1 / (2 * n)))
    sigsq <- vy
    sigb  <- rep(0.2, nk)

    temp <- y - mu
    if (!is.null(XF)) {
      B.mcmc <- matrix(0, nrow = nCum, ncol = nBet)
      temp <- temp - XF %*% Bet
    }
    temp <- temp - Reduce('+', u)
    nSel <- 0L
    vrb <- as.numeric(verbose)

    ## ---- hoist per-kernel structures out of the sampling loop -------------
    ## Repeated Ei[[j]][[1]]$U lookups inside 5000 x nk iterations are pure
    ## interpreter overhead. Resolve every invariant once, into flat vectors.
    ##
    ## NOTE: use single-bracket assignment for rng1. `rng1[[j]] <- NULL` would
    ## DELETE the element and shrink the list, causing a subscript-out-of-bounds
    ## on the next dense kernel. `rng1[j] <- list(NULL)` stores a real NULL.
    k_single <- vapply(Ei, function(e) length(e) == 1L, logical(1))
    U1   <- vector("list", nk)   # dense / single-block eigenvectors
    s1   <- vector("list", nk)
    dv1  <- vector("list", nk)
    nr1  <- integer(nk)
    rng1 <- vector("list", nk)   # NULL for dense, index range for single-block
    for (j in seq_len(nk)) {
      if (k_single[j]) {
        E1 <- Ei[[j]][[1]]
        U1[[j]] <- E1$U; s1[[j]] <- E1$s; dv1[[j]] <- E1$deltav
        nr1[j] <- E1$nr
        rng1[j] <- list(bd_meta[[j]]$single_rng)   # may legitimately be NULL
      } else {
        m <- bd_meta[[j]]
        s1[[j]] <- m$s; dv1[[j]] <- m$deltav; nr1[j] <- length(m$s)
        rng1[j] <- list(NULL)
      }
    }

    for (i in seq_len(ite)) {
      ## proc.time() every iteration costs ~1% of total runtime; only pay it
      ## when the user actually asked for iteration timing.
      if (vrb != 0) time.init <- proc.time()[3]

      ## mu
      temp <- temp + mu
      mu   <- rnorm(1, mean(temp), sqrt(sigsq / n))
      temp <- temp - mu

      ## fixed effects
      if (!is.null(XF)) {
        temp  <- temp + XF %*% Bet
        media <- tXX %*% crossprod(XF, temp)
        Bet   <- rmvnor_chol(nBet, media, sqrt(sigsq) * cholTXX)
        temp  <- temp - XF %*% Bet
      }

      ## random effects, kernel by kernel
      for (j in seq_len(nk)) {
        lambda <- sigb[j]
        uj <- u[[j]]
        temp <- temp + uj

        if (k_single[j]) {
          Uj  <- U1[[j]]
          rng <- rng1[[j]]
          d   <- if (is.null(rng)) crossprod(Uj, temp) else crossprod(Uj, temp[rng])

          sl   <- s1[[j]] * lambda
          vari <- sl / (1 + sl * tau)
          nr   <- nr1[j]
          b    <- rnorm(nr, tau * vari * d, sqrt(vari))

          if (is.null(rng)) {
            uj <- Uj %*% b
          } else {
            ## fresh allocation: a shared scratch buffer would alias across
            ## kernels, since u[[j]] <- uj does not deep-copy.
            uj <- numeric(n)
            uj[rng] <- Uj %*% b
          }
          deltav <- dv1[[j]]

        } else {                                    # multi-block BD
          m      <- bd_meta[[j]]
          blocks <- Ei[[j]]
          d <- numeric(length(m$s))
          for (k in seq_len(m$nsk)) {
            d[m$posi[k]:m$posf[k]] <-
              crossprod(blocks[[k]]$U, temp[m$pos[k, 1]:m$pos[k, 2]])
          }
          sl   <- m$s * lambda
          vari <- sl / (1 + sl * tau)
          nr   <- nr1[j]
          b    <- rnorm(nr, tau * vari * d, sqrt(vari))

          uj <- numeric(n)
          for (k in seq_len(m$nsk)) {
            uj[m$pos[k, 1]:m$pos[k, 2]] <-
              blocks[[k]]$U %*% b[m$posi[k]:m$posf[k]]
          }
          deltav <- dv1[[j]]
        }

        u[[j]] <- uj
        temp   <- temp - uj
        sigb[j] <- dcondsigb(b, deltav, nr, nu, Sc[j])
      }

      ## residual variance
      sigsq <- dcondsigsq(temp, n, nu, Sce)
      tau   <- 1 / sigsq

      ## impute missing phenotypes (only the NA entries are needed)
      if (nNa > 0) {
        uhat_na <- Reduce('+', lapply(u, function(x) x[whichNa]))
        aux <- if (!is.null(XF)) as.vector(XF[whichNa, , drop = FALSE] %*% Bet) else 0
        y[whichNa]    <- aux + mu + uhat_na + rnorm(nNa, sd = sqrt(sigsq))
        temp[whichNa] <- y[whichNa] - uhat_na - aux - mu
      }

      ## store thinned draw
      if (i %% thin == 0) {
        nSel <- nSel + 1L
        chain_varE[nSel] <- sigsq
        chain_mu[nSel]   <- mu
        chain_varU[nSel, ] <- sigb
        if (!is.null(XF)) B.mcmc[nSel, ] <- Bet
        if (store_chains) {
          for (j in seq_len(nk)) cpred[[j]][nSel, ] <- u[[j]]
        } else {
          for (j in seq_len(nk)) {
            uj <- u[[j]]
            u_sum[[j]] <- u_sum[[j]] + uj
            u_sq[[j]]  <- u_sq[[j]]  + uj * uj
            um <- mean(uj)
            vu_dr[nSel, j] <- sum((uj - um)^2) / (n - 1)
          }
        }
      }

      if (vrb != 0 && i %% vrb == 0)
        cat("Iter: ", i, " time: ", round(proc.time()[3] - time.init, 3), "\n")
    }

    ## ---- output ------------------------------------------------------------
    ndraw <- sum(draw)
    mu.est <- mean(chain_mu[draw])
    yHat <- mu.est
    if (!is.null(XF)) {
      B <- colMeans(B.mcmc[draw, , drop = FALSE])
      yHat <- yHat + as.vector(XF %*% B)
    }
    fixed_part <- if (!is.null(XF)) as.vector(XF %*% B) else 0

    out <- list()
    out$K <- vector("list", nk); names(out$K) <- names(K)

    if (store_chains) {
      u.est <- vapply(cpred, function(x) colMeans(x[draw, , drop = FALSE]), numeric(n))
      yHat <- yHat + rowSums(u.est)

      ## posterior prediction intervals -- recycled arithmetic, no sweep()
      U_sum <- Reduce('+', lapply(cpred, function(x) x[draw, , drop = FALSE]))
      if (isTRUE(compute_CI)) {
        yd <- U_sum + chain_mu[draw]                        # recycles down columns
        if (!is.null(XF)) yd <- yd + rep(fixed_part, each = ndraw)
        out$yHat.CI <- t(apply(yd, 2, quantile, probs = probs, na.rm = TRUE))
        colnames(out$yHat.CI) <- c("lower", "upper")
        rm(yd)
      }
      rm(U_sum)

      for (jj in seq_len(nk)) {
        xd <- cpred[[jj]][draw, , drop = FALSE]
        cm <- colMeans(xd)
        out$K[[jj]]$u    <- cm
        m2 <- colMeans(xd * xd)
        out$K[[jj]]$u.sd <- sqrt(pmax((m2 - cm^2) * ndraw / (ndraw - 1), 0))
        ## vectorised realized variance: (rowSums(x^2) - n*rowMean^2)/(n-1)
        rm_  <- rowMeans(xd)
        vdr  <- (rowSums(xd * xd) - n * rm_^2) / (n - 1)
        if (variance == "realized") {
          out$K[[jj]]$varu <- mean(vdr); out$K[[jj]]$varu.sd <- sd(vdr)
        } else {
          out$K[[jj]]$varu <- mean(chain_varU[draw, jj])
          out$K[[jj]]$varu.sd <- sd(chain_varU[draw, jj])
        }
        out$K[[jj]]$sigb <- mean(chain_varU[draw, jj])
      }
    } else {
      ## running-moment path (no per-draw matrices retained)
      for (jj in seq_len(nk)) {
        cm <- u_sum[[jj]] / nCum
        out$K[[jj]]$u    <- cm
        m2 <- u_sq[[jj]] / nCum
        out$K[[jj]]$u.sd <- sqrt(pmax((m2 - cm^2) * nCum / (nCum - 1), 0))
        vdr <- vu_dr[draw, jj]
        if (variance == "realized") {
          out$K[[jj]]$varu <- mean(vdr); out$K[[jj]]$varu.sd <- sd(vdr)
        } else {
          out$K[[jj]]$varu <- mean(chain_varU[draw, jj])
          out$K[[jj]]$varu.sd <- sd(chain_varU[draw, jj])
        }
        out$K[[jj]]$sigb <- mean(chain_varU[draw, jj])
      }
      yHat <- yHat + Reduce('+', lapply(out$K, function(z) z$u))
    }

    out$yHat  <- as.vector(yHat)
    out$varE  <- mean(chain_varE[draw])
    out$varE.sd <- sd(chain_varE[draw])
    out$chain <- list(mu = chain_mu, varE = chain_varE,
                      K = setNames(lapply(seq_len(nk),
                                          function(j) list(varU = chain_varU[, j])),
                                   names(K)))
    out$ite <- ite; out$burn <- burn; out$thin <- thin
    out$kernel_names <- names(K)
    if (isTRUE(keep_effects) && store_chains) out$Uchain <- cpred
    out$y <- y
    class(out) <- "BGGE"
    out
  }

  ## =========================================================================
  Vcomp.BGGE <- function(model, env, gid, digits = digits, alfa = .10) {
    t <- length(unique(gid)); e <- length(unique(env)); n <- t * e
    K <- model$K; size <- length(K)
    comps <- data.frame(matrix(NA, ncol = 3, nrow = size))
    VarE  <- data.frame(matrix(NA, ncol = 3, nrow = 1))
    names(comps) <- names(VarE) <- c("K", "Var", "SD.var")
    for (k in seq_len(size)) {
      comps[k, 1] <- names(K)[k]
      comps[k, 2] <- round(K[[k]]$varu, digits)
      comps[k, 3] <- round(K[[k]]$varu.sd, digits)
    }
    VarE[1, ] <- list("Residual", round(model$varE, digits),
                      round(model$varE.sd, digits))
    comps <- rbind(comps, VarE)

    comps$Type <- NA_character_
    ## ORDER MATTERS: 'KG_GE' matches BOTH '^KGE_'(no) and '^KG_'(yes), so the
    ## generic genotype rule must not overwrite a GxE assignment. Test every
    ## GxE spelling FIRST, then fill only the still-unassigned rows.
    is_gxe <- grepl("^KGE_", comps$K) | grepl("GE", comps$K)
    comps$Type[is_gxe] <- "GxE"
    comps$Type[is.na(comps$Type) & grepl("^KE_", comps$K)] <- "Environment (E)"
    comps$Type[is.na(comps$Type) & grepl("^KG_", comps$K)] <- "Genotype (G)"
    comps$Type[comps$K == "Residual"]   <- "Residual"

    comps$CI_upper <- NA; comps$CI_lower <- NA
    ENV <- which(comps$Type == 'Environment (E)')
    GID <- which(comps$Type == 'Genotype (G)')
    GE  <- which(comps$Type == 'GxE')
    R   <- which(comps$Type == 'Residual')

    comps$CI_upper[ENV] <- (n - e)     * comps$Var[ENV] / qchisq(alfa / 2, n - e)
    comps$CI_upper[GID] <- (n - t)     * comps$Var[GID] / qchisq(alfa / 2, n - t)
    comps$CI_upper[GE]  <- (n - t - e) * comps$Var[GE]  / qchisq(alfa / 2, n - t - e)
    comps$CI_upper[R]   <- (n - t - e) * comps$Var[R]   / qchisq(alfa / 2, n - t - e)
    comps$CI_lower[ENV] <- (n - e)     * comps$Var[ENV] / qchisq(1 - alfa / 2, n - e)
    comps$CI_lower[GID] <- (n - t)     * comps$Var[GID] / qchisq(1 - alfa / 2, n - t)
    comps$CI_lower[GE]  <- (n - t - e) * comps$Var[GE]  / qchisq(1 - alfa / 2, n - t - e)
    comps$CI_lower[R]   <- (n - t - e) * comps$Var[R]   / qchisq(1 - alfa / 2, n - t - e)

    comps$CI_upper <- round(comps$CI_upper, digits)
    comps$CI_lower <- round(comps$CI_lower, digits)
    comps <- comps[, c(4, 1:2, 6, 5, 3)]
    for (p in c('KG_', 'KE_', 'KGE_'))
      comps$K <- gsub(x = comps$K, pattern = p, replacement = '')
    comps
  }

  ## ---- fit -----------------------------------------------------------------
  .step(sprintf("[5/6] Gibbs sampling (%d iter)", iterations))
  start <- Sys.time()
  ## verbose = FALSE for the inner sampler: its per-iteration trace would break
  ## the step layout. Progress is summarised by the step tracker instead.
  fit <- BGGE(y = Y$y, K = random, XF = fixed,
              ite = iterations, burn = burnin, thin = thining,
              verbose = FALSE, R2 = R2, probs = probs,
              variance = variance, keep_effects = keep_effects,
              compute_CI = compute_CI, kscale = kscale)
  end <- Sys.time()
  .done(sprintf("%d draws kept", length(seq_len(iterations)[
         seq_len(iterations) %% thining == 0 & seq_len(iterations) > burnin])))

  ret <- list(yHat = fit$yHat,
              yHat.CI = fit$yHat.CI,
              varE = fit$varE,
              random = fit$K,
              ## Renamed from $BGGE. The sampler is still Granato's, credited in
              ## the roxygen @author block and the verbose banner; only the
              ## output slot name changes, so downstream code reads $fit.
              fit = fit,
              meta = Y,
              VarComp = Vcomp.BGGE(model = fit, env = Y$env, gid = Y$gid,
                                   digits = digits),
              runtime = as.numeric(difftime(end, start, units = "secs")))

  ## ---- breeder-facing reports (require keep_effects) -----------------------
  if (isTRUE(keep_effects) && !is.null(fit$Uchain)) {
    ite. <- fit$ite; burn. <- fit$burn; thin. <- fit$thin
    draw. <- which(seq_len(ite.) %% thin. == 0) > burn.
    U. <- fit$Uchain; nmU <- names(U.)

    eff_chain <- function(pattern = NULL) {
      sel <- if (is.null(pattern)) seq_along(U.) else grep(pattern, nmU)
      if (length(sel) == 0) return(NULL)
      Reduce(`+`, lapply(sel, function(j) U.[[j]][draw., , drop = FALSE]))
    }

    ei. <- as.character(Y$env); gi. <- as.character(Y$gid)
    envs. <- sort(unique(ei.)); gids. <- sort(unique(gi.))

    Gtot <- eff_chain(NULL)
    Ggen <- eff_chain("^KG_|^KGE_")
    Genv <- eff_chain("^KE_")

    ## (a) predictions with PEV / reliability
    preds <- data.frame(env = Y$env, gid = Y$gid,
                        yHat = as.numeric(fit$yHat), row.names = NULL)
    if (!is.null(fit$yHat.CI)) {
      preds$lower <- fit$yHat.CI[, "lower"]; preds$upper <- fit$yHat.CI[, "upper"]
    }
    if (!is.null(Gtot)) {
      nd <- nrow(Gtot); cm <- colMeans(Gtot)
      pev <- (colMeans(Gtot * Gtot) - cm^2) * nd / (nd - 1)
      vg  <- var(cm)
      preds$PEV <- pev
      preds$reliability <- pmax(0, pmin(1, 1 - pev / vg))
    }
    ret$predictions <- preds

    ## (b) genetic covariance / correlation among environments
    e. <- length(envs.); g. <- length(gids.)
    ie <- match(ei., envs.); ig <- match(gi., gids.)
    idx <- cbind(ig, ie)                       # hoisted: constant across draws
    balanced <- !anyDuplicated(idx) && nrow(idx) == g. * e.
    vc_sum <- matrix(0, e., e.); vc_sq <- matrix(0, e., e.); cr_sum <- matrix(0, e., e.)
    valid <- 0L
    for (d in seq_len(nrow(Ggen))) {
      M <- matrix(NA_real_, g., e.)
      M[idx] <- Ggen[d, ]
      V <- if (balanced) stats::cov(M) else stats::cov(M, use = "pairwise.complete.obs")
      if (anyNA(V)) next
      sdv <- sqrt(diag(V))
      if (any(!is.finite(sdv)) || any(sdv <= 0)) next
      R <- V / tcrossprod(sdv)                 # cov2cor without the checks
      vc_sum <- vc_sum + V; vc_sq <- vc_sq + V * V; cr_sum <- cr_sum + R
      valid <- valid + 1L
    }
    if (valid > 0L) {
      vcov. <- vc_sum / valid
      vsd.  <- sqrt(pmax(vc_sq / valid - vcov.^2, 0))
      cor.  <- cr_sum / valid
      dimnames(vcov.) <- dimnames(cor.) <- dimnames(vsd.) <- list(envs., envs.)
      ret$genetic_cor_env <- list(vcov = vcov., cor = cor., vcov_sd = vsd.,
                                  n_draws_used = valid)
    } else ret$genetic_cor_env <- NULL

    ## (c) genotype effects -- one split, not five tapply passes
    nd <- nrow(Ggen); gval <- colMeans(Ggen)
    gsd <- sqrt(pmax((colMeans(Ggen * Ggen) - gval^2) * nd / (nd - 1), 0))
    fg <- factor(gi.)
    sp_val <- split(gval, fg); sp_sd <- split(gsd, fg)
    gmean <- vapply(sp_val, mean, 0); gmin <- vapply(sp_val, min, 0)
    gmax  <- vapply(sp_val, max, 0)
    gsdev <- vapply(sp_val, function(z) if (length(z) > 1) sd(z) else NA_real_, 0)
    gunc  <- vapply(sp_sd, mean, 0)
    geff <- data.frame(gid = names(gmean),
                       mean_effect = as.numeric(gmean),
                       min_effect = as.numeric(gmin),
                       max_effect = as.numeric(gmax),
                       range = as.numeric(gmax - gmin),
                       stability_sd = as.numeric(gsdev),
                       post_sd = as.numeric(gunc), row.names = NULL)
    ret$genotype_effects <- geff[order(-geff$mean_effect), ]

    ## (d) environment summary + heritability
    env_effect <- if (!is.null(Genv)) tapply(colMeans(Genv), ei., mean)
                  else tapply(as.numeric(fit$yHat), ei., mean)
    varG_env <- tapply(gval, ei., var)
    vE. <- fit$varE
    h2_env <- varG_env / (varG_env + vE.)
    n_gid_env <- tapply(gi., ei., function(x) length(unique(x)))
    connect <- if (!is.null(ret$genetic_cor_env)) {
      Rc <- ret$genetic_cor_env$cor; diag(Rc) <- NA
      rowMeans(Rc, na.rm = TRUE)[envs.]
    } else rep(NA_real_, length(envs.))
    by_env <- data.frame(env = envs.,
                         mean_env_effect = as.numeric(env_effect[envs.]),
                         genetic_var = as.numeric(varG_env[envs.]),
                         genomic_h2 = as.numeric(h2_env[envs.]),
                         n_genotypes = as.integer(n_gid_env[envs.]),
                         mean_rg_with_others = as.numeric(connect),
                         row.names = NULL)

    ## variance proportions -- vectorised, no nested apply
    nobs <- ncol(U.[[1]])
    vu <- vapply(nmU, function(nm) {
      x <- U.[[nm]][draw., , drop = FALSE]
      rm_ <- rowMeans(x)
      (rowSums(x * x) - nobs * rm_^2) / (nobs - 1)
    }, numeric(sum(draw.)))
    if (is.null(dim(vu)))
      vu <- matrix(vu, ncol = length(nmU), dimnames = list(NULL, nmU))
    ve. <- fit$chain$varE[draw.]
    tot. <- rowSums(vu) + ve.
    propm <- cbind(vu, Residual = ve.) / tot.
    qs <- t(apply(propm, 2, quantile, probs = probs, na.rm = TRUE))
    var_prop <- data.frame(component = colnames(propm),
                           proportion = colMeans(propm),
                           lower = qs[, 1], upper = qs[, 2], row.names = NULL)
    g_idx <- grep("^KG_", nmU)
    H2 <- if (length(g_idx)) rowSums(vu[, g_idx, drop = FALSE]) / tot. else NULL
    overall <- c(genomic_h2_overall = mean(h2_env, na.rm = TRUE), residual_var = vE.)
    if (!is.null(H2)) overall["broad_sense_H2"] <- mean(H2)
    overall <- c(overall, setNames(var_prop$proportion,
                                   paste0("prop_", var_prop$component)))
    ret$env_summary <- list(by_env = by_env, overall = overall)
    ret$variance_proportions <- var_prop
    if (!is.null(H2))
      attr(ret$variance_proportions, "H2") <-
        c(mean = mean(H2), lower = quantile(H2, probs[1], names = FALSE),
          upper = quantile(H2, probs[2], names = FALSE))
  }

  .step("[6/6] Variance components and heritability")
  ## The fallback note is deferred to the .done() suffix rather than printed
  ## inline, so it does not split the step's timing across two lines.
  .vc_fallback <- FALSE
  ret$varcomp <- tryCatch({
    m <- varcomp_method
    if (identical(m, "postmean") && is.null(fit$Uchain)) {
      ## plain `<-`: the tryCatch body evaluates in this frame, so this
      ## assigns the local. `<<-` would skip it and write a global instead.
      .vc_fallback <- TRUE
      m <- "gibbs"
    }
    varcomp_summary(ret, method = m, digits = digits)
  }, error = function(e) {
    warning("varcomp_summary() failed: ", conditionMessage(e)); NULL
  })
  ret$varcomp_gibbs <- tryCatch(
    varcomp_summary(ret, method = "gibbs", digits = digits),
    error = function(e) NULL)
  .done(if (!is.null(ret$varcomp))
          sprintf("method '%s', h2 = %.3f%s", ret$varcomp$method,
                  ret$varcomp$h2_entry_mean,
                  if (.vc_fallback) "  [postmean needs keep_effects = TRUE]" else "")
        else "unavailable")

  class(ret) <- "kernel_model"

  if (.v) {
    cat("---------------------------------------------------------\n")
    cat(sprintf(" Total running time: %.2fs  (sampler %.2fs, %.0f%%)\n",
                proc.time()[3] - .t_all, ret$runtime,
                100 * ret$runtime / max(proc.time()[3] - .t_all, 1e-9)))
    cat(sprintf(" Residual variance  : %.4f\n", ret$varE))
    cat("=========================================================\n\n")
  }
  ret
}


## ===========================================================================
## S3 helpers
## ===========================================================================

#' @export
print.kernel_model <- function(x, ...) {
  f <- x$fit %||% x$BGGE
  cat("<kernel_model> Bayesian G + GxE fit (BGGE sampler, Granato et al. 2018)\n")
  cat("  observations :", length(x$yHat), "\n")
  cat("  kernels      :", paste(f$kernel_names, collapse = ", "), "\n")
  cat("  iterations   :", f$ite, "(burn-in", f$burn,
      ", thinning", f$thin, ")\n")
  cat("  residual var :", round(x$varE, 4), "\n")
  cat("  runtime (s)  :", round(x$runtime, 2), "\n")
  cat("  variance components:\n")
  print(x$VarComp, row.names = FALSE)
  if (!is.null(x$predictions)) {
    cat("  reports in list: predictions, genetic_cor_env, genotype_effects,",
        "env_summary, variance_proportions\n")
  } else {
    cat("  (re-fit with keep_effects = TRUE to include the report tables)\n")
  }
  if (!is.null(x$varcomp)) {
    cat("  validated varcomp (method = '", x$varcomp$method, "'):\n", sep = "")
    cat(sprintf("    h2 entry-mean %.4f | h2 plot %.4f | rg %.4f\n",
                x$varcomp$h2_entry_mean, x$varcomp$h2_plot, x$varcomp$rg))
    if (isTRUE(x$varcomp$few_environments))
      cat("    ! q < 20 environments: sigma2_GE only approximate\n")
  }
  invisible(x)
}

#' @export
summary.kernel_model <- function(object, ...)
  if (!is.null(object$varcomp)) object$varcomp else object$VarComp


## ===========================================================================
## Multi-fit utilities
## ===========================================================================

#' k-fold cross-validation predictive ability (out-of-sample)
#'
#' In-sample yHat overstates accuracy. This masks a fold of phenotypes as NA,
#' refits, and correlates held-out predictions against observed values.
#'
#' The kernels are decomposed ONCE and reused across all folds: they depend only
#' on genotypes and environments, never on the phenotype vector, so masking y
#' does not invalidate them. This is the single largest saving in the package.
#'
#' @param y character. Name of the response column in \code{data}.
#' @param data data.frame with the response, environment and genotype columns.
#' @param random named list of kernels passed to \code{kernel_model}.
#' @param env character. Name of the environment column in \code{data}.
#' @param gid character. Name of the genotype column in \code{data}.
#' @param folds integer. Number of cross-validation folds. Default 5.
#' @param scheme "cv1" (random cells) or "cv0" (leave-one-environment-out).
#' @param seed integer or NULL. Optional seed for reproducible fold assignment.
#' @param tol numeric. Eigenvalue tolerance passed to the kernel decomposition.
#' @param scale_kernels logical. If \code{TRUE} (default) kernels are scaled before fitting.
#' @param ... further arguments passed to \code{kernel_model}.
#' @export
kernel_cv <- function(y, data, random, env, gid, folds = 5,
                      scheme = c("cv1", "cv0"), seed = NULL,
                      tol = 1e-10, scale_kernels = TRUE, ...) {
  scheme <- match.arg(scheme)
  if (!is.null(seed)) set.seed(as.integer(seed))
  n <- nrow(data); obs <- data[[y]]

  ## decompose once for every fold
  if (!.kd_valid(random, tol = tol, scale_kernels = scale_kernels, n = n)) {
    ne <- as.vector(table(factor(as.character(data[[env]]),
                                 levels = unique(as.character(data[[env]])))))
    random <- decompose_kernels(random, ne = ne, tol = tol,
                                scale_kernels = scale_kernels)
  }

  fold_id <- if (scheme == "cv1") {
    sample(rep(seq_len(folds), length.out = n))
  } else {
    envs <- unique(data[[env]])
    efold <- sample(rep(seq_len(folds), length.out = length(envs)))
    efold[match(data[[env]], envs)]
  }

  pred <- rep(NA_real_, n)
  per_fold <- data.frame(fold = seq_len(folds), accuracy = NA_real_, rmse = NA_real_)
  for (k in seq_len(folds)) {
    test <- which(fold_id == k)
    dtrain <- data; dtrain[[y]][test] <- NA
    fit <- suppressMessages(
      kernel_model(y = y, data = dtrain, random = random, env = env, gid = gid,
                   tol = tol, scale_kernels = scale_kernels, ...))
    pred[test] <- as.numeric(fit$yHat)[test]
    ok <- is.finite(obs[test]) & is.finite(pred[test])
    per_fold$accuracy[k] <- if (sum(ok) > 2) cor(obs[test][ok], pred[test][ok]) else NA
    per_fold$rmse[k] <- if (any(ok)) sqrt(mean((obs[test][ok] - pred[test][ok])^2)) else NA
  }
  ok <- is.finite(obs) & is.finite(pred)
  list(scheme = scheme, folds = folds,
       accuracy = cor(obs[ok], pred[ok]),
       rmse = sqrt(mean((obs[ok] - pred[ok])^2)),
       per_fold = per_fold, fold_id = fold_id, predicted = pred, observed = obs)
}


#' Multi-chain fit with convergence diagnostics (Gelman-Rubin R-hat, ESS)
#'
#' Runs independent chains with different seeds and reports R-hat and a crude
#' ESS for the residual and each kernel-scale parameter. R-hat above 1.01 or a
#' tiny ESS means the chain has not mixed: increase `iterations`.
#'
#' The decomposition is shared across chains.
#'
#' @param ... arguments passed to \code{kernel_model} (e.g. \code{y}, \code{data},
#'   \code{random}, \code{env}, \code{gid}).
#' @param n_chains integer. Number of independent chains. Default 3.
#' @param seed integer. Base seed; chain \code{c} uses \code{seed + c}. Default 1.
#' @param tol numeric. Eigenvalue tolerance passed to the kernel decomposition.
#' @param scale_kernels logical. If \code{TRUE} (default) kernels are scaled before fitting.
#' @export
kernel_model_mc <- function(..., n_chains = 3, seed = 1,
                            tol = 1e-10, scale_kernels = TRUE) {
  dots <- list(...)
  if (!is.null(dots$random) && !is.null(dots$data)) {
    n <- nrow(dots$data)
    if (!.kd_valid(dots$random, tol, scale_kernels, n)) {
      ev <- as.character(dots$data[[dots$env]])
      dots$random <- decompose_kernels(dots$random,
                                       ne = as.vector(table(factor(ev, levels = unique(ev)))),
                                       tol = tol, scale_kernels = scale_kernels)
    }
  }
  fits <- vector("list", n_chains)
  for (c in seq_len(n_chains))
    fits[[c]] <- suppressMessages(do.call(kernel_model,
                                          c(dots, list(seed = seed + c,
                                                       tol = tol,
                                                       scale_kernels = scale_kernels))))

  get_scalar_chains <- function(f0) {
    b <- f0$fit %||% f0$BGGE
    ch <- b$chain
    draw <- which(seq_len(b$ite) %% b$thin == 0) > b$burn
    m <- cbind(varE = ch$varE[draw])
    for (nm in names(ch$K)) m <- cbind(m, ch$K[[nm]]$varU[draw])
    colnames(m) <- c("varE", names(ch$K))
    m
  }
  chains <- lapply(fits, get_scalar_chains)
  L <- min(vapply(chains, nrow, integer(1)))
  chains <- lapply(chains, function(m) m[seq_len(L), , drop = FALSE])
  params <- colnames(chains[[1]])

  rhat <- ess <- setNames(numeric(length(params)), params)
  for (p in params) {
    mat <- vapply(chains, function(m) m[, p], numeric(L))
    m_means <- colMeans(mat)
    B <- L * var(m_means); W <- mean(apply(mat, 2, var))
    varplus <- (1 - 1 / L) * W + B / L
    rhat[p] <- if (W > 0) sqrt(varplus / W) else NA
    ac1 <- mean(apply(mat, 2, function(x) {
      x <- x - mean(x); if (sum(x^2) == 0) return(0)
      sum(x[-1] * x[-length(x)]) / sum(x^2)
    }))
    ess[p] <- (n_chains * L) * (1 - ac1) / (1 + ac1)
  }

  out <- fits[[1]]
  out$diagnostics <- data.frame(parameter = params, Rhat = round(rhat, 4),
                                ESS = round(ess, 1), row.names = NULL)
  out$all_chains <- fits
  if (any(rhat > 1.05, na.rm = TRUE))
    warning("Some R-hat > 1.05: chains have not converged. Increase `iterations`.")
  out
}


## ===========================================================================
## Environment clustering
## ===========================================================================

#' Group environments into mega-environments (clusters)
#'
#' @description
#' Delineates mega-environments by clustering the \strong{genetic correlation}
#' among environments, i.e. how similarly genotypes rank from one environment to
#' another. Environments in which genotypes respond alike are placed in the same
#' cluster, which is the operational definition of a mega-environment in a
#' target population of environments (TPE).
#'
#' The input is either a fitted \code{\link{kernel_model}} (carrying the
#' environment genetic-correlation matrix in \code{$genetic_cor_env}, produced
#' with \code{keep_effects = TRUE}) or any square environment correlation matrix
#' you already have (for example from \code{\link{env_cor}} or an external GxE
#' analysis).
#'
#' @details
#' The genetic correlation \eqn{r_g} is turned into a distance and the chosen
#' algorithm is run on that distance:
#' \itemize{
#'   \item \code{distance = "one_minus"} uses \eqn{d = 1 - r_g} (default).
#'   \item \code{distance = "sqrt"} uses \eqn{d = \sqrt{1 - r_g}}, which spreads
#'         out highly correlated environments.
#' }
#' Correlations are clamped to \eqn{[-1, 1]} and the distance matrix is
#' symmetrised before clustering. Three algorithms are available through
#' \code{method}: agglomerative hierarchical clustering (\code{"hclust"},
#' controlled by \code{hclust_linkage}), partitioning around medoids
#' (\code{"pam"}), and \code{"kmeans"} run on a classical MDS embedding of the
#' distance matrix.
#'
#' When \code{k} is \code{NULL} the number of clusters is selected automatically
#' by maximising the \strong{average silhouette width} over the candidate values
#' in \code{k_range}; supplying \code{k} skips this search. At least three
#' environments are required.
#'
#' @param object a \code{\link{kernel_model}} fit carrying
#'   \code{$genetic_cor_env} (fit with \code{keep_effects = TRUE}), or a square
#'   environment-by-environment correlation matrix.
#' @param k integer or NULL. Number of clusters. When \code{NULL} (default) it is
#'   chosen by maximum average silhouette width over \code{k_range}.
#' @param method character. Clustering algorithm: \code{"hclust"} (default),
#'   \code{"pam"} or \code{"kmeans"}.
#' @param distance character. Distance built from the genetic correlation:
#'   \code{"one_minus"} (\eqn{1 - r_g}, default) or \code{"sqrt"}
#'   (\eqn{\sqrt{1 - r_g}}).
#' @param hclust_linkage character. Linkage passed to \code{\link[stats]{hclust}}
#'   when \code{method = "hclust"}. Default \code{"ward.D2"}.
#' @param k_range integer vector. Candidate numbers of clusters tested when
#'   \code{k} is \code{NULL}. Default \code{2:6}.
#' @param seed integer or NULL. Optional seed for reproducibility of the
#'   randomised starts in \code{"pam"} and \code{"kmeans"}.
#'
#' @return A list with:
#' \describe{
#'   \item{clusters}{named integer vector giving the cluster of each environment.}
#'   \item{k}{the number of clusters used (selected or supplied).}
#'   \item{method}{the clustering algorithm used.}
#'   \item{distance}{the symmetric distance matrix that was clustered.}
#'   \item{silhouette}{per-environment silhouette widths plus their
#'     \code{average} (higher is a tighter, better-separated solution).}
#'   \item{medoids}{representative environment of each cluster
#'     (\code{"pam"}/\code{"kmeans"}; \code{NULL} for \code{"hclust"}).}
#'   \item{hclust}{the \code{\link[stats]{hclust}} tree for \code{method = "hclust"},
#'     otherwise \code{NULL}.}
#'   \item{table}{a data.frame (\code{env}, \code{cluster}, \code{sil_width})
#'     ordered by cluster and decreasing silhouette width.}
#' }
#'
#' @seealso \code{\link{kernel_model}} (fit with \code{keep_effects = TRUE} to
#'   obtain \code{$genetic_cor_env}), \code{\link{kernel_model_clustered}} to fit
#'   a model informed by the clusters, \code{\link{env_cluster}} for clustering
#'   from an environmental covariable matrix, and \code{\link{env_cor}}.
#'
#' @examples
#' ## ------------------------------------------------------------------
#' ## 1. Cluster a genetic-correlation matrix directly (no model needed)
#' ## ------------------------------------------------------------------
#' ## Six environments forming two obvious groups
#' R <- matrix(0.2, 6, 6)
#' R[1:3, 1:3] <- 0.9
#' R[4:6, 4:6] <- 0.85
#' diag(R) <- 1
#' dimnames(R) <- list(paste0("E", 1:6), paste0("E", 1:6))
#'
#' ## Automatic number of clusters via silhouette width
#' cl <- cluster_environments(R)
#' cl$k
#' cl$table
#' cl$silhouette["average"]
#'
#' ## ------------------------------------------------------------------
#' ## 2. Force a specific number of clusters
#' ## ------------------------------------------------------------------
#' cluster_environments(R, k = 2)$clusters
#'
#' ## ------------------------------------------------------------------
#' ## 3. Try different algorithms and distances
#' ## ------------------------------------------------------------------
#' cluster_environments(R, k = 2, method = "pam", seed = 1)$medoids
#' cluster_environments(R, method = "kmeans", k_range = 2:4, seed = 1)$clusters
#' cluster_environments(R, distance = "sqrt",
#'                      hclust_linkage = "complete")$clusters
#'
#' \dontrun{
#' ## ------------------------------------------------------------------
#' ## 4. Full workflow from phenotypes, kernels and a fitted model
#' ## ------------------------------------------------------------------
#' data("maizeYield"); data("maizeG"); data("maizeWTH")
#'
#' ECs <- W_matrix(env.data = maizeWTH[maizeWTH$daysFromStart < 100, ],
#'                 var.id = c("FRUE", "PETP", "SRAD", "T2M_MAX"))
#' K <- get_kernel(K_G = list(G = maizeG),
#'                 K_E = list(W = env_kernel(env.data = ECs)[[2]]),
#'                 data = maizeYield, model = "RNMM")
#'
#' ## keep_effects = TRUE is required for $genetic_cor_env
#' fit <- kernel_model(y = "value", env = "env", gid = "gid",
#'                     data = maizeYield, random = K,
#'                     keep_effects = TRUE)
#'
#' megaenv <- cluster_environments(fit)
#' megaenv$table
#'
#' ## 5. Feed the clusters into a cluster-aware model
#' fit_cl <- kernel_model_clustered(y = "value", data = maizeYield,
#'                                  random = K, env = "env", gid = "gid",
#'                                  clusters = megaenv, use = "block")
#' }
#' @export
cluster_environments <- function(object, k = NULL,
                                 method = c("hclust", "pam", "kmeans"),
                                 distance = c("one_minus", "sqrt"),
                                 hclust_linkage = "ward.D2",
                                 k_range = 2:6, seed = NULL) {
  method <- match.arg(method); distance <- match.arg(distance)
  if (!is.null(seed)) set.seed(as.integer(seed))

  if (inherits(object, "kernel_model")) {
    if (is.null(object$genetic_cor_env))
      stop("No $genetic_cor_env in the fit. Re-fit kernel_model(..., keep_effects = TRUE).")
    R <- object$genetic_cor_env$cor
  } else if (is.matrix(object) && nrow(object) == ncol(object)) {
    R <- object
  } else stop("`object` must be a kernel_model fit or a square correlation matrix.")

  envs <- rownames(R); if (is.null(envs)) envs <- paste0("E", seq_len(nrow(R)))
  ne <- nrow(R)
  if (ne < 3) stop("Need at least 3 environments to cluster; got ", ne, ".")

  R <- pmin(pmax(R, -1), 1)
  D <- if (distance == "sqrt") sqrt(pmax(1 - R, 0)) else (1 - R)
  diag(D) <- 0; D <- (D + t(D)) / 2
  dobj <- stats::as.dist(D)

  silhouette_width <- function(labels, Dm) {
    n <- length(labels); s <- numeric(n)
    for (i in seq_len(n)) {
      ci <- labels[i]; same <- labels == ci & seq_len(n) != i
      a <- if (any(same)) mean(Dm[i, same]) else 0
      others <- setdiff(unique(labels), ci)
      b <- if (length(others))
        min(vapply(others, function(cc) mean(Dm[i, labels == cc]), numeric(1))) else 0
      s[i] <- if (max(a, b) > 0) (b - a) / max(a, b) else 0
    }
    s
  }

  pam_base <- function(Dm, k) {
    n <- nrow(Dm); medoids <- sample(seq_len(n), k)
    assign_clusters <- function(med) apply(Dm[, med, drop = FALSE], 1, which.min)
    lab <- assign_clusters(medoids)
    for (it in seq_len(100)) {
      new_med <- medoids
      for (c in seq_len(k)) {
        members <- which(lab == c)
        if (!length(members)) next
        costs <- vapply(members, function(m) sum(Dm[m, members]), numeric(1))
        new_med[c] <- members[which.min(costs)]
      }
      new_lab <- assign_clusters(new_med)
      if (identical(new_lab, lab) && identical(new_med, medoids)) break
      medoids <- new_med; lab <- new_lab
    }
    list(clustering = lab, medoids = medoids)
  }

  fit_k <- function(k) {
    if (method == "hclust") {
      hc <- stats::hclust(dobj, method = hclust_linkage)
      list(labels = stats::cutree(hc, k = k), hc = hc, medoids = NULL)
    } else if (method == "pam") {
      pr <- pam_base(D, k)
      list(labels = pr$clustering, hc = NULL, medoids = envs[pr$medoids])
    } else {
      mds <- stats::cmdscale(dobj, k = min(ne - 1, max(2, k)))
      km <- stats::kmeans(mds, centers = k, nstart = 25)
      med <- vapply(seq_len(k), function(c) {
        idx <- which(km$cluster == c); cen <- km$centers[c, , drop = FALSE]
        idx[which.min(rowSums((mds[idx, , drop = FALSE] -
              matrix(cen, length(idx), ncol(mds), byrow = TRUE))^2))]
      }, integer(1))
      list(labels = km$cluster, hc = NULL, medoids = envs[med])
    }
  }

  if (is.null(k)) {
    cand <- k_range[k_range >= 2 & k_range <= ne - 1]
    if (!length(cand)) cand <- 2
    avg_sil <- vapply(cand, function(kk) {
      lab <- fit_k(kk)$labels
      if (length(unique(lab)) < 2) return(-Inf)
      mean(silhouette_width(lab, D))
    }, numeric(1))
    k <- cand[which.max(avg_sil)]
  }

  res <- fit_k(k); labels <- res$labels; names(labels) <- envs
  sil <- silhouette_width(labels, D); names(sil) <- envs
  tab <- data.frame(env = envs, cluster = as.integer(labels),
                    sil_width = round(sil, 4), row.names = NULL)
  tab <- tab[order(tab$cluster, -tab$sil_width), ]

  list(clusters = labels, k = k, method = method, distance = D,
       silhouette = c(setNames(round(sil, 4), envs), average = round(mean(sil), 4)),
       medoids = res$medoids, hclust = res$hc, table = tab)
}


#' Fit kernel_model() informed by an environment-cluster table
#'
#' "fixed" adds mega-environment means as fixed effects; "block" zeroes
#' between-cluster entries of the GxE kernel so genotypes may rank differently
#' across mega-environments; "both" does each.
#'
#' Block-diagonalisation MODIFIES the GxE kernels, so any cached decomposition
#' for those kernels is dropped and recomputed.
#'
#' @param y character. Name of the response column in \code{data}.
#' @param data data.frame with the response, environment and genotype columns.
#' @param random named list of kernels passed to \code{kernel_model}.
#' @param env character. Name of the environment column in \code{data}.
#' @param gid character. Name of the genotype column in \code{data}.
#' @param clusters an environment-cluster mapping: a \code{cluster_environments} result,
#'   a data.frame with columns \code{env} and \code{cluster}, or a named vector.
#' @param use character. How clusters enter the model: \code{"block"} (default),
#'   \code{"fixed"} or \code{"both"}.
#' @param gxe_pattern character. Regular expression selecting the GxE kernels to
#'   block-diagonalise. Default \code{"^KGE_|GE"}.
#' @param fixed optional fixed-effects specification passed to \code{kernel_model}.
#' @param ... further arguments passed to \code{kernel_model}.
#' @export
kernel_model_clustered <- function(y, data, random, env, gid,
                                   clusters, use = c("block", "fixed", "both"),
                                   gxe_pattern = "^KGE_|GE",
                                   fixed = NULL, ...) {
  use <- match.arg(use)

  if (is.list(clusters) && !is.data.frame(clusters) && !is.null(clusters$clusters))
    clusters <- clusters$clusters
  if (is.data.frame(clusters)) {
    if (!all(c("env", "cluster") %in% names(clusters)))
      stop("`clusters` data.frame must have columns `env` and `cluster`.")
    cl_map <- setNames(as.character(clusters$cluster), as.character(clusters$env))
  } else if (!is.null(names(clusters))) {
    cl_map <- setNames(as.character(clusters), names(clusters))
  } else stop("`clusters` must be a data.frame (env, cluster), a named vector, ",
              "or a cluster_environments() result.")

  env_vec <- as.character(data[[env]])
  miss <- setdiff(unique(env_vec), names(cl_map))
  if (length(miss))
    stop("These environments have no cluster assignment: ", paste(miss, collapse = ", "))
  obs_cluster <- cl_map[env_vec]
  k_used <- length(unique(cl_map))

  if (use %in% c("fixed", "both")) {
    cl_f <- factor(obs_cluster)
    Xcl <- stats::model.matrix(~ cl_f)[, -1, drop = FALSE]
    colnames(Xcl) <- paste0("cluster_", levels(cl_f)[-1])
    fixed <- if (is.null(fixed)) Xcl else cbind(as.matrix(fixed), Xcl)
  }

  random2 <- random
  if (use %in% c("block", "both")) {
    same_cluster <- outer(obs_cluster, obs_cluster, `==`) * 1
    hit <- grep(gxe_pattern, names(random2))
    if (length(hit) == 0)
      warning("No kernels matched gxe_pattern = '", gxe_pattern,
              "'; nothing was block-diagonalised.")
    for (i in hit) {
      random2[[i]]$Kernel <- random2[[i]]$Kernel * same_cluster
      random2[[i]]$Type   <- "D"
      random2[[i]]$decomp <- NULL          # kernel changed: invalidate cache
      random2[[i]]$base <- random2[[i]]$idx <- NULL
    }
    attr(random2, "kd_meta") <- NULL       # force revalidation
  }

  fit <- kernel_model(y = y, data = data, random = random2, fixed = fixed,
                      env = env, gid = gid, ...)
  fit$cluster_info <- list(
    mapping = cl_map, k = k_used, mechanism = use, gxe_pattern = gxe_pattern,
    n_fixed_cols = if (use %in% c("fixed", "both")) k_used - 1L else 0L,
    blocked_kernels = if (use %in% c("block", "both"))
      names(random2)[grep(gxe_pattern, names(random2))] else character(0))
  fit
}


## ############################################################################
## PART 3 -- PREDICTION AT UNTESTED ENVIRONMENTS   scan_untested_envs()
## ############################################################################

## ============================================================================
## scan_untested_envs.R
## ============================================================================

#' @keywords internal
.pne_standardize <- function(W, center, scale) {
  sc <- scale; sc[!is.finite(sc) | sc == 0] <- 1
  out <- sweep(sweep(W, 2, center, "-"), 2, sc, "/")
  out[!is.finite(out)] <- 0
  out
}

#' Ridge-stabilised inverse, falling back to a pseudo-inverse.
#' @keywords internal
.pne_inv <- function(K, lambda = 1e-6) {
  n <- nrow(K)
  A <- K + diag(lambda * mean(diag(K)) + 1e-10, n)
  out <- tryCatch(solve(A), error = function(e) NULL)
  if (is.null(out)) {
    ei <- eigen(A, symmetric = TRUE)
    keep <- ei$values > max(ei$values) * 1e-10
    out <- ei$vectors[, keep, drop = FALSE] %*%
           diag(1 / ei$values[keep], sum(keep)) %*%
           t(ei$vectors[, keep, drop = FALSE])
  }
  out
}

#' @keywords internal
.pne_chain <- function(fit, pattern) {
  U <- fit$fit$Uchain
  if (is.null(U)) U <- fit$BGGE$Uchain
  if (is.null(U)) return(NULL)
  sel <- grep(pattern, names(U))
  if (!length(sel)) return(NULL)
  Reduce(`+`, U[sel])
}

#' @title Predict Trained Genotypes at Untested Environments Without Refitting
#'
#' @description
#' Transfers a fitted \code{\link{kernel_model}} to environments that were never
#' planted -- new locations, new planting dates, or both -- using only their
#' environmental covariates. No re-fitting and no phenotypes at the new sites
#' are required.
#'
#' The environmental main effect estimated at the training sites is
#' \emph{kriged} onto each new site through a linear environmental kernel, and
#' the whole posterior is carried draw-by-draw so that credible intervals,
#' win probabilities and an extrapolation diagnostic are all internally
#' consistent.
#'
#' \strong{The transfer.} Covariates are standardised on the training moments,
#' \deqn{\tilde{W} = (W - \bar{W}_{tr}) \; / \; s_{tr},}
#' and a linear (product) environmental kernel is formed from the \eqn{p}
#' covariates,
#' \deqn{K_{tt} = \tilde{W}_{tr}\tilde{W}_{tr}^{\top}/p, \quad
#'       K_{nt} = \tilde{W}_{n}\tilde{W}_{tr}^{\top}/p, \quad
#'       K_{nn} = \tilde{W}_{n}\tilde{W}_{n}^{\top}/p.}
#' The new-site environmental effects are a kriging (Gaussian-process
#' conditional-mean) map of the training effects,
#' \deqn{H = K_{nt}\,(K_{tt} + \lambda I)^{-1}, \qquad
#'       \hat{e}_{n}^{(d)} = H\,e_{tr}^{(d)},}
#' applied to every retained MCMC draw \eqn{d}. When a GxE kernel is present its
#' interaction effects are kriged the same way. The ridge \eqn{\lambda}
#' stabilises the inverse; \eqn{H} is \emph{not} constrained to be positive, and
#' large negative row weights are the signal that a site is being reached by
#' extrapolation (see \code{weight_negativity} below).
#'
#' \strong{The prediction.} For draw \eqn{d}, genotype \eqn{g} and new site
#' \eqn{j},
#' \deqn{\hat{y}^{(d)}_{gj} = \mu + u_g^{(d)} + \hat{e}^{(d)}_{j}
#'        + \widehat{(ge)}^{(d)}_{gj},}
#' with \eqn{\mu} the training grand mean. Point predictions are posterior means
#' over draws and the credible interval uses the empirical quantiles at
#' \eqn{\alpha=(1-\texttt{level})/2} and \eqn{1-\alpha}. The probability that a
#' genotype is the best at a site is the fraction of draws in which it attains
#' the row maximum, so it respects the genotype correlation \emph{within} a
#' draw. Per-site genomic heritability is reported as
#' \eqn{h^2 = \sigma^2_g / (\sigma^2_g + \sigma^2_e)}.
#'
#' @param object a fitted model from \code{\link{kernel_model}}, fitted with
#'   \code{keep_effects = TRUE}. The per-draw effect chains (not posterior
#'   means) are required; without them the intervals would be fabricated.
#' @param W_train training environmental covariate matrix,
#'   environments \eqn{\times} covariates (\eqn{e \times p}). If it has no
#'   \code{rownames} they are taken to be the sorted training environments.
#' @param W_new new-environment covariate matrix, sites \eqn{\times} covariates
#'   (\eqn{m \times p}). Must have the \emph{same columns} as \code{W_train}.
#'   Unnamed rows are labelled \code{NewEnv1..NewEnvm}.
#' @param K_G optional genomic relationship matrix. Used only to inform the GxE
#'   transfer; when omitted, GxE is kriged from environmental similarity alone
#'   (a warning is issued).
#' @param standardize logical. Standardise covariates on the training mean and
#'   SD before building the kernels. Default \code{TRUE} and strongly advised
#'   when covariates are on different scales.
#' @param lambda numeric ridge added to \eqn{K_{tt}} before inversion, as a
#'   fraction of its mean diagonal. Default \code{1e-6}. Raise it if the kernel
#'   is near-singular (highly collinear covariates).
#' @param level credible-interval mass. Default \code{0.95}.
#' @param max_draws integer. Cap on the number of MCMC draws used; draws are
#'   thinned by even spacing when the chain is longer. Default \code{500}.
#' @param envelope what to do about sites that fall outside the training
#'   covariate envelope: \code{"keep"} (default, return everything unchanged),
#'   \code{"flag"} (return everything, warn), \code{"mask"} (set flagged
#'   predictions to \code{NA} but keep the matrix structure so indexing by
#'   environment name still resolves), or \code{"drop"} (remove flagged
#'   rows/columns and record them in \code{$excluded}).
#' @param envelope_level which sites count as flagged: \code{"gross"} (default)
#'   or \code{"mild"}. Only \code{"gross"} separates cleanly in calibration;
#'   \code{"mild"} acts on a weak signal and is not recommended as a gate.
#'
#' @details
#' \strong{Extrapolation diagnostics.} The scan reports \emph{where} each new
#' site sits relative to the training covariate cloud. Several complementary
#' measures are returned because no single one is sufficient:
#' \itemize{
#'   \item \emph{Cosine similarity} to the nearest training site,
#'     \eqn{s_{ij} = K_{nt,ij}/\sqrt{K_{nn,ii}\,K_{tt,jj}}}, summarised as
#'     \code{max_similarity} and \code{mean_similarity}.
#'   \item A \emph{Mahalanobis-type distance} in the leading principal-component
#'     subspace (enough PCs to reach 95\% variance),
#'     \eqn{d_i = \sqrt{\sum_k z_{ik}^2}}, divided by the 95th percentile of the
#'     training distances to give \code{mahalanobis_ratio}.
#'   \item A \emph{marginal range test} per covariate (\code{prop_outside},
#'     \code{max_exceedance_sd}). This is marginal: a site inside every
#'     covariate's range individually can still occupy a \emph{combination} that
#'     never occurred, which only the distance and weight measures catch.
#'   \item \emph{Kriging-weight geometry}: \code{weight_negativity}
#'     \eqn{= \sum_j \min(H_{ij},0)} and the effective number of training sites
#'     leaned on, \eqn{n_{\mathrm{eff},i} = (\sum_j |H_{ij}|)^2 / \sum_j H_{ij}^2}.
#'   \item A combined \code{divergence_index}
#'     \eqn{= \min\!\big(1, \max(0,\; 0.5(1-s_{\max}) + 0.5\,t)\big)} with
#'     \eqn{t = \mathrm{clamp}((\texttt{mahalanobis\_ratio}-0.5)/2,\,0,\,1)}.
#' }
#' Sites are labelled \code{"interpolation"}, \code{"edge"},
#' \code{"mild_extrapolation"} or \code{"gross_extrapolation"} from these.
#'
#' \strong{Position is geometry, not a correctness guarantee.} Calibration
#' against 60 held-out environments found the envelope measures correlate only
#' \eqn{r \approx -0.37} to \eqn{-0.66} with realised accuracy. Treat the
#' position label as a description of where a site lies, never as a verdict that
#' the prediction is right or wrong. The model-internal \code{ci_width_ratio}
#' was the strongest single accuracy correlate (\eqn{r=-0.66}); to know your
#' true accuracy, run leave-environments-out cross-validation on your own data.
#'
#' @return
#' An object of class \code{"scan_untested_envs"}: a list with
#' \describe{
#'   \item{\code{yHat_matrix}, \code{yHat_lower}, \code{yHat_upper}}{genotype
#'     \eqn{\times} new-site matrices of the posterior mean and credible bounds.}
#'   \item{\code{yHat}}{the same in long form, with \code{ci_width}.}
#'   \item{\code{yHat_prob_win}}{probability each genotype is best at each site.}
#'   \item{\code{genetic_var}}{per-site genetic variance, its interval and
#'     \code{genomic_h2}.}
#'   \item{\code{vcov}, \code{cor}}{between-new-site genetic covariance and
#'     correlation matrices.}
#'   \item{\code{interpolation}, \code{divergence}, \code{diagnostics}}{the
#'     extrapolation measures described above, one row per new site.}
#'   \item{\code{exceedance_detail}}{one row per (site \eqn{\times} offending
#'     covariate) that fell outside the training range.}
#'   \item{\code{excluded}}{sites masked or dropped by \code{envelope}, else
#'     \code{NULL}.}
#'   \item{\code{meta}}{dimensions, settings, similarity and kriging weights.}
#' }
#'
#' @examples
#' \dontrun{
#' data("maizeWTH"); data("maizeYield"); data("maizeG")
#'
#' ## Training covariates (first 100 days after sowing)
#' vars <- c("FRUE", "PETP", "GDD", "T2M_MAX")
#' Wtr  <- W_matrix(env.data = maizeWTH[maizeWTH$daysFromStart < 100, ],
#'                  var.id = vars, statistic = "mean")
#'
#' ## Fit KEEPING the effect chains -- required by the scan
#' K   <- get_kernel(K_G = list(G = maizeG),
#'                   K_E = list(W = env_kernel(env.data = Wtr)[[2]]),
#'                   data = maizeYield, model = "RNMM",
#'                   env = "env", gid = "gid", y = "value")
#' fit <- kernel_model(y = "value", data = maizeYield, random = K,
#'                     env = "env", gid = "gid", keep_effects = TRUE)
#'
#' ## Predict the trained genotypes at new sites (reuse the training W here)
#' sc <- scan_untested_envs(fit, W_train = Wtr, W_new = Wtr, K_G = maizeG)
#' sc                                   # print method: per-site diagnostics
#' head(sc$yHat)                        # predictions + credible intervals
#' sc$interpolation[, c("env", "position", "mahalanobis_ratio")]
#' }
#'
#' @seealso \code{\link{kernel_model}}, \code{\link{get_kernel}},
#'   \code{\link{grid_scan}}, \code{\link{map_scan}},
#'   \code{\link{plot_planting_window}}, \code{\link{best_planting_date}}
#'
#' @references
#' Costa-Neto, G., Fritsche-Neto, R. & Crossa, J. (2021) Nonlinear kernels,
#' dominance, and envirotyping data increase the accuracy of genome-based
#' prediction in multi-environment trials. \emph{Heredity} 126, 92-106.
#'
#' Jarquin, D. et al. (2014) A reaction norm model for genomic selection using
#' high-dimensional genomic and environmental data. \emph{Theoretical and
#' Applied Genetics} 127, 595-607.
#'
#' Cressie, N. (1993) \emph{Statistics for Spatial Data}. Wiley. (kriging as a
#' Gaussian-process conditional mean.)
#'
#' @export
scan_untested_envs <- function(object, W_train, W_new, K_G = NULL,
                          standardize = TRUE, lambda = 1e-6, level = 0.95,
                          max_draws = 500,
                          envelope = c("keep", "flag", "mask", "drop"),
                          envelope_level = c("gross", "mild")) {

  ## envelope: what to DO about sites outside the training covariate envelope.
  ##   keep  - return everything unchanged (default; current behaviour)
  ##   flag  - return everything, warn prominently about flagged sites
  ##   mask  - set flagged predictions to NA, PRESERVING matrix structure so
  ##           downstream indexing (pred$yHat_matrix[, "SiteX"]) still works
  ##   drop  - remove flagged rows/columns; recorded in $excluded
  ## envelope_level: which sites count as flagged. "gross" (default) is the
  ##   only regime that separates unambiguously in calibration; "mild" acts on
  ##   a signal correlating only r ~ -0.37 with accuracy and is NOT recommended.
  envelope <- match.arg(envelope)
  envelope_level <- match.arg(envelope_level)

  .uch <- object$fit$Uchain
  if (is.null(.uch)) .uch <- object$BGGE$Uchain
  if (is.null(.uch))
    stop("No effect chains. Re-fit with kernel_model(..., keep_effects = TRUE).")
  if (is.null(object$meta))
    stop("Fit has no $meta. Re-fit with the current kernel_model().")

  W_train <- as.matrix(W_train); W_new <- as.matrix(W_new)
  if (ncol(W_train) != ncol(W_new))
    stop("W_train and W_new must have identical covariate columns (",
         ncol(W_train), " vs ", ncol(W_new), ").")
  if (is.null(rownames(W_new)))
    rownames(W_new) <- paste0("NewEnv", seq_len(nrow(W_new)))

  env_tr <- as.character(object$meta$env)
  gid_tr <- as.character(object$meta$gid)
  envs <- sort(unique(env_tr)); gids <- sort(unique(gid_tr))
  e <- length(envs); g <- length(gids); m <- nrow(W_new); p <- ncol(W_train)

  if (is.null(rownames(W_train))) {
    if (nrow(W_train) != e)
      stop("W_train has no rownames and ", nrow(W_train), " rows != ", e,
           " training environments.")
    rownames(W_train) <- envs
  }
  miss <- setdiff(envs, rownames(W_train))
  if (length(miss))
    stop("W_train missing training environments: ", paste(miss, collapse = ", "))
  W_train <- W_train[envs, , drop = FALSE]

  ctr <- colMeans(W_train); sdv <- apply(W_train, 2, stats::sd)
  Wtr <- if (standardize) .pne_standardize(W_train, ctr, sdv) else W_train
  Wnw <- if (standardize) .pne_standardize(W_new,   ctr, sdv) else W_new

  KE_tt <- tcrossprod(Wtr) / p
  KE_nt <- Wnw %*% t(Wtr) / p
  KE_nn <- tcrossprod(Wnw) / p

  Wgt <- KE_nt %*% .pne_inv(KE_tt, lambda)

  UE  <- .pne_chain(object, "^KE_")
  UG  <- .pne_chain(object, "^KG_")
  UGE <- .pne_chain(object, "^KGE_")
  if (is.null(UE)) stop("No environment kernel (^KE_) in the fit.")
  nd <- nrow(UE)
  use <- if (nd > max_draws) round(seq(1, nd, length.out = max_draws)) else seq_len(nd)
  nd_use <- length(use)

  ie <- match(env_tr, envs); ig <- match(gid_tr, gids)

  env_tr_eff <- t(apply(UE[use, , drop = FALSE], 1,
                        function(u) tapply(u, ie, mean)[as.character(seq_len(e))]))
  env_tr_eff[!is.finite(env_tr_eff)] <- 0
  env_new <- env_tr_eff %*% t(Wgt)

  if (!is.null(UG)) {
    gen_eff <- t(apply(UG[use, , drop = FALSE], 1,
                       function(u) tapply(u, ig, mean)[as.character(seq_len(g))]))
    gen_eff[!is.finite(gen_eff)] <- 0
  } else gen_eff <- matrix(0, nd_use, g)

  ge_new <- array(0, dim = c(nd_use, g, m))
  have_ge <- !is.null(UGE)
  if (have_ge) {
    if (is.null(K_G))
      warning("K_G not supplied: GxE krigged from environment similarity only.")
    cell <- (ie - 1L) * g + ig          ## column-major (fixed)
    for (d in seq_len(nd_use)) {
      M <- matrix(0, g, e)
      M[cell] <- UGE[use[d], ]
      ge_new[d, , ] <- M %*% t(Wgt)
    }
  }

  mu <- mean(object$meta$y, na.rm = TRUE)
  yd <- array(0, dim = c(nd_use, g, m))
  for (d in seq_len(nd_use))
    yd[d, , ] <- mu + outer(gen_eff[d, ], env_new[d, ], "+") + ge_new[d, , ]

  a <- (1 - level) / 2
  yhat <- apply(yd, c(2, 3), mean)
  ylo  <- apply(yd, c(2, 3), stats::quantile, probs = a,     na.rm = TRUE)
  yhi  <- apply(yd, c(2, 3), stats::quantile, probs = 1 - a, na.rm = TRUE)
  dimnames(yhat) <- dimnames(ylo) <- dimnames(yhi) <- list(gids, rownames(W_new))

  ## posterior probability each genotype is the best at each site, taken from
  ## the joint draws so it respects genotype correlation within a draw
  pwin <- matrix(0, g, m, dimnames = list(gids, rownames(W_new)))
  for (j in seq_len(m)) {
    wj <- max.col(matrix(yd[, , j], nrow = nd_use, ncol = g),
                  ties.method = "first")
    pwin[, j] <- tabulate(wj, nbins = g) / nd_use
  }

  yHat_long <- data.frame(
    gid = rep(gids, times = m), env = rep(rownames(W_new), each = g),
    yHat = as.vector(yhat), lower = as.vector(ylo), upper = as.vector(yhi),
    ci_width = as.vector(yhi - ylo), stringsAsFactors = FALSE)

  gv <- array(0, dim = c(nd_use, g, m))
  for (d in seq_len(nd_use))
    gv[d, , ] <- matrix(gen_eff[d, ], g, m) + ge_new[d, , ]

  gvar_draw <- matrix(NA_real_, nd_use, m)
  for (d in seq_len(nd_use))
    gvar_draw[d, ] <- apply(matrix(gv[d, , ], g, m), 2, stats::var)

  genetic_var <- data.frame(
    env = rownames(W_new),
    genetic_var = colMeans(gvar_draw, na.rm = TRUE),
    lower = apply(gvar_draw, 2, stats::quantile, probs = a,     na.rm = TRUE),
    upper = apply(gvar_draw, 2, stats::quantile, probs = 1 - a, na.rm = TRUE),
    stringsAsFactors = FALSE, row.names = NULL)

  s2e <- if (!is.null(object$varcomp)) object$varcomp$var_resid else object$varE
  genetic_var$genomic_h2 <- genetic_var$genetic_var / (genetic_var$genetic_var + s2e)

  vc_sum <- matrix(0, m, m); cr_sum <- matrix(0, m, m); nv <- 0
  for (d in seq_len(nd_use)) {
    S <- stats::cov(matrix(gv[d, , ], g, m))
    if (all(is.finite(S))) {
      vc_sum <- vc_sum + S
      dd <- sqrt(diag(S)); dd[dd <= 0] <- NA
      R <- S / outer(dd, dd); R[!is.finite(R)] <- 0
      cr_sum <- cr_sum + R; nv <- nv + 1
    }
  }
  vcov_new <- if (nv) vc_sum / nv else matrix(NA_real_, m, m)
  cor_new  <- if (nv) cr_sum / nv else matrix(NA_real_, m, m)
  diag(cor_new) <- 1
  dimnames(vcov_new) <- dimnames(cor_new) <- list(rownames(W_new), rownames(W_new))

  dt <- sqrt(diag(KE_tt)); dn <- sqrt(diag(KE_nn))
  simil <- KE_nt / outer(dn, dt); simil[!is.finite(simil)] <- 0
  max_sim <- apply(simil, 1, max); mean_sim <- rowMeans(simil)

  pca <- stats::prcomp(Wtr, center = TRUE, scale. = FALSE)
  cumv <- cumsum(pca$sdev^2) / sum(pca$sdev^2)
  k_pc <- max(2L, min(which(cumv >= 0.95)))
  k_pc <- min(k_pc, ncol(pca$rotation))
  R_pc <- pca$rotation[, seq_len(k_pc), drop = FALSE]

  S_tr <- Wtr %*% R_pc
  S_nw <- Wnw %*% R_pc
  ctr_s <- colMeans(S_tr)
  sd_s  <- apply(S_tr, 2, stats::sd); sd_s[!is.finite(sd_s) | sd_s <= 0] <- 1
  Z_tr <- sweep(sweep(S_tr, 2, ctr_s, "-"), 2, sd_s, "/")
  Z_nw <- sweep(sweep(S_nw, 2, ctr_s, "-"), 2, sd_s, "/")

  mahal    <- sqrt(rowSums(Z_nw^2))
  mahal_tr <- sqrt(rowSums(Z_tr^2))
  ref      <- stats::quantile(mahal_tr, 0.95, names = FALSE)
  if (!is.finite(ref) || ref <= 0) ref <- max(sqrt(k_pc), 1)
  mahal_ratio <- mahal / ref

  ## ---- covariate envelope -------------------------------------------------
  ## Marginal range test. NOTE this is per-covariate: a site can lie inside the
  ## range of EVERY covariate individually yet occupy a COMBINATION that never
  ## occurred (hot AND wet, when training saw hot-dry and cool-wet). The
  ## marginal test reports 0% outside in that case; mahalanobis_ratio and
  ## weight_negativity are what catch it. Report several, collapse none.
  rng <- apply(Wtr, 2, range)
  sd_cov <- apply(Wtr, 2, stats::sd); sd_cov[!is.finite(sd_cov) | sd_cov <= 0] <- 1
  below <- sweep(Wnw, 2, rng[1, ], "<")
  above <- sweep(Wnw, 2, rng[2, ], ">")
  outside <- below | above
  n_out <- rowSums(outside)
  prop_out <- n_out / p

  ## how far outside, in training SD units (0 when inside)
  exc <- matrix(0, m, p, dimnames = list(rownames(W_new), colnames(Wtr)))
  lo_gap <- sweep(sweep(-Wnw, 2, -rng[1, ], "+"), 2, sd_cov, "/")
  hi_gap <- sweep(sweep(Wnw,  2,  rng[2, ], "-"), 2, sd_cov, "/")
  exc[below] <- lo_gap[below]
  exc[above] <- hi_gap[above]
  exc[!is.finite(exc)] <- 0
  max_exc  <- apply(exc, 1, max)
  mean_exc <- sapply(seq_len(m), function(i)
    if (n_out[i] > 0) mean(exc[i, outside[i, ]]) else 0)

  ## ---- kriging-weight geometry -------------------------------------------
  ## H = K_*,train (K_E + lambda I)^-1 is NOT constrained positive. Genuine
  ## interpolation gives mostly-positive weights summing near 1. Extrapolation
  ## produces large NEGATIVE weights: the predictor subtracts one training site
  ## from another to reach past them. That is an extrapolation signature in the
  ## space the model actually uses, rather than in raw covariate space, and it
  ## costs nothing because Wgt is already computed.
  w_sum <- rowSums(Wgt)
  w_neg <- rowSums(pmin(Wgt, 0))
  w_abs <- rowSums(abs(Wgt)); w_abs[w_abs <= 0] <- NA
  ## effective number of training sites leaned on (inverse participation ratio)
  w_eff <- (rowSums(abs(Wgt))^2) / rowSums(Wgt^2)

  ## ---- CI inflation vs the training sites --------------------------------
  ## Model-internal: what the posterior says about its own uncertainty. Its
  ## close relative ci_width_ratio is the strongest accuracy correlate in the
  ## 60-environment calibration (r = -0.66), well ahead of any envelope test.
  mean_ci <- colMeans(yhi - ylo)
  spread  <- apply(yhat, 2, stats::sd)
  ci_ratio <- mean_ci / ifelse(spread > 0, spread, NA)
  ci_ref <- stats::median(mean_ci, na.rm = TRUE)
  ci_infl <- if (is.finite(ci_ref) && ci_ref > 0) mean_ci / ci_ref else NA_real_

  dist_term <- pmin(1, pmax(0, (mahal_ratio - 0.5) / 2))
  divergence_index <- pmin(1, pmax(0, 0.5 * (1 - max_sim) + 0.5 * dist_term))

  ## ---- position relative to the training envelope -------------------------
  ## Geometric statement ONLY. This says where the site sits relative to the
  ## covariate cloud; it does NOT say the prediction is wrong. Calibration
  ## against 60 held-out environments gave prop_covariates_outside r = -0.38
  ## and mahalanobis_ratio r = -0.36 against accuracy, and the best threshold
  ## split mean accuracy 0.870 vs 0.818 across an observed 0.604-0.922 range.
  ## Gross extrapolation is the one regime that separates unambiguously
  ## (measured 7.2 for a 6x-perturbed site vs 0.25-1.47 for real held-out sites).
  position <- rep("interpolation", m)
  edge <- (n_out > 0) | (mahal_ratio > 1)
  mild <- (mahal_ratio > 2) | (prop_out > 0.25) | (max_exc > 1)
  grossv <- (mahal_ratio > 3) | (prop_out > 0.5)
  position[edge]   <- "edge"
  position[mild]   <- "mild_extrapolation"
  position[grossv] <- "gross_extrapolation"

  interpolation <- data.frame(
    env = rownames(W_new),
    position = position,
    n_outside = n_out, prop_outside = prop_out,
    max_exceedance_sd = max_exc, mean_exceedance_sd = mean_exc,
    mahalanobis_ratio = mahal_ratio,
    max_similarity = max_sim,
    weight_sum = w_sum, weight_negativity = w_neg,
    effective_n_train = w_eff,
    ci_width_ratio = ci_ratio, ci_inflation = ci_infl,
    nearest_training_env = envs[apply(simil, 1, which.max)],
    stringsAsFactors = FALSE, row.names = NULL)

  ## per-(environment x offending covariate) detail, so a user can see WHICH
  ## covariate put a site outside rather than only that something did
  ed <- which(outside, arr.ind = TRUE)
  exceedance_detail <- if (nrow(ed)) data.frame(
    env = rownames(W_new)[ed[, 1]],
    covariate = colnames(Wtr)[ed[, 2]],
    value = Wnw[ed],
    train_min = rng[1, ][ed[, 2]],
    train_max = rng[2, ][ed[, 2]],
    exceedance_sd = exc[ed],
    direction = ifelse(below[ed], "below", "above"),
    stringsAsFactors = FALSE, row.names = NULL
  ) else data.frame(env = character(0), covariate = character(0),
                    value = numeric(0), train_min = numeric(0),
                    train_max = numeric(0), exceedance_sd = numeric(0),
                    direction = character(0), stringsAsFactors = FALSE)
  if (nrow(exceedance_detail))
    exceedance_detail <- exceedance_detail[
      order(exceedance_detail$env, -exceedance_detail$exceedance_sd), ]

  divergence <- data.frame(
    env = rownames(W_new),
    divergence_index = divergence_index,
    max_similarity = max_sim, mean_similarity = mean_sim,
    mahalanobis = mahal, mahalanobis_ratio = mahal_ratio,
    n_covariates_outside_range = n_out, prop_covariates_outside = prop_out,
    nearest_training_env = envs[apply(simil, 1, which.max)],
    stringsAsFactors = FALSE, row.names = NULL)

  div_tr <- 0.5 * (1 - apply(
    (KE_tt / outer(dt, dt)) - diag(Inf, e), 1, max)) +
    0.5 * pmin(1, pmax(0, (mahal_tr / ref - 0.5) / 2))
  div_pct <- sapply(divergence_index, function(z) mean(div_tr <= z, na.rm = TRUE))

  diagnostics <- data.frame(
    env = rownames(W_new),
    position = position,
    mean_prediction = colMeans(yhat),
    ci_lower = colMeans(ylo), ci_upper = colMeans(yhi),
    mean_ci_width = mean_ci,
    ci_width_ratio = ci_ratio,
    ci_inflation = ci_infl,
    divergence_index = divergence_index,
    divergence_pctile_vs_training = div_pct,
    max_similarity = max_sim,
    mahalanobis_ratio = mahal_ratio,
    prop_covariates_outside = prop_out,
    weight_negativity = w_neg,
    stringsAsFactors = FALSE, row.names = NULL)

  ## ---- act on the envelope ------------------------------------------------
  ## Default "keep" preserves prior behaviour. "mask" is the recommended safety
  ## setting: it NAs the prediction but KEEPS the column, so downstream
  ## indexing by environment name still resolves and the NA propagates
  ## visibly. "drop" removes rows/columns and records them in $excluded.
  flagged <- if (identical(envelope_level, "mild")) (mild | grossv) else grossv
  excluded <- NULL
  if (any(flagged) && envelope %in% c("mask", "drop")) {
    excluded <- interpolation[flagged, , drop = FALSE]
    if (identical(envelope, "mask")) {
      yhat[, flagged] <- NA_real_
      ylo[,  flagged] <- NA_real_
      yhi[,  flagged] <- NA_real_
      yHat_long$yHat[yHat_long$env %in% rownames(W_new)[flagged]] <- NA_real_
      yHat_long$lower[yHat_long$env %in% rownames(W_new)[flagged]] <- NA_real_
      yHat_long$upper[yHat_long$env %in% rownames(W_new)[flagged]] <- NA_real_
      yHat_long$ci_width[yHat_long$env %in% rownames(W_new)[flagged]] <- NA_real_
      genetic_var[flagged, c("genetic_var", "lower", "upper", "genomic_h2")] <- NA_real_
      vcov_new[flagged, ] <- NA_real_; vcov_new[, flagged] <- NA_real_
      cor_new[flagged, ]  <- NA_real_; cor_new[, flagged]  <- NA_real_
    } else {
      kp <- !flagged
      yhat <- yhat[, kp, drop = FALSE]
      ylo  <- ylo[,  kp, drop = FALSE]
      yhi  <- yhi[,  kp, drop = FALSE]
      yHat_long <- yHat_long[yHat_long$env %in% rownames(W_new)[kp], , drop = FALSE]
      genetic_var <- genetic_var[kp, , drop = FALSE]
      vcov_new <- vcov_new[kp, kp, drop = FALSE]
      cor_new  <- cor_new[kp,  kp, drop = FALSE]
      interpolation <- interpolation[kp, , drop = FALSE]
      divergence    <- divergence[kp, , drop = FALSE]
      diagnostics   <- diagnostics[kp, , drop = FALSE]
      simil <- simil[kp, , drop = FALSE]
      Wgt   <- Wgt[kp, , drop = FALSE]
    }
  }

  ## ---- warnings -----------------------------------------------------------
  if (e < 20)
    warning("Only ", e, " training environments. This tool is designed for a ",
            "LARGE training network; with few environments the covariate-to-",
            "performance mapping is poorly estimated.")
  if (any(grossv)) {
    act <- switch(envelope, keep = "returned unchanged", flag = "returned, flagged",
                  mask = "set to NA", drop = "removed")
    warning(sum(grossv), " of ", m, " new environments are GROSS extrapolations ",
            "(mahalanobis_ratio > 3 or >50% of covariates outside the training ",
            "range): ", paste(rownames(W_new)[grossv], collapse = ", "),
            ". Predictions there are unsupported by the training data (", act, ").")
  }
  if (identical(envelope, "flag") && any(mild & !grossv))
    warning(sum(mild & !grossv), " new environments sit at mild extrapolation. ",
            "This is a geometric statement, not a reliability verdict: in ",
            "calibration, envelope measures correlated only r ~ -0.37 with ",
            "accuracy. Rank sites on these numbers; do not treat them as a gate.")
  if (!have_ge)
    warning("No GxE kernel in the fit: every new environment will rank the ",
            "genotypes identically.")

  out <- list(
    yHat_matrix = yhat, yHat = yHat_long,
    yHat_lower = ylo, yHat_upper = yhi,
    yHat_prob_win = pwin,
    genetic_var = genetic_var,
    vcov = vcov_new, cor = cor_new,
    interpolation = interpolation,
    exceedance_detail = exceedance_detail,
    excluded = excluded,
    divergence = divergence, diagnostics = diagnostics,
    meta = list(n_train_env = e, n_new_env = m, n_genotypes = g,
                n_covariates = p, n_draws = nd_use, level = level,
                residual_var = s2e, has_gxe = have_ge,
                envelope = envelope, envelope_level = envelope_level,
                n_flagged = sum(flagged),
                similarity = simil, kriging_weights = Wgt),
    call = match.call())
  class(out) <- "scan_untested_envs"
  out
}


#' @title Print a Scan of Untested Environments
#'
#' @description
#' Compact console summary of a \code{\link{scan_untested_envs}} result: the
#' per-site prediction and its extrapolation diagnostics, a tally of position
#' labels, the covariates that fell outside the training range, and a reminder
#' of how the geometric measures relate to calibrated accuracy.
#'
#' @param x an object of class \code{"scan_untested_envs"}.
#' @param ... ignored; present for S3 compatibility.
#'
#' @return \code{x}, invisibly.
#'
#' @seealso \code{\link{scan_untested_envs}}
#'
#' @export
print.scan_untested_envs <- function(x, ...) {
  mm <- x$meta
  cat("<scan_untested_envs>", mm$n_genotypes, "genotypes x", mm$n_new_env,
      "new environments\n")
  cat("  trained on", mm$n_train_env, "environments |", mm$n_covariates,
      "covariates |", mm$n_draws, "draws\n")
  cat("  GxE included:", mm$has_gxe, "| credible level:", mm$level,
      "| envelope:", mm$envelope, "\n\n")

  ip <- x$interpolation
  tb <- data.frame(env = ip$env, position = ip$position,
                   pred = round(x$diagnostics$mean_prediction, 3),
                   ci_ratio = round(ip$ci_width_ratio, 2),
                   mahal = round(ip$mahalanobis_ratio, 2),
                   out = paste0(round(100 * ip$prop_outside), "%"),
                   w_neg = round(ip$weight_negativity, 2),
                   max_sim = round(ip$max_similarity, 3),
                   stringsAsFactors = FALSE)
  print(tb, row.names = FALSE)

  ct <- table(factor(ip$position, levels = c("interpolation", "edge",
                                             "mild_extrapolation",
                                             "gross_extrapolation")))
  cat("\n  position:", paste(sprintf("%s %d", names(ct), as.integer(ct)),
                             collapse = " | "), "\n")

  if (!is.null(x$excluded) && nrow(x$excluded))
    cat("  ", nrow(x$excluded), " environment(s) ",
        if (identical(mm$envelope, "mask")) "masked (NA)" else "dropped",
        ": ", paste(x$excluded$env, collapse = ", "), "\n", sep = "")

  if (nrow(x$exceedance_detail)) {
    cat("\n  covariates outside the training range (top 5):\n")
    hd <- utils::head(x$exceedance_detail, 5)
    for (i in seq_len(nrow(hd)))
      cat(sprintf("   - %s / %s: %.2f is %.2f SD %s [%.2f, %.2f]\n",
                  hd$env[i], hd$covariate[i], hd$value[i],
                  hd$exceedance_sd[i], hd$direction[i],
                  hd$train_min[i], hd$train_max[i]))
  }

  cat("\n  Position is GEOMETRIC: where the site sits relative to the training\n")
  cat("  covariate cloud. It is not a reliability verdict. Calibrated on 60\n")
  cat("  held-out environments:\n")
  cat("    ci_width_ratio   vs accuracy       r = -0.66  (best accuracy guide)\n")
  cat("    divergence_index vs accuracy       r = -0.41\n")
  cat("    max_similarity   vs accuracy       r = +0.12  (NOT predictive)\n")
  cat("    max_similarity   vs gain-over-mean r = +0.69  (best gain guide)\n")
  cat("  Lower ci_width_ratio = more precise. Higher max_similarity = the\n")
  cat("  covariates are more likely to beat a plain genotype mean. To know your\n")
  cat("  real accuracy, run leave-environments-out CV on your own data.\n")
  invisible(x)
}


## ############################################################################
## PART 4 -- (B) PLANTING DATES       plot_planting_window()
## ############################################################################

## ============================================================================
## plot_planting_window.R
##
## (B) Performance of genotypes across candidate PLANTING DATES, with dates
##     outside the training covariate envelope visibly marked.
##
## The input is a scan_untested_envs() result in which each "environment" is a
## candidate sowing date at one location. This function does no prediction of
## its own -- it reshapes and plots what the scan already computed, so the
## credible intervals are the scan's genuine posterior intervals.
##
## ---------------------------------------------------------------------------
## DESIGN DECISION: the envelope layer is VISUALLY DOMINANT, not a footnote.
## ---------------------------------------------------------------------------
## Rendering extrapolated predictions in the same style as interpolated ones
## invites exactly the overclaim that was stripped out of the scan diagnostics.
## Gross extrapolations are therefore desaturated AND shaded AND counted in the
## subtitle. A user who ignores all three has been told three times.
##
## Position is a GEOMETRIC statement (where the site sits relative to the
## training covariate cloud), NOT a reliability verdict. In the 60-environment
## calibration, envelope measures correlated only r ~ -0.37 with accuracy.
## ============================================================================

#' @keywords internal
.ppw_need <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE))
    stop("Package '", pkg, "' is required for this function. ",
         "Install it or use `engine = \"base\"`.", call. = FALSE)
}

#' @keywords internal
.ppw_pos_levels <- c("interpolation", "edge",
                     "mild_extrapolation", "gross_extrapolation")

#' @keywords internal
.ppw_pos_cols <- c(interpolation        = "#1B7837",
                   edge                 = "#7FBC41",
                   mild_extrapolation   = "#FDB863",
                   gross_extrapolation  = "#B2182B")

## ---------------------------------------------------------------------------
## assemble the long table
## ---------------------------------------------------------------------------

#' Tabulate scan predictions against planting date
#'
#' Reshapes a `scan_untested_envs()` result into one row per
#' genotype x planting date, carrying the credible bounds and the envelope
#' position. Useful on its own when you want the numbers rather than a plot.
#'
#' @param scan a `scan_untested_envs` object whose new "environments" are
#'   candidate planting dates.
#' @param dates either (a) a data.frame with a column `env` matching the scan's
#'   new-environment names and a column giving the date, or (b) a vector of
#'   dates in the same order as the scan's new environments.
#' @param date_col name of the date column when `dates` is a data.frame.
#' @param genotypes optional character vector to subset genotypes.
#' @param group optional named vector or data.frame (`gid`, `group`) assigning
#'   genotypes to groups, for group-mean curves.
#'
#' @return data.frame with `gid`, `env`, `date`, `yHat`, `lower`, `upper`,
#'   `ci_width`, `position`, `mahalanobis_ratio`, `prop_outside`,
#'   `weight_negativity`, and `group` when supplied.
#' @export
planting_window_table <- function(scan, dates, date_col = "date",
                                  genotypes = NULL, group = NULL) {

  if (!inherits(scan, "scan_untested_envs"))
    stop("`scan` must be a scan_untested_envs object.", call. = FALSE)

  envs <- colnames(scan$yHat_matrix)
  if (is.null(envs))
    stop("scan$yHat_matrix has no column names.", call. = FALSE)

  ## ---- resolve dates ------------------------------------------------------
  if (is.data.frame(dates)) {
    if (!"env" %in% names(dates))
      stop("`dates` data.frame needs an `env` column matching the scan's ",
           "new environments.", call. = FALSE)
    if (!date_col %in% names(dates))
      stop("`dates` has no column '", date_col, "'.", call. = FALSE)
    miss <- setdiff(envs, as.character(dates$env))
    if (length(miss))
      stop("`dates` is missing these scanned environments: ",
           paste(miss, collapse = ", "), call. = FALSE)
    dd <- dates[match(envs, as.character(dates$env)), , drop = FALSE]
    dvec <- dd[[date_col]]
  } else {
    if (length(dates) != length(envs))
      stop("`dates` has length ", length(dates), " but the scan has ",
           length(envs), " new environments. Supply a data.frame with an ",
           "`env` column if the order is not guaranteed.", call. = FALSE)
    dvec <- dates
  }

  ## coerce to Date when it looks like one; otherwise leave numeric/ordered
  if (inherits(dvec, "Date")) {
    dparsed <- dvec
  } else if (is.numeric(dvec)) {
    dparsed <- dvec                       # day-of-year or index
  } else {
    dparsed <- suppressWarnings(as.Date(as.character(dvec)))
    if (all(is.na(dparsed)))
      stop("Could not interpret the date column as Date or numeric. ",
           "Supply Date objects or a numeric day-of-year.", call. = FALSE)
  }
  if (anyNA(dparsed))
    warning(sum(is.na(dparsed)), " date(s) could not be parsed and are NA.")

  ## ---- long predictions ---------------------------------------------------
  yh <- scan$yHat_matrix
  lo <- scan$yHat_lower
  hi <- scan$yHat_upper
  gids <- rownames(yh)
  if (!is.null(genotypes)) {
    keep <- intersect(genotypes, gids)
    if (!length(keep))
      stop("None of `genotypes` are present in the scan. Available: ",
           paste(utils::head(gids, 5), collapse = ", "),
           if (length(gids) > 5) ", ..." else "", call. = FALSE)
    if (length(keep) < length(genotypes))
      warning(length(genotypes) - length(keep),
              " requested genotype(s) not in the scan and were ignored.")
    yh <- yh[keep, , drop = FALSE]
    lo <- lo[keep, , drop = FALSE]
    hi <- hi[keep, , drop = FALSE]
    gids <- keep
  }

  g <- length(gids); m <- length(envs)
  out <- data.frame(
    gid   = rep(gids, times = m),
    env   = rep(envs, each = g),
    date  = rep(dparsed, each = g),
    yHat  = as.vector(yh),
    lower = as.vector(lo),
    upper = as.vector(hi),
    stringsAsFactors = FALSE)
  out$ci_width <- out$upper - out$lower

  ## ---- attach envelope position ------------------------------------------
  ip <- scan$interpolation
  if (is.null(ip))
    stop("scan has no $interpolation table. Re-run with the current ",
         "scan_untested_envs().", call. = FALSE)
  idx <- match(out$env, ip$env)
  out$position          <- ip$position[idx]
  out$mahalanobis_ratio <- ip$mahalanobis_ratio[idx]
  out$prop_outside      <- ip$prop_outside[idx]
  out$weight_negativity <- ip$weight_negativity[idx]
  out$position <- factor(out$position, levels = .ppw_pos_levels)

  ## ---- optional grouping --------------------------------------------------
  if (!is.null(group)) {
    if (is.data.frame(group)) {
      if (!all(c("gid", "group") %in% names(group)))
        stop("`group` data.frame needs columns `gid` and `group`.", call. = FALSE)
      out$group <- group$group[match(out$gid, group$gid)]
    } else {
      if (is.null(names(group)))
        stop("`group` vector must be named by genotype.", call. = FALSE)
      out$group <- unname(group[out$gid])
    }
    if (anyNA(out$group))
      warning(sum(is.na(out$group)), " genotype-date rows have no group ",
              "assignment (NA).")
  }

  out <- out[order(out$gid, out$date), ]
  rownames(out) <- NULL
  attr(out, "n_gross") <- sum(ip$position == "gross_extrapolation")
  attr(out, "n_env")   <- m
  out
}


## ---------------------------------------------------------------------------
## plot
## ---------------------------------------------------------------------------

#' Plot genotype performance across candidate planting dates
#'
#' Curves of predicted performance against sowing date, with credible ribbons
#' and out-of-envelope dates marked. Dates classified
#' `gross_extrapolation` are desaturated and shaded, and counted in the
#' subtitle, because position is reported but must not be read as a
#' reliability guarantee.
#'
#' @param scan a `scan_untested_envs` object, or a table from
#'   [planting_window_table()].
#' @param dates,date_col,genotypes,group passed to [planting_window_table()]
#'   when `scan` is a scan object.
#' @param summarise one of `"genotype"` (one curve per genotype), `"group"`
#'   (mean curve per group, requires `group`), or `"overall"` (single mean
#'   curve across all selected genotypes).
#' @param ribbon draw credible ribbons. Switched off automatically when more
#'   than `ribbon_max` curves are drawn, since overlapping ribbons are
#'   unreadable.
#' @param ribbon_max curve count above which ribbons are suppressed.
#' @param mark_envelope shade and desaturate dates outside the envelope.
#' @param envelope_level which positions to mark: `"gross"` (default) or
#'   `"mild"` (also marks `mild_extrapolation`).
#' @param engine `"ggplot2"` (default) or `"base"`.
#' @param ylab,title axis label and title.
#'
#' @return a ggplot object (invisibly for `engine = "base"`).
#' @export
plot_planting_window <- function(scan, dates = NULL, date_col = "date",
                                 genotypes = NULL, group = NULL,
                                 summarise = c("genotype", "group", "overall"),
                                 ribbon = TRUE, ribbon_max = 6,
                                 mark_envelope = TRUE,
                                 envelope_level = c("gross", "mild"),
                                 engine = c("ggplot2", "base"),
                                 ylab = "Predicted performance",
                                 title = NULL) {

  summarise <- match.arg(summarise)
  envelope_level <- match.arg(envelope_level)
  engine <- match.arg(engine)

  tab <- if (inherits(scan, "scan_untested_envs")) {
    if (is.null(dates))
      stop("`dates` is required when passing a scan object.", call. = FALSE)
    planting_window_table(scan, dates, date_col, genotypes, group)
  } else if (is.data.frame(scan)) {
    scan
  } else stop("`scan` must be a scan_untested_envs object or a ",
              "planting_window_table().", call. = FALSE)

  if (summarise == "group" && !"group" %in% names(tab))
    stop("summarise = \"group\" needs a `group` assignment.", call. = FALSE)

  ## ---- which dates are marked --------------------------------------------
  marked <- if (envelope_level == "mild")
    c("mild_extrapolation", "gross_extrapolation") else "gross_extrapolation"
  tab$flagged <- tab$position %in% marked

  ## ---- aggregate ----------------------------------------------------------
  key <- switch(summarise,
                genotype = "gid",
                group    = "group",
                overall  = NULL)

  if (is.null(key)) {
    agg <- stats::aggregate(
      cbind(yHat, lower, upper) ~ date, data = tab, FUN = mean, na.rm = TRUE)
    agg$series <- "all genotypes"
  } else {
    f <- stats::as.formula(paste("cbind(yHat, lower, upper) ~ date +", key))
    agg <- stats::aggregate(f, data = tab, FUN = mean, na.rm = TRUE)
    agg$series <- as.character(agg[[key]])
  }
  ## envelope position is a property of the DATE, not the series
  pos_by_date <- tab[!duplicated(tab$date), c("date", "position", "flagged")]
  agg <- merge(agg, pos_by_date, by = "date", all.x = TRUE)
  agg <- agg[order(agg$series, agg$date), ]

  n_series <- length(unique(agg$series))
  if (ribbon && n_series > ribbon_max) {
    ribbon <- FALSE
    message("Ribbons suppressed: ", n_series, " curves exceeds ribbon_max = ",
            ribbon_max, ". Overlapping ribbons are unreadable. ",
            "Raise ribbon_max or subset `genotypes` to restore them.")
  }

  n_flag <- sum(pos_by_date$flagged, na.rm = TRUE)
  n_date <- nrow(pos_by_date)
  sub <- sprintf("%d of %d planting dates outside the training envelope (%s)",
                 n_flag, n_date, envelope_level)
  if (is.null(title))
    title <- "Predicted performance across planting dates"

  ## ---- base engine --------------------------------------------------------
  if (engine == "base") {
    op <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(op))
    xs <- sort(unique(agg$date))
    graphics::plot(range(xs), range(c(agg$lower, agg$upper), na.rm = TRUE),
                   type = "n", xlab = "Planting date", ylab = ylab,
                   main = title)
    graphics::mtext(sub, side = 3, line = 0.2, cex = 0.8)
    if (mark_envelope) {
      fl <- pos_by_date$date[pos_by_date$flagged]
      for (d in fl)
        graphics::abline(v = d, col = grDevices::adjustcolor("#B2182B", 0.18),
                         lwd = 8)
    }
    cols <- grDevices::hcl.colors(n_series, "Dark 3")
    for (i in seq_along(unique(agg$series))) {
      s <- unique(agg$series)[i]
      a <- agg[agg$series == s, ]
      graphics::lines(a$date, a$yHat, col = cols[i], lwd = 2)
    }
    graphics::legend("topright", legend = unique(agg$series), col = cols,
                     lwd = 2, bty = "n", cex = 0.8)
    return(invisible(agg))
  }

  ## ---- ggplot2 engine -----------------------------------------------------
  .ppw_need("ggplot2")

  p <- ggplot2::ggplot(agg, ggplot2::aes(x = .data$date, y = .data$yHat))

  ## envelope shading UNDER the curves
  if (mark_envelope && n_flag > 0) {
    fl <- pos_by_date[pos_by_date$flagged, , drop = FALSE]
    p <- p + ggplot2::geom_vline(
      data = fl, ggplot2::aes(xintercept = .data$date),
      colour = "#B2182B", alpha = 0.16, linewidth = 3, inherit.aes = FALSE)
  }

  if (ribbon)
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$lower, ymax = .data$upper,
                   fill = .data$series), alpha = 0.18, colour = NA)

  p <- p +
    ggplot2::geom_line(ggplot2::aes(colour = .data$series), linewidth = 0.9) +
    ggplot2::labs(x = "Planting date", y = ylab, title = title,
                  subtitle = sub, colour = NULL, fill = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = if (n_series > 12) "none" else "right",
                   panel.grid.minor = ggplot2::element_blank())

  ## desaturate the flagged region by overlaying points at flagged dates
  if (mark_envelope && n_flag > 0)
    p <- p + ggplot2::geom_point(
      data = agg[agg$flagged %in% TRUE, , drop = FALSE],
      ggplot2::aes(colour = .data$series), shape = 4, size = 1.8,
      show.legend = FALSE)

  p
}


## ---------------------------------------------------------------------------
## best-date summary
## ---------------------------------------------------------------------------

#' Best planting date per genotype
#'
#' @param scan,dates,date_col,genotypes,group as in [planting_window_table()].
#' @param exclude_flagged drop dates outside the envelope before ranking.
#'   Default TRUE: an optimum found in a region the training data does not
#'   cover is not an optimum you can act on.
#' @param envelope_level `"gross"` (default) or `"mild"`.
#' @return data.frame: `gid`, `best_date`, `yHat`, `lower`, `upper`,
#'   `position`, `n_dates_considered`.
#' @examples
#' ## A scan is normally produced by scan_untested_envs(); here we build a
#' ## minimal one by hand -- 3 genotypes across 6 candidate sowing windows --
#' ## so the example is self-contained.
#' gids <- paste0("G", 1:3)
#' envs <- paste0("sow", 1:6)
#' mu <- matrix(c(6, 5.5, 5), 3, 6) +
#'       matrix(rep(c(0, .4, .7, .6, .2, -.3), each = 3), 3, 6)  # peaks mid-window
#' dimnames(mu) <- list(gids, envs)
#' scan <- list(
#'   yHat_matrix = mu, yHat_lower = mu - 0.5, yHat_upper = mu + 0.5,
#'   interpolation = data.frame(
#'     env               = envs,
#'     position          = c("interpolation", "interpolation", "edge",
#'                           "interpolation", "mild_extrapolation",
#'                           "gross_extrapolation"),
#'     mahalanobis_ratio = c(0.6, 0.8, 1.4, 0.9, 2.1, 3.6),
#'     prop_outside      = c(0, 0, 0.1, 0, 0.3, 0.7),
#'     weight_negativity = 0, stringsAsFactors = FALSE))
#' class(scan) <- "scan_untested_envs"
#'
#' ## Map each scanned environment to a calendar sowing date
#' dates <- data.frame(env = envs,
#'                     date = as.Date("2024-04-01") + c(0, 14, 28, 42, 56, 70))
#'
#' ## Best *supported* sowing date per genotype (gross extrapolations dropped)
#' best_planting_date(scan, dates)
#'
#' ## Rank every date, including the out-of-envelope ones
#' best_planting_date(scan, dates, exclude_flagged = FALSE)
#'
#' ## Stricter: also drop mild extrapolations before ranking
#' best_planting_date(scan, dates, envelope_level = "mild")
#' @export
best_planting_date <- function(scan, dates = NULL, date_col = "date",
                               genotypes = NULL, group = NULL,
                               exclude_flagged = TRUE,
                               envelope_level = c("gross", "mild")) {
  envelope_level <- match.arg(envelope_level)
  tab <- if (inherits(scan, "scan_untested_envs"))
    planting_window_table(scan, dates, date_col, genotypes, group) else scan

  marked <- if (envelope_level == "mild")
    c("mild_extrapolation", "gross_extrapolation") else "gross_extrapolation"

  n_all <- length(unique(tab$date))
  if (exclude_flagged) {
    tab <- tab[!tab$position %in% marked, , drop = FALSE]
    if (!nrow(tab))
      stop("Every planting date is outside the training envelope at ",
           "envelope_level = \"", envelope_level, "\". Nothing to rank. ",
           "Either sample environmental conditions in this window or set ",
           "exclude_flagged = FALSE and treat the result as unsupported.",
           call. = FALSE)
  }

  sp <- split(tab, tab$gid)
  res <- do.call(rbind, lapply(sp, function(z) {
    i <- which.max(z$yHat)
    data.frame(gid = z$gid[i], best_date = z$date[i], yHat = z$yHat[i],
               lower = z$lower[i], upper = z$upper[i],
               position = as.character(z$position[i]),
               n_dates_considered = nrow(z),
               stringsAsFactors = FALSE)
  }))
  rownames(res) <- NULL
  attr(res, "n_dates_total") <- n_all
  attr(res, "excluded_flagged") <- exclude_flagged
  res
}


## ############################################################################
## PART 5 -- (A,C) MAPS               map_scan() / grid_scan()
## ############################################################################

## ============================================================================
## map_scan.R
##
## (A) Spatial maps of scanned predictions, and
## (C) the same faceted by planting date.
##
## Two map types:
##   type = "performance" -- predicted value per genotype, one panel per
##                           genotype (A) or per planting date (C)
##   type = "coverage"    -- WHERE the training network supports prediction and
##                           where you must go and sample. Genotype-free.
##
## The coverage map is arguably the more valuable output: it is a SAMPLING
## DESIGN tool, not a prediction. It says which regions your training network
## fails to cover.
##
## ---------------------------------------------------------------------------
## MEMORY: grid_scan() exists because of this
## ---------------------------------------------------------------------------
## scan_untested_envs() allocates an ndraws x g x m array. For a 10,000-cell
## raster with 25 genotypes and 500 draws that is 10,000 x 25 x 500 x 8 bytes
## = 10 GB. grid_scan() chunks over grid cells and reduces draws, keeping peak
## memory at roughly  chunk_size x g x max_draws x 8.
##
## Defaults (chunk 500, max_draws 200, 25 genotypes) -> ~20 MB per chunk.
##
## ---------------------------------------------------------------------------
## The envelope layer is VISUALLY DOMINANT, by the same argument as in
## plot_planting_window(): rendering extrapolated cells in the same colour
## scale as interpolated ones invites the overclaim the scan diagnostics were
## built to avoid. Gross-extrapolation cells are hatched or greyed, and the
## count appears in the subtitle.
##
## Position is GEOMETRIC -- where a site sits relative to the training
## covariate cloud. It is NOT a reliability verdict. Calibration gave envelope
## measures r ~ -0.37 against accuracy.
## ============================================================================

#' @importFrom ggplot2 .data
NULL

#' @keywords internal
.ms_need <- function(pkg, what) {
  if (!requireNamespace(pkg, quietly = TRUE))
    stop("Package '", pkg, "' is required for ", what, ".", call. = FALSE)
}

#' @keywords internal
.ms_pos_levels <- c("interpolation", "edge",
                    "mild_extrapolation", "gross_extrapolation")

#' @keywords internal
.ms_pos_cols <- c(interpolation       = "#1B7837",
                  edge                = "#A6DBA0",
                  mild_extrapolation  = "#FDB863",
                  gross_extrapolation = "#B2182B")


## ---------------------------------------------------------------------------
## chunked scanning over many sites
## ---------------------------------------------------------------------------

#' Scan many sites in memory-bounded chunks
#'
#' [scan_untested_envs()] allocates an `ndraws x g x m` array, which is
#' infeasible for raster-sized `m`. This wrapper splits `W_new` into chunks,
#' scans each, and stitches the per-site results. Quantities that are defined
#' ACROSS new sites -- `vcov` and `cor` -- cannot be recovered from chunks and
#' are returned as `NULL` with a warning if requested.
#'
#' @param object,W_train,K_G as in [scan_untested_envs()].
#' @param W_new covariates for every grid cell / site (m x p).
#' @param chunk_size sites per chunk.
#' @param max_draws MCMC draws retained per chunk. Lower than the scan default
#'   because mapping needs the posterior mean far more than a tight interval.
#' @param verbose report chunk progress.
#' @param ... passed to [scan_untested_envs()].
#'
#' @return list with `yHat_matrix`, `yHat_lower`, `yHat_upper`,
#'   `yHat_prob_win`, `genetic_var`, `interpolation`, `exceedance_detail`,
#'   `meta`.
#'   Class `scan_untested_envs` so downstream functions accept it, but with
#'   `vcov` and `cor` set to `NULL`.
#' @export
grid_scan <- function(object, W_train, W_new, K_G = NULL,
                      chunk_size = 500, max_draws = 200,
                      verbose = TRUE, ...) {

  W_new <- as.matrix(W_new)
  m <- nrow(W_new)
  if (is.null(rownames(W_new)))
    rownames(W_new) <- paste0("cell", seq_len(m))

  idx <- split(seq_len(m), ceiling(seq_len(m) / chunk_size))
  nch <- length(idx)
  if (verbose)
    message("grid_scan: ", m, " sites in ", nch, " chunk(s) of <= ",
            chunk_size, ", ", max_draws, " draws.")

  parts <- vector("list", nch)
  for (k in seq_len(nch)) {
    if (verbose && (nch > 1))
      message("  chunk ", k, "/", nch, " (", length(idx[[k]]), " sites)")
    parts[[k]] <- suppressWarnings(scan_untested_envs(
      object, W_train = W_train, W_new = W_new[idx[[k]], , drop = FALSE],
      K_G = K_G, max_draws = max_draws, ...))
  }

  bind_col <- function(f) do.call(cbind, lapply(parts, f))
  bind_row <- function(f) do.call(rbind, lapply(parts, f))

  out <- list(
    yHat_matrix       = bind_col(function(z) z$yHat_matrix),
    yHat_lower        = bind_col(function(z) z$yHat_lower),
    yHat_upper        = bind_col(function(z) z$yHat_upper),
    yHat_prob_win     = bind_col(function(z) z$yHat_prob_win),
    yHat              = bind_row(function(z) z$yHat),
    genetic_var       = bind_row(function(z) z$genetic_var),
    interpolation     = bind_row(function(z) z$interpolation),
    exceedance_detail = bind_row(function(z) z$exceedance_detail),
    divergence        = bind_row(function(z) z$divergence),
    diagnostics       = bind_row(function(z) z$diagnostics),
    vcov = NULL, cor = NULL, excluded = NULL)

  mt <- parts[[1]]$meta
  mt$n_new_env <- m
  mt$n_draws   <- max_draws
  mt$chunked   <- TRUE
  mt$n_chunks  <- nch
  mt$similarity <- NULL       # meaningless once stitched
  mt$kriging_weights <- NULL
  out$meta <- mt
  out$call <- match.call()
  class(out) <- "scan_untested_envs"

  if (verbose)
    message("grid_scan: done. vcov/cor are NULL (not recoverable from chunks).")
  out
}


## ---------------------------------------------------------------------------
## assemble a spatial table
## ---------------------------------------------------------------------------

#' Join scan predictions to coordinates
#'
#' @param scan a `scan_untested_envs` (or [grid_scan()]) result.
#' @param coords data.frame with `env` matching the scan's new environments,
#'   plus longitude and latitude columns.
#' @param lon_col,lat_col coordinate column names.
#' @param date_col optional column giving a planting date per site, for (C).
#' @param genotypes optional subset.
#' @return long data.frame: `gid`, `env`, `lon`, `lat`, [`date`], `yHat`,
#'   `lower`, `upper`, `position`, `mahalanobis_ratio`, `prop_outside`,
#'   `weight_negativity`, `prob_win` (posterior probability the genotype is the
#'   best at that site; `NA` for scans that do not carry it).
#' @export
scan_spatial_table <- function(scan, coords, lon_col = "lon", lat_col = "lat",
                               date_col = NULL, genotypes = NULL) {

  if (!inherits(scan, "scan_untested_envs"))
    stop("`scan` must be a scan_untested_envs object.", call. = FALSE)
  for (cc in c("env", lon_col, lat_col))
    if (!cc %in% names(coords))
      stop("`coords` has no column '", cc, "'.", call. = FALSE)

  envs <- colnames(scan$yHat_matrix)
  miss <- setdiff(envs, as.character(coords$env))
  if (length(miss))
    stop("`coords` is missing ", length(miss), " scanned environment(s): ",
         paste(utils::head(miss, 5), collapse = ", "),
         if (length(miss) > 5) ", ..." else "", call. = FALSE)
  cd <- coords[match(envs, as.character(coords$env)), , drop = FALSE]

  yh <- scan$yHat_matrix; lo <- scan$yHat_lower; hi <- scan$yHat_upper
  pw <- scan$yHat_prob_win
  gids <- rownames(yh)
  if (!is.null(genotypes)) {
    keep <- intersect(genotypes, gids)
    if (!length(keep))
      stop("None of `genotypes` are in the scan.", call. = FALSE)
    yh <- yh[keep, , drop = FALSE]; lo <- lo[keep, , drop = FALSE]
    hi <- hi[keep, , drop = FALSE]; gids <- keep
    if (!is.null(pw)) pw <- pw[keep, , drop = FALSE]
  }
  g <- length(gids); m <- length(envs)

  out <- data.frame(
    gid = rep(gids, times = m), env = rep(envs, each = g),
    lon = rep(cd[[lon_col]], each = g), lat = rep(cd[[lat_col]], each = g),
    yHat = as.vector(yh), lower = as.vector(lo), upper = as.vector(hi),
    stringsAsFactors = FALSE)
  if (!is.null(date_col)) {
    if (!date_col %in% names(coords))
      stop("`coords` has no column '", date_col, "'.", call. = FALSE)
    out$date <- rep(cd[[date_col]], each = g)
  }

  ip <- scan$interpolation
  i <- match(out$env, ip$env)
  out$position          <- factor(ip$position[i], levels = .ms_pos_levels)
  out$mahalanobis_ratio <- ip$mahalanobis_ratio[i]
  out$prop_outside      <- ip$prop_outside[i]
  out$weight_negativity <- ip$weight_negativity[i]
  out$ci_width          <- out$upper - out$lower
  out$prob_win          <- if (!is.null(pw)) as.vector(pw) else NA_real_
  rownames(out) <- NULL
  out
}


## ---------------------------------------------------------------------------
## the map
## ---------------------------------------------------------------------------

#' Map scanned predictions, or training-network coverage
#'
#' @details
#' Three kinds of map are available through `type`:
#' \itemize{
#'   \item `"performance"`: the predicted value (`yHat`) of one or more
#'     genotypes across the scanned sites, one panel per genotype (or a single
#'     panel of the mean with `summarise = "mean"`). Sites outside the training
#'     covariate envelope are drawn as grey crosses rather than on the value
#'     scale.
#'   \item `"coverage"`: where the training network actually supports
#'     prediction, colouring each site by its geometric position relative to
#'     the training covariate cloud (interpolation, edge, mild- or
#'     gross-extrapolation).
#'   \item `"which_won_where"`: a classic \emph{which-won-where} view for
#'     selection. It answers two questions with two maps (and two tables).
#'     \emph{Who won where} is the single top genotype per site, by posterior
#'     mean (`yHat`; ties broken by the first genotype). \emph{Winners} is how
#'     many genotypes are statistically tied for the win, taken from the scan's
#'     posterior draws: a genotype is a co-winner when its probability of being
#'     the best at that site (`prob_win`) is at least `win_threshold`. The count
#'     (`n_winners`) and the tied genotypes (`winners`, most probable first) are
#'     reported per site, so you can see where the decision is clear (one
#'     winner) versus where several genotypes are indistinguishable. This needs
#'     a scan from the current [scan_untested_envs()], which stores the per-site
#'     win probabilities.
#' }
#'
#' @param scan a `scan_untested_envs` / [grid_scan()] result.
#' @param coords,lon_col,lat_col,date_col as in [scan_spatial_table()].
#' @param type `"performance"` (predicted value), `"coverage"` (where the
#'   training network supports prediction), or `"which_won_where"` (the top
#'   genotype per site, and how many genotypes tie it).
#' @param genotypes genotypes to map. Required for `type = "performance"`
#'   unless `summarise = "mean"`. For `type = "which_won_where"` it restricts
#'   the pool of genotypes that compete for the win (default: all).
#' @param summarise for `type = "performance"`: `"none"` (panel per genotype)
#'   or `"mean"` (single panel, mean over selected genotypes).
#' @param facet_by `"genotype"` (A) or `"date"` (C). `"date"` needs `date_col`.
#' @param shp optional `sf` polygon layer (e.g. a state or country border from
#'   `geobr::read_state()`) drawn underneath the points; also used to clip when
#'   `clip = TRUE`.
#' @param clip if `TRUE` and `shp` is supplied, keep only the scanned sites that
#'   fall inside `shp` (points are matched to `shp` in its CRS).
#' @param mark_envelope hatch/grey cells outside the envelope.
#' @param envelope_level `"gross"` (default) or `"mild"`.
#' @param win_threshold for `type = "which_won_where"`: the minimum posterior
#'   probability of being the best genotype for another genotype to count as a
#'   co-winner at a site (default `0.25`).
#' @param point_size,palette cosmetic.
#' @param output `"plot"` (default) returns a ggplot; `"table"` returns the
#'   data.frame that would have been plotted (one row per site, or per
#'   site-by-genotype for `type = "performance"`), so you can build your own
#'   map. Clipping and `summarise` are applied first, and a logical `flagged`
#'   column marks sites outside the envelope. `"table"` needs no ggplot2.
#' @return a ggplot object, or a data.frame when `output = "table"`. For
#'   `type = "which_won_where"` a two-element list is returned instead:
#'   `who_won_where` (the single top genotype per site, by posterior mean) and
#'   `winners` (how many genotypes have posterior `prob_win >= win_threshold`),
#'   each a ggplot, or a data.frame when `output = "table"`.
#' @examples
#' ## A scan is normally produced by scan_untested_envs(); here we build a
#' ## minimal one by hand -- 2 genotypes predicted at 9 sites on a 3x3 grid.
#' sites <- paste0("S", 1:9)
#' grid  <- expand.grid(lon = 1:3, lat = 1:3)
#' yh <- rbind(G1 = 6 + grid$lon * 0.3 - grid$lat * 0.1,
#'             G2 = 5 + grid$lat * 0.4)
#' colnames(yh) <- sites
#' scan <- list(
#'   yHat_matrix = yh, yHat_lower = yh - 0.4, yHat_upper = yh + 0.4,
#'   ## posterior P(best) per genotype x site (illustrative softmax of yHat)
#'   yHat_prob_win = prop.table(exp(yh), margin = 2),
#'   interpolation = data.frame(
#'     env               = sites,
#'     position          = rep(c("interpolation", "edge",
#'                               "gross_extrapolation"), each = 3),
#'     mahalanobis_ratio = seq(0.5, 3.5, length.out = 9),
#'     prop_outside      = seq(0, 0.8, length.out = 9),
#'     weight_negativity = 0, stringsAsFactors = FALSE))
#' class(scan) <- "scan_untested_envs"
#'
#' ## Site coordinates, keyed by the scanned environment name
#' coords <- data.frame(env = sites, lon = grid$lon, lat = grid$lat)
#'
#' \donttest{
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   ## Where does the training network support prediction?
#'   map_scan(scan, coords, type = "coverage")
#'
#'   ## Predicted performance for a single genotype
#'   map_scan(scan, coords, type = "performance", genotypes = "G1")
#'
#'   ## Mean over several genotypes in one panel
#'   map_scan(scan, coords, type = "performance",
#'            genotypes = c("G1", "G2"), summarise = "mean")
#'
#'   ## Which genotype wins at each site, and how many genotypes tie it
#'   ww <- map_scan(scan, coords, type = "which_won_where")
#'   ww$who_won_where
#'   ww$winners
#' }
#' }
#'
#' ## Export the plotted data instead of a figure, to map it your own way
#' head(map_scan(scan, coords, type = "coverage", output = "table"))
#'
#' ## Both winner tables (no ggplot2 needed), one row per site
#' ww_tab <- map_scan(scan, coords, type = "which_won_where", output = "table")
#' head(ww_tab$who_won_where)
#' head(ww_tab$winners)
#'
#' \dontrun{
#' ## ---- Realistic scan: 50 new locations across Mato Grosso, Brazil --------
#' ## The bundled maize datasets are multi-environment trials from Mato Grosso:
#' ## maizeWTH (daily weather), maizeYield (grain yield) and maizeG (genomic
#' ## relationships). We train on them, then predict the trained genotypes at
#' ## 50 candidate sites spread across the state. Needs BGGE (fitting) and a
#' ## network connection (get_weather() queries NASA POWER).
#' data("maizeWTH"); data("maizeYield"); data("maizeG")
#'
#' ## 1. Training environmental covariates (first 100 days after sowing)
#' vars <- c("FRUE", "PETP", "GDD", "T2M_MAX")
#' Wtr  <- W_matrix(env.data = maizeWTH[maizeWTH$daysFromStart < 100, ],
#'                  var.id = vars, statistic = "mean")
#'
#' ## 2. Kernels + a fit that KEEPS the effect chains (required by any scan)
#' K   <- get_kernel(K_G = list(G = maizeG),
#'                   K_E = list(W = env_kernel(env.data = Wtr)[[2]]),
#'                   data = maizeYield, model = "RNMM",
#'                   env = "env", gid = "gid", y = "value")
#' fit <- kernel_model(y = "value", data = maizeYield, random = K,
#'                     env = "env", gid = "gid",
#'                     iterations = 5000, burnin = 1000, keep_effects = TRUE)
#'
#' ## 3. Fifty candidate locations across Mato Grosso (WGS84 decimal degrees).
#' ##    The state spans roughly 7.3-18.0 S and 50.2-61.6 W.
#' set.seed(2024)
#' coords <- data.frame(
#'   env = sprintf("MT%02d", 1:50),
#'   lon = runif(50, -61.6, -50.2),
#'   lat = runif(50, -18.0,  -7.3))
#'
#' ## 4. Daily weather for each site over one summer season, processed into the
#' ##    SAME covariates as the training matrix.
#' wth_new <- get_weather(env.id = coords$env, lat = coords$lat, lon = coords$lon,
#'                        start.day = "2023-11-01", end.day = "2024-02-15",
#'                        variables.names = c("T2M", "T2M_MAX", "T2M_MIN",
#'                                            "T2MDEW", "PRECTOT", "RH2M",
#'                                            "ALLSKY_SFC_SW_DWN"))
#' wth_new <- processWTH(wth_new)          # derives FRUE, PETP, GDD, ...
#' Wnew <- W_matrix(env.data = wth_new[wth_new$daysFromStart < 100, ],
#'                  var.id = vars, statistic = "mean")
#'
#' ## 5. Predict the trained genotypes at all 50 sites
#' sc <- scan_untested_envs(fit, W_train = Wtr, W_new = Wnew, K_G = maizeG)
#'
#' ## 6. Coverage first: where does the trial network actually support prediction?
#' map_scan(sc, coords, type = "coverage")
#'
#' ## Predicted performance of one genotype across Mato Grosso
#' g1 <- rownames(maizeG)[1]
#' map_scan(sc, coords, type = "performance", genotypes = g1)
#'
#' ## Mean predicted performance over the first ten genotypes
#' map_scan(sc, coords, type = "performance",
#'          genotypes = rownames(maizeG)[1:10], summarise = "mean")
#'
#' ## 7. Overlay the Mato Grosso state border for context. geobr serves the
#' ##    official IBGE boundaries for Brazil as sf objects; clip = TRUE also
#' ##    drops any scanned site that falls outside the polygon.
#' mt <- geobr::read_state(code_state = "MT", year = 2020)
#' map_scan(sc, coords, type = "coverage", shp = mt)
#' map_scan(sc, coords, type = "performance", genotypes = g1,
#'          shp = mt, clip = TRUE)
#' }
#' @export
map_scan <- function(scan, coords, lon_col = "lon", lat_col = "lat",
                     date_col = NULL,
                     type = c("performance", "coverage", "which_won_where"),
                     genotypes = NULL, summarise = c("none", "mean"),
                     facet_by = c("genotype", "date"),
                     shp = NULL, clip = FALSE,
                     mark_envelope = TRUE,
                     envelope_level = c("gross", "mild"),
                     win_threshold = 0.25,
                     point_size = 2, palette = "viridis",
                     output = c("plot", "table")) {

  type <- match.arg(type); summarise <- match.arg(summarise)
  facet_by <- match.arg(facet_by); envelope_level <- match.arg(envelope_level)
  output <- match.arg(output)
  if (output == "plot") .ms_need("ggplot2", "mapping")

  if (facet_by == "date" && is.null(date_col))
    stop("facet_by = \"date\" requires `date_col`.", call. = FALSE)

  tab <- scan_spatial_table(scan, coords, lon_col, lat_col, date_col,
                            genotypes = if (type == "coverage") NULL else genotypes)

  ## optionally clip the scanned points to the shapefile's polygons
  if (!is.null(shp) && isTRUE(clip)) {
    .ms_need("sf", "clipping points to a shapefile")
    pts <- sf::st_as_sf(tab, coords = c("lon", "lat"),
                        crs = sf::st_crs(shp), remove = FALSE)
    inside <- lengths(sf::st_intersects(pts, sf::st_union(shp))) > 0
    if (!any(inside))
      stop("clip = TRUE removed every site: none fall inside `shp`. ",
           "Check that `shp` and the coordinates share a CRS.", call. = FALSE)
    tab <- tab[inside, , drop = FALSE]
  }

  marked <- if (envelope_level == "mild")
    c("mild_extrapolation", "gross_extrapolation") else "gross_extrapolation"

  ## ---- which-won-where ----------------------------------------------------
  if (type == "which_won_where") {
    if (!is.numeric(win_threshold) || win_threshold < 0 || win_threshold > 1)
      stop("`win_threshold` must be a probability in [0, 1].", call. = FALSE)
    if (all(is.na(tab$prob_win)))
      stop("This scan carries no posterior win probabilities (`yHat_prob_win`); ",
           "re-run scan_untested_envs() to use type = \"which_won_where\".",
           call. = FALSE)
    idx_by_env <- split(seq_len(nrow(tab)), tab$env)
    res <- do.call(rbind, lapply(names(idx_by_env), function(en) {
      sub <- tab[idx_by_env[[en]], , drop = FALSE]
      o   <- which.max(sub$yHat)                 # top genotype by posterior mean
      co  <- sub$prob_win >= win_threshold       # statistical co-winners
      co[o] <- TRUE                              # the top genotype always counts
      ord <- order(sub$prob_win[co], decreasing = TRUE)
      data.frame(
        env = en, lon = sub$lon[1L], lat = sub$lat[1L],
        date = if (!is.null(date_col)) sub$date[1L] else NA,
        winner = sub$gid[o], winner_yHat = sub$yHat[o],
        winner_prob = sub$prob_win[o],
        n_winners = sum(co),
        winners = paste(sub$gid[co][ord], collapse = ", "),
        max_prob = max(sub$prob_win),
        position = as.character(sub$position[o]),
        stringsAsFactors = FALSE)
    }))
    res$position <- factor(res$position, levels = .ms_pos_levels)
    res$flagged  <- res$position %in% marked
    if (is.null(date_col)) res$date <- NULL

    keep_who <- c("env", "lon", "lat", if (!is.null(date_col)) "date",
                  "winner", "winner_yHat", "winner_prob",
                  "position", "flagged")
    keep_win <- c("env", "lon", "lat", if (!is.null(date_col)) "date",
                  "n_winners", "winners", "max_prob", "position", "flagged")
    tab_who <- res[, keep_who, drop = FALSE]; rownames(tab_who) <- NULL
    tab_win <- res[, keep_win, drop = FALSE]; rownames(tab_win) <- NULL
    if (output == "table")
      return(list(who_won_where = tab_who, winners = tab_win))

    facet_date <- facet_by == "date" && !is.null(date_col)

    ## map 1: the single top genotype per site, comparing yHat only
    p_who <- ggplot2::ggplot()
    if (!is.null(shp)) {
      .ms_need("sf", "shapefile overlays")
      p_who <- p_who + ggplot2::geom_sf(data = shp, fill = "grey96",
                                        colour = "grey70", linewidth = 0.3)
    }
    p_who <- p_who +
      ggplot2::geom_point(
        data = res, ggplot2::aes(.data$lon, .data$lat, colour = .data$winner),
        size = point_size) +
      ggplot2::labs(
        title = "Who won where",
        subtitle = sprintf("%d sites, %d distinct winning genotype(s) by yHat",
                           nrow(res), length(unique(res$winner))),
        colour = "winning genotype", x = NULL, y = NULL,
        caption = "Top genotype per site, by posterior mean (yHat).") +
      ggplot2::theme_minimal(base_size = 11)
    if (facet_date) p_who <- p_who + ggplot2::facet_wrap(~ date)

    ## map 2: how many genotypes tie the winner (overlapping CIs)
    p_win <- ggplot2::ggplot()
    if (!is.null(shp)) {
      .ms_need("sf", "shapefile overlays")
      p_win <- p_win + ggplot2::geom_sf(data = shp, fill = "grey96",
                                        colour = "grey70", linewidth = 0.3)
    }
    p_win <- p_win +
      ggplot2::geom_point(
        data = res,
        ggplot2::aes(.data$lon, .data$lat, colour = .data$n_winners),
        size = point_size)
    if (palette == "viridis")
      p_win <- p_win + ggplot2::scale_colour_viridis_c(name = "co-winners")
    else
      p_win <- p_win + ggplot2::scale_colour_gradientn(
        colours = grDevices::hcl.colors(11, palette), name = "co-winners")
    p_win <- p_win +
      ggplot2::labs(
        title = "How many genotypes tie for the win",
        subtitle = sprintf(
          "co-winners (P(win) >= %.2f) per site: median %g, range %d-%d",
          win_threshold, stats::median(res$n_winners),
          min(res$n_winners), max(res$n_winners)),
        x = NULL, y = NULL,
        caption = sprintf(paste("A genotype is a co-winner when its posterior",
                                "probability of being best at the site is >=",
                                "%.2f."), win_threshold)) +
      ggplot2::theme_minimal(base_size = 11)
    if (facet_date) p_win <- p_win + ggplot2::facet_wrap(~ date)

    return(list(who_won_where = p_who, winners = p_win))
  }

  ## ---- coverage map -------------------------------------------------------
  if (type == "coverage") {
    cov <- tab[!duplicated(tab$env), , drop = FALSE]
    cov$flagged <- cov$position %in% marked
    if (output == "table") { rownames(cov) <- NULL; return(cov) }
    n_bad <- sum(cov$flagged)
    sub <- sprintf(
      "%d of %d sites outside the training envelope (%s) -- sample here",
      n_bad, nrow(cov), envelope_level)

    p <- ggplot2::ggplot()
    if (!is.null(shp)) {
      .ms_need("sf", "shapefile overlays")
      p <- p + ggplot2::geom_sf(data = shp, fill = "grey96",
                                colour = "grey70", linewidth = 0.3)
    }
    p <- p +
      ggplot2::geom_point(
        data = cov,
        ggplot2::aes(.data$lon, .data$lat, colour = .data$position),
        size = point_size) +
      ggplot2::scale_colour_manual(values = .ms_pos_cols, drop = FALSE,
                                   name = "envelope position") +
      ggplot2::labs(
        title = "Training-network coverage",
        subtitle = sub, x = NULL, y = NULL,
        caption = paste("Geometric position relative to the training",
                        "covariate cloud -- not a reliability verdict.")) +
      ggplot2::theme_minimal(base_size = 11)
    if (facet_by == "date" && !is.null(date_col))
      p <- p + ggplot2::facet_wrap(~ date)
    return(p)
  }

  ## ---- performance map ----------------------------------------------------
  if (is.null(genotypes) && summarise == "none") {
    ng <- length(unique(tab$gid))
    if (ng > 12)
      stop("Mapping ", ng, " genotypes would produce ", ng, " panels. ",
           "Supply `genotypes` or set summarise = \"mean\".", call. = FALSE)
  }
  if (summarise == "mean") {
    by <- c("env", "lon", "lat", "position",
            if (!is.null(date_col)) "date")
    tab <- stats::aggregate(
      tab[, c("yHat", "lower", "upper")],
      by = as.list(tab[, by, drop = FALSE]), FUN = mean, na.rm = TRUE)
    tab$gid <- "mean of selected genotypes"
  }

  tab$flagged <- tab$position %in% marked
  if (output == "table") { rownames(tab) <- NULL; return(tab) }
  n_bad <- sum(tab$flagged[!duplicated(tab$env)])
  n_site <- length(unique(tab$env))
  sub <- sprintf("%d of %d sites outside the training envelope (%s)%s",
                 n_bad, n_site, envelope_level,
                 if (mark_envelope && n_bad) " -- marked" else "")

  p <- ggplot2::ggplot()
  if (!is.null(shp)) {
    .ms_need("sf", "shapefile overlays")
    p <- p + ggplot2::geom_sf(data = shp, fill = "grey96",
                              colour = "grey70", linewidth = 0.3)
  }

  ## in-envelope cells carry the colour scale
  good <- tab[!tab$flagged, , drop = FALSE]
  bad  <- tab[ tab$flagged, , drop = FALSE]

  p <- p + ggplot2::geom_point(
    data = good, ggplot2::aes(.data$lon, .data$lat, colour = .data$yHat),
    size = point_size)

  if (palette == "viridis")
    p <- p + ggplot2::scale_colour_viridis_c(name = "predicted")
  else
    p <- p + ggplot2::scale_colour_gradientn(
      colours = grDevices::hcl.colors(11, palette), name = "predicted")

  ## flagged cells: grey X, never on the value scale
  if (mark_envelope && nrow(bad))
    p <- p + ggplot2::geom_point(
      data = bad, ggplot2::aes(.data$lon, .data$lat),
      colour = "grey55", shape = 4, size = point_size * 0.9, stroke = 0.7)
  else if (!mark_envelope && nrow(bad))
    p <- p + ggplot2::geom_point(
      data = bad, ggplot2::aes(.data$lon, .data$lat, colour = .data$yHat),
      size = point_size)

  ## NB: vars() EVALUATES its arguments, unlike aes(), so ggplot2::vars(.data$gid)
  ## errors with "object '.data' not found". Using !! would work but pulls in
  ## rlang as a hard runtime dependency. facet_wrap() accepts a plain formula,
  ## which needs no tidy-eval at all -- use that.
  fw <- switch(facet_by,
               genotype = stats::as.formula("~ gid"),
               date     = stats::as.formula("~ date"))
  p <- p + ggplot2::facet_wrap(fw) +
    ggplot2::labs(
      title = if (facet_by == "date")
        "Predicted performance by planting date" else
        "Predicted performance by genotype",
      subtitle = sub, x = NULL, y = NULL,
      caption = paste("Grey x = outside the training covariate envelope.",
                      "Geometric position, not a reliability verdict.")) +
    ggplot2::theme_minimal(base_size = 11)
  p
}


#' Summarise coverage of a scanned region
#'
#' Tabulates how much of a scanned area falls in each envelope class, and
#' which covariates most often drive sites outside. The second table answers
#' "what should I go and measure".
#'
#' @param scan a `scan_untested_envs` / [grid_scan()] result.
#' @return list of `by_position` and `by_covariate` data.frames.
#' @export
coverage_summary <- function(scan) {
  ip <- scan$interpolation
  if (is.null(ip)) stop("scan has no $interpolation table.", call. = FALSE)
  pos <- factor(ip$position, levels = .ms_pos_levels)
  by_position <- data.frame(
    position = levels(pos),
    n = as.integer(table(pos)),
    prop = round(as.numeric(table(pos)) / nrow(ip), 4),
    stringsAsFactors = FALSE)

  ed <- scan$exceedance_detail
  by_covariate <- if (!is.null(ed) && nrow(ed)) {
    ag <- stats::aggregate(exceedance_sd ~ covariate, ed,
                           function(z) c(n = length(z), max = max(z),
                                         mean = mean(z)))
    res <- data.frame(covariate = ag$covariate,
                      n_sites_outside = ag$exceedance_sd[, "n"],
                      max_exceedance_sd = round(ag$exceedance_sd[, "max"], 3),
                      mean_exceedance_sd = round(ag$exceedance_sd[, "mean"], 3),
                      stringsAsFactors = FALSE)
    res[order(-res$n_sites_outside), ]
  } else data.frame(covariate = character(0), n_sites_outside = integer(0),
                    max_exceedance_sd = numeric(0),
                    mean_exceedance_sd = numeric(0))
  rownames(by_covariate) <- NULL
  list(by_position = by_position, by_covariate = by_covariate)
}


# =========================================================================
# PART 6 - SIMULATION OF MULTI-ENVIRONMENT TRIALS WITH A KNOWN ANSWER
# =========================================================================
#
# Two core simulators, plus accessories:
#
#   sim_met()  CORE      simulate a MET: genotype x environment phenotypes with
#                        a known genetic architecture (kinship K among lines,
#                        genetic correlation C among environments, per-env h2).
#   sim_W()    CORE      simulate an envirome (q x n_var covariable matrix) whose
#                        environmental kernel recovers a KNOWN share of a target
#                        environment correlation C_env.
#
#   sim_met_C()          accessory: build a valid q x q environment correlation.
#   env_cor()            accessory: extract the environment correlation from a
#                        sim_met object or a bare matrix, so hand-off never
#                        depends on internal structure.
#   sim_W_grid()         accessory: sweep sim_W() over size/quality for power
#                        analysis.
#
# THE CENTRAL DISTINCTION the whole design turns on:
#   sim_met(K = ...)    K     is n x n, a GENOMIC KINSHIP AMONG LINES.
#   sim_W(C_env = ...)  C_env is q x q, a CORRELATION AMONG ENVIRONMENTS.
# These are different objects with different dimensions. C_env is NOT a genomic
# relationship matrix.
#
# Reference for the EPA premise behind sim_W(): Costa-Neto et al. (2023) G3
# 13(2) jkac313.
# -------------------------------------------------------------------------

#' Symmetric PSD square root with eigenvalue clamping.
#' Empirical type-B correlation matrices from unbalanced METs are routinely
#' indefinite; that is the normal case, not an error.
#' @noRd
.sim_psd_sqrt <- function(C, tol = 1e-10, warn = TRUE) {
  ev <- eigen(C, symmetric = TRUE)
  neg <- ev$values < tol
  clamped <- sum(pmin(ev$values, 0))
  ev$values[neg] <- 0
  if (warn && any(neg) && abs(clamped) > 1e-8)
    warning(sprintf("C_env was not positive definite: %d eigenvalue(s) clamped, total mass %.3g.",
                    sum(neg), abs(clamped)), call. = FALSE)
  ev$vectors %*% diag(sqrt(ev$values), nrow(C)) %*% t(ev$vectors)
}

#' Make a kinship positive definite. Real marker-based K is routinely singular
#' (maizeG included) and a bare chol() fails on exactly the matrices users have.
#' @noRd
.sim_pd <- function(K, ridge = 1e-6) {
  K <- as.matrix(K)
  ev <- min(eigen(K, symmetric = TRUE, only.values = TRUE)$values)
  if (ev < ridge) K <- K + diag(ridge - ev + 1e-8, nrow(K))
  K
}

#' Orthonormalised Gaussian draw: t(U) %*% U = n I.
#' This is what makes the realised moments EXACT rather than merely centred on
#' the target. With a plain Gaussian draw, noise = 0 gives r ~ 0.97, not 1.
#' @noRd
.sim_orth <- function(n, k) {
  Z <- matrix(stats::rnorm(n * k), n, k)
  s <- svd(Z)
  (s$u %*% t(s$v)) * sqrt(n)
}

#' Environmental kernel, matching env_kernel()'s convention: WW'/(trace(WW')/q).
#' Do not re-derive this elsewhere or the numbers will not line up with what
#' users get downstream.
#' @noRd
.sim_GB <- function(W) {
  K <- tcrossprod(W)
  K / (sum(diag(K)) / nrow(K))
}

#' Off-diagonal Pearson correlation between two symmetric matrices.
#' Off-diagonals ONLY: including the diagonal inflates this badly, since both
#' matrices have large near-constant diagonals.
#' @noRd
.sim_offdiag_cor <- function(A, B) {
  lt <- lower.tri(A)
  stats::cor(A[lt], B[lt])
}

#' Validate a q x q correlation among ENVIRONMENTS.
#'
#' Tests the diagonal ELEMENTWISE, not on its mean. A genomic kinship is the
#' most likely wrong input and it is square, symmetric, PSD, with a diagonal
#' averaging ~1 -- every cheap check passes. Verified on a simulated 150-line
#' kinship: mean(diag(K)) = 0.996 but elementwise range 0.866-1.121.
#' @noRd
.sim_check_Cenv <- function(C, n_env = NULL, arg = "C_env") {
  if (!is.matrix(C) || !is.numeric(C))
    stop("'", arg, "' must be a numeric matrix.", call. = FALSE)
  if (nrow(C) != ncol(C))
    stop("'", arg, "' must be square; got ", nrow(C), " x ", ncol(C), ".",
         call. = FALSE)
  if (nrow(C) < 3L)
    stop("'", arg, "' needs at least 3 environments; got ", nrow(C), ".",
         call. = FALSE)
  if (max(abs(C - t(C))) > 1e-8)
    stop("'", arg, "' must be symmetric.", call. = FALSE)

  d <- unname(diag(C))                      # unname(): all.equal compares names
  if (max(abs(d - 1)) > 1e-6) {
    if (all(d > 0)) {
      stop("'", arg, "' must have unit diagonal (max deviation ",
           signif(max(abs(d - 1)), 3), ").\n",
           "  If you passed a GENOMIC KINSHIP among lines, that is the wrong ",
           "object: sim_W() takes a q x q correlation among ENVIRONMENTS.\n",
           "  To convert a covariance, use cov2cor().", call. = FALSE)
    }
    stop("'", arg, "' must have unit diagonal.", call. = FALSE)
  }
  if (!is.null(n_env) && nrow(C) != n_env)
    stop("'", arg, "' is ", nrow(C), " x ", nrow(C), " but n_env = ", n_env,
         ". Rows of ", arg, " are ENVIRONMENTS, not genotypes.", call. = FALSE)
  invisible(TRUE)
}

#' AR(1) correlation matrix of size q with parameter rho.
#' @noRd
.sim_ar1_cor <- function(q, rho) rho^abs(outer(seq_len(q), seq_len(q), "-"))

#' C-vine partial-correlation sampler: always PSD, no rejection needed.
#' Partial correlations are drawn in [min.cor, max.cor]; the realised marginal
#' correlations then follow from the vine recursion (range only approximate).
#' @noRd
.sim_C_vine <- function(q, min.cor, max.cor) {
  P <- matrix(0, q, q); C <- diag(q)
  for (k in seq_len(q - 1L)) for (i in (k + 1L):q) {
    P[k, i] <- stats::runif(1, min.cor, max.cor)
    p <- P[k, i]
    if (k > 1L) for (l in (k - 1L):1L)
      p <- p * sqrt((1 - P[l, i]^2) * (1 - P[l, k]^2)) + P[l, i] * P[l, k]
    C[k, i] <- C[i, k] <- p
  }
  C
}

#' Structured environment-correlation generators. All are PSD by construction
#' except 'block' with a hostile between-correlation, which the caller repairs.
#' @noRd
.sim_C_structured <- function(q, structure, min.cor, max.cor, rho,
                              n_blocks, between.cor, n_factor) {
  if (structure == "equi") {
    r <- if (!is.null(rho)) rho else (min.cor + max.cor) / 2
    C <- matrix(r, q, q); diag(C) <- 1
  } else if (structure == "ar1") {
    r <- if (!is.null(rho)) rho else max.cor
    C <- .sim_ar1_cor(q, r)
  } else if (structure == "block") {
    wr  <- if (!is.null(rho)) rho else max.cor
    grp <- rep_len(seq_len(max(1L, n_blocks)), q)
    C <- matrix(between.cor, q, q)
    for (b in unique(grp)) C[grp == b, grp == b] <- wr
    diag(C) <- 1
    attr(C, "blocks") <- grp
  } else if (structure == "factor") {
    m  <- max(1L, n_factor)
    sd <- sqrt(max(1e-6, (min.cor + max.cor) / 2))
    L  <- matrix(stats::rnorm(q * m, sd = sd), q, m)
    C  <- stats::cov2cor(tcrossprod(L) + diag(stats::runif(q, 0.1, 0.5), q))
  } else {  # "vine"
    C <- .sim_C_vine(q, min.cor, max.cor)
  }
  C
}

#' Sample an AR(1) x AR(1) spatial field, returned as the first 'nplots' cells.
#' Field covariance is Ac (x) Ar; the diagonal is 1, so variance equals 'varc'.
#' @noRd
.sim_field_sample <- function(nplots, dim = NULL, rho = 0.6, varc = 1) {
  if (is.null(dim)) { nr <- ceiling(sqrt(nplots)); nc <- ceiling(nplots / nr) }
  else { nr <- dim[1L]; nc <- dim[2L] }
  if (nr * nc < nplots)
    stop("'field.dim' has fewer cells (", nr * nc, ") than plots in an ",
         "environment (", nplots, ").", call. = FALSE)
  Lr <- t(chol(.sim_ar1_cor(nr, rho)))
  Uc <- chol(.sim_ar1_cor(nc, rho))
  f  <- Lr %*% matrix(stats::rnorm(nr * nc), nr, nc) %*% Uc
  as.numeric(f)[seq_len(nplots)] * sqrt(varc)
}

#' Map Gaussian columns to a target marginal while preserving rank order
#' (a simple nonparanormal / copula transform) so sim_W covariables are not
#' forced to be Gaussian.
#' @noRd
.sim_marginal <- function(W, marginal = c("gaussian", "right-skewed",
                                          "heavy-tailed", "bounded")) {
  marginal <- match.arg(marginal)
  if (marginal == "gaussian") return(W)
  q <- nrow(W)
  trans <- function(z) {
    u <- rank(z, ties.method = "average") / (q + 1)
    switch(marginal,
      "right-skewed" = stats::qgamma(u, shape = 1.5, rate = 1),
      "heavy-tailed" = stats::qt(u, df = 3),
      "bounded"      = stats::qbeta(u, 2, 2))
  }
  matrix(apply(W, 2, trans), nrow = q, dimnames = dimnames(W))
}

#' Gaussian RBF kernel on covariables with a median-heuristic bandwidth,
#' normalised to the env_kernel() trace convention.
#' @noRd
.sim_rbf <- function(W, bandwidth = NULL) {
  D <- as.matrix(stats::dist(W))
  h <- if (is.null(bandwidth)) stats::median(D[lower.tri(D)]) else bandwidth
  if (!is.finite(h) || h <= 0) h <- 1
  K <- exp(-(D^2) / (2 * h^2))
  K / (sum(diag(K)) / nrow(K))
}


#' @title Simulate a Multi-Environment Trial with a Known Genetic Architecture
#'
#' @description
#' Core simulator. Generates genotype-by-environment phenotypes whose truth is
#' known exactly: a genomic kinship \eqn{K} among lines, a genetic (type-B)
#' correlation \eqn{C} among environments, and a per-environment heritability.
#' It is the reference-answer generator against which \code{\link{kernel_model}}
#' and \code{\link{scan_untested_envs}} can be checked, and it pairs with
#' \code{\link{sim_W}} to study how well an envirome recovers \eqn{C}.
#'
#' @details
#' \strong{Model.} For \eqn{n} lines and \eqn{q} environments, genetic values are
#' drawn from a separable (Kronecker) Gaussian
#' \deqn{\mathrm{vec}(G) \sim \mathcal{N}\!\left(0,\; K \otimes \Sigma_g\right),
#'       \qquad \Sigma_g = D^{1/2}\, C\, D^{1/2},\; D = \mathrm{diag}(\sigma^2_g),}
#' realised as \eqn{G = L_K\, U\, L_C^{\top}} with \eqn{L_K L_K^{\top}=K},
#' \eqn{L_C L_C^{\top}=\Sigma_g}. When \code{exact = TRUE} the driver \eqn{U} is
#' orthonormalised so that \eqn{U^{\top}U = nI}, making the realised moments
#' match their targets instead of merely scattering around them.
#'
#' \strong{Phenotypes.} For replicate \eqn{r},
#' \deqn{y_{ijr} = \mu_j + g_{ij} + \varepsilon_{ijr}, \qquad
#'       \varepsilon_{ijr}\sim\mathcal{N}(0,\sigma^2_{e,j}),\quad
#'       \sigma^2_{e,j} = \sigma^2_{g,j}\,\frac{1-h^2_j}{h^2_j}.}
#'
#' \strong{Heritability, two bases.} The plot-basis value
#' \eqn{h^2 = V_g/(V_g+V_e)} is what \code{min.h2}/\code{max.h2} target; the
#' line-mean basis \eqn{h^2/(h^2 + (1-h^2)/n_{rep})} is also reported.
#'
#' \strong{Realised environment correlation is measured on the K-whitened scale.}
#' Because lines are related through \eqn{K}, \code{cor(G)} across lines is not
#' \eqn{C}; forcing it to be is wrong. The realised correlation is defined as
#' \eqn{\mathrm{cor}(L_K^{-1} G)} and stored for validation.
#'
#' \strong{The two matrices are different objects.} \code{K} is \eqn{n\times n}
#' among lines; the environment correlation is \eqn{q\times q}. A genomic kinship
#' passed where an environment correlation is expected is the most common error
#' and is rejected by an elementwise unit-diagonal check.
#'
#' \strong{Alternative genetic mechanisms.} Besides the default separable draw,
#' three opt-in mechanisms generate the genetic values:
#' \itemize{
#'   \item \emph{Reaction norm} (supply \code{W}): genetic values are built
#'     causally from an envirome, \eqn{g_{ij} = m_i + \sum_k w_{jk} b_{ik}} with a
#'     main effect \eqn{m} and line sensitivities \eqn{b_{\cdot k}} both
#'     \eqn{K}-structured. The environment correlation is then \emph{induced} by
#'     \code{W} (and reported), not imposed. \code{rn.main} splits variance between
#'     the main effect and the reaction-norm interaction.
#'   \item \emph{Marker-based} (supply \code{markers}): QTL effects are sampled per
#'     environment with the target correlation \code{C}, giving large-effect
#'     architectures instead of the infinitesimal draw.
#'   \item \emph{Multiple kinships} (supply a \emph{list} for \code{K} with
#'     \code{K.weights}): independent additive, dominance or epistatic components
#'     are summed.
#' }
#' \code{gxe.specific} adds a line-independent, environment-specific deviation that
#' makes the covariance \emph{non-separable}. \code{rep.var} and \code{spatial}
#' (an AR1\eqn{\times}AR1 field) add trial realism to the phenotypes.
#'
#' @param K n x n genomic kinship AMONG LINES (e.g. \code{maizeG}), or a named
#'   \emph{list} of such matrices for multiple variance components (combined with
#'   \code{K.weights}). Made positive definite internally, so a singular
#'   marker-based kinship is fine. May be \code{NULL} when \code{markers} is given.
#' @param n_env integer. Number of environments; inferred from \code{C} if given.
#' @param min.h2,max.h2 numeric. Plot-basis heritability range, used when
#'   \code{h2} is not supplied. Must satisfy \eqn{0 < min.h2 \le max.h2 < 1}.
#' @param min.cor,max.cor numeric. Genetic correlation range among environments,
#'   used when \code{C} is not supplied. Mutually exclusive with \code{C}.
#' @param C q x q correlation AMONG ENVIRONMENTS (see \code{\link{sim_met_C}}).
#'   Mutually exclusive with \code{min.cor}/\code{max.cor}.
#' @param h2 numeric vector of per-environment heritabilities (length 1 or
#'   \code{n_env}). Mutually exclusive with \code{min.h2}/\code{max.h2}.
#' @param n_rep integer. Replicates per genotype x environment cell.
#' @param sg2 numeric. Genetic variance per environment (length 1 or \code{n_env}).
#' @param exact logical. Orthonormalise the draw so realised moments match the
#'   targets. Default \code{TRUE}.
#' @param env_effects numeric. Optional fixed environment means \eqn{\mu_j}.
#' @param subset numeric in (0, 1]. Fraction of cells retained, to induce
#'   unbalancedness.
#' @param subset.by "cell" or "gid_env". Whether missingness is drawn per
#'   observation or per genotype-by-environment cell.
#' @param K.weights numeric. Non-negative weights (rescaled to sum 1) for a list
#'   of kinships in \code{K}; ignored for a single matrix. Default equal weights.
#' @param W q x k envirome matrix. When supplied, genetic values are generated as
#'   reaction norms and the environment correlation is induced from \code{W}.
#' @param reaction.norm logical. Use the reaction-norm mechanism. Defaults to
#'   \code{TRUE} when \code{W} is supplied.
#' @param rn.main numeric in [0, 1]. Share of genetic variance assigned to the
#'   main (across-environment) effect under the reaction-norm mechanism; the rest
#'   is reaction-norm interaction. Default 0.5.
#' @param gxe.specific numeric in [0, 1). Share of genetic variance made
#'   line-independent and environment-specific, breaking separability. Default 0.
#' @param markers n x m marker matrix. When supplied, a marker/QTL mechanism
#'   replaces the infinitesimal draw and \code{K} is derived from the markers if
#'   not given.
#' @param n.qtl integer. Number of markers acting as QTL (sampled) when
#'   \code{markers} is used; \code{NULL} uses all markers.
#' @param rep.var numeric. Variance of a random replicate/block effect added to
#'   phenotypes. Default 0.
#' @param spatial logical. Add an AR1\eqn{\times}AR1 spatial field to each
#'   environment. Default \code{FALSE}.
#' @param field.dim integer length-2 \code{c(nrow, ncol)}. Field layout per
#'   environment; defaults to a near-square grid sized to the plots.
#' @param spatial.rho numeric. Lag-1 spatial autocorrelation. Default 0.6.
#' @param spatial.var numeric. Variance of the spatial field. Default 1.
#' @param seed integer. RNG seed. The caller's RNG stream is restored on exit.
#' @param verbose logical. Print a summary. Default \code{TRUE}.
#'
#' @return
#' An object of class \code{"sim_met"}: a list with
#' \describe{
#'   \item{\code{data}}{long \code{data.frame} with \code{env}, \code{gid},
#'     \code{rep}, \code{value}.}
#'   \item{\code{C_env}}{the target environment correlation, with the realised
#'     correlation attached as attribute \code{"realised"}.}
#'   \item{\code{truth}}{list of everything known: genetic values \code{g},
#'     target/realised correlations, target/realised heritabilities on both
#'     bases, variances, means.}
#'   \item{\code{K}}{the positive-definite kinship actually used.}
#' }
#'
#' @examples
#' \dontrun{
#' data(maizeG)
#'
#' ## 1. Baseline: random C and random h2 within a range
#' met <- sim_met(maizeG, n_env = 10, min.cor = 0.2, max.cor = 0.8,
#'                min.h2 = 0.3, max.h2 = 0.7, seed = 1)
#' head(met$data)
#' env_cor(met, "realised")
#'
#' ## 2. Supply an explicit environment correlation (e.g. a year-like series)
#' C   <- sim_met_C(q = 8, structure = "ar1", rho = 0.6)
#' met <- sim_met(maizeG, C = C, h2 = 0.5, seed = 1)
#'
#' ## 3. Fixed per-environment heritabilities and environment means
#' met <- sim_met(maizeG, n_env = 4, h2 = c(0.3, 0.5, 0.6, 0.8),
#'                env_effects = c(2, 4, 6, 8), seed = 1)
#'
#' ## 4. Replicates and an unbalanced design (retain 70% of g x e cells)
#' met <- sim_met(maizeG, n_env = 6, n_rep = 3,
#'                subset = 0.7, subset.by = "gid_env", seed = 1)
#'
#' ## 5. Several genetic variance components (additive + dominance + epistasis)
#' met <- sim_met(K = list(add = maizeG, dom = maizeG, epi = maizeG),
#'                K.weights = c(0.6, 0.3, 0.1), n_env = 8, seed = 1)
#'
#' ## 6. Causal reaction-norm: C is INDUCED by the envirome, not imposed
#' Wenv <- scale(matrix(rnorm(10 * 6), 10, 6))
#' rn   <- sim_met(maizeG, W = Wenv, rn.main = 0.4, seed = 1)
#' env_cor(rn, "target")
#'
#' ## 7. Marker / QTL architecture instead of the infinitesimal draw
#' X   <- matrix(rbinom(200 * 500, 2, 0.3), 200, 500)
#' qtl <- sim_met(markers = X, n_env = 6, n.qtl = 30, seed = 1)
#'
#' ## 8. Non-separable GxE: 25% of genetic variance is environment-specific
#' met <- sim_met(maizeG, n_env = 6, gxe.specific = 0.25, seed = 1)
#'
#' ## 9. Field realism: replicate/block noise plus an AR1 x AR1 spatial trend
#' met <- sim_met(maizeG, n_env = 4, n_rep = 2, rep.var = 0.5,
#'                spatial = TRUE, field.dim = c(20, 15),
#'                spatial.rho = 0.7, spatial.var = 1.5, seed = 1)
#'
#' ## 10. Pair with a simulated envirome and diagnose the recovery
#' W <- sim_W(met$C_env, noise = 0.3, n_var = 50, seed = 1)
#' plot(met)
#' }
#'
#' @seealso \code{\link{sim_W}}, \code{\link{sim_met_C}}, \code{\link{env_cor}},
#'   \code{\link{kernel_model}}, \code{\link{scan_untested_envs}}
#'
#' @references
#' Costa-Neto, G., et al. (2023). Envirome-wide associations enhance
#' multi-environment prediction. \emph{G3} 13(2), jkac313.
#' @export
sim_met <- function(K = NULL, n_env = NULL, min.h2 = 0.2, max.h2 = 0.8,
                    min.cor = 0.0, max.cor = 0.8,
                    C = NULL, h2 = NULL,
                    n_rep = 1, sg2 = 1, exact = TRUE,
                    env_effects = NULL, subset = NULL,
                    subset.by = c("cell", "gid_env"),
                    K.weights = NULL,
                    W = NULL, reaction.norm = !is.null(W), rn.main = 0.5,
                    gxe.specific = 0,
                    markers = NULL, n.qtl = NULL,
                    rep.var = 0, spatial = FALSE, field.dim = NULL,
                    spatial.rho = 0.6, spatial.var = 1,
                    seed = NULL, verbose = TRUE) {

  subset.by <- match.arg(subset.by)

  # RNG hygiene: seed locally and restore the caller's stream on exit.
  if (!is.null(seed)) {
    if (exists(".Random.seed", envir = .GlobalEnv)) {
      .oldseed <- get(".Random.seed", envir = .GlobalEnv)
      on.exit(assign(".Random.seed", .oldseed, envir = .GlobalEnv), add = TRUE)
    }
    set.seed(seed)
  }

  # ---- mutually exclusive specifications: error, do not silently prefer ----
  cor_given <- !missing(min.cor) || !missing(max.cor)
  if (!is.null(C) && cor_given)
    stop("Supply either 'C' or ('min.cor', 'max.cor'), not both.", call. = FALSE)
  h2_given <- !missing(min.h2) || !missing(max.h2)
  if (!is.null(h2) && h2_given)
    stop("Supply either 'h2' or ('min.h2', 'max.h2'), not both.", call. = FALSE)
  if (gxe.specific < 0 || gxe.specific >= 1)
    stop("'gxe.specific' must lie in [0, 1).", call. = FALSE)
  if (rn.main < 0 || rn.main > 1) stop("'rn.main' must lie in [0, 1].", call. = FALSE)

  use_markers <- !is.null(markers)
  use_rn      <- isTRUE(reaction.norm) && !use_markers && !is.null(W)

  # ---- lines, kinship and its Cholesky factor(s) ---------------------------
  if (use_markers) {
    markers <- as.matrix(markers)
    n <- nrow(markers)
    gid_names <- rownames(markers); if (is.null(gid_names)) gid_names <- paste0("G", seq_len(n))
    Zall <- scale(markers); Zall[is.na(Zall)] <- 0
    Kpd  <- if (!is.null(K)) .sim_pd(as.matrix(K)) else .sim_pd(tcrossprod(Zall) / ncol(Zall))
    L_K_list <- list(t(chol(Kpd))); K.weights <- 1
  } else if (is.list(K)) {
    Ks <- lapply(K, function(k) .sim_pd(as.matrix(k)))
    if (length(unique(vapply(Ks, nrow, 1L))) != 1L)
      stop("All kinships in 'K' must share the same dimension.", call. = FALSE)
    if (is.null(K.weights)) K.weights <- rep(1, length(Ks))
    if (length(K.weights) != length(Ks))
      stop("'K.weights' must have one weight per matrix in 'K'.", call. = FALSE)
    K.weights <- K.weights / sum(K.weights)
    Kpd <- Reduce(`+`, Map(function(k, w) w * k, Ks, K.weights))
    L_K_list <- lapply(Ks, function(k) t(chol(k)))
    gid_names <- rownames(Ks[[1]])
    n <- nrow(Kpd); if (is.null(gid_names)) gid_names <- paste0("G", seq_len(n))
  } else {
    if (is.null(K)) stop("Supply 'K' (or 'markers').", call. = FALSE)
    K <- as.matrix(K)
    if (nrow(K) != ncol(K)) stop("'K' must be square.", call. = FALSE)
    Kpd <- .sim_pd(K)
    L_K_list <- list(t(chol(Kpd))); K.weights <- 1
    n <- nrow(Kpd); gid_names <- rownames(K)
    if (is.null(gid_names)) gid_names <- paste0("G", seq_len(n))
  }
  L_K <- t(chol(Kpd))

  # ---- environment correlation (may be induced by the envirome) ------------
  if (!is.null(C)) {
    .sim_check_Cenv(C, n_env = n_env, arg = "C")
    q <- nrow(C); C_target <- unclass(C)
  } else if (use_rn) {
    q <- nrow(as.matrix(W)); C_target <- NULL
  } else {
    if (is.null(n_env)) stop("Supply 'n_env' when 'C' is not given.", call. = FALSE)
    q <- n_env
    C_target <- unclass(sim_met_C(q, min.cor, max.cor))
  }
  env_names <- if (!is.null(C_target)) rownames(C_target) else rownames(as.matrix(W))
  if (is.null(env_names)) env_names <- paste0("E", seq_len(q))
  if (!is.null(C_target)) dimnames(C_target) <- list(env_names, env_names)

  # ---- heritabilities ------------------------------------------------------
  if (!is.null(h2)) {
    if (length(h2) == 1L) h2 <- rep(h2, q)
    if (length(h2) != q) stop("'h2' must be length 1 or n_env.", call. = FALSE)
    if (any(h2 <= 0 | h2 >= 1)) stop("'h2' must lie in (0, 1).", call. = FALSE)
  } else {
    if (min.h2 <= 0 || max.h2 >= 1 || min.h2 > max.h2)
      stop("Require 0 < min.h2 <= max.h2 < 1.", call. = FALSE)
    h2 <- stats::runif(q, min.h2, max.h2)
  }
  names(h2) <- env_names
  if (length(sg2) == 1L) sg2 <- rep(sg2, q)
  if (length(sg2) != q) stop("'sg2' must be length 1 or n_env.", call. = FALSE)

  # ---- genetic values by mechanism -----------------------------------------
  if (use_markers) {
    if (is.null(C_target)) C_target <- unclass(sim_met_C(q, min.cor, max.cor))
    dimnames(C_target) <- list(env_names, env_names)
    L_C <- t(chol(.sim_pd(diag(sqrt(sg2), q) %*% C_target %*% diag(sqrt(sg2), q), 1e-10)))
    Z <- Zall
    if (!is.null(n.qtl) && n.qtl < ncol(Z))
      Z <- Z[, sort(sample.int(ncol(Z), n.qtl)), drop = FALSE]
    A  <- matrix(stats::rnorm(ncol(Z) * q), ncol(Z), q) %*% t(L_C)
    Gv <- Z %*% A
    Gv <- sweep(Gv, 2, sqrt(sg2 / pmax(apply(Gv, 2, stats::var), 1e-12)), "*")
    mech <- "marker"
  } else if (use_rn) {
    Wc <- scale(as.matrix(W)); Wc[is.na(Wc)] <- 0
    p  <- ncol(Wc); sgbar <- mean(sg2)
    m0 <- as.numeric(L_K %*% stats::rnorm(n))
    m0 <- m0 / stats::sd(m0) * sqrt(sgbar * rn.main)
    RN <- (L_K %*% matrix(stats::rnorm(n * p), n, p)) %*% t(Wc)
    RN <- RN / sqrt(mean(apply(RN, 2, stats::var))) * sqrt(sgbar * (1 - rn.main))
    Gv <- matrix(m0, n, q) + RN
    Sig <- rn.main * matrix(1, q, q) + (1 - rn.main) * (tcrossprod(Wc) / p)
    C_target <- stats::cov2cor(Sig + diag(1e-8, q))
    dimnames(C_target) <- list(env_names, env_names)
    mech <- "reaction-norm"
  } else {
    Sg  <- diag(sqrt(sg2), q) %*% C_target %*% diag(sqrt(sg2), q)
    L_C <- t(chol(.sim_pd(Sg, 1e-10)))
    Gv  <- matrix(0, n, q)
    for (cc in seq_along(L_K_list)) {
      U <- if (isTRUE(exact)) .sim_orth(n, q) else matrix(stats::rnorm(n * q), n, q)
      Gv <- Gv + sqrt(K.weights[cc]) * (L_K_list[[cc]] %*% U %*% t(L_C))
    }
    mech <- "kronecker"
  }
  dimnames(Gv) <- list(gid_names, env_names)

  # ---- optional non-separable, line-independent GxE ------------------------
  if (gxe.specific > 0) {
    vg0 <- apply(Gv, 2, stats::var)
    Gv  <- sqrt(1 - gxe.specific) * Gv +
           sweep(matrix(stats::rnorm(n * q), n, q), 2, sqrt(gxe.specific * vg0), "*")
    dimnames(Gv) <- list(gid_names, env_names)
  }

  # ---- realised genetic correlation, on the K-WHITENED scale ---------------
  # cor(Gv) across lines is NOT C, because the lines are related through K;
  # the parameter is defined as cor(L_K^{-1} G).
  C_realised <- stats::cor(solve(L_K) %*% Gv)
  dimnames(C_realised) <- list(env_names, env_names)

  # ---- phenotypes ----------------------------------------------------------
  vg  <- apply(Gv, 2, stats::var)
  se2 <- if (mech == "kronecker") sg2 * (1 - h2) / h2 else vg * (1 - h2) / h2
  names(se2) <- env_names
  mu  <- if (is.null(env_effects)) rep(0, q) else {
    if (length(env_effects) != q)
      stop("'env_effects' must be length n_env.", call. = FALSE)
    env_effects
  }

  d <- expand.grid(gid = gid_names, env = env_names, rep = seq_len(n_rep),
                   KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  gi <- match(d$gid, gid_names); ei <- match(d$env, env_names)
  d$value <- mu[ei] + Gv[cbind(gi, ei)] + stats::rnorm(nrow(d), 0, sqrt(se2[ei]))

  if (rep.var > 0)
    d$value <- d$value + stats::rnorm(q * n_rep, 0, sqrt(rep.var))[(ei - 1L) * n_rep + d$rep]
  if (isTRUE(spatial))
    for (e in seq_len(q)) {
      rows <- which(ei == e)
      d$value[rows] <- d$value[rows] +
        .sim_field_sample(length(rows), field.dim, spatial.rho, spatial.var)
    }

  d$env <- factor(d$env, levels = env_names)
  d$gid <- factor(d$gid, levels = gid_names)
  d <- d[, c("env", "gid", "rep", "value")]

  # ---- unbalancedness ------------------------------------------------------
  if (!is.null(subset)) {
    if (subset <= 0 || subset > 1) stop("'subset' must be in (0, 1].", call. = FALSE)
    if (subset.by == "cell") {
      keep <- stats::runif(nrow(d)) < subset
    } else {
      cells <- unique(d[, c("gid", "env")])
      kc <- cells[stats::runif(nrow(cells)) < subset, , drop = FALSE]
      keep <- paste(d$gid, d$env) %in% paste(kc$gid, kc$env)
    }
    d <- d[keep, , drop = FALSE]
  }
  rownames(d) <- NULL

  # ---- realised heritabilities, BOTH bases ---------------------------------
  extra   <- rep.var + if (isTRUE(spatial)) spatial.var else 0
  h2_plot <- vg / (vg + se2 + extra)
  h2_line <- h2_plot / (h2_plot + (1 - h2_plot) / n_rep)

  C_out <- C_target
  attr(C_out, "realised") <- C_realised
  class(C_out) <- c("C_env", "matrix", "array")

  out <- list(
    data  = d,
    C_env = C_out,
    truth = list(g = Gv, mechanism = mech,
                 C_target = C_target, C_realised = C_realised,
                 h2_target = h2,
                 h2_realised_plot = h2_plot,
                 h2_realised_linemean = h2_line,
                 sg2 = sg2, se2 = se2, vg = vg, mu = mu,
                 n_rep = n_rep, exact = exact, K.weights = K.weights,
                 rn.main = if (use_rn) rn.main else NA_real_,
                 gxe.specific = gxe.specific, rep.var = rep.var,
                 spatial.var = if (isTRUE(spatial)) spatial.var else 0),
    K = Kpd)
  class(out) <- "sim_met"

  if (verbose) print(out)
  out
}

#' @title Print a Simulated Multi-Environment Trial
#' @param x a \code{sim_met} object.
#' @param ... ignored.
#' @return \code{x}, invisibly.
#' @seealso \code{\link{sim_met}}
#' @export
print.sim_met <- function(x, ...) {
  q <- nrow(x$C_env); n <- nrow(x$K)
  ct <- x$truth$C_target[lower.tri(x$truth$C_target)]
  cr <- x$truth$C_realised[lower.tri(x$truth$C_realised)]
  cat("<sim_met>\n")
  cat(sprintf("  mechanism .............. %s\n", x$truth$mechanism))
  cat(sprintf("  lines x environments ... %d x %d   (K is %d x %d, C_env is %d x %d)\n",
              n, q, n, n, q, q))
  cat(sprintf("  observations ........... %d  (n_rep = %d%s)\n", nrow(x$data),
              x$truth$n_rep,
              if (nrow(x$data) < n * q * x$truth$n_rep) ", unbalanced" else ""))
  cat(sprintf("  h2 target .............. %.3f - %.3f\n",
              min(x$truth$h2_target), max(x$truth$h2_target)))
  cat(sprintf("  h2 realised (plot) ..... %.3f - %.3f   max err %.3f\n",
              min(x$truth$h2_realised_plot), max(x$truth$h2_realised_plot),
              max(abs(x$truth$h2_realised_plot - x$truth$h2_target))))
  cat(sprintf("  h2 realised (line-mean). %.3f - %.3f\n",
              min(x$truth$h2_realised_linemean), max(x$truth$h2_realised_linemean)))
  cat(sprintf("  gcor target ............ %.3f - %.3f\n", min(ct), max(ct)))
  cat(sprintf("  gcor realised (whitened) %.3f - %.3f   max err %.3g\n",
              min(cr), max(cr), max(abs(cr - ct))))
  cat(sprintf("  exact .................. %s\n", x$truth$exact))
  invisible(x)
}

#' @title Coerce a Simulated MET to a data.frame
#' @description Returns the long phenotype table, so a \code{sim_met} object can
#'   be handed straight to modelling functions or inspected as a plain frame.
#' @param x a \code{sim_met} object.
#' @param row.names,optional passed to the default method for consistency; unused.
#' @param ... ignored.
#' @return A \code{data.frame} with one row per plot and columns \code{env},
#'   \code{gid}, \code{rep} and \code{value}.
#' @seealso \code{\link{sim_met}}
#' @export
as.data.frame.sim_met <- function(x, row.names = NULL, optional = FALSE, ...) {
  d <- x$data
  if (!is.null(row.names)) rownames(d) <- row.names
  d
}

#' @title Diagnostic Plot for a Simulated MET
#' @description Scatter of the target against the realised (K-whitened)
#'   off-diagonal environment correlations, with the 1:1 line.
#' @param x a \code{sim_met} object.
#' @param ... passed to \code{\link[graphics]{plot}}.
#' @return \code{x}, invisibly.
#' @seealso \code{\link{sim_met}}
#' @export
plot.sim_met <- function(x, ...) {
  ct <- x$truth$C_target[lower.tri(x$truth$C_target)]
  cr <- x$truth$C_realised[lower.tri(x$truth$C_realised)]
  graphics::plot(ct, cr, xlab = "target genetic correlation",
                 ylab = "realised (K-whitened) correlation",
                 main = sprintf("sim_met (%s)", x$truth$mechanism),
                 pch = 19, col = grDevices::adjustcolor("steelblue", 0.6), ...)
  graphics::abline(0, 1, lty = 2, col = "grey40")
  invisible(x)
}


#' @title Simulate an Envirome with a Known Relationship to an Environment Correlation
#'
#' @description
#' Core simulator. Given a target correlation \eqn{C_{env}} among environments,
#' it builds a \eqn{q \times n_{var}} covariable matrix \eqn{W} whose
#' environmental kernel recovers a \emph{known} share of \eqn{C_{env}}. It answers:
#' if the true similarity among my environments is \eqn{C_{env}}, what would an
#' envirome of this size and quality look like, and how much of that similarity
#' would it recover? Pairs with \code{\link{sim_met}}.
#'
#' @details
#' \strong{Construction.} With \eqn{C_{env}^{1/2}} the symmetric PSD square root
#' of the target, an orthonormal driver \eqn{O} (\eqn{q\times n_{var}}) and iid
#' noise \eqn{E},
#' \deqn{W = \sqrt{1-a}\; C_{env}^{1/2} O \; + \; \sqrt{a}\; E,}
#' and the environmental kernel follows the \code{\link{env_kernel}} convention
#' \deqn{K_W = \frac{W W^{\top}}{\mathrm{tr}(W W^{\top})/q}.}
#' Recovery is measured by the off-diagonal Pearson correlation
#' \eqn{r = \mathrm{cor}\big(\mathrm{offdiag}(K_W),\, \mathrm{offdiag}(C_{env})\big)}.
#'
#' \strong{Calibration.} When \code{calibrate = TRUE}, \code{noise} is on the
#' \emph{outcome} scale: the internal variance share \eqn{a} is solved by
#' \code{\link[stats]{uniroot}} so that \eqn{E[r] = 1 - \mathrm{noise}}. Because
#' \eqn{r} decreases in \eqn{a}, \code{noise = 0} gives \eqn{a = 0} (full
#' recovery) and \code{noise = 1} gives \eqn{a = 1} (pure noise).
#'
#' \strong{Collinearity.} With \code{collinearity > 0}, a share of columns is
#' collapsed onto \code{n_blocks} driver families, mimicking the redundancy of
#' real enviromes and lowering the effective rank.
#'
#' \strong{\code{C_env} is an environment correlation, not a kinship.} It must be
#' \eqn{q\times q} with unit diagonal. Passing a genomic kinship is rejected.
#'
#' \strong{Realism knobs.} \code{marginal} maps each covariable to a non-Gaussian
#' shape (rank-preserving), \code{kernel = "rbf"} builds a Gaussian instead of a
#' linear kernel, \code{family.sizes} groups columns into named covariable
#' families (as in real enviromes), and \code{na.frac} injects missing values
#' into the returned \eqn{W} \emph{after} the truth is computed, so downstream
#' cleaning (e.g. \code{\link{W_matrix}} imputation) can be tested against a known
#' answer.
#'
#' @param C_env q x q correlation AMONG ENVIRONMENTS, from any source
#'   (e.g. \code{\link{sim_met}}\code{$C_env}, \code{\link{env_cor}},
#'   \code{cov2cor(env_kernel(W)$envCov)}, or hand-specified). NOT a kinship.
#' @param noise numeric in [0, 1]. With \code{calibrate = TRUE}, the expected
#'   off-diagonal correlation between the simulated kernel and \code{C_env} is
#'   \code{1 - noise}.
#' @param n_var integer. Number of covariables (columns of \eqn{W}).
#' @param collinearity numeric in [0, 1). Share of columns collapsed onto drivers
#'   (also the within-family strength when \code{family.sizes} is used).
#' @param n_blocks integer. Collinear columns grouped into this many families.
#' @param calibrate logical. Solve for the internal variance share \eqn{a}
#'   numerically so \code{noise} is on the outcome scale. Default \code{TRUE}.
#' @param cal.nrep integer. Monte-Carlo replicates used during calibration.
#' @param marginal character. Marginal distribution of the covariables:
#'   \code{"gaussian"} (default), \code{"right-skewed"}, \code{"heavy-tailed"}
#'   or \code{"bounded"}. Applied rank-preserving, so recovery is barely changed.
#' @param kernel character. Environmental kernel: \code{"linear"} (default,
#'   \eqn{WW^\top}) or \code{"rbf"} (Gaussian).
#' @param bandwidth numeric or \code{NULL}. RBF bandwidth; \code{NULL} uses the
#'   median heuristic.
#' @param na.frac numeric in [0, 1). Fraction of entries set to \code{NA} in the
#'   returned \eqn{W} (truth is computed on the complete matrix first).
#' @param family.sizes integer vector (optionally named). Column counts per
#'   covariable family; must sum to \code{n_var}. Overrides \code{n_blocks}
#'   redundancy and names the columns by family.
#' @param var_names,env_names character. Optional names for columns / rows.
#' @param seed integer. RNG seed. The caller's RNG stream is restored on exit.
#' @param verbose logical. Print a summary. Default \code{TRUE}.
#'
#' @return
#' A \eqn{q \times n_{var}} matrix of class \code{"sim_W"}. Attributes:
#' \code{"C_env"} (the target), \code{"K_W"} (the realised kernel) and
#' \code{"explained"} (a list with \code{r}, \code{r2}, per-PC variance shares,
#' effective rank \code{eff_rank}, the internal share \code{a_internal}, the
#' \code{marginal}/\code{kernel}/\code{na.frac} used, and the requested/target
#' values).
#'
#' @examples
#' \dontrun{
#' C <- sim_met_C(q = 10, min.cor = 0.2, max.cor = 0.8, seed = 1)
#'
#' ## 1. An envirome that recovers ~70% of C (noise is on the outcome scale)
#' W <- sim_W(C, noise = 0.3, n_var = 50, seed = 1)
#' attr(W, "explained")$r
#'
#' ## 2. Perfect vs pure-noise enviromes
#' W_good  <- sim_W(C, noise = 0.0, n_var = 50)
#' W_noise <- sim_W(C, noise = 1.0, n_var = 50)
#'
#' ## 3. Redundant envirome: many collinear columns, low effective rank
#' W <- sim_W(C, noise = 0.3, n_var = 200, collinearity = 0.8, n_blocks = 5)
#'
#' ## 4. Non-Gaussian covariables (rank-preserving marginal transform)
#' W <- sim_W(C, noise = 0.3, n_var = 60, marginal = "right-skewed")
#' W <- sim_W(C, noise = 0.3, n_var = 60, marginal = "heavy-tailed")
#' W <- sim_W(C, noise = 0.3, n_var = 60, marginal = "bounded")
#'
#' ## 5. Gaussian RBF kernel instead of the linear WW'
#' W <- sim_W(C, noise = 0.3, n_var = 60, kernel = "rbf")
#' W <- sim_W(C, noise = 0.3, n_var = 60, kernel = "rbf", bandwidth = 2)
#'
#' ## 6. Inject missing values AFTER the truth is computed
#' W <- sim_W(C, noise = 0.3, n_var = 60, na.frac = 0.1)
#' sum(is.na(W))
#'
#' ## 7. Named covariable families (sizes must sum to n_var)
#' W <- sim_W(C, n_var = 12, family.sizes = c(temp = 4, rain = 5, rad = 3))
#' colnames(W)
#'
#' ## 8. Custom names and a faster/looser calibration
#' W <- sim_W(C, noise = 0.4, n_var = 30, cal.nrep = 5,
#'            env_names = paste0("Loc", 1:10),
#'            var_names = paste0("bio", 1:30), seed = 1)
#'
#' ## 9. Raw (uncalibrated) internal variance share
#' W <- sim_W(C, noise = 0.5, n_var = 50, calibrate = FALSE)
#' plot(W)
#' }
#'
#' @seealso \code{\link{sim_met}}, \code{\link{sim_met_C}}, \code{\link{sim_W_grid}},
#'   \code{\link{env_cor}}, \code{\link{env_kernel}}
#'
#' @references
#' Costa-Neto, G., et al. (2023). Envirome-wide associations enhance
#' multi-environment prediction. \emph{G3} 13(2), jkac313.
#' @export
sim_W <- function(C_env, noise = 0, n_var = 50, collinearity = 0,
                  n_blocks = 1, calibrate = TRUE, cal.nrep = 15L,
                  marginal = c("gaussian", "right-skewed", "heavy-tailed",
                               "bounded"),
                  kernel = c("linear", "rbf"), bandwidth = NULL, na.frac = 0,
                  family.sizes = NULL, var_names = NULL, env_names = NULL,
                  seed = NULL, verbose = TRUE) {

  marginal <- match.arg(marginal)
  kernel   <- match.arg(kernel)

  # RNG hygiene: seed locally and restore the caller's stream on exit.
  if (!is.null(seed)) {
    if (exists(".Random.seed", envir = .GlobalEnv)) {
      .oldseed <- get(".Random.seed", envir = .GlobalEnv)
      on.exit(assign(".Random.seed", .oldseed, envir = .GlobalEnv), add = TRUE)
    }
    set.seed(seed)
  }
  if (inherits(C_env, "sim_met"))
    stop("Pass the correlation matrix, not the sim_met object: use ",
         "sim_W(met$C_env, ...) or sim_W(env_cor(met), ...).", call. = FALSE)

  C_env <- unclass(C_env)
  .sim_check_Cenv(C_env, arg = "C_env")
  q <- nrow(C_env); k <- n_var
  if (noise < 0 || noise > 1) stop("'noise' must lie in [0, 1].", call. = FALSE)
  if (collinearity < 0 || collinearity >= 1)
    stop("'collinearity' must lie in [0, 1).", call. = FALSE)
  if (na.frac < 0 || na.frac >= 1)
    stop("'na.frac' must lie in [0, 1).", call. = FALSE)
  if (!is.null(family.sizes)) {
    if (sum(family.sizes) != k)
      stop("'family.sizes' must sum to n_var (", k, ").", call. = FALSE)
    fam.grp <- rep(seq_along(family.sizes), times = family.sizes)
    fam.cor <- if (collinearity > 0) collinearity else 0.8
  }

  Csq  <- .sim_psd_sqrt(C_env)
  Kfun <- if (kernel == "rbf") function(W) .sim_rbf(W, bandwidth) else .sim_GB

  # ---- one draw at internal variance share a ------------------------------
  draw <- function(a) {
    S <- Csq %*% .sim_orth(q, k)
    if (!is.null(family.sizes)) {
      for (f in unique(fam.grp)) {
        cols <- which(fam.grp == f)
        drv  <- S[, cols[1]]
        S[, cols] <- sqrt(fam.cor) * drv +
                     sqrt(1 - fam.cor) * S[, cols, drop = FALSE]
      }
    } else if (collinearity > 0) {
      nb <- round(k * collinearity)
      if (nb >= 2L) {
        idx <- sample.int(k, nb)
        grp <- split(idx, rep_len(seq_len(n_blocks), length(idx)))
        for (g in grp) {
          drv <- S[, g[1]]
          S[, g] <- sqrt(collinearity) * drv +
                    sqrt(1 - collinearity) * S[, g, drop = FALSE]
        }
      }
    }
    E  <- matrix(stats::rnorm(q * k), q, k)
    Wd <- sqrt(1 - a) * S + sqrt(a) * E
    .sim_marginal(Wd, marginal)
  }

  # ---- calibration: r is DECREASING in a, so bracketing is inverted -------
  r_at <- function(a, nrep = cal.nrep)
    mean(vapply(seq_len(nrep),
                function(i) .sim_offdiag_cor(Kfun(draw(a)), C_env),
                numeric(1)))

  if (isTRUE(calibrate)) {
    if (noise <= 0)      a <- 0
    else if (noise >= 1) a <- 1            # uniroot cannot bracket a boundary
    else {
      f <- function(a) r_at(a) - (1 - noise)
      a <- tryCatch(stats::uniroot(f, c(1e-6, 1 - 1e-6), tol = 1e-3)$root,
                    error = function(e) {
                      warning("Calibration failed; using noise as the raw ",
                              "variance share.", call. = FALSE)
                      noise
                    })
    }
  } else a <- noise

  W <- draw(a)
  rn <- if (!is.null(env_names)) env_names else rownames(C_env)
  if (is.null(rn)) rn <- paste0("E", seq_len(q))
  if (!is.null(var_names)) cn <- var_names
  else if (!is.null(family.sizes)) {
    fnm <- names(family.sizes); if (is.null(fnm)) fnm <- paste0("F", seq_along(family.sizes))
    cn  <- unlist(Map(function(nm, s) paste0(nm, "_", seq_len(s)), fnm, family.sizes),
                  use.names = FALSE)
  } else cn <- paste0("V", seq_len(k))
  dimnames(W) <- list(rn, cn)

  # ---- truth computed on the COMPLETE matrix, before injecting NA ----------
  K_W <- Kfun(W)
  r   <- .sim_offdiag_cor(K_W, C_env)
  ev  <- svd(scale(W, TRUE, FALSE))$d^2
  eff <- sum(ev)^2 / sum(ev^2)

  if (na.frac > 0) {
    n_na <- round(na.frac * length(W))
    if (n_na > 0) W[sample.int(length(W), n_na)] <- NA_real_
  }

  attr(W, "C_env") <- C_env
  attr(W, "K_W")   <- K_W
  attr(W, "explained") <- list(
    r = r, r2 = r^2,
    per_pc = ev / sum(ev), eff_rank = eff,
    a_internal = a, noise_requested = noise, r_target = 1 - noise,
    marginal = marginal, kernel = kernel, na.frac = na.frac)
  class(W) <- c("sim_W", "matrix", "array")

  if (verbose) print(W)
  W
}

#' @title Print a Simulated Envirome
#' @param x a \code{sim_W} object.
#' @param ... ignored.
#' @return \code{x}, invisibly.
#' @seealso \code{\link{sim_W}}
#' @export
print.sim_W <- function(x, ...) {
  e <- attr(x, "explained")
  cat("<sim_W>\n")
  cat(sprintf("  environments x covariables . %d x %d\n", nrow(x), ncol(x)))
  cat(sprintf("  r(K_W, C_env) off-diagonal . %.4f   (target %.4f)\n",
              e$r, e$r_target))
  cat(sprintf("  r2  %% of env relatedness ... %.1f%%\n", 100 * e$r2))
  cat(sprintf("  effective rank ............. %.2f of %d covariables\n",
              e$eff_rank, ncol(x)))
  cat(sprintf("  PC1 share .................. %.1f%%\n", 100 * e$per_pc[1]))
  cat(sprintf("  internal variance share a .. %.4f\n", e$a_internal))
  cat(sprintf("  marginal / kernel .......... %s / %s\n",
              e$marginal, e$kernel))
  if (!is.null(e$na.frac) && e$na.frac > 0)
    cat(sprintf("  missing values injected .... %.1f%%\n", 100 * e$na.frac))
  invisible(x)
}

#' @title Coerce a Simulated Envirome to a Plain Matrix
#' @description Strips the \code{sim_W} class and its diagnostic attributes,
#'   returning the bare \eqn{q \times n_{var}} covariable matrix for use with
#'   \code{\link{env_kernel}}, \code{\link{W_matrix}} consumers, or base matrix ops.
#' @param x a \code{sim_W} object.
#' @param ... ignored.
#' @return A numeric matrix (environments in rows, covariables in columns) with
#'   no \code{sim_W} class or attributes.
#' @seealso \code{\link{sim_W}}
#' @export
as.matrix.sim_W <- function(x, ...) {
  m <- unclass(x)
  attr(m, "C_env") <- NULL
  attr(m, "K_W") <- NULL
  attr(m, "explained") <- NULL
  m
}

#' @title Diagnostic Plot for a Simulated Envirome
#' @description Two panels: the realised kernel off-diagonals against the target
#'   environment correlation, and the scree of per-PC variance shares.
#' @param x a \code{sim_W} object.
#' @param ... passed to \code{\link[graphics]{plot}}.
#' @return \code{x}, invisibly.
#' @seealso \code{\link{sim_W}}
#' @export
plot.sim_W <- function(x, ...) {
  e  <- attr(x, "explained")
  kw <- attr(x, "K_W"); ce <- attr(x, "C_env")
  op <- graphics::par(mfrow = c(1, 2)); on.exit(graphics::par(op), add = TRUE)
  graphics::plot(ce[lower.tri(ce)], kw[lower.tri(kw)],
                 xlab = "target C_env (off-diagonal)",
                 ylab = "realised kernel (off-diagonal)",
                 main = sprintf("recovery r = %.3f", e$r),
                 pch = 19, col = grDevices::adjustcolor("darkgreen", 0.6), ...)
  graphics::abline(stats::lm(kw[lower.tri(kw)] ~ ce[lower.tri(ce)]),
                   lty = 2, col = "grey40")
  graphics::plot(seq_along(e$per_pc), e$per_pc, type = "b", pch = 19,
                 xlab = "principal component", ylab = "variance share",
                 main = sprintf("effective rank = %.1f", e$eff_rank))
  invisible(x)
}


#' @title Simulate a q x q Genetic Correlation Among Environments
#'
#' @description
#' Accessory to \code{\link{sim_met}} and \code{\link{sim_W}}. Builds a valid
#' \eqn{q \times q} correlation matrix among environments, so the same \code{C}
#' can seed both simulators explicitly.
#'
#' @details
#' Sampling off-diagonals uniformly from \code{[min.cor, max.cor]} does not
#' generally give a positive semi-definite matrix. The equicorrelated bound is
#' \eqn{\rho \ge -1/(q-1)}, so for \eqn{q = 6} nothing below \eqn{-0.2} is
#' attainable. A naive eigenvalue repair silently returns something else, so this
#' function errors up front when the request is infeasible and rejection-samples
#' otherwise, only falling back to eigenvalue clamping (with a warning) if
#' sampling fails.
#'
#' @param q integer. Number of environments (>= 3).
#' @param min.cor,max.cor numeric. Range for off-diagonal correlations. For the
#'   structured generators these bound the random inputs (vine partial correlations,
#'   factor loadings) rather than the realised off-diagonals exactly.
#' @param structure character. How the matrix is built: \code{"random"} (default,
#'   rejection-sampled uniform off-diagonals, as before), \code{"vine"} (C-vine
#'   partial correlations, always PSD without rejection), \code{"equi"}
#'   (equicorrelated), \code{"ar1"} (first-order autoregressive, for ordered
#'   environments such as a year series), \code{"block"} (block-diagonal
#'   mega-environments) or \code{"factor"} (a random factor model
#'   \eqn{\Lambda\Lambda^\top + \Psi}).
#' @param rho numeric or \code{NULL}. The correlation used by \code{"equi"}
#'   (off-diagonal), \code{"ar1"} (lag-1) and \code{"block"} (within-block). If
#'   \code{NULL}, a sensible value is derived from \code{min.cor}/\code{max.cor}.
#' @param n_blocks integer. Number of mega-environments when \code{structure = "block"}.
#' @param between.cor numeric. Between-block correlation when \code{structure = "block"}.
#' @param n_factor integer. Number of latent factors when \code{structure = "factor"}.
#' @param env_names character. Optional environment names.
#' @param max.tries integer. Rejection-sampling attempts before repair.
#' @param seed integer. RNG seed.
#'
#' @return A \eqn{q \times q} correlation matrix of class \code{"C_env"} with a
#'   logical attribute \code{"repaired"} and an attribute \code{"requested"}
#'   recording the target range.
#'
#' @examples
#' C <- sim_met_C(q = 8, min.cor = 0.2, max.cor = 0.8, seed = 1)
#' range(C[lower.tri(C)])
#'
#' ## Structured mega-environments and an autoregressive year series
#' Cb <- sim_met_C(q = 9, structure = "block", n_blocks = 3, rho = 0.7,
#'                 between.cor = 0.1)
#' Car <- sim_met_C(q = 8, structure = "ar1", rho = 0.6)
#'
#' @seealso \code{\link{sim_met}}, \code{\link{sim_W}}, \code{\link{env_cor}}
#' @export
sim_met_C <- function(q, min.cor = 0.0, max.cor = 0.8,
                      structure = c("random", "vine", "equi", "ar1",
                                    "block", "factor"),
                      rho = NULL, n_blocks = 2L, between.cor = 0,
                      n_factor = 2L, env_names = NULL,
                      max.tries = 200L, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  structure <- match.arg(structure)
  if (q < 3L) stop("'q' must be at least 3.", call. = FALSE)
  if (min.cor > max.cor)
    stop("'min.cor' (", min.cor, ") exceeds 'max.cor' (", max.cor, ").",
         call. = FALSE)
  if (max.cor > 1 || min.cor < -1)
    stop("Correlations must lie in [-1, 1].", call. = FALSE)

  nm <- if (!is.null(env_names)) env_names else paste0("E", seq_len(q))

  # ---- constructive (structured) generators: PSD by construction ----------
  if (structure != "random") {
    C <- .sim_C_structured(q, structure, min.cor, max.cor, rho,
                           n_blocks, between.cor, n_factor)
    blocks <- attr(C, "blocks")
    ev <- eigen(C, symmetric = TRUE, only.values = TRUE)$values
    repaired <- FALSE
    if (min(ev) < 1e-8) {
      e <- eigen(C, symmetric = TRUE); e$values[e$values < 1e-8] <- 1e-8
      C <- stats::cov2cor(e$vectors %*% diag(e$values) %*% t(e$vectors))
      repaired <- TRUE
    }
    dimnames(C) <- list(nm, nm)
    attr(C, "requested") <- c(min.cor = min.cor, max.cor = max.cor)
    attr(C, "structure") <- structure
    attr(C, "repaired")  <- repaired
    if (!is.null(blocks)) attr(C, "blocks") <- blocks
    class(C) <- c("C_env", "matrix", "array")
    return(C)
  }

  bound <- -1 / (q - 1)
  if (max.cor < bound)
    stop(sprintf(paste0("Infeasible correlation range for q = %d.\n",
                        "  The equicorrelated lower bound is rho >= %.4f, ",
                        "but max.cor = %.4f.\n",
                        "  No PSD correlation matrix of this size has all ",
                        "off-diagonals that negative."),
                 q, bound, max.cor), call. = FALSE)

  nod <- q * (q - 1) / 2
  C <- NULL
  for (i in seq_len(max.tries)) {
    v <- stats::runif(nod, min.cor, max.cor)
    Ctry <- diag(q)
    Ctry[lower.tri(Ctry)] <- v
    Ctry[upper.tri(Ctry)] <- t(Ctry)[upper.tri(Ctry)]
    if (min(eigen(Ctry, symmetric = TRUE, only.values = TRUE)$values) > 1e-8) {
      C <- Ctry
      attr(C, "tries") <- i
      attr(C, "repaired") <- FALSE
      break
    }
  }

  if (is.null(C)) {
    v <- stats::runif(nod, min.cor, max.cor)
    C <- diag(q); C[lower.tri(C)] <- v; C[upper.tri(C)] <- t(C)[upper.tri(C)]
    ev <- eigen(C, symmetric = TRUE)
    ev$values[ev$values < 1e-8] <- 1e-8
    C <- ev$vectors %*% diag(ev$values) %*% t(ev$vectors)
    C <- stats::cov2cor(C)
    got <- range(C[lower.tri(C)])
    warning(sprintf(paste0("Rejection sampling failed in %d tries; C was ",
                           "repaired by eigenvalue clamping.\n",
                           "  Requested off-diagonal range [%.3f, %.3f]; ",
                           "realised [%.3f, %.3f]."),
                    max.tries, min.cor, max.cor, got[1], got[2]), call. = FALSE)
    attr(C, "tries") <- max.tries
    attr(C, "repaired") <- TRUE
  }

  nm <- if (!is.null(env_names)) env_names else paste0("E", seq_len(q))
  dimnames(C) <- list(nm, nm)
  attr(C, "requested") <- c(min.cor = min.cor, max.cor = max.cor)
  attr(C, "structure") <- "random"
  class(C) <- c("C_env", "matrix", "array")
  C
}


#' @title Extract the Correlation Among Environments
#'
#' @description
#' Accessory generic that pulls the environment correlation out of a
#' \code{\link{sim_met}} object or a bare matrix, so a hand-off to
#' \code{\link{sim_W}} never depends on internal structure.
#'
#' @param x an object (\code{sim_met}, matrix, ...).
#' @param type "target" (the specified correlation) or "realised" (what the
#'   simulation actually produced, on the K-whitened scale).
#' @param ... unused.
#'
#' @return A \eqn{q \times q} correlation matrix.
#'
#' @examples
#' \dontrun{
#' met <- sim_met(maizeG, n_env = 6, seed = 1)
#' env_cor(met, "target")
#' env_cor(met, "realised")
#' }
#'
#' @seealso \code{\link{sim_met}}, \code{\link{sim_W}}
#' @export
env_cor <- function(x, type = c("target", "realised"), ...) UseMethod("env_cor")

#' @rdname env_cor
#' @export
env_cor.sim_met <- function(x, type = c("target", "realised"), ...) {
  type <- match.arg(type)
  if (type == "target") x$C_env else attr(x$C_env, "realised")
}

#' @rdname env_cor
#' @export
env_cor.matrix <- function(x, type = c("target", "realised"), ...) {
  .sim_check_Cenv(x, arg = "x"); x
}

#' @rdname env_cor
#' @export
env_cor.sim_W <- function(x, type = c("target", "realised"), ...)
  attr(x, "C_env")

#' @rdname env_cor
#' @export
env_cor.default <- function(x, type = c("target", "realised"), ...)
  stop("No env_cor() method for class <", paste(class(x), collapse = "/"), ">.",
       call. = FALSE)


#' @title Sweep Envirome Size and Quality (Power Analysis)
#'
#' @description
#' Accessory to \code{\link{sim_W}}. Repeatedly simulates enviromes across a grid
#' of \code{noise}, \code{n_var} and \code{collinearity}, returning the mean
#' recovery and effective rank at each design point.
#'
#' @param C_env q x q correlation AMONG ENVIRONMENTS (see \code{\link{sim_W}}).
#' @param noise,n_var,collinearity numeric vectors defining the grid.
#' @param n_rep integer. Replicate simulations per design point.
#' @param calibrate logical. Passed to \code{\link{sim_W}}.
#' @param cal.nrep integer. Calibration replicates per \code{\link{sim_W}} call.
#' @param seed integer. RNG seed; each design point gets a reproducible sub-seed.
#' @param verbose logical. Print per-design-point progress. Default \code{TRUE}.
#'
#' @return A \code{data.frame} with one row per design point: \code{noise},
#'   \code{n_var}, \code{collinearity}, mean recovery \code{r} and its SD
#'   \code{r_sd}, \code{r2}, and mean \code{eff_rank}.
#'
#' @examples
#' \dontrun{
#' C <- sim_met_C(q = 10, min.cor = 0.2, max.cor = 0.8, seed = 1)
#' grid <- sim_W_grid(C, noise = c(0, 0.3, 0.6), n_var = c(20, 100),
#'                    collinearity = c(0, 0.8), n_rep = 5)
#' }
#'
#' @seealso \code{\link{sim_W}}, \code{\link{sim_met_C}}
#' @export
sim_W_grid <- function(C_env, noise = seq(0, 1, 0.25), n_var = c(10, 50, 200),
                       collinearity = c(0, 0.5, 0.9), n_rep = 10,
                       calibrate = TRUE, cal.nrep = 15L,
                       seed = NULL, verbose = TRUE) {
  if (!is.null(seed)) set.seed(seed)
  g <- expand.grid(noise = noise, n_var = n_var, collinearity = collinearity,
                   KEEP.OUT.ATTRS = FALSE)
  np <- nrow(g)
  cell_seeds <- if (!is.null(seed)) seed + seq_len(np) else rep(list(NULL), np)
  res <- do.call(rbind, lapply(seq_len(np), function(i) {
    if (isTRUE(verbose))
      message(sprintf("  [%d/%d] noise=%.2f  n_var=%d  collinearity=%.2f",
                      i, np, g$noise[i], g$n_var[i], g$collinearity[i]))
    cs <- if (is.list(cell_seeds)) cell_seeds[[i]] else cell_seeds[i]
    rr <- vapply(seq_len(n_rep), function(j) {
      W <- sim_W(C_env, noise = g$noise[i], n_var = g$n_var[i],
                 collinearity = g$collinearity[i], calibrate = calibrate,
                 cal.nrep = cal.nrep,
                 seed = if (is.null(cs)) NULL else cs * 1000L + j,
                 verbose = FALSE)
      e <- attr(W, "explained")
      c(e$r, e$r2, e$eff_rank)
    }, numeric(3))
    data.frame(g[i, , drop = FALSE],
               r = mean(rr[1, ]), r_sd = stats::sd(rr[1, ]),
               r2 = mean(rr[2, ]), eff_rank = mean(rr[3, ]))
  }))
  rownames(res) <- NULL
  res
}
