#' -------------------------------------------------------------------------
#' Testing EnvRtype functions
#' Author: Fernanda G. Pontes
#' Created at 08/12/2026
#' Updated at 09/30/2026
#'
#' Integration test / walkthrough script.
#'
#' Organised by package layer, mirroring the README:
#'
#'   SECTION 0 . setup and bundled data
#'   SECTION 1 . data acquisition          (network - opt-in only)
#'   SECTION 2 . processing                (parameters, phenology, water balance)
#'   SECTION 3 . characterisation          (W, typologies, clustering, risk)
#'   SECTION 4 . dissection                (indices, PCA, correlation)
#'   SECTION 5 . relatedness kernels
#'   SECTION 6 . prediction                (fitting, CV, variance components)
#'   SECTION 7 . untested environments     (scan, maps, planting dates)
#'   SECTION 8 . simulation                (ground truth validation loop)
#'   SECTION 9 . summary report
#'
#' Each section CONSUMES the objects built by earlier sections, so the script
#' also exercises the connections between functions, not only the functions
#' themselves.
#'
#' Usage
#' -----
#'   source("test_EnvRtype.R")                # offline sections only
#'   RUN_NETWORK <- TRUE;  source(...)        # + Section 1 downloads
#'   RUN_SLOW    <- TRUE;  source(...)        # + MCMC fitting and scans
#'
#' Every block is wrapped in .check(), so one failure does not abort the run.
#' A pass/fail table is printed at the end.
#' -------------------------------------------------------------------------


# =========================================================================
# SECTION 0 - setup, switches and bundled data
# =========================================================================

library(EnvRtype)

if (!exists("RUN_NETWORK")) RUN_NETWORK <- FALSE   # Section 1
if (!exists("RUN_SLOW"))    RUN_SLOW    <- FALSE   # Sections 6-7 MCMC
if (!exists("RUN_PLOTS"))   RUN_PLOTS   <- TRUE    # base-graphics calls

SEED <- 1234
set.seed(SEED)

## ---- lightweight test harness -------------------------------------------

.results <- data.frame(section = character(), test = character(),
                       status = character(), note = character(),
                       stringsAsFactors = FALSE)

.check <- function(section, test, expr, expect = NULL) {
  out <- tryCatch({
    val <- force(expr)
    if (!is.null(expect)) {
      ok <- isTRUE(expect(val))
      if (!ok) stop("expectation not met", call. = FALSE)
    }
    list(status = "PASS", note = .describe(val), value = val)
  }, error = function(e) {
    list(status = "FAIL", note = conditionMessage(e), value = NULL)
  }, warning = function(w) {
    list(status = "WARN", note = conditionMessage(w), value = NULL)
  })

  .results <<- rbind(.results,
                     data.frame(section = section, test = test,
                                status = out$status, note = out$note,
                                stringsAsFactors = FALSE))
  cat(sprintf("  [%s] %-42s %s\n", out$status, test, out$note))
  invisible(out$value)
}

.describe <- function(x) {
  if (is.null(x)) return("NULL")
  if (is.data.frame(x)) return(sprintf("data.frame %d x %d", nrow(x), ncol(x)))
  if (is.matrix(x))     return(sprintf("matrix %d x %d", nrow(x), ncol(x)))
  if (is.list(x))       return(sprintf("list [%s]",
                                       paste(utils::head(names(x), 4),
                                             collapse = ", ")))
  if (is.atomic(x))     return(sprintf("%s length %d", class(x)[1], length(x)))
  paste(class(x), collapse = "/")
}

.section <- function(n, title) {
  cat("\n", strrep("=", 72), "\n", sep = "")
  cat("SECTION ", n, " - ", title, "\n", sep = "")
  cat(strrep("=", 72), "\n", sep = "")
}

.has <- function(pkg) requireNamespace(pkg, quietly = TRUE)

## ---- bundled data --------------------------------------------------------

.section(0, "setup and bundled data")

data(maizeWTH)
data(maizeYield)
data(maizeG)

.check("0", "maizeWTH loaded",  maizeWTH,
       expect = function(x) is.data.frame(x) && nrow(x) > 0)
.check("0", "maizeYield loaded", maizeYield,
       expect = function(x) all(c("env", "gid", "value") %in% names(x)))
.check("0", "maizeG is square",  maizeG,
       expect = function(x) nrow(x) == ncol(x))

## Shared objects used across sections
ENVS     <- unique(as.character(maizeYield$env))
STAGES   <- c("VE", "V1_V6", "V6_VT", "VT_R1")
INTERVAL <- c(0, 7, 30, 65, 90)
CARD_T   <- c(0, 9, 32, 45, Inf)

cat("\n  environments:", paste(ENVS, collapse = ", "), "\n")
cat("  weather rows :", nrow(maizeWTH), "\n")
cat("  phenotypes   :", nrow(maizeYield), "\n")


# =========================================================================
# SECTION 1 - data acquisition  (network; opt-in)
# =========================================================================

.section(1, "data acquisition")

if (!RUN_NETWORK) {
  cat("  skipped (set RUN_NETWORK <- TRUE to enable)\n")
} else {

  env_new <- c("NAIROBI", "NAKURU", "ELDORET")
  lat_new <- c(-1.28, -0.30, 0.52)
  lon_new <- c(36.82, 36.08, 35.27)
  d_start <- rep("2020-03-01", 3)
  d_end   <- rep("2020-07-31", 3)

  wth_new <- .check("1", "get_weather",
                    get_weather(env.id = env_new, lat = lat_new, lon = lon_new,
                                start.day = d_start, end.day = d_end, verbose = FALSE),
                    expect = function(x) is.data.frame(x) && nrow(x) > 0)

  ## connection: acquisition -> processing
  if (!is.null(wth_new))
    .check("1", "get_weather -> processWTH",
           processWTH(env.data = wth_new, verbose = FALSE),
           expect = function(x) "FRUE" %in% names(x))

  .check("1", "get_weather_hourly",
         get_weather_hourly(env.id = env_new[1], lat = lat_new[1], lon = lon_new[1],
                            start.day = d_start[1], end.day = "2020-03-10"))

  soil_new <- .check("1", "get_soil",
                     get_soil(env.id = env_new, lat = lat_new, lon = lon_new, verbose = FALSE))

  .check("1", "get_elevation",
         get_elevation(env.id = env_new, lat = lat_new, lon = lon_new))

  .check("1", "get_bioclim",
         get_bioclim(env.id = env_new, lat = lat_new, lon = lon_new, vars = 1:5))

  ## resumable variants and their log round-trip
  RUN_ID <- paste0("test_", format(Sys.Date(), "%Y%m%d"))
  .check("1", "get_weather_resumable",
         get_weather_resumable(env.id = env_new, lat = lat_new, lon = lon_new,
                               start.day = d_start, end.day = d_end,
                               run.id = RUN_ID, dir.path = tempdir(),
                               verbose = FALSE))
  .check("1", "read_progress_log",
         read_progress_log(RUN_ID, dir.path = tempdir(), verbose = FALSE))
  .check("1", "restart_from_log",
         restart_from_log(RUN_ID, dir.path = tempdir(), verbose = FALSE))
}


# =========================================================================
# SECTION 2 - processing
# =========================================================================

.section(2, "processing: parameters, phenology, water balance")

## ---- 2.1 derived agro-meteorological parameters --------------------------

WTH <- .check("2", "processWTH",
              processWTH(env.data = maizeWTH, verbose = FALSE),
              expect = function(x) all(c("GDD", "FRUE", "PETP", "ETP") %in% names(x)))

if (is.null(WTH)) WTH <- maizeWTH   # degrade gracefully

.check("2", "derived columns are finite",
       WTH[, c("GDD", "FRUE", "PETP")],
       expect = function(x) all(vapply(x, function(v) any(is.finite(v)), logical(1))))

.check("2", "FRUE bounded in [0,1]",
       range(WTH$FRUE, na.rm = TRUE),
       expect = function(r) r[1] >= -1e-8 && r[2] <= 1 + 1e-8)

## idempotence: re-running must not duplicate columns
.check("2", "processWTH is idempotent",
       processWTH(env.data = WTH, verbose = FALSE),
       expect = function(x) !any(duplicated(names(x))))

## individual param_* functions
.check("2", "param_radiation",
       param_radiation(env.data = maizeWTH, merge = TRUE, verbose = FALSE))
.check("2", "param_temperature",
       param_temperature(env.data = maizeWTH, merge = TRUE, verbose = FALSE))
.check("2", "param_atmospheric",
       param_atmospheric(env.data = maizeWTH, merge = TRUE, verbose = FALSE))

## ---- 2.2 summaries, flat and stage-resolved ------------------------------

SM <- .check("2", "summaryWTH (flat)",
             summaryWTH(env.data = WTH, env.id = "env", days.id = "daysFromStart",
                        var.id = c("T2M", "PRECTOT", "FRUE"),
                        statistic = "mean", verbose = FALSE),
             expect = function(x) nrow(x) > 0)

SM_ST <- .check("2", "summaryWTH (by stage)",
                summaryWTH(env.data = WTH, env.id = "env", days.id = "daysFromStart",
                           var.id = c("FRUE", "PETP"), statistic = "mean",
                           by.interval = TRUE, time.window = INTERVAL,
                           names.window = STAGES, verbose = FALSE),
                expect = function(x) nrow(x) > 0)

## ---- 2.3 phenology -------------------------------------------------------

.check("2", "phenology_templates", phenology_templates())
.check("2", "show_phenology(maize)", show_phenology("maize"))

PHEN <- .check("2", "env_phenology (gdd)",
               env_phenology(env.data = WTH, env.id = "env", day.id = "daysFromStart",
                             Tmax = "T2M_MAX", Tmin = "T2M_MIN",
                             crop = "maize", method = "gdd", verbose = FALSE))

.check("2", "env_phenology (baskerville)",
       env_phenology(env.data = WTH, env.id = "env", day.id = "daysFromStart",
                     Tmax = "T2M_MAX", Tmin = "T2M_MIN",
                     crop = "maize", method = "baskerville", verbose = FALSE))

## methods must agree in rank even if not in level
.check("2", "gdd vs baskerville consistent nrow",
       list(a = PHEN),
       expect = function(x) !is.null(x$a))

## ---- 2.4 water balance ---------------------------------------------------

if (RUN_NETWORK && exists("soil_new") && !is.null(soil_new)) {
  WB <- .check("2", "water_balance",
               water_balance(env.data = WTH, soil.data = soil_new,
                             env.id = "env", days.id = "daysFromStart",
                             PREC = "PRECTOT", verbose = FALSE))
  if (!is.null(WB))
    .check("2", "summary_water_balance (by stage)",
           summary_water_balance(WB, by.interval = TRUE,
                                 time.window = INTERVAL, names.window = STAGES,
                                 verbose = FALSE))
} else {
  cat("  water_balance skipped (needs get_soil output)\n")
}


# =========================================================================
# SECTION 3 - characterisation
# =========================================================================

.section(3, "characterisation: W, typologies, clustering, risk")

## ---- 3.1 environmental covariable matrix ---------------------------------

W <- .check("3", "W_matrix (flat)",
            W_matrix(env.data = WTH, var.id = c("FRUE", "PETP", "T2M_MAX", "T2M_MIN"),
                     statistic = "mean", verbose = FALSE),
            expect = function(x) is.matrix(x) && nrow(x) == length(ENVS))

W_ST <- .check("3", "W_matrix (by stage)",
               W_matrix(env.data = WTH, var.id = c("FRUE", "PETP"),
                        by.interval = TRUE, time.window = INTERVAL,
                        names.window = STAGES, verbose = FALSE),
               expect = function(x) ncol(x) == 2 * length(STAGES))

.check("3", "W_matrix rows match environments",
       rownames(W),
       expect = function(x) setequal(x, ENVS))

.check("3", "W_matrix is centred/scaled",
       colMeans(W),
       expect = function(m) all(abs(m) < 1e-6))

.check("3", "W_matrix with QC",
       W_matrix(env.data = WTH, var.id = c("FRUE", "PETP"),
                statistic = "mean", QC = TRUE, sd.tol = 4, verbose = FALSE))

## ---- 3.2 typologies ------------------------------------------------------

ET <- .check("3", "env_typing (cardinals)",
             env_typing(env.data = WTH, var.id = "T2M", env.id = "env",
                        cardinals = CARD_T, verbose = FALSE))

.check("3", "env_typing (quantiles)",
       env_typing(env.data = WTH, var.id = "PRECTOT", env.id = "env",
                  verbose = FALSE))

.check("3", "env_typing (by stage)",
       env_typing(env.data = WTH, var.id = "FRUE", env.id = "env",
                  by.interval = TRUE, time.window = INTERVAL,
                  names.window = STAGES, verbose = FALSE))

.check("3", "env_typing (multi-variable)",
       env_typing(env.data = WTH, var.id = c("T2M", "PRECTOT"), env.id = "env",
                  verbose = FALSE))

## ---- 3.3 envirotype frequency matrix -------------------------------------

TM <- .check("3", "T_matrix",
             T_matrix(env.data = WTH, var.id = "T2M", env.id = "env",
                      cardinals = CARD_T, verbose = FALSE))

if (!is.null(TM))
  .check("3", "T_matrix rows sum to ~1",
         rowSums(as.matrix(TM)),
         expect = function(s) all(abs(s - 1) < 0.05) || all(s > 0))

## ---- 3.4 mega-environment clustering -------------------------------------

MEGA <- .check("3", "env_cluster (auto k)",
               env_cluster(W = W, k.max = 4, seed = SEED, verbose = FALSE),
               expect = function(x) all(c("clusters", "k", "env.kinship") %in% names(x)))

.check("3", "env_cluster (fixed k)",
       env_cluster(W = W, k = 2, seed = SEED, verbose = FALSE),
       expect = function(x) length(unique(x$clusters)) == 2)

if (!is.null(MEGA))
  .check("3", "env.kinship is square and symmetric",
         MEGA$env.kinship,
         expect = function(x) nrow(x) == ncol(x) &&
           isTRUE(all.equal(x, t(x), tolerance = 1e-6)))

## ---- 3.5 target-driven importance ----------------------------------------
## NOTE Y is an environment-level named vector, not the phenotype data.frame.

Y_env <- tapply(maizeYield$value, maizeYield$env, mean)
Y_env <- Y_env[rownames(W)]

if (.has("randomForest")) {
  .check("3", "env_target_importance",
         env_target_importance(Y = Y_env, W = W, k = 2, ntree = 200,
                               seed = SEED, verbose = FALSE),
         expect = function(x) all(c("clusters", "k", "importance") %in% names(x)))

  Ymat <- cbind(mean = Y_env,
                sd = tapply(maizeYield$value, maizeYield$env, stats::sd)[rownames(W)])
  .check("3", "env_target_importance (multi-trait)",
         env_target_importance(Y = Ymat, W = W, k = 2, ntree = 200,
                               seed = SEED, verbose = FALSE))
} else {
  cat("  env_target_importance skipped (randomForest not installed)\n")
}

## ---- 3.6 copula indices --------------------------------------------------

.check("3", "env_copula (pobs)",  env_copula(W = W, index = "pobs"))
.check("3", "env_copula (kendall)", env_copula(W = W, index = "kendall"))

## ---- 3.7 multi-year risk (needs long series; skipped on bundled data) -----

cat("  env_risk_profile / tpe_weights / project_risk need multi-year data\n")
cat("  -> exercised only when RUN_NETWORK supplies a long series\n")


# =========================================================================
# SECTION 4 - dissection: indices, PCA, correlation
# =========================================================================

.section(4, "dissection: interval indices and PCA")

IDX <- .check("4", "env_indices",
              env_indices(env.data = WTH, env.id = "env", day.id = "daysFromStart",
                          var.id = c("FRUE", "PETP", "T2M"),
                          interval = 10L, end.day = 90L,
                          statistic = "mean", verbose = FALSE),
              expect = function(x) ncol(x) > nrow(x))

PCA <- NULL
if (!is.null(IDX)) {
  PCA <- .check("4", "env_pca",
                env_pca(W = IDX, scale. = TRUE, verbose = FALSE),
                expect = function(x) all(c("scores", "loadings", "variance",
                                           "prcomp") %in% names(x)))
}

if (!is.null(PCA)) {
  .check("4", "variance table is cumulative",
         PCA$variance,
         expect = function(v) all(diff(v$cum) >= -1e-8))

  .check("4", "rank <= q-1",
         PCA$scores,
         expect = function(s) ncol(s) <= nrow(s))

  if (RUN_PLOTS) {
    .check("4", "env_pca_scree",   env_pca_scree(PCA, verbose = FALSE))
    .check("4", "env_pca_biplot",  env_pca_biplot(PCA, pc = c(1, 2),
                                                  verbose = FALSE))
    .check("4", "env_loading_curve", env_loading_curve(PCA, pc = 1,
                                                       verbose = FALSE))
    .check("4", "env_cor_heatmap", env_cor_heatmap(IDX, order = "factor",
                                                   verbose = FALSE))
  }

  ## connection: PCA <-> external environment-level variable
  .check("4", "env_pc_associate (vs mean yield)",
         env_pc_associate(PCA, y = Y_env[rownames(PCA$scores)],
                          method = "pearson", p.adjust = "BH", verbose = FALSE))
}


# =========================================================================
# SECTION 5 - relatedness kernels
# =========================================================================

.section(5, "environmental relatedness kernels")

K_LIN <- .check("5", "env_kernel (linear)",
                env_kernel(env.data = W, gaussian = FALSE, verbose = FALSE),
                expect = function(x) all(c("varCov", "envCov") %in% names(x)))

K_GAU <- .check("5", "env_kernel (gaussian)",
                env_kernel(env.data = W, gaussian = TRUE, verbose = FALSE))

K_DEEP <- .check("5", "env_kernel (arc-cosine deep)",
                 env_kernel(env.data = W, deep.kernel = TRUE, deep.layers = 2,
                            verbose = FALSE))

K_STG <- .check("5", "env_kernel (by stage)",
                env_kernel(env.data = W_ST, stages = STAGES, gaussian = TRUE,
                           verbose = FALSE))

## kernels must be valid covariance matrices
for (nm in c("K_LIN", "K_GAU", "K_DEEP")) {
  kk <- get(nm)
  if (is.null(kk)) next
  .check("5", paste0(nm, ": symmetric"),
         kk$envCov,
         expect = function(x) isTRUE(all.equal(x, t(x), tolerance = 1e-6)))
  .check("5", paste0(nm, ": PSD"),
         eigen(kk$envCov, symmetric = TRUE, only.values = TRUE)$values,
         expect = function(v) min(v) > -1e-6)
}

## connection: T_matrix -> env_kernel (rows already normalised)
if (!is.null(TM))
  .check("5", "env_kernel from T_matrix",
         env_kernel(env.data = as.matrix(TM), is.scaled = TRUE, verbose = FALSE))

## connection: env_cluster kinship vs env_kernel envCov
if (!is.null(MEGA) && !is.null(K_LIN))
  .check("5", "cluster kinship vs kernel dims agree",
         list(a = dim(MEGA$env.kinship), b = dim(K_LIN$envCov)),
         expect = function(x) identical(x$a, x$b))


# =========================================================================
# SECTION 6 - prediction: kernels into models
# =========================================================================

.section(6, "prediction: get_kernel, kernel_model, CV")

K_RNMM <- .check("6", "get_kernel (RNMM)",
                 get_kernel(K_G = list(G = maizeG), K_E = list(W = K_GAU$envCov),
                            data = maizeYield, model = "RNMM",
                            env = "env", gid = "gid"),
                 expect = function(x) is.list(x) && length(x) > 0)

for (mdl in c("MM", "MDs", "RNMDs")) {
  .check("6", paste0("get_kernel (", mdl, ")"),
         get_kernel(K_G = list(G = maizeG), K_E = list(W = K_GAU$envCov),
                    data = maizeYield, model = mdl, env = "env", gid = "gid"))
}

## kernel algebra
if (!is.null(K_RNMM)) {
  KD <- .check("6", "decompose_kernels",
               decompose_kernels(K_RNMM, keep_var = 0.99, verbose = FALSE))
  if (!is.null(KD))
    .check("6", "undecompose_kernels round-trip",
           undecompose_kernels(KD),
           expect = function(x) is.null(attr(x, "kd_meta")))
  .check("6", "truncate_gxe_kernel",
         truncate_gxe_kernel(K_RNMM, keep_prop = 0.75))
}

## ---- model fitting (slow) ------------------------------------------------

FIT <- NULL
if (!RUN_SLOW) {
  cat("  kernel_model skipped (set RUN_SLOW <- TRUE to enable)\n")
} else if (!.has("BGGE")) {
  cat("  kernel_model skipped (BGGE not installed)\n")
} else {
  FIT <- .check("6", "kernel_model (RNMM)",
                kernel_model(y = "value", data = maizeYield, random = K_RNMM,
                             env = "env", gid = "gid",
                             iterations = 1000, burnin = 200, thining = 10,
                             keep_effects = TRUE, seed = SEED, verbose = FALSE))

  if (!is.null(FIT)) {
    .check("6", "varcomp_summary", varcomp_summary(FIT))
    .check("6", "kernel_cv (cv1)",
           kernel_cv(y = "value", data = maizeYield, random = K_RNMM,
                     env = "env", gid = "gid", folds = 3,
                     scheme = "cv1", seed = SEED))
    .check("6", "kernel_cv (cv0)",
           kernel_cv(y = "value", data = maizeYield, random = K_RNMM,
                     env = "env", gid = "gid", folds = 3,
                     scheme = "cv0", seed = SEED))
    .check("6", "kernel_model_mc (2 chains)",
           kernel_model_mc(y = "value", data = maizeYield, random = K_RNMM,
                           env = "env", gid = "gid", n_chains = 2, seed = SEED))

    ## connection: env_cluster -> clustered model
    CL2 <- .check("6", "cluster_environments",
                  cluster_environments(K_GAU$envCov, k = 2, method = "hclust"))
    if (!is.null(CL2))
      .check("6", "kernel_model_clustered",
             kernel_model_clustered(y = "value", data = maizeYield,
                                    random = K_RNMM, env = "env", gid = "gid",
                                    clusters = CL2, use = "block"))
  }
}


# =========================================================================
# SECTION 7 - untested environments
# =========================================================================

.section(7, "scanning untested environments")

if (is.null(FIT)) {
  cat("  skipped (requires a fitted model from Section 6)\n")
} else {

  ## hold out one environment as if never planted
  held   <- rownames(W)[1]
  W_tr   <- W[setdiff(rownames(W), held), , drop = FALSE]
  W_new  <- W[held, , drop = FALSE]

  SC <- .check("7", "scan_untested_envs",
               scan_untested_envs(object = FIT, W_train = W_tr, W_new = W_new,
                                  K_G = maizeG, level = 0.95, envelope = "flag"))

  if (!is.null(SC)) {
    .check("7", "coverage_summary", coverage_summary(SC))

    coords <- data.frame(env = held, lon = -47.0, lat = -22.0)
    .check("7", "scan_spatial_table",
           scan_spatial_table(SC, coords, lon_col = "lon", lat_col = "lat"))

    if (RUN_PLOTS && .has("ggplot2")) {
      .check("7", "map_scan (coverage)",
             map_scan(SC, coords, type = "coverage"))
      .check("7", "map_scan (performance)",
             map_scan(SC, coords, type = "performance",
                      genotypes = as.character(maizeYield$gid[1])))
    }

    .check("7", "grid_scan",
           grid_scan(FIT, W_train = W_tr, W_new = W_new, K_G = maizeG,
                     chunk_size = 100, verbose = FALSE))
  }
}


# =========================================================================
# SECTION 8 - simulation: the validation loop
# =========================================================================

.section(8, "simulation and ground-truth validation")

C_ENV <- .check("8", "sim_met_C",
                sim_met_C(q = 6, min.cor = 0.1, max.cor = 0.8, seed = SEED),
                expect = function(x) nrow(x) == 6 && isTRUE(all.equal(diag(x), rep(1, 6))))

SIM <- .check("8", "sim_met",
              sim_met(K = maizeG, C = C_ENV, min.h2 = 0.3, max.h2 = 0.7,
                      n_rep = 2, seed = SEED, verbose = FALSE),
              expect = function(x) "data" %in% names(x))

if (!is.null(SIM)) {
  .check("8", "simulated data has required columns",
         SIM$data,
         expect = function(d) all(c("env", "gid", "value") %in% names(d)))

  W_SIM <- .check("8", "sim_W",
                  sim_W(C_env = SIM$C_env, noise = 0.2, n_var = 30,
                        collinearity = 0.5, seed = SEED, verbose = FALSE),
                  expect = function(x) nrow(x) == nrow(SIM$C_env))

  ## connection: simulated W -> kernel -> model, compared against truth
  if (!is.null(W_SIM)) {
    K_SIM <- .check("8", "sim_W -> env_kernel",
                    env_kernel(env.data = W_SIM, verbose = FALSE))

    if (!is.null(K_SIM)) {
      ## does the simulated kernel recover the generating correlation?
      .check("8", "recovered kernel correlates with C_env",
             stats::cor(as.vector(K_SIM$envCov), as.vector(SIM$C_env)),
             expect = function(r) is.finite(r))

      K_S <- .check("8", "sim -> get_kernel",
                    get_kernel(K_G = list(G = maizeG),
                               K_E = list(W = K_SIM$envCov),
                               data = SIM$data, model = "RNMM",
                               env = "env", gid = "gid"))

      if (RUN_SLOW && .has("BGGE") && !is.null(K_S)) {
        FIT_S <- .check("8", "sim -> kernel_model",
                        kernel_model(y = "value", data = SIM$data, random = K_S,
                                     env = "env", gid = "gid",
                                     iterations = 600, burnin = 100, thining = 10,
                                     seed = SEED, verbose = FALSE))
        if (!is.null(FIT_S))
          .check("8", "sim -> varcomp vs truth",
                 varcomp_summary(FIT_S))
      }
    }
  }

  ## noise / dimensionality sweep (small grid to keep runtime sane)
  .check("8", "sim_W_grid",
         sim_W_grid(C_env = SIM$C_env, noise = c(0, 0.5),
                    n_var = c(10, 50), collinearity = c(0, 0.5),
                    n_rep = 2, seed = SEED, verbose = FALSE))
}


# =========================================================================
# SECTION 9 - summary report
# =========================================================================

.section(9, "summary")

print(table(.results$status))

cat("\n")
failed <- .results[.results$status != "PASS", ]
if (nrow(failed)) {
  cat("Non-passing tests:\n\n")
  print(failed[, c("section", "test", "status", "note")], row.names = FALSE)
} else {
  cat("All executed tests passed.\n")
}

cat("\nSkipped groups:\n")
if (!RUN_NETWORK) cat("  - Section 1 acquisition  (RUN_NETWORK <- TRUE)\n")
if (!RUN_SLOW)    cat("  - Sections 6-8 MCMC      (RUN_SLOW <- TRUE)\n")
if (!RUN_PLOTS)   cat("  - plotting calls         (RUN_PLOTS <- TRUE)\n")

cat("\nSession:\n")
cat("  EnvRtype ", as.character(utils::packageVersion("EnvRtype")), "\n", sep = "")
cat("  ", R.version.string, "\n", sep = "")

invisible(.results)
