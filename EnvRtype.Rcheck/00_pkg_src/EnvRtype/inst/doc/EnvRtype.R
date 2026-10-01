## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse   = TRUE,
  comment    = "#>",
  fig.width  = 6.5,
  fig.height = 4.2,
  dpi        = 96
)

## ----setup--------------------------------------------------------------------
library(EnvRtype)

## ----echo = FALSE, results = "asis"-------------------------------------------
cat("
| Layer | Purpose | Key functions |
|---|---|---|
| 1 · Acquire | Pull weather, soil and geodata | `get_weather()`, `get_soil()`, `get_elevation()` |
| 2 · Process | Derive agro-meteorological parameters | `processWTH()`, `water_balance()`, `env_phenology()` |
| 3 · Characterise | Build covariables, typologies, indices | `W_matrix()`, `env_typing()`, `env_pca()` |
| 4 · Model & scan | Kernels, prediction, untested environments | `get_kernel()`, `kernel_model()`, `scan_untested_envs()` |
| Simulate | Ground truth to validate every layer | `sim_met()`, `sim_W()` |
")

## -----------------------------------------------------------------------------
data(maizeWTH)
data(maizeYield)
data(maizeG)

dim(maizeWTH)
dim(maizeYield)
dim(maizeG)

## ----eval = FALSE-------------------------------------------------------------
# env   <- c("NAIROBI", "NAKURU", "ELDORET")
# lat   <- c(-1.28, -0.30, 0.52)
# lon   <- c(36.82, 36.08, 35.27)
# start <- rep("2020-03-01", 3)
# end   <- rep("2020-07-31", 3)
# 
# wth <- get_weather(env.id = env, lat = lat, lon = lon,
#                    start.day = start, end.day = end)

## ----eval = FALSE-------------------------------------------------------------
# wth_h <- get_weather_hourly(env.id = env, lat = lat, lon = lon,
#                             start.day = start, end.day = end)

## ----eval = FALSE-------------------------------------------------------------
# wth <- get_weather_resumable(env.id = env, lat = lat, lon = lon,
#                              start.day = start, end.day = end,
#                              run.id = "kenya2020")
# 
# read_progress_log("kenya2020")           # what succeeded, what failed
# wth <- restart_from_log("kenya2020")     # resume, retrying failures

## ----eval = FALSE-------------------------------------------------------------
# soil <- get_soil(env.id = env, lat = lat, lon = lon)
# elev <- get_elevation(env.id = env, lat = lat, lon = lon)
# bioc <- get_bioclim(env.id = env, lat = lat, lon = lon, vars = 1:19)
# aez  <- get_AEZ(env.id = env, lat = lat, lon = lon)

## ----eval = FALSE-------------------------------------------------------------
# fut <- get_climate_scenario(env.id = env, lat = lat, lon = lon,
#                             scenario = "ssp245")

## -----------------------------------------------------------------------------
maizeWTH <- processWTH(env.data = maizeWTH, verbose = FALSE)
names(maizeWTH)

## -----------------------------------------------------------------------------
head(maizeWTH[, c("env", "daysFromStart", "GDD", "FRUE", "PETP")])

## ----eval = FALSE-------------------------------------------------------------
# d <- param_radiation(env.data = maizeWTH, merge = TRUE)
# d <- param_atmospheric(env.data = d, Alt = 600)
# d <- param_temperature(env.data = d, Tbase1 = 9, Tbase2 = 45,
#                        Topt1 = 26, Topt2 = 32)

## -----------------------------------------------------------------------------
sm <- summaryWTH(env.data = maizeWTH, env.id = "env",
                 days.id = "daysFromStart",
                 var.id = c("T2M", "PRECTOT", "FRUE"),
                 verbose = FALSE)
head(sm)

## -----------------------------------------------------------------------------
stages   <- c("VE", "V1_V6", "V6_VT", "VT_R1")
interval <- c(0, 7, 30, 65, 90)

sm_stage <- summaryWTH(env.data = maizeWTH, env.id = "env",
                       days.id = "daysFromStart",
                       var.id = c("FRUE", "PETP"),
                       by.interval = TRUE,
                       time.window = interval,
                       names.window = stages,
                       verbose = FALSE)
head(sm_stage)

## ----eval = FALSE-------------------------------------------------------------
# show_phenology("maize")        # inspect the bundled template
# phenology_templates()          # all available crops
# 
# ph <- env_phenology(env.data = maizeWTH, env.id = "env",
#                     day.id = "daysFromStart",
#                     crop = "maize", method = "gdd")

## ----eval = FALSE-------------------------------------------------------------
# wb <- water_balance(env.data = maizeWTH, soil.data = soil,
#                     env.id = "env", days.id = "daysFromStart",
#                     PREC = "PRECTOT",
#                     root.depth = 1000, root.init = 300,
#                     dap.root.max = 60)
# 
# summary_water_balance(wb, by.interval = TRUE,
#                       time.window = interval,
#                       names.window = stages)

## -----------------------------------------------------------------------------
W <- W_matrix(env.data = maizeWTH,
              var.id = c("FRUE", "PETP", "T2M_MAX", "T2M_MIN"),
              statistic = "mean",
              verbose = FALSE)
dim(W)
round(W[, 1:4], 2)

## -----------------------------------------------------------------------------
W_stage <- W_matrix(env.data = maizeWTH,
                    var.id = c("FRUE", "PETP"),
                    by.interval = TRUE,
                    time.window = interval,
                    names.window = stages,
                    verbose = FALSE)
dim(W_stage)
colnames(W_stage)

## -----------------------------------------------------------------------------
card_temp <- c( 9, 32, 45,Inf)

ET <- env_typing(env.data = maizeWTH, var.id = "T2M", env.id = "env",
                 cardinals = card_temp, verbose = FALSE)
head(ET)

## -----------------------------------------------------------------------------
ET_rain <- env_typing(env.data = maizeWTH, var.id = "PRECTOT",
                      env.id = "env", verbose = FALSE)
head(ET_rain)

## -----------------------------------------------------------------------------
Tm <- T_matrix(env.data = maizeWTH, var.id = "T2M", env.id = "env",
               cardinals = card_temp, verbose = FALSE)
dim(Tm)

## -----------------------------------------------------------------------------
cl <- env_cluster(W = W, k.max = 4, seed = 1234, verbose = FALSE)
cl$k
cl$clusters

## -----------------------------------------------------------------------------
cl$silhouette

## -----------------------------------------------------------------------------
round(cl$env.kinship, 2)

## -----------------------------------------------------------------------------
cl2 <- env_cluster(W = W, k = 2, seed = 1234, verbose = FALSE)
table(cl2$clusters)

## -----------------------------------------------------------------------------
Y_env <- tapply(maizeYield$value, maizeYield$env, mean)
Y_env <- Y_env[rownames(W)]
round(Y_env, 3)

## -----------------------------------------------------------------------------
if (requireNamespace("randomForest", quietly = TRUE)) {
  imp <- env_target_importance(Y = Y_env, W = W, k = 2, ntree = 200,
                               seed = 1234, verbose = FALSE)
  print(imp$k)
  print(imp$clusters)
  head(imp$importance)
}

## -----------------------------------------------------------------------------
if (requireNamespace("randomForest", quietly = TRUE)) {
  Ymat <- cbind(
    mean = Y_env,
    sd   = tapply(maizeYield$value, maizeYield$env, stats::sd)[rownames(W)]
  )
  imp2 <- env_target_importance(Y = Ymat, W = W, k = 2, ntree = 200,
                                seed = 1234, verbose = FALSE)
  head(imp2$importance)
}

## ----eval = FALSE-------------------------------------------------------------
# risk <- env_risk_profile(env.data = long_term_wth,
#                          site.id = "env", date.id = "YYYYMMDD",
#                          var.id = c("T2M", "PETP"),
#                          season.length = 150,
#                          crop = "maize", by.stage = TRUE)
# 
# tpe_weights(risk, stage = "flowering")

## ----eval = FALSE-------------------------------------------------------------
# project_risk(baseline = risk, scenario = fut)

## -----------------------------------------------------------------------------
cop <- env_copula(W = W, index = "pobs")
dim(cop)

## -----------------------------------------------------------------------------
idx <- env_indices(env.data = maizeWTH, env.id = "env",
                   day.id = "daysFromStart",
                   var.id = c("FRUE", "PETP", "T2M"),
                   interval = 10L, end.day = 90L,
                   statistic = "mean", verbose = FALSE)
dim(idx)

## -----------------------------------------------------------------------------
pca <- env_pca(W = idx, scale. = TRUE, verbose = FALSE)
head(pca$variance)

## -----------------------------------------------------------------------------
head(pca$meta)
dim(pca$scores)
dim(pca$loadings)

## ----fig.alt = "Scree plot of variance explained by environmental principal components"----
env_pca_scree(pca, verbose = FALSE)

## ----fig.alt = "Environments plotted in the space of the first two principal components"----
env_pca_biplot(pca, pc = c(1, 2), verbose = FALSE)

## ----fig.alt = "Loess-smoothed PC1 loadings across the growing season by weather factor"----
env_loading_curve(pca, pc = 1, span = 0.4, verbose = FALSE)

## ----fig.alt = "Heatmap of correlations among interval-resolved environmental indices"----
env_cor_heatmap(idx, order = "factor", verbose = FALSE)

## -----------------------------------------------------------------------------
env_mean <- tapply(maizeYield$value, maizeYield$env, mean)
env_mean <- env_mean[rownames(pca$scores)]

env_pc_associate(pca, y = env_mean, method = "pearson",
                 p.adjust = "BH", verbose = FALSE)

## ----eval = FALSE-------------------------------------------------------------
# zones <- soil_classification(soil.data = soil)

## -----------------------------------------------------------------------------
K_lin <- env_kernel(env.data = W, gaussian = FALSE, verbose = FALSE)
round(K_lin$envCov, 2)

## -----------------------------------------------------------------------------
K_lin <- env_kernel(env.data = W, gaussian = FALSE, verbose = FALSE)
round(K_lin$envCov, 2)

## -----------------------------------------------------------------------------
K_gau <- env_kernel(env.data = W, gaussian = TRUE, verbose = FALSE)
round(K_gau$envCov, 2)

## -----------------------------------------------------------------------------
K_deep <- env_kernel(env.data = W, deep.kernel = TRUE, deep.layers = 2,
                     verbose = FALSE)
round(K_deep$envCov, 2)

## -----------------------------------------------------------------------------
K_stage <- env_kernel(env.data = W_stage, stages = stages,
                      gaussian = TRUE, verbose = FALSE)
names(K_stage$envCov)

## -----------------------------------------------------------------------------
K <- get_kernel(K_G = list(G = maizeG),
                K_E = list(W = K_gau$envCov),
                data = maizeYield,
                model = "RNMM",
                env = "env", gid = "gid")
names(K)

## ----eval = FALSE-------------------------------------------------------------
# fit <- kernel_model(y = "value", data = maizeYield, random = K,
#                     env = "env", gid = "gid",
#                     iterations = 1000, burnin = 200, thining = 10,
#                     keep_effects = TRUE)
# 
# varcomp_summary(fit)

## ----eval = FALSE-------------------------------------------------------------
# Kd  <- decompose_kernels(K, keep_var = 0.99)
# fit <- kernel_model(y = "value", data = maizeYield, random = Kd,
#                     env = "env", gid = "gid", keep_effects = TRUE)

## ----eval = FALSE-------------------------------------------------------------
# cv <- kernel_cv(y = "value", data = maizeYield, random = K,
#                 env = "env", gid = "gid",
#                 folds = 5, scheme = "cv1", seed = 1234)

## ----eval = FALSE-------------------------------------------------------------
# fit_mc <- kernel_model_mc(y = "value", data = maizeYield, random = K,
#                           env = "env", gid = "gid", n_chains = 3, seed = 1)
# 
# cl2 <- cluster_environments(K_gau$envCov, k = 2, method = "hclust")
# fit_cl <- kernel_model_clustered(y = "value", data = maizeYield, random = K,
#                                  env = "env", gid = "gid",
#                                  clusters = cl2, use = "block")

## ----eval = FALSE-------------------------------------------------------------
# sc <- scan_untested_envs(object = fit,
#                          W_train = W, W_new = W_new,
#                          K_G = maizeG,
#                          level = 0.95,
#                          envelope = "flag")
# 
# coverage_summary(sc)

## ----eval = FALSE-------------------------------------------------------------
# map_scan(sc, coords, type = "coverage")
# map_scan(sc, coords, type = "performance", genotypes = "G01")
# map_scan(sc, coords, type = "which_won_where")
# 
# scan_spatial_table(sc, coords, lon_col = "lon", lat_col = "lat")

## ----eval = FALSE-------------------------------------------------------------
# gs <- grid_scan(fit, W_train = W, W_new = W_grid, K_G = maizeG,
#                 chunk_size = 500)

## ----eval = FALSE-------------------------------------------------------------
# planting_window_table(sc, dates = sowing_dates)
# plot_planting_window(sc, dates = sowing_dates, summarise = "genotype")
# best_planting_date(sc, dates = sowing_dates, exclude_flagged = TRUE)

## ----eval = FALSE-------------------------------------------------------------
# map_scan(sc, coords, date_col = "date", facet_by = "date",
#          genotypes = "G01")

## -----------------------------------------------------------------------------
C_env <- sim_met_C(q = 6, min.cor = 0.1, max.cor = 0.8, seed = 42)
round(C_env, 2)

## -----------------------------------------------------------------------------
set.seed(42)
sim <- sim_met(K = maizeG, C = C_env, min.h2 = 0.3, max.h2 = 0.7,
               n_rep = 2, seed = 42)
str(sim$data)

## -----------------------------------------------------------------------------
W_sim <- sim_W(C_env = sim$C_env, noise = 0.2, n_var = 30,
               collinearity = 0.5, seed = 42, verbose = FALSE)
dim(W_sim)

## ----eval = FALSE-------------------------------------------------------------
# K_sim <- env_kernel(env.data = W_sim, verbose = FALSE)
# Ks    <- get_kernel(K_G = list(G = maizeG),
#                     K_E = list(W = K_sim$envCov),
#                     data = sim$data, model = "RNMM",
#                     env = "env", gid = "gid")
# fit_s <- kernel_model(y = "value", data = sim$data, random = Ks,
#                       env = "env", gid = "gid")
# varcomp_summary(fit_s)      # compare against sim$truth

## ----eval = FALSE-------------------------------------------------------------
# sim_W_grid(C_env = sim$C_env,
#            noise = seq(0, 1, 0.25),
#            n_var = c(10, 50, 200),
#            collinearity = c(0, 0.5, 0.9),
#            n_rep = 10, seed = 42)

## ----eval = FALSE-------------------------------------------------------------
# citation("EnvRtype")

