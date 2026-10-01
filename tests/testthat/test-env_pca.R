test_that("env_pca runs on an environmental covariable matrix", {
  data("maizeWTH", package = "EnvRtype")

  W <- W_matrix(env.data = maizeWTH,
                var.id = c("T2M", "T2M_MAX", "PRECTOT", "SRAD"),
                statistic = "mean",
                verbose = FALSE)

  pca <- env_pca(W, verbose = FALSE)

  expect_s3_class(pca, "env_pca")
  expect_true(is.list(pca))
})
