test_that("env_kernel builds a symmetric environmental relatedness kernel", {
  data("maizeWTH", package = "EnvRtype")

  W <- W_matrix(env.data = maizeWTH,
                var.id = c("T2M", "PRECTOT", "SRAD"),
                statistic = "mean",
                verbose = FALSE)

  K <- env_kernel(env.data = W, verbose = FALSE)

  expect_type(K, "list")
  expect_true("envCov" %in% names(K))

  KE <- K$envCov
  expect_true(is.matrix(KE))
  expect_equal(nrow(KE), ncol(KE))
  expect_equal(nrow(KE), nrow(W))
  expect_true(isSymmetric(unname(KE)))
})
