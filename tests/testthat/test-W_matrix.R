test_that("W_matrix builds an environment-by-covariable matrix", {
  data("maizeWTH", package = "EnvRtype")

  W <- W_matrix(env.data = maizeWTH,
                var.id = c("T2M", "PRECTOT"),
                statistic = "mean",
                verbose = FALSE)

  expect_true(is.matrix(W))
  expect_equal(nrow(W), length(unique(maizeWTH$env)))
  expect_setequal(rownames(W), unique(as.character(maizeWTH$env)))
  expect_true(is.numeric(W))
  expect_false(anyNA(W))
})

test_that("W_matrix respects the requested covariables", {
  data("maizeWTH", package = "EnvRtype")

  W <- W_matrix(env.data = maizeWTH,
                var.id = c("T2M"),
                statistic = "mean",
                center = FALSE, scale = FALSE,
                verbose = FALSE)

  expect_true(all(grepl("T2M", colnames(W))))
})
