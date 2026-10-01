test_that("summaryWTH returns a long summary of weather variables", {
  data("maizeWTH", package = "EnvRtype")

  s <- summaryWTH(env.data = maizeWTH,
                  env.id = "env",
                  var.id = c("T2M", "PRECTOT"),
                  statistic = "mean",
                  verbose = FALSE)

  expect_s3_class(s, "data.frame")
  expect_true("env" %in% names(s))
  expect_gt(nrow(s), 0L)
})
