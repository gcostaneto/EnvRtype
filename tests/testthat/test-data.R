test_that("bundled maize data sets have the documented structure", {
  data("maizeYield", package = "EnvRtype")
  data("maizeWTH",   package = "EnvRtype")
  data("maizeG",     package = "EnvRtype")

  expect_s3_class(maizeYield, "data.frame")
  expect_identical(names(maizeYield), c("env", "gid", "value"))
  expect_equal(nrow(maizeYield), 750L)
  expect_type(maizeYield$value, "double")

  expect_s3_class(maizeWTH, "data.frame")
  expect_true(all(c("env", "daysFromStart", "T2M", "PRECTOT") %in% names(maizeWTH)))
  expect_gt(length(unique(maizeWTH$env)), 1L)

  expect_true(is.matrix(maizeG))
  expect_equal(nrow(maizeG), ncol(maizeG))
  expect_equal(nrow(maizeG), 150L)
  expect_true(isSymmetric(unname(maizeG)))
})
