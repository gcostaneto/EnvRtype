test_that("env_pc_associate examples run without external workspace objects", {
  example_file <- tempfile(fileext = ".R")
  plot_file <- tempfile(fileext = ".pdf")
  on.exit(unlink(c(example_file, plot_file)), add = TRUE)
  grDevices::pdf(plot_file)
  on.exit(grDevices::dev.off(), add = TRUE)

  rd_path <- test_path("..", "..", "man", "env_pc_associate.Rd")
  rd <- if (file.exists(rd_path)) rd_path else
    tools::Rd_db("EnvRtype")[["env_pc_associate.Rd"]]
  tools::Rd2ex(rd, out = example_file, commentDontrun = TRUE,
               commentDonttest = FALSE)
  example_env <- new.env(parent = asNamespace("EnvRtype"))
  suppressMessages(sys.source(example_file, envir = example_env))

  expect_equal(dim(example_env$Y), c(160L, 20L))
  expect_equal(dim(example_env$M), c(160L, 90L))
  expect_gt(ncol(example_env$EM), 0L)
  expect_setequal(rownames(example_env$EM), rownames(example_env$pc$scores))
  expect_identical(names(example_env$assoc),
                   c("column", "PC", "r", "p", "p.adj", "n"))
  expect_true(all(is.finite(example_env$assoc$r)))
  expect_true(all(example_env$assoc$n == 20L))
  top <- example_env$assoc[1, ]
  scores <- example_env$pc$scores
  expect_equal(top$r, unname(cor(scores[, top$PC],
                                example_env$EM[rownames(scores), top$column])))
  expect_length(example_env$perm.max, 500L)
  expect_true(all(is.finite(example_env$perm.max)))
})
