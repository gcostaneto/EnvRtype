clustered_fixture <- function(engine = "factored") {
  dat <- expand.grid(gid = paste0("G", 1:4), env = paste0("E", 1:3))
  dat$value <- c(2, 4, 3, 6, 5, 3, 7, 4, 3, 6, 4, 8)
  G <- diag(4) + matrix(0.2, 4, 4)
  E <- diag(3) + matrix(0.3, 3, 3)
  list(data = dat, clusters = c(E1 = 1, E2 = 1, E3 = 2),
       kernels = get_kernel(K_G = list(G = G), K_E = list(W = E),
                            data = dat, engine = engine))
}

test_that("cluster blocking handles dense and factored kernels consistently", {
  dense <- clustered_fixture("dense")
  mask <- outer(dense$clusters[as.character(dense$data$env)],
                dense$clusters[as.character(dense$data$env)], `==`)
  local_mocked_bindings(kernel_model = function(random, fixed, ...) {
    list(random = random, fixed = fixed)
  })

  for (engine in c("factored", "dense")) {
    fixture <- clustered_fixture(engine)
    original <- fixture$kernels
    for (use in c("block", "both", "fixed")) {
      fit <- kernel_model_clustered(
        "value", fixture$data, fixture$kernels, "env", "gid",
        fixture$clusters, use = use
      )
      if (use == "fixed") {
        expect_identical(fit$random, original)
      } else {
        blocked <- fit$random$KGE_GW
        expect_equal(blocked$Kernel, dense$kernels$KGE_GW$Kernel * mask,
                     ignore_attr = TRUE)
        expect_equal(dim(blocked$Kernel), c(12L, 12L))
        expect_identical(blocked$Type, "D")
        expect_null(blocked$decomp)
        expect_null(attr(fit$random, "kd_meta"))
        expect_identical(fit$random$KG_G, original$KG_G)
        expect_identical(fit$random$KE_W, original$KE_W)
      }
      expect_identical(fixture$kernels, original)
    }
  }
})

test_that("cluster blocking invalidates dense kernel decomposition shortcuts", {
  fixture <- clustered_fixture("dense")
  fixture$kernels$KGE_GW$kron <- list(stale = TRUE)
  fixture$kernels$KGE_GW$base <- diag(3)
  fixture$kernels$KGE_GW$idx <- rep(1:3, each = 4)
  local_mocked_bindings(kernel_model = function(random, ...) list(random = random))

  fit <- kernel_model_clustered(
    "value", fixture$data, fixture$kernels, "env", "gid", fixture$clusters
  )
  expect_null(fit$random$KGE_GW$kron)
  expect_null(fit$random$KGE_GW$base)
  expect_null(fit$random$KGE_GW$idx)
})

test_that("clustered models fit with the default factored kernel engine", {
  fixture <- clustered_fixture()
  for (use in c("block", "both", "fixed")) {
    fit <- kernel_model_clustered(
      "value", fixture$data, fixture$kernels, "env", "gid",
      fixture$clusters, use = use, iterations = 30, burnin = 10,
      thining = 2, seed = 1
    )
    expect_s3_class(fit, "kernel_model")
    expect_identical(fit$cluster_info$mechanism, use)
    expect_equal(fit$cluster_info$k, 2L)
  }
})
