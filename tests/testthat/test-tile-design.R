test_that("exp_distance validates probabilities and rates", {
  exp_distance <- getFromNamespace("exp_distance", "ConScapeRtools")

  expect_equal(exp_distance(0.5, lambda = 2), -log(0.5) / 2)
  expect_error(exp_distance(0, lambda = 2), "between 0 and 1")
  expect_error(exp_distance(1, lambda = 2), "between 0 and 1")
  expect_error(exp_distance(0.5, lambda = 0), "positive")
})

test_that("exp_opt penalizes decay values below the threshold", {
  exp_opt <- getFromNamespace("exp_opt", "ConScapeRtools")
  dists <- matrix(c(1, 10, 100), ncol = 1)

  expect_equal(exp_opt(1, dists = dists, cell = 3, threshold = 0.9), -9999)
  expect_gt(exp_opt(1, dists = dists, cell = 1, threshold = 0.001), -1)
})

test_that("tile_design uses target as source when source is omitted", {
  r <- make_test_raster(n = 5, vals = 1)

  testthat::local_mocked_bindings(
    juliaSetupOk = function() TRUE,
    ConScapeR_setup = function(...) invisible(TRUE),
    Grid = function(...) "grid",
    GridRSP = function(...) "rsp",
    expected_cost = function(...) matrix(seq_len(400), ncol = 1),
    juliaEval = function(...) TRUE,
    stopJulia = function() invisible(TRUE),
    .env = asNamespace("ConScapeRtools")
  )

  expect_message(
    out <- tile_design(
      r_mov = r,
      r_target = r,
      max_d = 2,
      theta = 0.2,
      jl_home = "C:/Julia/bin"
    ),
    "r_source not provided"
  )

  expect_s3_class(out, "ConScapeRtools_design")
  expect_named(out, c("exp_d", "tile_d", "tile_trim", "theta"))
  expect_equal(out$theta, 0.2)
  expect_gt(out$tile_d, 0)
  expect_gt(out$tile_trim, 0)
})

test_that("tile_design uses source as target when target is omitted", {
  r <- make_test_raster(n = 5, vals = 1)

  testthat::local_mocked_bindings(
    juliaSetupOk = function() TRUE,
    ConScapeR_setup = function(...) invisible(TRUE),
    Grid = function(...) "grid",
    GridRSP = function(...) "rsp",
    expected_cost = function(...) matrix(seq_len(400), ncol = 1),
    juliaEval = function(...) TRUE,
    stopJulia = function() invisible(TRUE),
    .env = asNamespace("ConScapeRtools")
  )

  expect_message(
    out <- tile_design(
      r_mov = r,
      r_source = r,
      max_d = 2,
      theta = 0.3,
      jl_home = "C:/Julia/bin"
    ),
    "r_target not provided"
  )
  expect_equal(out$theta, 0.3)
})

test_that("tile_design rejects invalid Julia setup", {
  r <- make_test_raster(n = 5, vals = 1)

  testthat::local_mocked_bindings(
    juliaSetupOk = function() FALSE,
    .env = asNamespace("ConScapeRtools")
  )

  expect_error(
    tile_design(
      r_mov = r,
      r_source = r,
      max_d = 2,
      jl_home = "bad-bin"
    ),
    "path to the Julia"
  )
})
