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
  expect_true(all(c(
    "distance_scale", "exp_d", "tile_d", "tile_trim", "theta",
    "landmark", "trim_threshold", "overlap_area_factor", "diagnostics"
  ) %in% names(out)))
  expect_equal(out$distance_scale, out$exp_d)
  expect_equal(out$theta, 0.2)
  expect_equal(out$trim_threshold, 0.01)
  expect_gt(out$tile_d, 0)
  expect_gt(out$tile_trim, 0)
  expect_true(all(c(
    "requested_tile_trim", "effective_tile_trim", "expected_tile_count",
    "overlap_area_factor", "proximity_at_effective_trim"
  ) %in% names(out$diagnostics)))
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

test_that("tile_design calibration raster is fully connected and center-referenced", {
  r <- make_test_raster(n = 5, vals = 1)
  captured <- NULL

  testthat::local_mocked_bindings(
    juliaSetupOk = function() TRUE,
    ConScapeR_setup = function(...) invisible(TRUE),
    Grid = function(affinities, sources, targets, costs) {
      captured <<- list(affinities = affinities, sources = sources, targets = targets)
      "grid"
    },
    GridRSP = function(...) "rsp",
    expected_cost = function(...) {
      n <- terra::ncell(captured$affinities)
      matrix(seq_len(n * n), nrow = n, ncol = n)
    },
    juliaEval = function(...) TRUE,
    stopJulia = function() invisible(TRUE),
    .env = asNamespace("ConScapeRtools")
  )

  out <- tile_design(
    r_mov = r,
    r_target = r,
    max_d = 2,
    theta = 0.2,
    jl_home = "C:/Julia/bin",
    trim_threshold = 0.02,
    overlap_area_factor = 1.5,
    landmark = 5L
  )

  expect_false(any(is.nan(terra::values(captured$affinities, mat = FALSE))))
  expect_equal(unique(terra::values(captured$affinities, mat = FALSE)), 1)
  expect_gt(out$diagnostics$source_cell, 1)
  expect_lte(out$diagnostics$proximity_at_effective_trim, 0.02)
  expect_equal(out$diagnostics$tile_width_source, "overlap_area_factor")
  expect_equal(out$landmark, 5L)
})

test_that("tile_design supports explicit tile width and max tile cells", {
  r <- make_test_raster(n = 20, vals = 1)

  testthat::local_mocked_bindings(
    juliaSetupOk = function() TRUE,
    ConScapeR_setup = function(...) invisible(TRUE),
    Grid = function(...) "grid",
    GridRSP = function(...) "rsp",
    expected_cost = function(...) matrix(seq_len(6400), nrow = 80),
    juliaEval = function(...) TRUE,
    stopJulia = function() invisible(TRUE),
    .env = asNamespace("ConScapeRtools")
  )

  by_width <- tile_design(
    r_mov = r,
    r_target = r,
    max_d = 2,
    jl_home = "C:/Julia/bin",
    tile_width = 8,
    landmark = 2L
  )
  expect_equal(by_width$diagnostics$tile_width_source, "tile_width")
  expect_gte(by_width$tile_d, 8)

  by_cells <- tile_design(
    r_mov = r,
    r_target = r,
    max_d = 2,
    jl_home = "C:/Julia/bin",
    max_tile_cells = 16,
    landmark = 2L
  )
  expect_equal(by_cells$diagnostics$tile_width_source, "max_tile_cells")
  expect_gte(by_cells$tile_d, 4)

  expect_error(
    tile_design(
      r_mov = r,
      r_target = r,
      max_d = 2,
      jl_home = "C:/Julia/bin",
      tile_width = 8,
      max_tile_cells = 16
    ),
    "either tile_width or max_tile_cells"
  )
})

test_that("tile_design can keep Julia open when requested", {
  r <- make_test_raster(n = 5, vals = 1)
  stop_calls <- 0L

  testthat::local_mocked_bindings(
    juliaSetupOk = function() TRUE,
    ConScapeR_setup = function(...) invisible(TRUE),
    Grid = function(...) "grid",
    GridRSP = function(...) "rsp",
    expected_cost = function(...) matrix(seq_len(400), ncol = 1),
    juliaEval = function(...) TRUE,
    stopJulia = function() {
      stop_calls <<- stop_calls + 1L
      invisible(TRUE)
    },
    .env = asNamespace("ConScapeRtools")
  )

  out <- tile_design(
    r_mov = r,
    r_target = r,
    max_d = 2,
    jl_home = "C:/Julia/bin",
    stop_julia = FALSE
  )

  expect_s3_class(out, "ConScapeRtools_design")
  expect_equal(stop_calls, 0L)
})
