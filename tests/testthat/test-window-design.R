test_that("window_design translates map-unit tile design to cells", {
  r <- make_test_raster(n = 40, vals = 1)

  out <- window_design(
    r = r,
    tile_d = 9,
    tile_trim = 3,
    landmark = 5L
  )

  expect_s3_class(out, "ConScapeRtools_window_design")
  expect_equal(out$centersize, 10L)
  expect_equal(out$buffer, 5L)
  expect_equal(out$tile_d, 10)
  expect_equal(out$tile_trim, 5)
  expect_equal(out$run_args$centersize, 10L)
  expect_equal(out$run_args$buffer, 5L)
  expect_equal(out$prep_args$target_mode, "center")
  expect_equal(out$classic_prep_args$tile_d, 10)
  expect_equal(out$diagnostics$total_window_width_cells, 20L)
  expect_equal(out$diagnostics$center_cells, 100L)
  expect_equal(out$diagnostics$source_window_cells, 400L)
  expect_equal(out$diagnostics$overlap_area_factor, 4)
  expect_equal(out$diagnostics$expected_windows, 9L)
})

test_that("window_design accepts explicit cell and map windows", {
  r <- make_test_raster(n = 20, vals = 1)

  cells <- window_design(r, centersize = 6, buffer = 2)
  expect_equal(cells$centersize, 6L)
  expect_equal(cells$buffer, 2L)
  expect_equal(cells$tile_d, 6)
  expect_equal(cells$tile_trim, 2)

  map <- window_design(
    r,
    centersize = 6,
    buffer = 2,
    window_units = "map"
  )
  expect_equal(map$centersize, 6L)
  expect_equal(map$buffer, 2L)
})

test_that("window_design consumes tile_design-style objects", {
  r <- make_test_raster(n = 40, vals = 1)
  design <- list(
    tile_d = 11,
    tile_trim = 4,
    landmark = 5L,
    theta = 0.2,
    distance_scale = 123
  )
  class(design) <- "ConScapeRtools_design"

  out <- window_design(r, design = design)

  expect_equal(out$centersize, 15L)
  expect_equal(out$buffer, 5L)
  expect_equal(out$theta, 0.2)
  expect_equal(out$distance_scale, 123)
})

test_that("window_design validates inputs", {
  r <- make_test_raster(n = 20, vals = 1)

  expect_error(window_design(matrix(1), tile_d = 1, tile_trim = 1), "SpatRaster")
  expect_error(window_design(r, tile_d = 1), "Supply design")
  expect_error(window_design(r, centersize = 1), "Supply both centersize and buffer")
  expect_error(window_design(r, centersize = 0, buffer = 1), "centersize")
  expect_error(window_design(r, centersize = 1, buffer = -1), "buffer")
  expect_error(window_design(r, tile_d = -1, tile_trim = 1), "tile_d")
  expect_error(window_design(r, tile_d = 1, tile_trim = -1), "tile_trim")
})
