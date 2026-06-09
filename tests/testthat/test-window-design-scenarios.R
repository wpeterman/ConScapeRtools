test_that("window_design_scenarios returns a tidy data frame from a window design", {
  r <- terra::rast(
    nrows = 60, ncols = 60, resolution = 1, xmin = 0, xmax = 60,
    ymin = 0, ymax = 60, crs = "EPSG:4326", vals = 1
  )
  wd <- window_design(r = r, tile_d = 20, tile_trim = 5, landmark = 5L)

  sc <- window_design_scenarios(r = r, design = wd)

  expect_s3_class(sc, "ConScapeRtools_scenarios")
  expect_s3_class(sc, "data.frame")
  expect_true(all(c("centersize", "buffer", "landmark",
                    "source_window_cells", "center_cells",
                    "target_cells_after_landmark",
                    "expected_windows", "ram_per_window_gb",
                    "total_ram_gb", "overlap_area_factor",
                    "workload_proxy") %in% names(sc)))
  expect_gte(nrow(sc), 1L)
  # centersize sorted descending
  expect_true(all(diff(sc$centersize) <= 0))
  # buffer/landmark constant across scenarios
  expect_equal(length(unique(sc$buffer)), 1L)
  expect_equal(length(unique(sc$landmark)), 1L)
  # arithmetic sanity
  expect_equal(sc$source_window_cells, (sc$centersize + 2 * sc$buffer)^2)
  expect_equal(sc$center_cells, sc$centersize^2)
  expect_equal(sc$target_cells_after_landmark,
               sc$center_cells / sc$landmark^2)
  expect_equal(sc$overlap_area_factor,
               sc$source_window_cells / sc$center_cells)
  # ram_per_window matches the documented worst-case envelope
  expect_equal(sc$ram_per_window_gb,
               8 * sc$source_window_cells *
                 sc$target_cells_after_landmark / 1e9)
})

test_that("default centersize sweep halves wd$centersize geometrically, snapped to landmark", {
  r <- terra::rast(
    nrows = 200, ncols = 200, resolution = 1, xmin = 0, xmax = 200,
    ymin = 0, ymax = 200, crs = "EPSG:4326", vals = 1
  )
  wd <- window_design(r = r, tile_d = 80, tile_trim = 10, landmark = 10L)

  sc <- window_design_scenarios(r = r, design = wd)

  # All candidates are multiples of landmark, at least landmark, and at
  # most the reference centersize.
  expect_true(all(sc$centersize %% sc$landmark == 0L))
  expect_true(all(sc$centersize >= sc$landmark))
  expect_true(all(sc$centersize <= wd$centersize))
})

test_that("custom centersize candidates are respected, dropped if too large", {
  r <- terra::rast(
    nrows = 50, ncols = 50, resolution = 1, xmin = 0, xmax = 50,
    ymin = 0, ymax = 50, crs = "EPSG:4326", vals = 1
  )
  wd <- window_design(r = r, tile_d = 20, tile_trim = 5, landmark = 5L)

  sc <- window_design_scenarios(
    r          = r,
    design     = wd,
    centersize = c(20, 10, 5)
  )
  expect_equal(sc$centersize, c(20L, 10L, 5L))

  # A centersize too large for the raster is dropped (expected_windows = 0).
  expect_error(
    window_design_scenarios(
      r = r, design = wd, centersize = c(200, 400)
    ),
    "No candidate"
  )
})

test_that("workers multiplier scales total_ram_gb but not per-window RAM", {
  r <- terra::rast(
    nrows = 60, ncols = 60, resolution = 1, xmin = 0, xmax = 60,
    ymin = 0, ymax = 60, crs = "EPSG:4326", vals = 1
  )
  wd <- window_design(r = r, tile_d = 20, tile_trim = 5, landmark = 5L)
  sc1 <- window_design_scenarios(r = r, design = wd, workers = 1L)
  sc8 <- window_design_scenarios(r = r, design = wd, workers = 8L)

  expect_equal(sc1$ram_per_window_gb, sc8$ram_per_window_gb)
  expect_equal(sc8$total_ram_gb, sc1$ram_per_window_gb * 8)
  expect_equal(attr(sc8, "workers"), 8L)
})

test_that("target_ram_gb flags the largest centersize fitting the budget", {
  r <- terra::rast(
    nrows = 60, ncols = 60, resolution = 1, xmin = 0, xmax = 60,
    ymin = 0, ymax = 60, crs = "EPSG:4326", vals = 1
  )
  wd <- window_design(r = r, tile_d = 20, tile_trim = 5, landmark = 5L)

  # Pick a budget such that the largest candidate exceeds it but at least
  # one smaller candidate does not.
  sc_unfiltered <- window_design_scenarios(r = r, design = wd)
  mid <- median(sc_unfiltered$total_ram_gb)
  sc <- window_design_scenarios(r = r, design = wd, target_ram_gb = mid)

  expect_true("recommended" %in% names(sc))
  expect_equal(sum(sc$recommended), 1L)
  rec_idx <- which(sc$recommended)
  # The recommended row's total_ram_gb is <= target.
  expect_lte(sc$total_ram_gb[rec_idx], mid)
  # And every larger centersize fails the budget.
  larger <- sc$centersize > sc$centersize[rec_idx]
  if (any(larger)) {
    expect_true(all(sc$total_ram_gb[larger] > mid))
  }
  expect_equal(attr(sc, "target_ram_gb"), mid)
})

test_that("target_ram_gb that no scenario fits returns no recommendation", {
  r <- terra::rast(
    nrows = 60, ncols = 60, resolution = 1, xmin = 0, xmax = 60,
    ymin = 0, ymax = 60, crs = "EPSG:4326", vals = 1
  )
  wd <- window_design(r = r, tile_d = 20, tile_trim = 5, landmark = 5L)

  sc <- window_design_scenarios(
    r = r, design = wd, target_ram_gb = 1e-20
  )
  expect_true("recommended" %in% names(sc))
  expect_equal(sum(sc$recommended), 0L)
})

test_that("window_design_scenarios errors helpfully without raster or design", {
  expect_error(
    window_design_scenarios(),
    "at least one of `r` or `design`"
  )
})

test_that("window_design_scenarios errors when buffer is missing", {
  r <- terra::rast(
    nrows = 60, ncols = 60, resolution = 1, xmin = 0, xmax = 60,
    ymin = 0, ymax = 60, crs = "EPSG:4326", vals = 1
  )
  expect_error(
    window_design_scenarios(r = r, centersize = c(20, 10)),
    "`buffer` is required"
  )
})

test_that("print method runs without error and includes the recommendation header", {
  r <- terra::rast(
    nrows = 60, ncols = 60, resolution = 1, xmin = 0, xmax = 60,
    ymin = 0, ymax = 60, crs = "EPSG:4326", vals = 1
  )
  wd <- window_design(r = r, tile_d = 20, tile_trim = 5, landmark = 5L)
  sc <- window_design_scenarios(r = r, design = wd,
                                target_ram_gb = 0.5, workers = 2L)
  expect_output(print(sc), "ConScapeRtools centersize scenarios")
  expect_output(print(sc), "workers modeled: 2")
  expect_output(print(sc), "target system RAM \\(GB\\): 0.5")
})
