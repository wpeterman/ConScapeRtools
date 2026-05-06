test_that("mosaic_conscape rebuilds a prepared tile set on the mask grid", {
  r <- make_test_raster()
  prep <- conscape_prep(
    tile_d = 8,
    tile_trim = 3,
    asc_dir = file.path(tempdir(), "mosaic-prep-test"),
    r_target = r,
    r_mov = r,
    r_src = r,
    clear_dir = TRUE,
    landmark = 5L,
    progress = FALSE
  )

  mask <- terra::rast(file.path(prep$asc_dir, "mask", "mask.asc"))
  out <- mosaic_conscape(prep$target, mask = mask, tile_trim = prep$tile_trim)

  expect_equal(dim(out), dim(r))
  expect_equal(as.vector(terra::ext(out)), as.vector(terra::ext(r)))
  expect_equal(terra::res(out), terra::res(r))
  expect_true(all(terra::values(out, mat = FALSE) == 1))
})

test_that("mosaic_conscape applies zero and NA mask cells", {
  root <- file.path(tempdir(), "mosaic-mask-test")
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  r <- make_test_raster(n = 5, vals = 5)
  mask <- make_test_raster(n = 5, vals = 1)
  mask_vals <- terra::values(mask, mat = FALSE)
  mask_vals[1] <- 0
  mask_vals[2] <- NA_real_
  terra::values(mask) <- mask_vals

  write_test_asc(r, file.path(root, "tile.asc"))
  out <- mosaic_conscape(root, mask = mask, tile_trim = 0)
  out_vals <- terra::values(out, mat = FALSE)

  expect_equal(out_vals[1], 0)
  expect_true(is.na(out_vals[2]))
  expect_equal(out_vals[3], 5)
})

test_that("mosaic_conscape validates inputs", {
  expect_error(
    mosaic_conscape(file.path(tempdir(), "does-not-exist"), tile_trim = -1),
    "tile_trim"
  )

  empty_dir <- file.path(tempdir(), "empty-mosaic-test")
  dir.create(empty_dir, recursive = TRUE, showWarnings = FALSE)
  expect_error(
    mosaic_conscape(empty_dir, tile_trim = 0),
    "No raster files"
  )
})
