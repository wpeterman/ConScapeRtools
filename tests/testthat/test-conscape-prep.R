test_that("conscape_prep writes aligned tile directories and a valid mask", {
  target <- make_test_raster()
  vals <- terra::values(target, mat = FALSE)
  vals[1] <- NA_real_
  vals[2] <- 0.05
  terra::values(target) <- vals

  mov <- make_test_raster()
  src <- make_test_raster()
  out_dir <- file.path(tempdir(), "conscape-prep-test")

  prep <- conscape_prep(
    tile_d = 8,
    tile_trim = 3,
    asc_dir = out_dir,
    r_target = target,
    r_mov = mov,
    r_src = src,
    target_threshold = 0.1,
    clear_dir = TRUE,
    landmark = 5L,
    progress = FALSE
  )

  expect_s3_class(prep, "ConScapeRtools_prep")
  expect_true(dir.exists(prep$target))
  expect_true(dir.exists(prep$src))
  expect_true(dir.exists(prep$mov))
  expect_equal(prep$tile_trim, 5)
  expect_equal(prep$landmark, 5L)
  expect_equal(length(prep$tile_num), nrow(prep$cs_tiles))
  expect_equal(length(list.files(prep$target, pattern = "\\.asc$")), nrow(prep$cs_tiles))

  mask <- terra::rast(file.path(prep$asc_dir, "mask", "mask.asc"))
  mask_vals <- terra::values(mask, mat = FALSE)
  expect_equal(mask_vals[1], 0)
  expect_equal(mask_vals[2], 0)
  expect_true(all(mask_vals[-c(1, 2)] == 1))
})

test_that("conscape_prep rejects incompatible raster geometry", {
  r <- make_test_raster()
  different_res <- terra::rast(
    nrows = 20, ncols = 10,
    xmin = 0, xmax = 20,
    ymin = 0, ymax = 20,
    vals = 1
  )

  expect_error(
    conscape_prep(
      tile_d = 8,
      tile_trim = 3,
      asc_dir = file.path(tempdir(), "bad-prep-test"),
      r_target = r,
      r_mov = different_res,
      r_src = r,
      clear_dir = TRUE,
      landmark = 5L,
      progress = FALSE
    ),
    "same resolution"
  )
})

test_that("conscape_prep protects non-empty output directories", {
  r <- make_test_raster()
  out_dir <- file.path(tempdir(), "non-empty-prep-test")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  writeLines("keep", file.path(out_dir, "existing.txt"))

  expect_error(
    conscape_prep(
      tile_d = 8,
      tile_trim = 3,
      asc_dir = out_dir,
      r_target = r,
      r_mov = r,
      r_src = r,
      clear_dir = FALSE,
      landmark = 5L,
      progress = FALSE
    ),
    "asc_dir is not empty"
  )
})
