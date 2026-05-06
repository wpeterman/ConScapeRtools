test_that("make_tiles aligns tile dimensions and overlap to landmark cells", {
  make_tiles <- getFromNamespace("make_tiles", "ConScapeRtools")
  r <- make_test_raster()

  tiles <- make_tiles(r, tile_d = 8, tile_trim = 3, landmark = 5L)

  expect_equal(tiles$overlap_cells, 5)
  expect_equal(tiles$tile_trim, 5)
  expect_equal(tiles$landmark, 5L)
  expect_equal(nrow(tiles$cs_tiles), 4)
  expect_equal(as.vector(terra::ext(tiles$cs_tiles)), as.vector(terra::ext(r)))
})

test_that("tile_rast writes overlapped ASCII tiles from a shared design", {
  make_tiles <- getFromNamespace("make_tiles", "ConScapeRtools")
  tile_rast <- getFromNamespace("tile_rast", "ConScapeRtools")
  r <- make_test_raster()
  tiles <- make_tiles(r, tile_d = 8, tile_trim = 3, landmark = 5L)
  out_dir <- file.path(tempdir(), "tile-rast-test")

  out <- tile_rast(r, tiles, out_dir = out_dir, clear_dir = TRUE, progress = FALSE)

  files <- list.files(out$asc_dir, pattern = "\\.asc$", full.names = TRUE)
  expect_equal(length(files), nrow(tiles$cs_tiles))
  first_tile <- terra::rast(files[1])
  expect_equal(dim(first_tile)[1:2], c(20, 20))
  expect_true(any(is.na(terra::values(first_tile, mat = FALSE))))
})

test_that("validate_conscape_tiles_idx removes empty tile triplets", {
  validate_tiles <- getFromNamespace("validate_conscape_tiles_idx", "ConScapeRtools")
  root <- file.path(tempdir(), "validate-tiles-test")
  target_dir <- file.path(root, "target")
  source_dir <- file.path(root, "source")
  move_dir <- file.path(root, "move")
  dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(move_dir, recursive = TRUE, showWarnings = FALSE)

  good <- make_test_raster(n = 5, vals = 1)
  bad <- make_test_raster(n = 5, vals = NA_real_)

  write_test_asc(good, file.path(target_dir, "r_1.asc"))
  write_test_asc(good, file.path(source_dir, "r_1.asc"))
  write_test_asc(good, file.path(move_dir, "r_1.asc"))
  write_test_asc(bad, file.path(target_dir, "r_2.asc"))
  write_test_asc(good, file.path(source_dir, "r_2.asc"))
  write_test_asc(good, file.path(move_dir, "r_2.asc"))

  out <- validate_tiles(target_dir, source_dir, move_dir, delete_bad = TRUE, quiet = TRUE)

  expect_equal(out$ok_idx, 1L)
  expect_equal(out$bad_idx, 2L)
  expect_false(file.exists(file.path(target_dir, "r_2.asc")))
  expect_false(file.exists(file.path(source_dir, "r_2.asc")))
  expect_false(file.exists(file.path(move_dir, "r_2.asc")))
})
