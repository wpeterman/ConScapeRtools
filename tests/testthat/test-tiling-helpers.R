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

test_that("make_tiles warns for non-square cells", {
  make_tiles <- getFromNamespace("make_tiles", "ConScapeRtools")
  r <- terra::rast(
    nrows = 10, ncols = 10,
    xmin = 0, xmax = 20,
    ymin = 0, ymax = 10,
    vals = 1
  )

  expect_warning(
    tiles <- make_tiles(r, tile_d = 4, tile_trim = 1, landmark = 2L),
    "Non-square cells"
  )
  expect_equal(tiles$landmark, 2L)
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

test_that("tile_rast rejects non-empty output directories unless cleared", {
  make_tiles <- getFromNamespace("make_tiles", "ConScapeRtools")
  tile_rast <- getFromNamespace("tile_rast", "ConScapeRtools")
  r <- make_test_raster()
  tiles <- make_tiles(r, tile_d = 8, tile_trim = 3, landmark = 5L)
  out_dir <- file.path(tempdir(), "tile-rast-non-empty")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  writeLines("existing", file.path(out_dir, "existing.txt"))

  expect_error(
    tile_rast(r, tiles, out_dir = out_dir, clear_dir = FALSE, progress = FALSE),
    "Files currently exist"
  )
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

test_that("validate_conscape_tiles_idx validates missing and mismatched tile sets", {
  validate_tiles <- getFromNamespace("validate_conscape_tiles_idx", "ConScapeRtools")
  root <- file.path(tempdir(), "validate-errors-test")
  target_dir <- file.path(root, "target")
  source_dir <- file.path(root, "source")
  move_dir <- file.path(root, "move")
  dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(move_dir, recursive = TRUE, showWarnings = FALSE)

  expect_error(validate_tiles(target_dir, source_dir, move_dir), "No target tiles")

  r <- make_test_raster(n = 5, vals = 1)
  write_test_asc(r, file.path(target_dir, "r_1.asc"))
  write_test_asc(r, file.path(source_dir, "r_1.asc"))
  write_test_asc(r, file.path(move_dir, "r_1.asc"))
  write_test_asc(r, file.path(source_dir, "r_2.asc"))

  expect_error(validate_tiles(target_dir, source_dir, move_dir), "Unequal number")
})

test_that("ensure_tile_id honors existing and supplied identifiers", {
  ensure_tile_id <- getFromNamespace("ensure_tile_id", "ConScapeRtools")
  make_tiles <- getFromNamespace("make_tiles", "ConScapeRtools")
  tiles <- make_tiles(make_test_raster(), tile_d = 8, tile_trim = 3, landmark = 5L)$cs_tiles

  tiles$tile <- seq_len(nrow(tiles))
  tiles$tile_id <- NULL
  out <- ensure_tile_id(tiles)
  expect_true(out$id_col %in% c("tile", "tile_id"))

  tiles$custom_id <- seq_len(nrow(tiles)) + 10L
  out <- ensure_tile_id(tiles, id_col = "custom_id")
  expect_equal(out$id_col, "custom_id")

  tiles$tile <- NULL
  tiles$custom_id <- NULL
  out <- ensure_tile_id(tiles)
  expect_equal(out$id_col, "tile_id")
  expect_equal(out$cs_tiles$tile_id, seq_len(nrow(tiles)))

  expect_error(ensure_tile_id(tiles, id_col = "missing"), "id_col not found")
})
