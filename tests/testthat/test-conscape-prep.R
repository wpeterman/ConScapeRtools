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

test_that("conscape_prep supports centersize and buffer with center-restricted backend targets", {
  r <- make_test_raster(n = 10, vals = 1)
  prep <- conscape_prep(
    centersize = 5,
    buffer = 2,
    asc_dir = file.path(tempdir(), "center-window-prep-test"),
    r_target = r,
    r_mov = r,
    r_src = r,
    clear_dir = TRUE,
    landmark = 1L,
    progress = FALSE
  )

  expect_equal(prep$target_mode, "center")
  expect_equal(prep$centersize, 5L)
  expect_equal(prep$buffer, 2L)
  expect_true("tile_diagnostics" %in% names(prep))

  # In center mode, the target raster written to disk for ConScape is the
  # center-restricted raster. The on-disk finite count should therefore match
  # the per-tile n_finite_output_target diagnostic, not n_finite_target (which
  # records the full buffered raster as available context).
  target_file <- list.files(prep$target, pattern = "\\.asc$", full.names = TRUE)[1]
  target_tile <- terra::rast(target_file)
  target_vals <- terra::values(target_tile, mat = FALSE)
  expect_true(any(is.na(target_vals)))
  expect_gt(sum(is.finite(target_vals)), 0)
  expect_equal(
    sum(is.finite(target_vals)),
    prep$tile_diagnostics$n_finite_output_target[1]
  )
  expect_lt(
    sum(is.finite(target_vals)),
    prep$tile_diagnostics$n_finite_target[1]
  )
  expect_true(all(
    prep$tile_diagnostics$n_positive_output_target <=
      prep$tile_diagnostics$n_positive_target
  ))
  expect_lt(
    sum(prep$tile_diagnostics$n_positive_output_target),
    sum(prep$tile_diagnostics$n_positive_target)
  )
})

test_that("center-target prep distinguishes output cells from backend target cells", {
  r <- make_test_raster(n = 10, vals = 1)
  classic <- conscape_prep(
    tile_d = 5,
    tile_trim = 2,
    asc_dir = file.path(tempdir(), "classic-efficiency-prep-test"),
    r_target = r,
    r_mov = r,
    r_src = r,
    clear_dir = TRUE,
    landmark = 1L,
    progress = FALSE,
    target_mode = "full"
  )
  centered <- conscape_prep(
    centersize = 5,
    buffer = 2,
    asc_dir = file.path(tempdir(), "center-efficiency-prep-test"),
    r_target = r,
    r_mov = r,
    r_src = r,
    clear_dir = TRUE,
    landmark = 1L,
    progress = FALSE
  )

  assessment <- conscape_efficiency_assessment(
    classic = classic,
    centered = centered,
    reference = "classic"
  )

  expect_s3_class(assessment, "ConScapeEfficiencyAssessment")
  expect_equal(assessment$scenario, c("classic", "centered"))
  expect_equal(assessment$target_mode, c("full", "center"))
  # Both preps report the same potential (buffered) target volume per tile:
  # the difference is in what each backend actually solves.
  expect_equal(assessment$positive_target_cells[2], assessment$positive_target_cells[1])
  # Center mode writes a smaller target raster than full mode, so the backend-
  # consumed target count and the corresponding work proxy are both reduced.
  expect_lt(assessment$positive_output_target_cells[2], assessment$positive_target_cells[2])
  expect_lt(
    assessment$estimated_source_target_cells[2],
    assessment$estimated_source_target_cells[1]
  )
  expect_gt(assessment$source_target_work_reduction[2], 0)
  expect_equal(assessment$target_cell_reduction[2], 0)
  expect_gt(assessment$output_target_cell_reduction[2], 0)
  expect_gt(assessment$output_source_target_work_reduction[2], 0)
})

test_that("conscape_prep rejects incompatible raster extent", {
  r <- make_test_raster()
  shifted <- terra::rast(
    nrows = 20, ncols = 20,
    xmin = 1, xmax = 21,
    ymin = 0, ymax = 20,
    vals = 1
  )

  expect_error(
    conscape_prep(
      tile_d = 8,
      tile_trim = 3,
      asc_dir = file.path(tempdir(), "bad-extent-prep-test"),
      r_target = r,
      r_mov = shifted,
      r_src = r,
      clear_dir = TRUE,
      landmark = 5L,
      progress = FALSE
    ),
    "same extent"
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
