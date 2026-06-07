integration_jl_home <- function() {
  candidates <- c(
    Sys.getenv("CONSCAPERTOOLS_JL_HOME", ""),
    Sys.getenv("JULIA_BINDIR", "")
  )
  candidates <- candidates[nzchar(candidates)]
  if (length(candidates)) candidates[[1]] else ""
}

crop_cells_from_top_left <- function(r, n = 40L) {
  x0 <- terra::xmin(r)
  y1 <- terra::ymax(r)
  terra::crop(
    r,
    terra::ext(
      x0,
      x0 + n * terra::res(r)[1],
      y1 - n * terra::res(r)[2],
      y1
    ),
    snap = "near"
  )
}

surface_difference_summary <- function(tiled, untiled) {
  diff <- terra::values(tiled - untiled, mat = FALSE)
  diff <- diff[is.finite(diff)]
  c(
    max_abs = max(abs(diff)),
    mean_abs = mean(abs(diff)),
    rmse = sqrt(mean(diff^2))
  )
}

surface_agreement_summary <- function(tiled, untiled, mode) {
  data.frame(
    mode = mode,
    metric = c("btwn", "fcon"),
    rbind(
      surface_difference_summary(tiled$btwn, untiled[["btwn"]]),
      surface_difference_summary(tiled$fcon, untiled[["fcon"]])
    ),
    row.names = NULL
  )
}

expect_tiled_matches_untiled <- function(tiled, untiled, mode) {
  summaries <- surface_agreement_summary(tiled, untiled, mode = mode)
  # In this degenerate-coverage test the tiled and untiled answers differ only
  # by floating-point solver noise from ConScape's dense LU decomposition of
  # (I - W). Absolute bounds vary by metric magnitude (btwn ~ 1e4 in package
  # data, fcon ~ 1e1), so we use a *relative* tolerance against each
  # untiled-layer's max-abs value.
  ref_scale <- vapply(c("btwn", "fcon"), function(layer) {
    vals <- abs(terra::values(untiled[[layer]], mat = FALSE))
    vals <- vals[is.finite(vals)]
    if (length(vals)) max(vals) else 1
  }, numeric(1))
  rel_tol <- 1e-6
  expect_lt(max(summaries$max_abs / pmax(ref_scale, 1)), rel_tol,
            label = paste0(mode, " max_abs relative to untiled max_abs"))
  expect_lt(max(summaries$rmse / pmax(ref_scale, 1)), rel_tol,
            label = paste0(mode, " rmse relative to untiled max_abs"))
  invisible(summaries)
}

## NOTE on what this test actually exercises:
##
## The 40x40 demo crop with centersize = 20 and buffer = 20 cells is a
## *degenerate* tiling case: every tile's buffered window (60x60 cropped to
## 40x40) already covers the entire raster. Both tiled modes therefore see
## the same W matrix as the untiled run, and per-target Z columns are
## identical to untiled. Math:
##
##   - Full mode: 4 identical full-target surfaces -> mean = untiled.
##   - Center mode (after the center-target fix): each tile solves Z over
##     center-only targets but with the same full-landscape W; centers tile
##     disjointly and sum-mosaic = untiled.
##
## So this test is a smoke check of the plumbing (file IO, mosaic, mask),
## not a test of buffer-truncation accuracy. The convergence test in
## test-tiled-untiled-convergence.R exercises non-degenerate buffers and
## reports actual error decay vs untiled.

test_that("full and center-buffer tiled surfaces match untiled surfaces after mosaicing", {
  skip_on_cran()
  skip_if_not(
    identical(Sys.getenv("RUN_CONSCAPERTOOLS_INTEGRATION", ""), "true"),
    "Set RUN_CONSCAPERTOOLS_INTEGRATION=true to run Julia-backed ConScape checks."
  )

  jl_home <- integration_jl_home()
  skip_if_not(
    nzchar(jl_home),
    "Set CONSCAPERTOOLS_JL_HOME or JULIA_BINDIR to the Julia bin directory."
  )

  habitat_file <- system.file("extdata", "suitability.asc", package = "ConScapeRtools")
  affinity_file <- system.file("extdata", "affinity.asc", package = "ConScapeRtools")
  habitat <- crop_cells_from_top_left(terra::rast(habitat_file), n = 40L)
  affinity <- crop_cells_from_top_left(terra::rast(affinity_file), n = 40L)

  expect_equal(dim(habitat)[1:2], c(40L, 40L))
  expect_equal(dim(affinity)[1:2], c(40L, 40L))

  cell_size <- terra::res(habitat)[1]
  tile_width <- 20 * cell_size
  tile_trim <- 20 * cell_size

  full_prep <- conscape_prep(
    tile_d = tile_width,
    tile_trim = tile_trim,
    asc_dir = file.path(tempdir(), "agreement-full-prep"),
    r_target = habitat,
    r_mov = affinity,
    r_src = habitat,
    clear_dir = TRUE,
    landmark = 5L,
    progress = FALSE,
    target_mode = "full"
  )

  center_prep <- conscape_prep(
    centersize = 20L,
    buffer = 20L,
    window_units = "cells",
    asc_dir = file.path(tempdir(), "agreement-center-prep"),
    r_target = habitat,
    r_mov = affinity,
    r_src = habitat,
    clear_dir = TRUE,
    landmark = 5L,
    progress = FALSE
  )

  expect_equal(length(full_prep$tile_num), 4L)
  expect_equal(length(center_prep$tile_num), 4L)
  expect_equal(center_prep$target_mode, "center")

  efficiency <- conscape_efficiency_assessment(
    full = full_prep,
    center = center_prep,
    reference = "full"
  )
  # With the corrected center-target backend, center mode actually consumes
  # fewer targets than full mode: estimated_source_target_cells reflects the
  # backend's real work and is lower for center mode.
  expect_lt(
    efficiency$estimated_source_target_cells[efficiency$scenario == "center"],
    efficiency$estimated_source_target_cells[efficiency$scenario == "full"]
  )
  expect_lt(
    efficiency$estimated_output_source_target_cells[efficiency$scenario == "center"],
    efficiency$estimated_output_source_target_cells[efficiency$scenario == "full"]
  )

  full_tiled <- run_conscape(
    conscape_prep = full_prep,
    out_dir = file.path(tempdir(), "agreement-full-tiled"),
    theta = 0.15,
    distance_scale = 150,
    jl_home = jl_home,
    progress = FALSE
  )

  center_tiled <- run_conscape(
    conscape_prep = center_prep,
    out_dir = file.path(tempdir(), "agreement-center-tiled"),
    theta = 0.15,
    distance_scale = 150,
    jl_home = jl_home,
    progress = FALSE
  )

  untiled <- run_conscape(
    target_qualities = habitat,
    source_qualities = habitat,
    affinities = affinity,
    out_dir = file.path(tempdir(), "agreement-untiled"),
    theta = 0.15,
    distance_scale = 150,
    landmark = 5L,
    tile_trim = full_prep$tile_trim,
    jl_home = jl_home,
    progress = FALSE
  )

  summaries <- rbind(
    expect_tiled_matches_untiled(full_tiled, untiled, mode = "full"),
    expect_tiled_matches_untiled(center_tiled, untiled, mode = "center")
  )
  message(sprintf(
    "Tiled versus untiled max_abs: %s",
    paste(
      paste(summaries$mode, summaries$metric, signif(summaries$max_abs, 3),
            sep = "="),
      collapse = ", "
    )
  ))
})
