## Tiled-vs-untiled convergence diagnostics for SENSITIVITY outputs.
##
## Companion test to test-tiled-untiled-convergence.R (which exercises
## the additive fcon / btwn outputs). This test exercises ConScape's
## landscape-summary derivatives:
##
##   - elasticity_quality           (sensitivity = "Q", unitless)
##   - elasticity_affinity_cost_linked  (sensitivity = "A&C=f(A)", unitless)
##
## Two things make sensitivity fundamentally different from fcon / btwn:
##
##   1. `ConScape.sensitivity` requires `target_equal_source = TRUE`. The
##      corrected `target_mode = "center"` workflow sets target = center
##      cells only and source = full buffered window, which violates that
##      precondition inside every tile. Therefore center-mode sensitivity
##      is rejected at the R API layer; this test only exercises
##      `target_mode = "full"`.
##
##   2. Sensitivity surfaces are tile-local landscape-summary derivatives,
##      not pairwise sums over target cells. There is no exact tile->whole
##      mosaic identity even in the degenerate-coverage limit (because each
##      tile's per-tile landscape summary differs from the global one).
##      `mosaic_method` is therefore mean (the "mosaic" reduction in
##      mosaic_conscape()), and the only claim we test is that the tiled
##      surface CONVERGES toward the untiled surface as buffer grows.
##
## The runtime is dominated by the untiled sensitivity reference (one
## landmark = 1 solve on the whole crop) plus the per-buffer tile runs.
## Keep n small (~40x40) so the test completes in a few minutes on a
## developer laptop.

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

surface_diff_stats <- function(tiled, untiled) {
  diff <- terra::values(tiled - untiled, mat = FALSE)
  ok <- is.finite(diff)
  diff <- diff[ok]
  vt <- terra::values(tiled, mat = FALSE)[ok]
  vu <- terra::values(untiled, mat = FALSE)[ok]
  cor_val <- if (length(diff) > 1L &&
                 stats::sd(vt) > 0 &&
                 stats::sd(vu) > 0) {
    suppressWarnings(stats::cor(vt, vu))
  } else {
    NA_real_
  }
  ref_scale <- max(abs(vu), 0)
  rel_max_abs <- max(abs(diff)) / pmax(ref_scale, 1)
  c(
    n           = length(diff),
    max_abs     = max(abs(diff)),
    mean_abs    = mean(abs(diff)),
    rmse        = sqrt(mean(diff^2)),
    cor         = cor_val,
    rel_max_abs = rel_max_abs
  )
}

stack_diff_stats <- function(tiled, untiled, layers, buffer) {
  rows <- lapply(layers, function(lyr) {
    if (!(lyr %in% names(tiled)) || !(lyr %in% names(untiled))) return(NULL)
    s <- surface_diff_stats(tiled[[lyr]], untiled[[lyr]])
    data.frame(
      layer       = lyr,
      buffer      = buffer,
      n           = s["n"],
      max_abs     = s["max_abs"],
      mean_abs    = s["mean_abs"],
      rmse        = s["rmse"],
      cor         = s["cor"],
      rel_max_abs = s["rel_max_abs"],
      row.names   = NULL,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, Filter(Negate(is.null), rows))
}

result_as_stack <- function(x, layers) {
  if (inherits(x, "SpatRaster")) {
    keep <- intersect(names(x), layers)
    if (length(keep) == 0L) return(NULL)
    return(x[[keep]])
  }
  present <- intersect(layers, names(x))
  if (length(present) == 0L) return(NULL)
  do.call(c, lapply(present, function(n) {
    r <- x[[n]]
    names(r) <- n
    r
  }))
}

test_that("Fix 1: sensitivity rejects center-mode prep at the R layer", {
  # Pure-R sanity check that does not require Julia. Confirms that the
  # combination of (a) target_mode = "center" prep and (b) a sensitivity
  # request errors before any Julia work is attempted, with a message that
  # points to the documented workaround.
  habitat_file  <- system.file("extdata", "suitability.asc", package = "ConScapeRtools")
  affinity_file <- system.file("extdata", "affinity.asc",  package = "ConScapeRtools")
  habitat       <- crop_cells_from_top_left(terra::rast(habitat_file),  n = 40L)
  affinity      <- crop_cells_from_top_left(terra::rast(affinity_file), n = 40L)

  center_prep <- conscape_prep(
    centersize    = 20L,
    buffer        = 5L,
    window_units  = "cells",
    asc_dir       = file.path(tempdir(), "sens-center-prep"),
    r_target      = habitat,
    r_mov         = affinity,
    r_src         = habitat,
    clear_dir     = TRUE,
    landmark      = 1L,
    progress      = FALSE
  )
  expect_identical(center_prep$target_mode, "center")

  expect_error(
    run_conscape(
      conscape_prep   = center_prep,
      out_dir         = file.path(tempdir(), "sens-center-run-should-error"),
      theta           = 0.1,
      distance_scale  = 3.8,
      jl_home         = "fake_path_never_reached",
      sensitivity     = conscape_sensitivity(wrt = "Q"),
      progress        = FALSE
    ),
    regexp = "target_mode = \"center\""
  )
})

test_that("tiled full-mean sensitivity converges to untiled as buffer grows", {
  skip_on_cran()
  skip_if_not(
    identical(Sys.getenv("RUN_CONSCAPERTOOLS_INTEGRATION", ""), "true"),
    "Set RUN_CONSCAPERTOOLS_INTEGRATION=true to run Julia-backed ConScape checks."
  )
  skip_if_not(
    identical(Sys.getenv("RUN_CONSCAPERTOOLS_SENSITIVITY", ""), "true"),
    paste0("Set RUN_CONSCAPERTOOLS_SENSITIVITY=true to run sensitivity ",
           "convergence checks (requires a ConScape build with ",
           "ConScape.sensitivity available; slow).")
  )

  jl_home <- integration_jl_home()
  skip_if_not(
    nzchar(jl_home),
    "Set CONSCAPERTOOLS_JL_HOME or JULIA_BINDIR to the Julia bin directory."
  )

  habitat_file  <- system.file("extdata", "suitability.asc", package = "ConScapeRtools")
  affinity_file <- system.file("extdata", "affinity.asc",  package = "ConScapeRtools")
  habitat       <- crop_cells_from_top_left(terra::rast(habitat_file),  n = 40L)
  affinity      <- crop_cells_from_top_left(terra::rast(affinity_file), n = 40L)

  expect_equal(dim(habitat)[1:2], c(40L, 40L))

  # Calibrated parameters from tile_design(max_d = 1200, theta = 0.1,
  # trim_threshold = 0.05, landmark = 5L) -- the same values used in the
  # main validation vignette. landmark forced to 1L because sensitivity
  # requires q_tgt = q_src per cell, and ConScape.sensitivity errors
  # otherwise.
  theta          <- 0.1
  distance_scale <- 3.8
  landmark       <- 1L
  centersize_c   <- 20L
  buffer_sweep   <- c(5L, 10L, 20L)   # 20 covers the 40-cell crop -> dgnrate
  sens_wrt       <- c("Q", "A&C=f(A)")
  sens_layers    <- c("elasticity_quality",
                      "elasticity_affinity_cost_linked")

  sens_spec <- conscape_sensitivity(
    wrt               = sens_wrt,
    landscape_measure = "sum",
    unitless          = TRUE
  )

  # ----- untiled reference -------------------------------------------------
  ref_out_dir <- file.path(tempdir(), "sens-convergence-untiled")
  untiled <- run_conscape(
    target_qualities = habitat,
    source_qualities = habitat,
    affinities       = affinity,
    out_dir          = ref_out_dir,
    theta            = theta,
    distance_scale   = distance_scale,
    landmark         = landmark,
    tile_trim        = 0,
    jl_home          = jl_home,
    sensitivity      = sens_spec,
    progress         = FALSE
  )
  untiled_stack <- result_as_stack(untiled, sens_layers)
  expect_false(is.null(untiled_stack))

  # ----- buffer sweep ------------------------------------------------------
  all_rows <- list()
  for (buf in buffer_sweep) {
    full_prep <- conscape_prep(
      centersize    = centersize_c,
      buffer        = buf,
      window_units  = "cells",
      asc_dir       = file.path(tempdir(),
                                paste0("sens-convergence-full-prep-", buf)),
      r_target      = habitat,
      r_mov         = affinity,
      r_src         = habitat,
      clear_dir     = TRUE,
      landmark      = landmark,
      progress      = FALSE,
      target_mode   = "full"   # sensitivity requires q_tgt == q_src per tile
    )
    full_run <- run_conscape(
      conscape_prep  = full_prep,
      out_dir        = file.path(tempdir(),
                                 paste0("sens-convergence-full-run-", buf)),
      theta          = theta,
      distance_scale = distance_scale,
      jl_home        = jl_home,
      sensitivity    = sens_spec,
      progress       = FALSE
    )
    full_stack <- result_as_stack(full_run, sens_layers)
    expect_false(is.null(full_stack))
    all_rows[[length(all_rows) + 1L]] <-
      stack_diff_stats(full_stack, untiled_stack,
                       layers = sens_layers, buffer = buf)
  }

  summary_df <- do.call(rbind, all_rows)
  attr(summary_df, "centersize")     <- centersize_c
  attr(summary_df, "distance_scale") <- distance_scale
  attr(summary_df, "theta")          <- theta

  message("\nSensitivity tiled (full mean) vs untiled (centersize = ",
          centersize_c, ", distance_scale = ", distance_scale,
          ", theta = ", theta, ")")
  message(paste(
    capture.output(print(summary_df, row.names = FALSE, digits = 4)),
    collapse = "\n"))

  # ----- assertions --------------------------------------------------------
  #
  # The assertions are deliberately MODEST. Unlike fcon / btwn under
  # center+sum, tiled sensitivity is NOT exact even in degenerate-coverage
  # (each tile computes its own landscape-summary derivative, which differs
  # from the global derivative because the per-tile landscape summary
  # itself differs). The contract this test enforces is monotone (or
  # near-monotone) convergence of the relative error as buffer grows.
  for (lyr in sens_layers) {
    rows <- subset(summary_df, layer == lyr)
    rows <- rows[order(rows$buffer), ]
    if (nrow(rows) < 2L) next

    # (1) relative max_abs should not get markedly worse as buffer grows
    # (allow 10% slack: sensitivity surfaces are noisier than fcon/btwn
    # and small landscapes have less room for monotone behavior).
    expect_lte(
      rows$rel_max_abs[nrow(rows)],
      rows$rel_max_abs[1L] * 1.10,
      label = paste0("rel_max_abs at max buffer for ", lyr,
                     " should not markedly exceed rel_max_abs at min buffer")
    )

    # (2) correlation with untiled at max buffer should be at least
    # moderately positive. Tile-local sensitivity surfaces can be quite
    # different from the global surface in shape, so we don't require
    # near-1 correlation -- only that the trend is positive enough to
    # confirm tile-local approximation is informative.
    expect_gt(
      rows$cor[nrow(rows)], 0.5,
      label = paste0("correlation with untiled at max buffer for ", lyr,
                     " should be moderately positive")
    )
  }

  invisible(summary_df)
})
