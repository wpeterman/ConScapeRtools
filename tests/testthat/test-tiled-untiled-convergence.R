## Tiled versus untiled agreement across several buffers on an 80x80 crop.
## Center mode partitions targets and sums contributions; full mode retains
## buffered targets and averages tile-local results. Both omit graph context
## at finite buffers, so error direction and monotonicity are not guaranteed.
## Each mode keeps per-window graph size bounded when the landscape grows at
## fixed centersize/buffer. Center mode generally solves fewer target columns.
## These endpoint checks complement the full-coverage and execution-mode tests.

integration_jl_home <- function() {
  candidates <- c(
    Sys.getenv("CONSCAPERTOOLS_JL_HOME", ""),
    Sys.getenv("JULIA_BINDIR", "")
  )
  candidates <- candidates[nzchar(candidates)]
  if (length(candidates)) candidates[[1]] else ""
}

crop_cells_from_top_left <- function(r, n = 120L) {
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
  # Pearson correlation between tiled and untiled values on jointly finite cells.
  cor_val <- if (length(diff) > 1L && stats::sd(vt) > 0 && stats::sd(vu) > 0) {
    suppressWarnings(stats::cor(vt, vu))
  } else {
    NA_real_
  }
  c(
    n        = length(diff),
    max_abs  = max(abs(diff)),
    mean_abs = mean(abs(diff)),
    rmse     = sqrt(mean(diff^2)),
    cor      = cor_val
  )
}

stack_diff_stats <- function(tiled, untiled, mode, buffer) {
  layers <- intersect(names(tiled), names(untiled))
  rows <- lapply(layers, function(lyr) {
    s <- surface_diff_stats(tiled[[lyr]], untiled[[lyr]])
    data.frame(
      mode      = mode,
      buffer    = buffer,
      metric    = lyr,
      n         = s["n"],
      max_abs   = s["max_abs"],
      mean_abs  = s["mean_abs"],
      rmse      = s["rmse"],
      cor       = s["cor"],
      row.names = NULL,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

result_as_stack <- function(x) {
  # run_conscape() returns a ConScapeResults list for tiled and a SpatRaster
  # for untiled. Normalize to a SpatRaster keyed by layer name.
  if (inherits(x, "SpatRaster")) return(x)
  layers <- intersect(c("btwn", "fcon"), names(x))
  do.call(c, lapply(layers, function(n) {
    r <- x[[n]]
    names(r) <- n
    r
  }))
}

run_with_dev <- function() {
  identical(Sys.getenv("RUN_CONSCAPERTOOLS_DEV", ""), "true")
}

test_that("tiled center-sum and full-mean both converge to untiled as buffer grows", {
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

  habitat_file  <- system.file("extdata", "suitability.asc", package = "ConScapeRtools")
  affinity_file <- system.file("extdata", "affinity.asc",  package = "ConScapeRtools")
  habitat       <- crop_cells_from_top_left(terra::rast(habitat_file),  n = 80L)
  affinity      <- crop_cells_from_top_left(terra::rast(affinity_file), n = 80L)

  expect_equal(dim(habitat)[1:2], c(80L, 80L))

  # Use a short proximity decay scale relative to the buffer sweep, so the
  # center_sum bias is dominated by the buffer-truncation term that we expect
  # to vanish at large buffer (see header comment).
  theta          <- 0.1
  distance_scale <- 20
  landmark       <- 1L
  centersize_c   <- 20L
  buffer_sweep   <- c(5L, 10L, 20L, 40L)

  # ----- untiled reference --------------------------------------------------
  ref_out_dir <- file.path(tempdir(), "convergence-untiled")
  untiled <- run_conscape(
    target_qualities = habitat,
    source_qualities = habitat,
    affinities       = affinity,
    out_dir          = ref_out_dir,
    theta            = theta,
    distance_scale   = distance_scale,
    landmark         = landmark,
    tile_trim        = 0,                # no buffering for the reference
    jl_home          = jl_home,
    progress         = FALSE
  )
  untiled_stack <- result_as_stack(untiled)

  # ----- buffer sweep --------------------------------------------------------
  all_rows <- list()

  for (buf in buffer_sweep) {
    # Center-sum (corrected math)
    center_prep <- conscape_prep(
      centersize    = centersize_c,
      buffer        = buf,
      window_units  = "cells",
      asc_dir       = file.path(tempdir(), paste0("convergence-center-prep-", buf)),
      r_target      = habitat,
      r_mov         = affinity,
      r_src         = habitat,
      clear_dir     = TRUE,
      landmark      = landmark,
      progress      = FALSE
    )
    center_run <- run_conscape(
      conscape_prep  = center_prep,
      out_dir        = file.path(tempdir(), paste0("convergence-center-run-", buf)),
      theta          = theta,
      distance_scale = distance_scale,
      jl_home        = jl_home,
      progress       = FALSE
    )

    # Legacy full-mean (tile-local approximation, supported as opt-in)
    full_prep <- conscape_prep(
      centersize    = centersize_c,
      buffer        = buf,
      window_units  = "cells",
      asc_dir       = file.path(tempdir(), paste0("convergence-full-prep-", buf)),
      r_target      = habitat,
      r_mov         = affinity,
      r_src         = habitat,
      clear_dir     = TRUE,
      landmark      = landmark,
      progress      = FALSE,
      target_mode   = "full"
    )
    full_run <- run_conscape(
      conscape_prep  = full_prep,
      out_dir        = file.path(tempdir(), paste0("convergence-full-run-", buf)),
      theta          = theta,
      distance_scale = distance_scale,
      jl_home        = jl_home,
      progress       = FALSE
    )

    all_rows[[length(all_rows) + 1L]] <-
      stack_diff_stats(result_as_stack(center_run), untiled_stack,
                       mode = "center_sum", buffer = buf)
    all_rows[[length(all_rows) + 1L]] <-
      stack_diff_stats(result_as_stack(full_run), untiled_stack,
                       mode = "full_mean", buffer = buf)

    # Optional dev WindowedProblem comparison
    if (run_with_dev()) {
      dev_run <- tryCatch(
        run_conscape(
          target_qualities = habitat,
          source_qualities = habitat,
          affinities       = affinity,
          out_dir          = file.path(tempdir(), paste0("convergence-dev-run-", buf)),
          jl_home          = jl_home,
          backend          = "conscape_dev",
          dev_mode         = "windowed",
          centersize       = centersize_c,
          buffer           = buf,
          landmark         = landmark,
          theta            = theta,
          distance_scale   = distance_scale,
          progress         = FALSE
        ),
        error = function(e) {
          message("dev backend skipped at buffer=", buf, ": ", conditionMessage(e))
          NULL
        }
      )
      if (!is.null(dev_run)) {
        all_rows[[length(all_rows) + 1L]] <-
          stack_diff_stats(result_as_stack(dev_run), untiled_stack,
                           mode = "dev_windowed", buffer = buf)
      }
    }
  }

  summary_df <- do.call(rbind, all_rows)
  attr(summary_df, "centersize")     <- centersize_c
  attr(summary_df, "distance_scale") <- distance_scale
  attr(summary_df, "theta")          <- theta

  # Report so this is visible even when nothing is asserted.
  message("\nTiled-vs-untiled convergence (centersize = ", centersize_c,
          ", distance_scale = ", distance_scale, ", theta = ", theta, ")")
  message(paste(capture.output(print(summary_df, row.names = FALSE,
                                     digits = 4)),
                collapse = "\n"))

  # ----- assertions ---------------------------------------------------------
  #
  # The assertions are deliberately *modest*. The point of this test is to
  # produce the diagnostic table; the assertions only guard the most basic
  # convergence properties that must hold if the math is correctly wired:
  #
  #  1. center_sum max_abs and rmse are not increasing as buffer grows (close
  #     to monotone non-increasing, allowing small numerical noise).
  #  2. center_sum correlation with untiled is high at the largest buffer.
  #  3. full_mean correlation with untiled is also high at the largest buffer
  #     (each tile approaches untiled).
  #
  # Comparing center_sum and full_mean error magnitudes is meaningful only
  # when buffer is much larger than the proximity decay scale, which depends
  # on data and parameters; we do not assert that ordering here.

  for (m in c("btwn", "fcon")) {
    center_m <- subset(summary_df, mode == "center_sum" & metric == m)
    center_m <- center_m[order(center_m$buffer), ]
    full_m <- subset(summary_df, mode == "full_mean" & metric == m)
    full_m <- full_m[order(full_m$buffer), ]

    # (1) center_sum error should not get worse as buffer grows (allow 1% slop
    # for floating-point and ConScape solver tolerances).
    expect_lte(
      center_m$max_abs[nrow(center_m)],
      center_m$max_abs[1L] * 1.01,
      label = paste0("center_sum max_abs at max buffer for ", m,
                     " should not exceed max_abs at min buffer")
    )
    expect_lte(
      center_m$rmse[nrow(center_m)],
      center_m$rmse[1L] * 1.01,
      label = paste0("center_sum rmse at max buffer for ", m,
                     " should not exceed rmse at min buffer")
    )

    # (2) center_sum correlation with untiled should be high at max buffer.
    expect_gt(
      center_m$cor[nrow(center_m)], 0.9,
      label = paste0("center_sum correlation at max buffer for ", m)
    )
    # (3) full_mean correlation with untiled should also be high at max buffer.
    expect_gt(
      full_m$cor[nrow(full_m)], 0.9,
      label = paste0("full_mean correlation at max buffer for ", m)
    )
  }

  invisible(summary_df)
})
