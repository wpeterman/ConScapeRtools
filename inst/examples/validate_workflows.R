## ============================================================================
## Identity validation: tiled vs untiled ConScape
## ============================================================================
##
## This script demonstrates that the corrected center-buffer tiled workflow in
## ConScapeRtools produces a numerically identical surface to a single untiled
## ConScape run, when the tile buffer is large enough to cover the relevant
## landscape context. It also exercises the legacy full-target tiled workflow
## and (optionally) the experimental ConScape dev WindowedProblem backend.
##
## Workflows compared on the same demo crop:
##
##   1. untiled           - run_conscape() on the whole raster as a single graph
##   2. classic_tiled     - conscape_prep(tile_d, tile_trim) + run_conscape()
##                          (target_mode = "center", sum mosaic)
##   3. center_buffer     - conscape_prep(centersize, buffer) + run_conscape()
##                          (cell-unit center-buffer prep, otherwise identical
##                          to classic_tiled)
##   4. full_legacy       - conscape_prep(target_mode = "full") + run_conscape()
##                          (biased mean-mosaic estimator; opt-in for legacy)
##   5. dev_windowed      - run_conscape(backend = "conscape_dev",
##                                       dev_mode = "windowed") [optional]
##   6. dev_batch         - run_conscape(backend = "conscape_dev",
##                                       dev_mode = "batch")    [optional]
##
## The "identity" claim is precise: when every tile's buffered window covers
## the full landscape (degenerate-coverage), the per-tile W matrix and Z
## columns are mathematically identical to the untiled solution, and the
## center-sum reduction sums disjoint center contributions back into the full
## untiled fcon/btwn surfaces. Differences in this regime are floating-point
## noise from ConScape's dense LU solver, bounded at ~1e-6 relative to the
## untiled surface magnitude.
##
## When the buffered window does NOT cover the full landscape, tiled outputs
## CONVERGE to untiled as buffer grows but do not equal it exactly. That
## convergence is exercised by tests/testthat/test-tiled-untiled-convergence.R.
##
## Usage:
##   Sys.setenv(CONSCAPERTOOLS_JL_HOME = "C:/path/to/julia/bin")
##   source("inst/examples/validate_workflows.R")
##   results <- validate_workflows()
##   print(results$summary)   # identity table (max_abs, rmse, cor, identical_fp)
##   print(results$timings)   # wall-clock timing table (seconds, speedup)
##
## Optional dev-backend run (also runs dev_windowed and dev_batch):
##   results <- validate_workflows(include_dev = TRUE)
##
## Re-using an existing ConScape dev project (skips setup timing):
##   dev_project <- conscape_dev_backend_setup(jl_home = jl_home,
##                                             rev = "alg_efficiency")
##   results <- validate_workflows(include_dev = TRUE,
##                                 dev_project = dev_project)
## ============================================================================

library(ConScapeRtools)
library(terra)

## -- helpers -----------------------------------------------------------------

crop_top_left <- function(r, n) {
  x0 <- terra::xmin(r); y1 <- terra::ymax(r)
  terra::crop(
    r,
    terra::ext(
      x0, x0 + n * terra::res(r)[1],
      y1 - n * terra::res(r)[2], y1
    ),
    snap = "near"
  )
}

as_stack <- function(x) {
  if (inherits(x, "SpatRaster")) return(x)
  layers <- intersect(c("btwn", "fcon"), names(x))
  do.call(c, lapply(layers, function(n) {
    r <- x[[n]]; names(r) <- n; r
  }))
}

agreement <- function(candidate, reference, label) {
  layers <- intersect(names(candidate), names(reference))
  do.call(rbind, lapply(layers, function(lyr) {
    diff <- terra::values(candidate[[lyr]] - reference[[lyr]], mat = FALSE)
    cand <- terra::values(candidate[[lyr]], mat = FALSE)
    refv <- terra::values(reference[[lyr]], mat = FALSE)
    ok <- is.finite(diff) & is.finite(cand) & is.finite(refv)
    diff <- diff[ok]; cand <- cand[ok]; refv <- refv[ok]
    ref_scale <- max(abs(refv))
    cor_val <- if (length(diff) > 1L && sd(cand) > 0 && sd(refv) > 0) {
      suppressWarnings(cor(cand, refv))
    } else NA_real_
    data.frame(
      workflow      = label,
      layer         = lyr,
      n             = length(diff),
      max_abs       = max(abs(diff)),
      mean_abs      = mean(abs(diff)),
      rmse          = sqrt(mean(diff^2)),
      cor           = cor_val,
      rel_max_abs   = max(abs(diff)) / pmax(ref_scale, 1),
      identical_fp  = max(abs(diff)) / pmax(ref_scale, 1) < 1e-6,
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }))
}

## -- main --------------------------------------------------------------------

## Time a single workflow invocation. Returns list(result, seconds, error).
## Errors are caught so one failing workflow does not abort the harness.
## NB: do not attach attributes to a NULL result -- R disallows this and the
## attempt itself raises. Carry the error message in the returned list slot.
time_run <- function(label, expr) {
  message("  [time] ", label, " ...")
  t0 <- proc.time()
  err <- NULL
  res <- tryCatch(force(expr),
                  error = function(e) {
                    err <<- conditionMessage(e)
                    message("    ", label, " failed: ", err)
                    NULL
                  })
  dt <- as.numeric((proc.time() - t0)["elapsed"])
  if (is.null(err)) {
    message("    ", label, " elapsed: ", sprintf("%.2f s", dt))
  } else {
    message("    ", label, " elapsed before error: ", sprintf("%.2f s", dt))
  }
  list(result = res, seconds = dt, error = err)
}

validate_workflows <- function(jl_home    = Sys.getenv("CONSCAPERTOOLS_JL_HOME"),
                               n_cells    = 40L,
                               centersize = 20L,
                               buffer     = 20L,
                               theta      = 0.15,
                               distance_scale = 150,
                               landmark   = 5L,
                               include_dev = FALSE,
                               dev_project = NULL,
                               dev_conscape_rev = "alg_efficiency",
                               work_dir   = tempfile("validate_workflows_")) {
  stopifnot(nzchar(jl_home))
  dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)

  hab_file <- system.file("extdata", "suitability.asc", package = "ConScapeRtools")
  aff_file <- system.file("extdata", "affinity.asc",  package = "ConScapeRtools")
  habitat  <- crop_top_left(terra::rast(hab_file),  n_cells)
  affinity <- crop_top_left(terra::rast(aff_file), n_cells)

  cell_size  <- terra::res(habitat)[1]
  tile_d_map <- centersize * cell_size
  tile_trim_map <- buffer    * cell_size

  message("Demo crop: ", n_cells, " x ", n_cells, " cells (",
          n_cells * cell_size, " x ", n_cells * cell_size, " map units)")
  message("centersize = ", centersize, " cells = ", tile_d_map, " map units")
  message("buffer     = ", buffer,    " cells = ", tile_trim_map, " map units")

  # Timing convention: each workflow's "seconds" is wall-clock elapsed time for
  # everything that workflow does end-to-end on the same demo crop, measured
  # with proc.time(). For workflows that use conscape_prep(), the prep step is
  # included so the timing reflects the user-visible cost (tile_design() is
  # excluded - it is shared design work, not per-workflow). For the dev
  # backend, conscape_dev_backend_setup() is excluded if a dev_project is
  # supplied, otherwise it is included once (the first dev call pays for it).
  timings <- list()
  workflow_errors <- list()

  ## (1) Untiled reference -----------------------------------------------------
  message("\n[1/5] Running untiled reference ...")
  t_untiled <- time_run("untiled",
    run_conscape(
      target_qualities = habitat,
      source_qualities = habitat,
      affinities       = affinity,
      out_dir          = file.path(work_dir, "untiled"),
      theta            = theta,
      distance_scale   = distance_scale,
      landmark         = landmark,
      tile_trim        = tile_trim_map,
      jl_home          = jl_home,
      progress         = FALSE
    )
  )
  untiled <- t_untiled$result
  untiled_s <- as_stack(untiled)
  timings$untiled <- t_untiled$seconds

  ## (2) Classic tiled (map-unit) ---------------------------------------------
  message("[2/5] Running classic tiled (tile_d / tile_trim) ...")
  t_classic <- time_run("classic_tiled", {
    classic_prep <- conscape_prep(
      tile_d    = tile_d_map,
      tile_trim = tile_trim_map,
      asc_dir   = file.path(work_dir, "classic-prep"),
      r_target  = habitat, r_src = habitat, r_mov = affinity,
      clear_dir = TRUE, landmark = landmark, progress = FALSE
    )
    run_conscape(
      conscape_prep  = classic_prep,
      out_dir        = file.path(work_dir, "classic-run"),
      theta          = theta, distance_scale = distance_scale,
      jl_home        = jl_home, progress = FALSE
    )
  })
  classic <- t_classic$result
  timings$classic_tiled <- t_classic$seconds

  ## (3) Center-buffer (cell-unit) --------------------------------------------
  message("[3/5] Running center-buffer (centersize / buffer) ...")
  t_center <- time_run("center_buffer", {
    center_prep <- conscape_prep(
      centersize    = centersize, buffer = buffer, window_units = "cells",
      asc_dir       = file.path(work_dir, "center-prep"),
      r_target      = habitat, r_src = habitat, r_mov = affinity,
      clear_dir     = TRUE, landmark = landmark, progress = FALSE
    )
    run_conscape(
      conscape_prep  = center_prep,
      out_dir        = file.path(work_dir, "center-run"),
      theta          = theta, distance_scale = distance_scale,
      jl_home        = jl_home, progress = FALSE
    )
  })
  center <- t_center$result
  timings$center_buffer <- t_center$seconds

  ## (4) Legacy full-target (mean mosaic) -------------------------------------
  message("[4/5] Running legacy full-target (target_mode = 'full') ...")
  t_full <- time_run("full_legacy", {
    full_prep <- conscape_prep(
      centersize    = centersize, buffer = buffer, window_units = "cells",
      asc_dir       = file.path(work_dir, "full-prep"),
      r_target      = habitat, r_src = habitat, r_mov = affinity,
      clear_dir     = TRUE, landmark = landmark, progress = FALSE,
      target_mode   = "full"
    )
    run_conscape(
      conscape_prep  = full_prep,
      out_dir        = file.path(work_dir, "full-run"),
      theta          = theta, distance_scale = distance_scale,
      jl_home        = jl_home, progress = FALSE
    )
  })
  full <- t_full$result
  timings$full_legacy <- t_full$seconds

  ## (5) Dev windowed/batch (optional) ---------------------------------------
  dev_w <- NULL; dev_b <- NULL
  dev_setup_seconds <- NA_real_
  if (isTRUE(include_dev)) {
    if (is.null(dev_project)) {
      message("[5/5] Setting up ConScape dev project (one-off) ...")
      t_setup <- time_run("conscape_dev_setup",
        conscape_dev_backend_setup(
          jl_home = jl_home,
          rev = dev_conscape_rev
        )
      )
      dev_project <- t_setup$result
      dev_setup_seconds <- t_setup$seconds
    }

    if (!is.null(dev_project)) {
      message("[5/5] Running ConScape dev WindowedProblem ...")
      t_dev_w <- time_run("dev_windowed",
        run_conscape(
          target_qualities = habitat, source_qualities = habitat,
          affinities = affinity,
          out_dir = file.path(work_dir, "dev-windowed"),
          jl_home = jl_home, backend = "conscape_dev",
          dev_mode = "windowed",
          centersize = centersize, buffer = buffer,
          dev_project = dev_project,
          landmark = landmark, theta = theta, distance_scale = distance_scale,
          progress = FALSE
        )
      )
      dev_w <- t_dev_w$result
      timings$dev_windowed <- t_dev_w$seconds
      if (!is.null(t_dev_w$error)) workflow_errors$dev_windowed <- t_dev_w$error

      message("    Running ConScape dev BatchProblem ...")
      t_dev_b <- time_run("dev_batch",
        run_conscape(
          target_qualities = habitat, source_qualities = habitat,
          affinities = affinity,
          out_dir = file.path(work_dir, "dev-batch"),
          jl_home = jl_home, backend = "conscape_dev",
          dev_mode = "batch",
          centersize = centersize, buffer = buffer,
          dev_project = dev_project,
          landmark = landmark, theta = theta, distance_scale = distance_scale,
          progress = FALSE
        )
      )
      dev_b <- t_dev_b$result
      timings$dev_batch <- t_dev_b$seconds
      if (!is.null(t_dev_b$error)) workflow_errors$dev_batch <- t_dev_b$error
    }
  } else {
    message("[5/5] Skipping dev backend (set include_dev = TRUE to run)")
  }

  ## -- compare ---------------------------------------------------------------
  rows <- list(
    agreement(as_stack(classic),  untiled_s, "classic_tiled"),
    agreement(as_stack(center),   untiled_s, "center_buffer"),
    agreement(as_stack(full),     untiled_s, "full_legacy")
  )
  if (!is.null(dev_w)) rows[[length(rows) + 1L]] <-
    agreement(as_stack(dev_w), untiled_s, "dev_windowed")
  if (!is.null(dev_b)) rows[[length(rows) + 1L]] <-
    agreement(as_stack(dev_b), untiled_s, "dev_batch")

  summary <- do.call(rbind, rows)
  attr(summary, "n_cells")        <- n_cells
  attr(summary, "centersize")     <- centersize
  attr(summary, "buffer")         <- buffer
  attr(summary, "theta")          <- theta
  attr(summary, "distance_scale") <- distance_scale

  ## Timing summary (wall-clock seconds). Untiled is the reference: speedup
  ## = untiled_seconds / workflow_seconds; >1 means the workflow is faster.
  untiled_seconds <- timings$untiled
  timing_workflows <- names(timings)
  timing_seconds <- as.numeric(unlist(timings))
  # Status: "ok" for normal runs, "failed" if the timed call returned NULL.
  status_for <- function(name) {
    if (name %in% names(workflow_errors) && !is.null(workflow_errors[[name]])) {
      "failed"
    } else {
      "ok"
    }
  }
  timing_df <- data.frame(
    workflow = timing_workflows,
    seconds  = timing_seconds,
    speedup_vs_untiled = round(untiled_seconds / timing_seconds, 2),
    status = vapply(timing_workflows, status_for, character(1)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  attr(timing_df, "dev_setup_seconds") <- dev_setup_seconds
  attr(timing_df, "errors") <- workflow_errors

  cat("\n----- Identity validation: tiled vs untiled -----\n")
  cat(sprintf("Demo: %dx%d cells, centersize=%d, buffer=%d, theta=%g, distance_scale=%g\n",
              n_cells, n_cells, centersize, buffer, theta, distance_scale))
  cat("identical_fp = TRUE when max_abs / max(|untiled|) < 1e-6 (floating-point tolerance)\n\n")
  print(summary, row.names = FALSE, digits = 4)

  cat("\n----- Wall-clock timing (proc.time elapsed seconds) -----\n")
  cat("speedup_vs_untiled = untiled_seconds / workflow_seconds\n")
  cat("(> 1 means workflow is faster than untiled on this crop)\n")
  cat("status = 'failed' means the workflow errored; see attr(timings,\"errors\")\n\n")
  print(timing_df, row.names = FALSE, digits = 3)
  if (is.finite(dev_setup_seconds)) {
    cat(sprintf(
      "\n(dev backend one-off setup: %.2f s, not included in dev workflow timings)\n",
      dev_setup_seconds))
  }
  if (length(workflow_errors)) {
    cat("\nWorkflow errors:\n")
    for (nm in names(workflow_errors)) {
      cat("  ", nm, ": ", workflow_errors[[nm]], "\n", sep = "")
    }
  }

  invisible(list(
    summary  = summary,
    timings  = timing_df,
    untiled  = untiled_s,
    classic  = as_stack(classic),
    center   = as_stack(center),
    full     = as_stack(full),
    dev_windowed = if (!is.null(dev_w)) as_stack(dev_w),
    dev_batch    = if (!is.null(dev_b)) as_stack(dev_b),
    dev_project  = dev_project,
    work_dir = work_dir
  ))
}
