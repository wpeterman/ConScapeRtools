#' Compare RAM and throughput across alternative tile sizes
#'
#' @description
#' After you have a calibrated tile or window design from [tile_design()] or
#' [window_design()], you usually want to know what changing `centersize`
#' would buy you: how matrix storage changes and how many more
#' (or fewer) windows you would solve. These proxies identify candidates
#' for benchmarking; they do not establish runtime, peak RAM, or accuracy.
#' `window_design_scenarios()` sweeps a vector of candidate `centersize`
#' values while holding `buffer`, `landmark`, and the landscape fixed, and
#' returns a one-row-per-scenario summary of those tradeoffs.
#'
#' @param r A `SpatRaster` describing the landscape's cell grid. Only its
#'   dimensions, resolution, and `centersize`/`buffer` window arithmetic
#'   are used. Optional when `design` carries enough geometry, but you
#'   typically want to pass the analysis raster so window counts match
#'   what [conscape_prep()] would produce.
#' @param design Optional object returned by [tile_design()] or
#'   [window_design()]. Used to default the fixed-across-scenarios
#'   parameters (`buffer`, `landmark`) and the reference `centersize` from
#'   which the default candidate sweep is built.
#' @param centersize Optional numeric vector of candidate center window
#'   sizes in cells. When `NULL` (default), a geometric sweep
#'   `c(1, 1/2, 1/4, 1/8) * design$centersize` is built, snapped to
#'   multiples of `landmark` so coarse-grained target cells stay aligned.
#'   Explicit candidates are rounded up to multiples of `landmark`.
#'   Candidates smaller than `landmark` or larger than the raster are
#'   dropped.
#' @param buffer Fixed buffer width in cells (the biological-calibration
#'   parameter). When `NULL` (default), taken from `design`. Rounded up to
#'   a multiple of `landmark`, with a minimum of `landmark` when it exceeds 1.
#' @param landmark Fixed coarse-graining window size in cells. When `NULL`
#'   (default), taken from `design`, then `10L`.
#' @param workers Integer number of parallel workers to model.
#'   `total_ram_gb` in the output is the combined single-matrix proxy,
#'   `workers * ram_per_window_gb`. Default `1L`.
#' @param target_ram_gb Optional numeric. If supplied, the output gains a
#'   logical column `recommended` flagging the **largest** `centersize`
#'   whose `total_ram_gb` is at or below the target. Larger centers are
#'   preferred at a given RAM budget because they minimize dispatch
#'   overhead. A flagged row is a candidate for benchmarking, not a guarantee
#'   that the full run fits in RAM.
#'
#' @details
#' The memory proxy counts one dense double-precision matrix of size
#' `source_window_cells * target_cells_after_landmark`. It is not a peak-RAM
#' estimate or upper bound. ConScape also needs sparse factorizations and
#' several dense work arrays; parallel processes have separate runtime costs.
#' Benchmark a representative tile and leave memory for these allocations.
#' The default expected-cost metrics now solve target columns in blocks of
#' at most 128. The reported full-matrix proxy is retained for comparison and
#' for other metrics, not a prediction of the default streaming solver's RAM.
#'
#' Smaller centers reduce the number of target columns per solve but increase
#' the number of graph factorizations. The workload proxy does not model this
#' tradeoff fully and cannot predict wall-clock time.
#'
#' The target metric is unchanged by centersize, but its tiled approximation
#' can change: different centers define different finite graph boundaries.
#' Check convergence after changing centersize as well as buffer. With
#' `landmark > 1L`, coarse graining is another approximation, especially near
#' missing cells and partial edge blocks. Landmark-aligned layouts preserve
#' nominal grid positions, not guaranteed agreement with an uncoarsened run.
#'
#' Counts describe classic, externally padded tiles created by
#' [conscape_prep()], before empty-target tiles are removed. Experimental dev
#' windows use different edge handling; use [window_design()] for those counts.
#'
#' @return
#' A data frame of class `"ConScapeRtools_scenarios"` with one row per
#' candidate `centersize`. Columns:
#'
#' * `centersize`, `buffer`, `landmark`: design parameters in cells.
#' * `source_window_cells`: `(centersize + 2*buffer)^2`, the per-window
#'   graph size.
#' * `center_cells`: `centersize^2`, retained as RSP targets.
#' * `target_cells_after_landmark`: `center_cells / landmark^2`, the
#'   potential number of `Z` columns before masking and edge effects.
#' * `expected_windows`: approximate tile count over `r`.
#' * `ram_per_window_gb`: storage for one dense source-by-target matrix in
#'   decimal GB, assuming a fully populated window. Excludes solver workspaces.
#' * `total_ram_gb`: `workers * ram_per_window_gb`, the system-wide
#'   single-matrix storage when running `workers` solves in parallel.
#' * `overlap_area_factor`: `source_window_cells / center_cells`. A
#'   measure of buffer overhead per retained center cell. Values near
#'   2 are efficient; large values mean most per-tile work is repeated
#'   buffer context.
#' * `workload_proxy`: `expected_windows * source_window_cells *
#'   center_cells`. A relative cell-pair count, not elapsed time.
#' * `recommended`: logical, only when `target_ram_gb` is supplied.
#'   `TRUE` for the largest `centersize` whose `total_ram_gb` is at or
#'   below the target.
#'
#' The companion `print` method adds a short header explaining the
#' columns and highlights the recommended scenario when one is flagged.
#'
#' @examples
#' library(ConScapeRtools)
#' r <- terra::rast(system.file("extdata", "suitability.asc",
#'                              package = "ConScapeRtools"))
#'
#' # Start from a known window design (here: a small synthetic example so
#' # the call does not need Julia).
#' wd <- window_design(
#'   r        = r,
#'   tile_d   = 5400,
#'   tile_trim = 1350,
#'   landmark = 5L
#' )
#'
#' # Default geometric sweep around wd$centersize.
#' scenarios <- window_design_scenarios(r = r, design = wd)
#' scenarios
#'
#' # Custom candidates and a target RAM budget of 2 GB per system.
#' window_design_scenarios(
#'   r            = r,
#'   design       = wd,
#'   centersize   = c(40, 20, 10, 5),
#'   target_ram_gb = 2,
#'   workers      = 4L
#' )
#'
#' @seealso [tile_design()], [window_design()], [conscape_prep()],
#'   [run_conscape()]
#' @author Bill Peterman
#' @export
window_design_scenarios <- function(r = NULL,
                                    design = NULL,
                                    centersize = NULL,
                                    buffer = NULL,
                                    landmark = NULL,
                                    workers = 1L,
                                    target_ram_gb = NULL) {

  if (is.null(r) && is.null(design)) {
    stop("Supply at least one of `r` or `design`.", call. = FALSE)
  }
  if (!is.null(r) && !inherits(r, "SpatRaster")) {
    stop("`r` must be a SpatRaster.", call. = FALSE)
  }

  ## ---- defaults from design -----------------------------------------------
  if (!is.null(design)) {
    if (is.null(buffer)     && !is.null(design$buffer))     buffer     <- design$buffer
    if (is.null(landmark)   && !is.null(design$landmark))   landmark   <- design$landmark
    if (is.null(centersize) && !is.null(design$centersize)) {
      ref_centersize <- design$centersize
    } else {
      ref_centersize <- NULL
    }
  } else {
    ref_centersize <- NULL
  }
  if (is.null(landmark)) landmark <- 10L
  validate_positive_scalar(landmark, "landmark")
  if (landmark != floor(landmark) || landmark > .Machine$integer.max) {
    stop("landmark must be a positive integer.", call. = FALSE)
  }
  landmark <- as.integer(landmark)

  if (is.null(buffer)) {
    stop("`buffer` is required when `design` does not provide it.",
         call. = FALSE)
  }
  validate_nonnegative_scalar(buffer, "buffer")
  buffer <- as.integer(ceiling(max(buffer, floor(landmark / 2)) / landmark) * landmark)
  if (is.na(buffer)) stop("buffer is too large.", call. = FALSE)

  if (!is.numeric(workers) || length(workers) != 1L ||
      is.na(workers) || !is.finite(workers) || workers < 1 ||
      workers != floor(workers) || workers > .Machine$integer.max) {
    stop("`workers` must be a positive integer.", call. = FALSE)
  }
  workers <- as.integer(workers)

  if (!is.null(target_ram_gb)) {
    if (!is.numeric(target_ram_gb) || length(target_ram_gb) != 1L ||
        is.na(target_ram_gb) || !is.finite(target_ram_gb) || target_ram_gb <= 0) {
      stop("`target_ram_gb` must be NULL or a single positive numeric.",
           call. = FALSE)
    }
  }

  ## ---- candidate centersizes ----------------------------------------------
  if (is.null(centersize)) {
    if (is.null(ref_centersize)) {
      stop("Supply `centersize` (vector of candidates) or pass a `design` ",
           "object whose `centersize` defines the sweep.", call. = FALSE)
    }
    raw <- ref_centersize * c(1, 1 / 2, 1 / 4, 1 / 8)
    # Snap to landmark multiples so coarse-graining target placement stays
    # aligned, deduplicate, drop anything smaller than the landmark, and
    # keep at most 8 candidates sorted descending.
    snapped <- pmax(landmark, round(raw / landmark) * landmark)
    centersize <- sort(unique(as.integer(snapped)), decreasing = TRUE)
  } else {
    if (!is.numeric(centersize) || any(is.na(centersize)) ||
        !length(centersize) || any(!is.finite(centersize)) || any(centersize <= 0) ||
        any(centersize > .Machine$integer.max)) {
      stop("`centersize` must be a positive numeric vector.", call. = FALSE)
    }
    centersize <- sort(unique(as.integer(ceiling(centersize / landmark) * landmark)), decreasing = TRUE)
  }

  ## ---- per-candidate diagnostics ------------------------------------------
  rows <- lapply(centersize, function(cs) {
    if (!is.null(r)) {
      ncols <- terra::ncol(r)
      nrows <- terra::nrow(r)
      expected_windows <- ceiling(ncols / cs) * ceiling(nrows / cs)
    } else if (!is.null(design$diagnostics$expected_windows) &&
               !is.null(design$centersize)) {
      # Fall back to scaling the reference window count by area ratio.
      ref_n <- design$diagnostics$expected_windows
      expected_windows <- as.integer(round(ref_n * (design$centersize / cs)^2))
    } else {
      expected_windows <- NA_integer_
    }

    source_window_cells <- (as.numeric(cs) + 2 * buffer)^2
    center_cells <- as.numeric(cs)^2
    target_after_lm <- center_cells / as.numeric(landmark)^2
    overlap_area_factor <- source_window_cells / center_cells

    # 8 bytes per double; dominant cost is the Z block of size S x K.
    ram_bytes <- 8 * source_window_cells * target_after_lm
    ram_gb <- ram_bytes / 1e9

    workload_proxy <- if (is.na(expected_windows)) {
      NA_real_
    } else {
      as.numeric(expected_windows) * source_window_cells * center_cells
    }

    data.frame(
      centersize                  = as.integer(cs),
      buffer                      = buffer,
      landmark                    = landmark,
      source_window_cells         = source_window_cells,
      center_cells                = center_cells,
      target_cells_after_landmark = target_after_lm,
      expected_windows            = expected_windows,
      ram_per_window_gb           = ram_gb,
      total_ram_gb                = ram_gb * workers,
      overlap_area_factor         = overlap_area_factor,
      workload_proxy              = workload_proxy,
      stringsAsFactors            = FALSE,
      row.names                   = NULL
    )
  })
  out <- do.call(rbind, rows)
  # Drop scenarios where the candidate centersize is larger than the raster
  # in either dimension (the window arithmetic still reports 1 tile, but
  # it represents an untiled run, not a tiling alternative) or where no
  # window fits at all.
  if (!is.null(r)) {
    max_dim <- min(terra::ncol(r), terra::nrow(r))
    too_big <- out$centersize > max_dim
    out <- out[!too_big & out$expected_windows >= 1L, , drop = FALSE]
  }
  if (nrow(out) == 0L) {
    stop("No candidate `centersize` fits on `r`. Candidates either exceed ",
         "the raster's smallest dimension or leave no room for the buffer. ",
         "Try smaller candidates or a smaller `buffer`.", call. = FALSE)
  }

  ## ---- target RAM recommendation ------------------------------------------
  if (!is.null(target_ram_gb)) {
    fits <- out$total_ram_gb <= target_ram_gb
    recommended <- logical(nrow(out))
    if (any(fits)) {
      # Pick the LARGEST centersize that fits (minimizes dispatch overhead).
      recommended[which(fits)[which.max(out$centersize[fits])]] <- TRUE
    }
    out$recommended <- recommended
    attr(out, "target_ram_gb") <- target_ram_gb
  }

  attr(out, "workers") <- workers
  rownames(out) <- NULL
  class(out) <- c("ConScapeRtools_scenarios", class(out))
  out
}

#' @rdname conscape_s3_methods
#' @method print ConScapeRtools_scenarios
#' @export
print.ConScapeRtools_scenarios <- function(x, ...) {
  workers <- attr(x, "workers")
  target  <- attr(x, "target_ram_gb")
  cat("ConScapeRtools centersize scenarios\n")
  cat("  workers modeled: ", workers, "\n", sep = "")
  if (!is.null(target)) {
    cat("  target system RAM (GB): ", target, "\n", sep = "")
  }
  cat("  ram_per_window_gb counts one matrix, not peak RAM or an upper bound ",
      "(8 * source_window_cells * center_cells / landmark^2 / 1e9)\n", sep = "")
  cat("  workload_proxy counts cell pairs; benchmark actual runtime and peak RAM\n\n")

  display <- as.data.frame(unclass(x))
  display$source_window_cells <- format(display$source_window_cells,
                                        big.mark = ",", scientific = FALSE)
  display$center_cells <- format(display$center_cells,
                                 big.mark = ",", scientific = FALSE)
  display$target_cells_after_landmark <-
    format(round(display$target_cells_after_landmark),
           big.mark = ",", scientific = FALSE)
  display$expected_windows <- format(display$expected_windows,
                                     big.mark = ",", scientific = FALSE)
  display$ram_per_window_gb <- signif(display$ram_per_window_gb, 3)
  display$total_ram_gb      <- signif(display$total_ram_gb, 3)
  display$overlap_area_factor <- round(display$overlap_area_factor, 2)
  display$workload_proxy    <- signif(display$workload_proxy, 3)
  print.data.frame(display, row.names = FALSE)
  if (!is.null(target) && any(x$recommended)) {
    rec <- x$centersize[x$recommended]
    cat("\n  centersize candidate under the matrix-storage budget ", target, " GB: ",
        rec, "; benchmark before committing to a full run\n", sep = "")
  } else if (!is.null(target)) {
    cat("\n  no candidate centersize fits within ", target, " GB on ",
        workers, " workers; consider smaller centers, fewer workers, ",
        "or a larger budget. Increasing landmark changes the approximation.\n", sep = "")
  }
  invisible(x)
}
