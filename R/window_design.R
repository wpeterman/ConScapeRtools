#' Translate tile design into WindowedProblem and BatchProblem parameters
#'
#' @description
#' Converts a map-unit ConScapeRtools tile design into the cell-unit
#' `centersize` and `buffer` parameters used by the experimental
#' `backend = "conscape_dev"` workflows. Use [tile_design()] first when the
#' goal is to choose biologically meaningful distances from a species'
#' dispersal ability. Use `window_design()` after that to make the same design
#' usable by ConScape development `WindowedProblem` and `BatchProblem` objects.
#'
#' @param r `SpatRaster` used to determine cell size, raster dimensions, and
#'   approximate job counts. This should be the analysis raster or a raster
#'   with the same resolution and extent.
#' @param design Optional object returned by [tile_design()]. If supplied,
#'   `tile_d`, `tile_trim`, `landmark`, `theta`, and `distance_scale` are read
#'   from this object unless explicitly overridden.
#' @param tile_d Interior tile or center-window width in map units. Optional
#'   when `design` is supplied or when `centersize` and `buffer` are supplied.
#' @param tile_trim Buffer or trim distance in map units. Optional when
#'   `design` is supplied or when `centersize` and `buffer` are supplied.
#' @param landmark Integer coarse-graining window size used for landmark
#'   alignment. If `NULL`, the value is taken from `design` when available,
#'   otherwise `10L`.
#' @param centersize Optional center window size. If supplied, `buffer` must
#'   also be supplied. By default this is interpreted in cells.
#' @param buffer Optional buffer width around each center window. If supplied,
#'   `centersize` must also be supplied. By default this is interpreted in
#'   cells.
#' @param window_units Units for explicitly supplied `centersize` and `buffer`:
#'   `"cells"` or `"map"`.
#'
#' @details
#' ConScape development windows separate two ideas that are combined in the
#' classic tiled workflow:
#'
#' * `centersize` is the retained target window. These are the cells whose
#'   outputs are intended to be kept from a given window.
#' * `buffer` is the contextual margin around the center. Source qualities,
#'   target qualities, and affinities in this margin are used in the
#'   calculation so the center cells see nearby landscape context.
#'
#' The total source/context window has width `centersize + 2 * buffer` cells.
#' For randomized shortest-path connectivity, a useful first-order workload
#' proxy is:
#'
#' `source_window_cells * center_cells * expected_windows`.
#'
#' This is not an exact runtime formula, but it tracks the key scaling pressure:
#' more source cells and more target cells make each window more expensive, and
#' smaller center windows increase the number of windows.
#'
#' When `design` or map-unit `tile_d` and `tile_trim` are supplied,
#' `window_design()` calls the same internal landmark-alignment logic used by
#' [conscape_prep()]. This makes the returned `centersize` and `buffer`
#' consistent with the classic tiled workflow:
#'
#' * `centersize = tile_cells`
#' * `buffer = overlap_cells`
#' * `tile_d = centersize * cell size`
#' * `tile_trim = buffer * cell size`
#'
#' Practical guidance:
#'
#' * Choose `buffer` from biology first. In practice this means using
#'   [tile_design()] with a defensible `max_d` and `trim_threshold`.
#' * Choose `centersize` from compute constraints second. Larger centers
#'   reduce repeated overlap work but increase per-window memory and solve
#'   size. Smaller centers are safer for memory but repeat more buffered
#'   context.
#' * Inspect `overlap_area_factor`. Values near 2 to 4 are often reasonable
#'   starting points. Very large values mean the buffer dominates the center,
#'   so most compute is repeated context.
#' * Use `dev_mode = "windowed"` first to test correctness and API behavior.
#'   Use `dev_mode = "batch"` when the same design should be scheduled as
#'   resumable on-disk jobs.
#'
#' @return
#' A named list of class `"ConScapeRtools_window_design"` with:
#'
#' * `centersize` and `buffer` in cells for `run_conscape()`.
#' * `tile_d` and `tile_trim` in map units for [conscape_prep()].
#' * `run_args`, a small list suitable for `run_conscape(...,
#'   backend = "conscape_dev")`.
#' * `prep_args`, a small list suitable for [conscape_prep()] center-buffer
#'   preparation.
#' * `diagnostics`, including expected window count, source/context window
#'   cells, center cells, overlap area factor, and a workload proxy.
#'
#' @export
#' @examples
#' library(ConScapeRtools)
#' r <- terra::rast(system.file("extdata", "suitability.asc",
#'                             package = "ConScapeRtools"))
#'
#' wd <- window_design(
#'   r = r,
#'   tile_d = 5400,
#'   tile_trim = 2700,
#'   landmark = 10L
#' )
#'
#' wd$run_args
#' wd$diagnostics
#'
#' @seealso [tile_design()], [conscape_prep()], [run_conscape()]
window_design <- function(r,
                          design = NULL,
                          tile_d = NULL,
                          tile_trim = NULL,
                          landmark = NULL,
                          centersize = NULL,
                          buffer = NULL,
                          window_units = c("cells", "map")) {
  window_units <- match.arg(window_units)

  if (!inherits(r, "SpatRaster")) {
    stop("r must be a SpatRaster.", call. = FALSE)
  }

  if (!is.null(design)) {
    if (is.null(tile_d) && !is.null(design$tile_d)) tile_d <- design$tile_d
    if (is.null(tile_trim) && !is.null(design$tile_trim)) tile_trim <- design$tile_trim
    if (is.null(landmark) && !is.null(design$landmark)) landmark <- design$landmark
  }
  if (is.null(landmark)) landmark <- 10L
  if (!is.numeric(landmark) || length(landmark) != 1L ||
      is.na(landmark) || landmark < 1) {
    stop("landmark must be a positive integer.", call. = FALSE)
  }
  landmark <- as.integer(landmark)

  resx <- terra::res(r)[1]
  resy <- terra::res(r)[2]
  if (!isTRUE(all.equal(resx, resy))) {
    warning("Non-square cells; using x-resolution for window design.")
  }

  using_window_args <- !is.null(centersize) || !is.null(buffer)
  if (using_window_args) {
    if (is.null(centersize) || is.null(buffer)) {
      stop("Supply both centersize and buffer, or supply neither.", call. = FALSE)
    }
    validate_positive_scalar(centersize, "centersize")
    validate_nonnegative_scalar(buffer, "buffer")

    if (identical(window_units, "cells")) {
      centersize_cells <- as.integer(ceiling(centersize))
      buffer_cells <- as.integer(ceiling(buffer))
      tile_d <- centersize_cells * resx
      tile_trim <- buffer_cells * resx
    } else {
      tile_d <- centersize
      tile_trim <- buffer
      centersize_cells <- as.integer(ceiling(tile_d / resx))
      buffer_cells <- as.integer(ceiling(tile_trim / resx))
    }
  } else {
    if (is.null(tile_d) || is.null(tile_trim)) {
      stop(
        "Supply design, tile_d and tile_trim, or centersize and buffer.",
        call. = FALSE
      )
    }
    validate_positive_scalar(tile_d, "tile_d")
    validate_nonnegative_scalar(tile_trim, "tile_trim")
    layout <- make_tiles(r, tile_d = tile_d, tile_trim = tile_trim,
                         landmark = landmark)
    centersize_cells <- as.integer(layout$tile_cells)
    buffer_cells <- as.integer(layout$overlap_cells)
    tile_d <- as.numeric(layout$tile_width)
    tile_trim <- as.numeric(layout$tile_trim)
  }

  if (centersize_cells < 1L) {
    stop("centersize must resolve to at least one raster cell.", call. = FALSE)
  }
  if (buffer_cells < 0L) {
    stop("buffer must resolve to a non-negative number of cells.", call. = FALSE)
  }

  ncols <- terra::ncol(r)
  nrows <- terra::nrow(r)
  starts_x <- window_starts(ncols, centersize_cells, buffer_cells)
  starts_y <- window_starts(nrows, centersize_cells, buffer_cells)
  expected_windows <- length(starts_x) * length(starts_y)
  source_window_cells <- (centersize_cells + 2L * buffer_cells)^2
  center_cells <- centersize_cells^2
  overlap_area_factor <- source_window_cells / center_cells
  workload_proxy <- expected_windows * source_window_cells * center_cells

  diagnostics <- list(
    nrow = nrows,
    ncol = ncols,
    resolution = c(x = resx, y = resy),
    starts_x = starts_x,
    starts_y = starts_y,
    expected_windows = expected_windows,
    center_cells = center_cells,
    source_window_cells = source_window_cells,
    total_window_width_cells = centersize_cells + 2L * buffer_cells,
    overlap_area_factor = overlap_area_factor,
    workload_proxy = workload_proxy,
    landmark_aligned = all(c(centersize_cells, max(buffer_cells, 1L)) %% landmark == 0L)
  )

  out <- list(
    centersize = centersize_cells,
    buffer = buffer_cells,
    tile_d = as.numeric(tile_d),
    tile_trim = as.numeric(tile_trim),
    centersize_map = centersize_cells * resx,
    buffer_map = buffer_cells * resx,
    landmark = landmark,
    theta = if (!is.null(design) && !is.null(design$theta)) design$theta else NA_real_,
    distance_scale = if (!is.null(design) && !is.null(design$distance_scale)) {
      design$distance_scale
    } else {
      NA_real_
    },
    run_args = list(
      centersize = centersize_cells,
      buffer = buffer_cells
    ),
    prep_args = list(
      centersize = centersize_cells,
      buffer = buffer_cells,
      window_units = "cells",
      target_mode = "center",
      landmark = landmark
    ),
    classic_prep_args = list(
      tile_d = as.numeric(tile_d),
      tile_trim = as.numeric(tile_trim),
      landmark = landmark
    ),
    diagnostics = diagnostics
  )
  class(out) <- "ConScapeRtools_window_design"
  out
}

validate_positive_scalar <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x <= 0) {
    stop(name, " must be a single positive numeric value.", call. = FALSE)
  }
}

validate_nonnegative_scalar <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x < 0) {
    stop(name, " must be a single non-negative numeric value.", call. = FALSE)
  }
}

window_starts <- function(n, centersize, buffer) {
  limit <- n - 2L * buffer
  if (limit < 1L) {
    return(integer())
  }
  as.integer(seq.int(1L, limit, by = centersize))
}
