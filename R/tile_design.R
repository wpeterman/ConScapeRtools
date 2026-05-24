#' Design tiling and decay parameters for ConScape
#'
#' @description
#' Uses the resolution and values of a movement affinity raster to derive a
#' calibrated exponential decay parameter (`distance_scale`) and suggested tile
#' width (`tile_d`) and overlap (`tile_trim`) for ConScape analyses, given a
#' presumed maximum dispersal distance in ideal habitat. Internally, this builds
#' a fully connected synthetic ideal-habitat landscape, runs ConScape's
#' randomized shortest-path model from the center cell, and calibrates the decay
#' so that proximity at `max_d` drops to `calibration_threshold`.
#'
#' @param r_mov `SpatRaster` representing movement affinities/permeability.
#'   Values should be in a projected coordinate reference system (map units in
#'   meters) and should be higher where movement is easier.
#' @param r_source Optional `SpatRaster` representing source quality
#'   (default `NULL`). If only `r_target` is provided, it is used for both
#'   source and target.
#' @param r_target Optional `SpatRaster` representing target quality
#'   (default `NULL`). If only `r_source` is provided, it is used for both
#'   source and target.
#' @param max_d Presumed maximum movement distance through ideal habitat,
#'   in the same map units as the raster CRS (typically meters). Used to
#'   calibrate the exponential distance–decay parameter.
#' @param theta Parameter controlling path randomness in the randomized
#'   shortest–path model (`GridRSP`). As `theta` approaches 0, movement
#'   becomes increasingly random; large values approach deterministic
#'   least–cost paths. Default is `0.1`.
#' @param jl_home Path to the `bin` directory where the Julia executable
#'   resides. This is passed to `JuliaConnectoR` and used to initialize
#'   the ConScape / ConScapeR environment.
#' @param calibration_threshold Proximity target used to calibrate
#'   `distance_scale` at `max_d`. Default is `0.025`.
#' @param trim_threshold Maximum proximity expected at the retained tile
#'   boundary after trimming. Smaller values create larger tile overlaps.
#'   Default is `0.01`.
#' @param tile_width Optional target interior tile width in map units. Use this
#'   when memory or runtime constraints require a specific tile size. If `NULL`,
#'   `max_tile_cells` or `overlap_area_factor` determines the tile width.
#' @param max_tile_cells Optional approximate maximum number of interior cells
#'   per tile. Used only when `tile_width` is `NULL`.
#' @param overlap_area_factor Desired ratio of extended tile area to retained
#'   interior tile area when `tile_width` and `max_tile_cells` are both `NULL`.
#'   The default `2` chooses the smallest interior tile width that keeps the
#'   per-tile computational area near twice the retained area before landmark
#'   rounding.
#' @param landmark Integer coarse-graining window size used to round tile width
#'   and trim to the same landmark-aligned grid used by [conscape_prep()].
#' @param stop_julia Logical. If `TRUE` (default), close the Julia session when
#'   `tile_design()` exits. Set to `FALSE` to keep Julia open for faster
#'   repeated calls with the same Julia path and process settings. Use
#'   [conscape_julia_stop()] when finished.
#'
#' @details
#' The input `r_mov` must be in a projected CRS so that cell resolution and
#' `max_d` are expressed in compatible distance units. Only the cell
#' resolution and maximum values of `r_mov`, `r_source`, and `r_target`
#' are used: the function constructs a fully connected synthetic landscape
#' whose size is chosen to span roughly four times `max_d`, fills the movement
#' grid with the maximum movement probability, and assigns constant maximum
#' source and target quality.
#'
#' A ConScape `Grid` and `GridRSP` object are then created, and expected
#' costs are computed from the center cell to all others. The decay parameter
#' `distance_scale` is chosen by 1D optimization so that the exponential
#' proximity at the cell closest to `max_d` drops to a specified threshold.
#' This value is then used to derive:
#'
#' * `tile_trim` from `trim_threshold`, so tile-edge influence is explicit,
#' * `tile_d` from `tile_width`, `max_tile_cells`, or `overlap_area_factor`,
#' * diagnostics showing requested and effective trim, expected tile count,
#'   overlap area factor, and proximity at the final trim distance.
#'
#' Pass `tile_d` and `tile_trim` to [conscape_prep()], and `distance_scale`
#' and `theta` to [run_conscape()].
#'
#' @return
#' A named list of class `"ConScapeRtools_design"` with elements:
#'
#' * `distance_scale` – calibrated exponential decay parameter (numerator of
#'   `exp(-dist / distance_scale)`) to pass to [run_conscape()]. Also
#'   accessible as `exp_d` for backwards compatibility.
#' * `exp_d` – same as `distance_scale`; retained for backwards compatibility.
#' * `tile_d` – suggested minimum interior tile width (map units) for
#'   [conscape_prep()].
#' * `tile_trim` – suggested minimum tile overlap / trim width (map units) for
#'   [conscape_prep()] and [mosaic_conscape()].
#' * `theta` – the `theta` value used for calibration, to pass to
#'   [run_conscape()].
#' * `landmark` – the landmark value used for tile and trim rounding.
#' * `trim_threshold` – the requested maximum proximity at the trim distance.
#' * `overlap_area_factor` – realized extended-area to retained-area ratio.
#' * `diagnostics` – a list describing requested trim, effective trim after
#'   landmark rounding, expected tile count, overlap area factor, and proximity
#'   at the effective trim.
#'
#' The function also prints a short, colorized summary of these design
#' parameters to the console.
#'
#' @export
#' @examples
#' \dontrun{
#' library(ConScapeRtools)
#'
#' ## Import example data
#' s <- system.file("extdata", "suitability.asc", package = "ConScapeRtools")
#' habitat <- terra::rast(s)
#'
#' a <- system.file("extdata", "affinity.asc", package = "ConScapeRtools")
#' affinity <- terra::rast(a)
#'
#' jl_home <- "/path/to/julia/bin"
#'
#' td <- tile_design(r_mov    = affinity,
#'                   r_target = habitat,
#'                   max_d    = 7000,
#'                   theta    = 0.1,
#'                   jl_home  = jl_home,
#'                   trim_threshold = 0.01,
#'                   overlap_area_factor = 2,
#'                   landmark = 5L)
#'
#' td$diagnostics
#'
#' ## Pass results directly to conscape_prep() and run_conscape()
#' # td$tile_d, td$tile_trim, td$landmark  --> conscape_prep()
#' # td$distance_scale, td$theta --> run_conscape()
#' }
#' @seealso [conscape_prep()], [run_conscape()]
#' @author Bill Peterman
#' @details
#' `r_mov` must be in a projected coordinate reference system.
#' The function relies on helper routines adapted from the `ConScapeR`
#' package (https://github.com/ConScape/ConScapeR) and uses Julia via
#' `JuliaConnectoR`.
#'
#' @rdname tile_design
#' @importFrom crayon %+% green red bold cyan
#' @importFrom JuliaConnectoR juliaEval juliaImport juliaSetupOk juliaFun stopJulia
#' @importFrom stats optimize


tile_design <- function(r_mov,
                        r_source = NULL,
                        r_target = NULL,
                        max_d,
                        theta = 0.1,
                        jl_home,
                        calibration_threshold = 0.025,
                        trim_threshold = 0.01,
                        tile_width = NULL,
                        max_tile_cells = NULL,
                        overlap_area_factor = 2,
                        landmark = 10L,
                        stop_julia = TRUE) {
  if (is.null(r_source) && is.null(r_target)) {
    stop("At least one of r_source or r_target must be provided.")
  }
  if (!is.numeric(max_d) || length(max_d) != 1L || is.na(max_d) || max_d <= 0) {
    stop("max_d must be a single positive numeric value.", call. = FALSE)
  }
  if (!is.numeric(theta) || length(theta) != 1L || is.na(theta) || theta <= 0) {
    stop("theta must be a single positive numeric value.", call. = FALSE)
  }
  if (!is.numeric(calibration_threshold) || length(calibration_threshold) != 1L ||
      is.na(calibration_threshold) || calibration_threshold <= 0 || calibration_threshold >= 1) {
    stop("calibration_threshold must be a single numeric value between 0 and 1.", call. = FALSE)
  }
  if (!is.numeric(trim_threshold) || length(trim_threshold) != 1L ||
      is.na(trim_threshold) || trim_threshold <= 0 || trim_threshold >= 1) {
    stop("trim_threshold must be a single numeric value between 0 and 1.", call. = FALSE)
  }
  if (!is.numeric(overlap_area_factor) || length(overlap_area_factor) != 1L ||
      is.na(overlap_area_factor) || overlap_area_factor <= 1) {
    stop("overlap_area_factor must be a single numeric value greater than 1.", call. = FALSE)
  }
  if (!is.null(tile_width) &&
      (!is.numeric(tile_width) || length(tile_width) != 1L || is.na(tile_width) || tile_width <= 0)) {
    stop("tile_width must be NULL or a single positive numeric value.", call. = FALSE)
  }
  if (!is.null(max_tile_cells) &&
      (!is.numeric(max_tile_cells) || length(max_tile_cells) != 1L ||
       is.na(max_tile_cells) || max_tile_cells < 1)) {
    stop("max_tile_cells must be NULL or a single numeric value of at least 1.", call. = FALSE)
  }
  if (!is.null(tile_width) && !is.null(max_tile_cells)) {
    stop("Use either tile_width or max_tile_cells, not both.", call. = FALSE)
  }
  if (!is.numeric(landmark) || length(landmark) != 1L || is.na(landmark) || landmark < 1) {
    stop("landmark must be a positive integer.", call. = FALSE)
  }
  if (!is.logical(stop_julia) || length(stop_julia) != 1L || is.na(stop_julia)) {
    stop("stop_julia must be TRUE or FALSE.", call. = FALSE)
  }
  landmark <- as.integer(landmark)
  if (isTRUE(stop_julia)) {
    on.exit(stop_conscape_julia(), add = TRUE)
  }

  conscape_julia_start(jl_home, quiet = TRUE)

  if (!is.null(r_source) && is.null(r_target)) {
    r_target <- r_source
    message("r_target not provided. Using r_source as r_target.")
  } else if (is.null(r_source) && !is.null(r_target)) {
    r_source <- r_target
    message("r_source not provided. Using r_target as r_source.")
  }

  r_res <- terra::res(r_mov)
  resx <- r_res[1]
  resy <- r_res[2]
  cal_ncol <- max(3L, ceiling((4 * max_d) / resx))
  cal_nrow <- max(3L, ceiling((4 * max_d) / resy))
  r_crs <- terra::crs(r_mov)

  ## Max Values
  mx_p <- terra::global(r_mov, 'max', na.rm = TRUE)[[1]]
  mx_src <- terra::global(r_source, 'max', na.rm = TRUE)[[1]]
  mx_target <- terra::global(r_target, 'max', na.rm = TRUE)[[1]]
  if (!is.finite(mx_p) || !is.finite(mx_src) || !is.finite(mx_target)) {
    stop("r_mov, r_source, and r_target must contain at least one finite value.", call. = FALSE)
  }

  ## Create a fully connected ideal-habitat calibration raster
  mov <- terra::rast(nrows = cal_nrow, ncols = cal_ncol,
                     vals = mx_p,
                     resolution = r_res, crs = r_crs,
                     xmin = 0, xmax = resx * cal_ncol,
                     ymin = 0, ymax = resy * cal_nrow)
  src <- terra::rast(nrows = cal_nrow, ncols = cal_ncol,
                     vals = mx_src,
                     resolution = r_res, crs = r_crs,
                     xmin = 0, xmax = resx * cal_ncol,
                     ymin = 0, ymax = resy * cal_nrow)
  target <- terra::rast(nrows = cal_nrow, ncols = cal_ncol,
                        vals = mx_target,
                        resolution = r_res, crs = r_crs,
                        xmin = 0, xmax = resx * cal_ncol,
                        ymin = 0, ymax = resy * cal_nrow)

  center_cell <- terra::cellFromRowCol(
    mov,
    row = ceiling(terra::nrow(mov) / 2),
    col = ceiling(terra::ncol(mov) / 2)
  )
  cell_xy <- terra::crds(target)
  center_coord <- cell_xy[center_cell, ]
  e_dist <- sqrt((cell_xy[, 1] - center_coord[1])^2 +
                   (cell_xy[, 2] - center_coord[2])^2)

  max_cell <- which.min(abs(e_dist - max_d))

  # Create ConScape Grid
  g <- Grid(affinities = mov,
            sources = src,
            targets = target,
            costs = "minuslog")

  # Create a ConScape GridRSP by providing the randomness parameter theta
  h <- GridRSP(g, theta = theta)

  dists <- expected_cost(h)
  source_col <- if (is.matrix(dists) && ncol(dists) >= center_cell) center_cell else 1L

  exp_d <- optimize(exp_opt,
                    dists = matrix(dists[, source_col], ncol = 1),
                    cell = max_cell,
                    threshold = calibration_threshold,
                    interval = c(1,2500),
                    tol = 0.001,
                    maximum = T)

  ## Trim and tile size
  requested_trim <- exp_distance(p = 1 - trim_threshold, lambda = 1/exp_d$maximum)
  cell_aligned_trim <- ceiling(requested_trim / resx) * resx

  if (!is.null(tile_width)) {
    raw_tile_width <- tile_width
    tile_width_source <- "tile_width"
  } else if (!is.null(max_tile_cells)) {
    raw_tile_width <- sqrt(max_tile_cells) * resx
    tile_width_source <- "max_tile_cells"
  } else {
    raw_tile_width <- 2 * cell_aligned_trim / (sqrt(overlap_area_factor) - 1)
    tile_width_source <- "overlap_area_factor"
  }
  requested_tile_width <- ceiling(raw_tile_width / resx) * resx
  tile_layout <- make_tiles(
    r = r_mov,
    tile_d = requested_tile_width,
    tile_trim = cell_aligned_trim,
    landmark = landmark
  )

  exp_d_val <- as.numeric(round(exp_d$maximum, 1))
  effective_tile_width <- as.numeric(tile_layout$tile_width)
  effective_tile_trim <- as.numeric(tile_layout$tile_trim)
  realized_overlap_area_factor <- ((effective_tile_width + 2 * effective_tile_trim) /
                                     effective_tile_width)^2
  proximity_at_effective_trim <- exp(-effective_tile_trim / exp_d$maximum)
  diagnostics <- list(
    calibration_threshold = calibration_threshold,
    trim_threshold = trim_threshold,
    requested_tile_trim = as.numeric(requested_trim),
    cell_aligned_tile_trim = as.numeric(cell_aligned_trim),
    effective_tile_trim = effective_tile_trim,
    requested_tile_width = as.numeric(requested_tile_width),
    effective_tile_width = effective_tile_width,
    tile_width_source = tile_width_source,
    max_tile_cells = max_tile_cells,
    landmark = landmark,
    expected_tile_count = length(tile_layout$tile_num),
    tile_cells = tile_layout$tile_cells,
    overlap_cells = tile_layout$overlap_cells,
    overlap_area_factor = as.numeric(realized_overlap_area_factor),
    proximity_at_requested_trim = exp(-requested_trim / exp_d$maximum),
    proximity_at_effective_trim = proximity_at_effective_trim,
    calibration_nrow = cal_nrow,
    calibration_ncol = cal_ncol,
    source_cell = center_cell,
    max_d_cell = max_cell,
    max_d_cell_distance = as.numeric(e_dist[max_cell])
  )
  out <- list(distance_scale = exp_d_val,
              exp_d = exp_d_val,
              tile_d = effective_tile_width,
              tile_trim = effective_tile_trim,
              theta = as.numeric(theta),
              landmark = landmark,
              trim_threshold = trim_threshold,
              overlap_area_factor = as.numeric(realized_overlap_area_factor),
              diagnostics = diagnostics)

  cat(green("\n" %+% cyan("     *** Tile Design Parameters ***") %+%"\n",
            "Given a maximum dispersal of " %+% red$bold(max_d) %+%" meters,\n",
            "`distance_scale` should be set to: " %+% red$bold(out$distance_scale) %+% "\n\n",
            "`tile_d` should be at least: " %+% red$bold(out$tile_d) %+%",\n\n",
            "`tile_trim` should be at least: " %+% red$bold(out$tile_trim) %+%",\n\n",
            "Expected tiles: " %+% red$bold(diagnostics$expected_tile_count) %+% "\n",
            "Overlap area factor: " %+% red$bold(round(out$overlap_area_factor, 2)) %+% "\n",
            "Proximity at trim: " %+% red$bold(signif(diagnostics$proximity_at_effective_trim, 3)) %+% "\n\n"))
  class(out) <- 'ConScapeRtools_design'
  invisible(suppressMessages(juliaEval('1+1')))
  return(out)
}

# Optimization function ------------------------------------------------

exp_opt <- function(x, dists, cell, threshold){
  prox <- exp(-dists[,1] * (1/x))
  y <- threshold - prox[cell]
  if(y > 0){
    y <- -9999
  }
  return(y)
}

#' Calculate the distance for a given cumulative density in a negative exponential distribution
#' @noRd

exp_distance <- function(p, lambda) {
  if (p <= 0 || p >= 1) {
    stop("Cumulative density `p` must be between 0 and 1 (exclusive).")
  }
  if (lambda <= 0) {
    stop("Rate parameter `lambda` must be positive.")
  }
  x <- -log(1 - p) / lambda
  return(x)
}

