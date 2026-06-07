#' Prepare rasters and tiles for `run_conscape()`
#'
#' @description
#' Convenience wrapper that (i) checks that all input rasters are compatible,
#' (ii) designs a landmark-aligned tiling scheme with `make_tiles()`, (iii)
#' applies the same tiling to the target, source, and movement rasters, and
#' (iv) writes a binary mask that is used later by
#' [run_conscape()] and [mosaic_conscape()].
#'
#' @param tile_d Target interior tile width in map units (e.g., meters).
#'   Internally, tile dimensions are converted to a number of cells and rounded
#'   up to the nearest multiple of `landmark` to ensure that coarse–graining
#'   windows in ConScape align across tiles. Optional when `centersize` and
#'   `buffer` are supplied.
#' @param tile_trim Minimum overlap width between neighbouring tiles, in map
#'   units. The actual overlap used may be increased so that the overlap in
#'   cells is a multiple of `landmark`. This overlap is removed again in
#'   [mosaic_conscape()]. Optional when `centersize` and `buffer` are supplied.
#' @param asc_dir Directory where `.asc` tiles and ancillary files will be
#'   written. If `NULL` (default), a new directory is created under the current
#'   R session's temporary directory.
#' @param r_target `SpatRaster` representing ConScape `target_qualities`. May
#'   be the same object as `r_src`. Values below `target_threshold` are treated
#'   as unsuitable and are masked out.
#' @param r_mov `SpatRaster` representing ConScape cell affinities /
#'   permeability. These values are converted to graph `affinities` with
#'   `ConScape.graph_matrix_from_raster()`; they should not be resistance
#'   values unless they have first been transformed to affinities.
#' @param r_src `SpatRaster` representing ConScape `source_qualities`. May be
#'   the same object as `r_target`.
#' @param target_threshold Numeric threshold applied to `r_target` to define
#'   suitable habitat (default `0`). Cells with values `< target_threshold` are
#'   coded as 0 in the mask and set to `NA` in the tiled target raster.
#' @param clear_dir Logical (default `FALSE`). If `TRUE`, any existing contents
#'   of `asc_dir` are deleted before tiles are written. If `FALSE` and
#'   `asc_dir` is not empty, the function stops with an error.
#' @param landmark Integer giving the number of pixels in the coarse–graining
#'   window used by ConScape (the `npix` argument in `coarse_graining()`;
#'   default `10L`). This value is used to ensure that tile dimensions and
#'   overlaps are aligned with the coarse–graining grid.
#' @param progress Logical. If `TRUE`, simple text messages are printed to
#'   report tiling progress.
#' @param centersize Optional center window size used by the buffered-window
#'   mode. By default this is interpreted as a number of raster cells. When
#'   supplied with `buffer`, `tile_d` is derived from `centersize` and the
#'   final mosaicked output is retained for the corresponding center tile.
#' @param buffer Optional buffer width used by the buffered-window mode. By
#'   default this is interpreted as a number of raster cells. When supplied
#'   with `centersize`, `tile_trim` is derived from `buffer`.
#' @param window_units Units for `centersize` and `buffer`: `"cells"`
#'   (default) or `"map"`.
#' @param target_mode Controls which cells in a buffered tile act as RSP
#'   absorbing targets when ConScape solves the tile.
#'
#'   * `"center"` (the default under `"auto"`) restricts the target raster
#'     written for each tile to its center cells; buffered source qualities
#'     and affinities are kept in full. Adjacent tiles' centers are disjoint,
#'     so every cell is a target in exactly one tile, and [mosaic_conscape()]
#'     (called by [run_conscape()]) combines outputs with `method = "sum"` to
#'     reconstruct the global source-to-target sums up to buffer truncation.
#'     This matches the per-window math used by ConScape's dev
#'     `WindowedProblem` and converges to the untiled solution as buffer
#'     grows.
#'   * `"full"` preserves the legacy ConScapeRtools behavior: every cell in
#'     the buffered tile is a target, every tile produces a full surface, and
#'     overlapping tiles are averaged with `method = "mosaic"` (mean). This is
#'     a smoothed approximation; individual tiles undercount targets that
#'     fall in other tiles, so the mean is biased low. Use only for backward
#'     compatibility.
#'   * `"auto"` (default) chooses `"center"` whenever `centersize`/`buffer`
#'     or `tile_d`/`tile_trim` are supplied. The center-target geometry is
#'     derived from `tile_d` and `tile_trim` (rounded to whole cells) when
#'     those arguments are used.
#'
#' @details
#' All three rasters (`r_target`, `r_mov`, `r_src`) must have identical extent
#' and resolution; this is checked at the start of the function. The CRS is not
#' modified but should be the same for all layers.
#'
#' The argument names mirror the original ConScapeRtools workflow, while the
#' returned `input_slots` element records the equivalent ConScape slots:
#' `r_target` becomes `target_qualities`, `r_src` becomes `source_qualities`,
#' and `r_mov` becomes the affinity/permeability layer used to construct graph
#' `affinities`. If your input layer is a resistance or cost surface, transform
#' it to an affinity surface before calling `conscape_prep()`.
#'
#' A binary mask (`mask.asc`) is written to `file.path(asc_dir, "mask")` and
#' contains 1 for cells with `r_target >= target_threshold` and 0 otherwise.
#' This mask is used later to restrict the mosaicked ConScape outputs to the
#' potentially occupiable landscape.
#'
#' Tiles are designed once internally using `make_tiles()` based on the thresholded target
#' raster, and the same tile geometry is then applied internally to all layers via the
#' `tile_rast()`. Tile sizes and overlaps are expressed in cell units and
#' rounded so that both tile width and overlap are multiples of `landmark`,
#' which ensures that ConScape's coarse–grained "landmarks" fall on a common
#' grid when tiles are reassembled.
#'
#' @return
#' A named list of class `"ConScapeRtools_prep"` with elements:
#'
#' * `cs_tiles` – `SpatVector` of interior tile polygons.
#' * `tile_num` – integer vector of tile IDs that contain usable data.
#' * `asc_dir` – path to the root directory where `.asc` tiles and the mask
#'   were written.
#' * `src`, `target`, `mov` – subdirectories containing the source, target,
#'   and movement tiles (as `.asc` files).
#' * `tile_trim` – effective overlap width in map units (may be larger than
#'   the user–supplied `tile_trim` if increased for landmark alignment).
#' * `landmark` – the coarse–graining window size passed in via `landmark`.
#' * `input_slots` – named paths showing how prepared rasters map to ConScape's
#'   `target_qualities`, `source_qualities`, and `affinities` slots.
#'
#' This object is intended to be passed directly to [run_conscape()].
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
#' ## Calibrate tile and decay parameters
#' td <- tile_design(r_mov    = affinity,
#'                   r_target = habitat,
#'                   max_d    = 7000,
#'                   theta    = 0.1,
#'                   jl_home  = jl_home,
#'                   landmark = 5L)
#'
#' ## Prepare tiled rasters
#' prep <- conscape_prep(tile_d    = td$tile_d,
#'                       tile_trim = td$tile_trim,
#'                       r_target  = habitat,
#'                       r_mov     = affinity,
#'                       r_src     = habitat,
#'                       clear_dir = TRUE,
#'                       landmark  = td$landmark)
#'
#' ## Run ConScape on tiles
#' cs_res <- run_conscape(conscape_prep  = prep,
#'                        out_dir        = "conscape_out",
#'                        theta          = td$theta,
#'                        distance_scale = td$distance_scale,
#'                        jl_home        = jl_home)
#' }
#'
#' @seealso [run_conscape()], [mosaic_conscape()]
#' @author Bill Peterman

#' @importFrom terra writeVector setGDALconfig
#' @importFrom utils txtProgressBar setTxtProgressBar

conscape_prep <- function(tile_d = NULL,
                          tile_trim = NULL,
                          asc_dir = NULL,
                          r_target,
                          r_mov,
                          r_src,
                          target_threshold = 0,
                          clear_dir = FALSE,
                          landmark = 10L,
                          progress = TRUE,
                          centersize = NULL,
                          buffer = NULL,
                          window_units = c("cells", "map"),
                          target_mode = c("auto", "full", "center")) {

  terra::setGDALconfig("GDAL_PAM_ENABLED", "FALSE")
  on.exit(terra::setGDALconfig("GDAL_PAM_ENABLED", "TRUE"), add = TRUE)

  window_units <- match.arg(window_units)
  target_mode <- match.arg(target_mode)

  ## ----- basic checks: extents and resolutions -----
  rasters <- list(r_target = r_target, r_mov = r_mov, r_src = r_src)

  # extents
  if (!all(vapply(rasters, function(x) isTRUE(all.equal(terra::ext(x), terra::ext(r_target))), logical(1)))) {
    stop("All rasters (r_target, r_mov, r_src) must have the same extent.")
  }
  # resolution
  if (!all(vapply(rasters, function(x) isTRUE(all.equal(terra::res(x), terra::res(r_target))), logical(1)))) {
    stop("All rasters (r_target, r_mov, r_src) must have the same resolution.")
  }

  resx <- terra::res(r_target)[1]

  using_window_args <- !is.null(centersize) || !is.null(buffer)
  if (using_window_args) {
    if (is.null(centersize) || is.null(buffer)) {
      stop("Supply both centersize and buffer, or supply neither.", call. = FALSE)
    }
    if (!is.numeric(centersize) || length(centersize) != 1L ||
        is.na(centersize) || centersize <= 0) {
      stop("centersize must be a single positive numeric value.", call. = FALSE)
    }
    if (!is.numeric(buffer) || length(buffer) != 1L ||
        is.na(buffer) || buffer < 0) {
      stop("buffer must be a single non-negative numeric value.", call. = FALSE)
    }
    if (identical(window_units, "cells")) {
      tile_d <- centersize * resx
      tile_trim <- buffer * resx
      centersize_cells <- as.integer(ceiling(centersize))
      buffer_cells <- as.integer(ceiling(buffer))
    } else {
      tile_d <- centersize
      tile_trim <- buffer
      centersize_cells <- as.integer(ceiling(centersize / resx))
      buffer_cells <- as.integer(ceiling(buffer / resx))
    }
    if (identical(target_mode, "auto")) target_mode <- "center"
  } else {
    if (is.null(tile_d) || is.null(tile_trim)) {
      stop("Supply tile_d and tile_trim, or supply centersize and buffer.", call. = FALSE)
    }
    # In tile_d/tile_trim mode, the user has not explicitly chosen a windowing
    # scheme. We derive centersize/buffer from tile_d/tile_trim so center mode
    # has well-defined geometry, then default to "center". Users who depend on
    # the legacy biased-mean tiled behavior can opt in with target_mode = "full".
    resx_cells <- resx
    centersize_cells <- as.integer(ceiling(tile_d / resx_cells))
    buffer_cells <- as.integer(ceiling(tile_trim / resx_cells))
    if (identical(target_mode, "auto")) target_mode <- "center"
  }

  if (!is.numeric(tile_d) || length(tile_d) != 1L || is.na(tile_d) || tile_d <= 0) {
    stop("tile_d must be a single positive numeric value.", call. = FALSE)
  }
  if (!is.numeric(tile_trim) || length(tile_trim) != 1L || is.na(tile_trim) || tile_trim < 0) {
    stop("tile_trim must be a single non-negative numeric value.", call. = FALSE)
  }

  ## ----- set up asc_dir and (optionally) clear it -----
  if (is.null(asc_dir)) {
    asc_dir <- file.path(tempdir(), "asc")
  }
  if (!dir.exists(asc_dir)) dir.create(asc_dir, recursive = TRUE)
  asc_dir <- normalizePath(asc_dir)

  # if asc_dir already has content and we don't want to clear, stop
  if (length(list.files(asc_dir, recursive = TRUE)) > 0 && !clear_dir) {
    stop("asc_dir is not empty. Use clear_dir = TRUE to overwrite existing contents.")
  }

  if (clear_dir) {
    unlink(file.path(asc_dir, "*"), recursive = TRUE, force = TRUE)
  }

  ## ----- create and write mask (this is what run_conscape expects) -----

  # binary mask based on target_threshold
  t_mask <- terra::ifel(is.na(r_target) | r_target < target_threshold, 0, 1)

  mask_dir <- file.path(asc_dir, "mask")
  if (!dir.exists(mask_dir)) dir.create(mask_dir, recursive = TRUE)

  terra::writeRaster(
    t_mask,
    filename = file.path(mask_dir, "mask.asc"),
    overwrite = TRUE,
    NAflag   = -9999
  )

  ## also apply threshold to r_target for tiling
  r_target_thr <- r_target
  r_target_thr[r_target_thr < target_threshold] <- NA

  ## ----- design tiles (npix / landmark aligned) -----
  if (progress) cat("Designing tiles...\n")
  tile_design <- make_tiles(r_target_thr,
                            tile_d    = tile_d,
                            tile_trim = tile_trim,
                            landmark  = landmark)

  tmp <- ensure_tile_id(tile_design$cs_tiles, id_col = NULL)
  tile_design$cs_tiles <- tmp$cs_tiles
  id_col <- tmp$id_col

  ## ----- tile all rasters with the same design -----
  if (progress) cat("Writing matched tile triplets...\n")
  tile_triplets <- tile_rast_triplets(
    r_target = r_target_thr,
    r_mov = r_mov,
    r_src = r_src,
    tiles = tile_design,
    root_dir = asc_dir,
    clear_dir = TRUE,
    progress = progress,
    target_mode = target_mode
  )

  if (progress) cat("\nValidating tiles...\n")
  # Ensure a stable sequential id exists (optional but helpful for tracking)
  if (!("tile_id" %in% names(tile_design$cs_tiles))) {
    tile_design$cs_tiles$tile_id <- seq_len(nrow(tile_design$cs_tiles))
  }

  val <- list(
    ok_idx = tile_triplets$ok_idx,
    bad_idx = tile_triplets$bad_idx,
    bad_files = tile_triplets$bad_files
  )
  if (length(val$ok_idx) == 0L) {
    stop("No valid ConScape tile triplets were written.", call. = FALSE)
  }
  if (progress) {
    message(sprintf("Tile validation: %d ok, %d bad (deleted=TRUE)",
                    length(val$ok_idx), length(val$bad_idx)))
  }

  # Subset polygons by row index, not filename-derived IDs
  tile_design$cs_tiles <- tile_design$cs_tiles[val$ok_idx, , drop = FALSE]
  tile_design$tile_num <- tile_design$tile_num[val$ok_idx]

  # Save bookkeeping
  tile_design$tile_idx <- val$ok_idx
  tile_design$tile_validation <- list(
    n_ok  = length(val$ok_idx),
    n_bad = length(val$bad_idx),
    bad   = val$bad_files,
    diagnostics = tile_triplets$diagnostics
  )

  ## ----- write tile polygons for reference -----
  terra::writeVector(tile_design$cs_tiles,
                     file.path(asc_dir, "tiles.shp"),
                     overwrite = TRUE)

  ## ----- return prep object -----
  terra::setGDALconfig("GDAL_PAM_ENABLED", "TRUE")

  out <- list(
    cs_tiles  = tile_design$cs_tiles,
    tile_num  = tile_design$tile_num,
    asc_dir   = asc_dir,
    src       = tile_triplets$src_dir,
    target    = tile_triplets$target_dir,
    mov       = tile_triplets$mov_dir,
    tile_trim = tile_design$tile_trim,   # map units, for mosaic_conscape
    landmark  = tile_design$landmark,
    target_mode = target_mode,
    centersize = centersize_cells,
    buffer = buffer_cells,
    window_units = window_units,
    tile_d = tile_design$tile_width,
    target_threshold = target_threshold,
    tile_diagnostics = tile_triplets$diagnostics,
    input_slots = list(
      target_qualities = tile_triplets$target_dir,
      source_qualities = tile_triplets$src_dir,
      affinities = tile_triplets$mov_dir
    )
  )
  class(out) <- "ConScapeRtools_prep"
  out
}


