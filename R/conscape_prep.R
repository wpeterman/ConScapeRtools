#' Prepare rasters and tiles for `run_conscape()`
#'
#' @description
#' Convenience wrapper that (i) checks that all input rasters are compatible,
#' (ii) designs a landmark–aligned tiling scheme with `make_tiles()`, (iii)
#' applies the same tiling to the target, source, and movement rasters via
#' `tile_rast()`, and (iv) writes a binary mask that is used later by
#' [run_conscape()] and [mosaic_conscape()].
#'
#' @param tile_d Target interior tile width in map units (e.g., meters).
#'   Internally, tile dimensions are converted to a number of cells and rounded
#'   up to the nearest multiple of `landmark` to ensure that coarse–graining
#'   windows in ConScape align across tiles.
#' @param tile_trim Minimum overlap width between neighbouring tiles, in map
#'   units. The actual overlap used may be increased so that the overlap in
#'   cells is a multiple of `landmark`. This overlap is removed again in
#'   [mosaic_conscape()].
#' @param asc_dir Directory where `.asc` tiles and ancillary files will be
#'   written. If `NULL` (default), a new directory is created under the current
#'   R session's temporary directory.
#' @param r_target `SpatRaster` representing target qualities. May be the same
#'   object as `r_src`. Values below `target_threshold` are treated as
#'   unsuitable and are masked out.
#' @param r_mov `SpatRaster` representing movement probabilities.
#' @param r_src `SpatRaster` representing source qualities. May be the same
#'   object as `r_target`.
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
#'
#' @details
#' All three rasters (`r_target`, `r_mov`, `r_src`) must have identical extent
#' and resolution; this is checked at the start of the function. The CRS is not
#' modified but should be the same for all layers.
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
#'
#' This object is intended to be passed directly to [run_conscape()].
#'
#' @export
#' @examples
#' \dontrun{
#' library(ConScapeRtools)
#'
#' ## Import data
#' s <- system.file("data/suitability.asc", package = "ConScapeRtools")
#' source <- terra::rast(s)
#'
#' a <- system.file("data/affinity.asc", package = "ConScapeRtools")
#' resist <- terra::rast(a)
#'
#' jl_home <- "/path/to/julia/bin"
#'
#' td <- tile_design(r_mov = resist,
#'                   r_source = source,
#'                   max_d = 7000,
#'                   theta = 0.1,
#'                   jl_home = jl_home)
#'
#' ## Tile dimension
#' tile_d <- td$tile_d
#'
#' # How much to trim tiles
#' tile_trim <- td$tile_trim
#'
#' # Makes computation more efficient
#' landmark <- 5L # Must be an integer, not numeric
#'
#' # Controls level of randomness of paths
#' theta <- td$theta
#'
#' # Controls rate of decay with distance
#' exp_d <- td$exp_d
#'
#' ## Prepare data for analysis
#' prep <- conscape_prep(tile_d = tile_d,
#'                       tile_trim = tile_trim,
#'                       r_target = source,
#'                       r_mov = resist,
#'                       r_src = source,
#'                       clear_dir = T,
#'                       landmark = landmark)
#'
#' cs_res <- run_conscape(conscape_prep = prep,
#'                        out_dir = "conscape_out",
#'                        jl_home = jl_home)
#' }
#'
#' @seealso [run_conscape()], [mosaic_conscape()]
#' @author Bill Peterman

#' @importFrom terra writeVector setGDALconfig
#' @importFrom utils txtProgressBar setTxtProgressBar

conscape_prep <- function(tile_d,
                          tile_trim,
                          asc_dir = NULL,
                          r_target,
                          r_mov,
                          r_src,
                          target_threshold = 0,
                          clear_dir = FALSE,
                          landmark = 10L,
                          progress = TRUE) {

  setGDALconfig("GDAL_PAM_ENABLED", "FALSE")

  ## ----- basic checks: extents and resolutions -----
  rasters <- list(r_target = r_target, r_mov = r_mov, r_src = r_src)

  # extents
  if (!all(sapply(rasters, function(x) all.equal(ext(x), ext(r_target))))) {
    stop("All rasters (r_target, r_mov, r_src) must have the same extent.")
  }
  # resolution
  if (!all(sapply(rasters, function(x) all.equal(res(x), res(r_target))))) {
    stop("All rasters (r_target, r_mov, r_src) must have the same resolution.")
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
  t_mask <- r_target
  t_mask[is.na(t_mask)] <- 0
  t_mask[t_mask <  target_threshold] <- 0
  t_mask[t_mask >= target_threshold] <- 1

  mask_dir <- file.path(asc_dir, "mask")
  if (!dir.exists(mask_dir)) dir.create(mask_dir, recursive = TRUE)

  writeRaster(
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

  ## ----- tile each raster with same design -----
  if (progress) cat("Tiling target raster...\n")
  target_tiles <- tile_rast(
    r          = r_target_thr,
    tiles = tile_design,
    out_dir    = file.path(asc_dir, "target"),
    clear_dir  = TRUE,
    progress   = progress
  )

  if (progress) cat("Tiling movement raster...\n")
  mov_tiles <- tile_rast(
    r          = r_mov,
    tiles = tile_design,
    out_dir    = file.path(asc_dir, "mov"),
    clear_dir  = TRUE,
    progress   = progress
  )

  if (progress) cat("Tiling source raster...\n")
  src_tiles <- tile_rast(
    r          = r_src,
    tiles = tile_design,
    out_dir    = file.path(asc_dir, "src"),
    clear_dir  = TRUE,
    progress   = progress
  )

  cat("\nValidating tiles...\n")
  # Ensure a stable sequential id exists (optional but helpful for tracking)
  if (!("tile_id" %in% names(tile_design$cs_tiles))) {
    tile_design$cs_tiles$tile_id <- seq_len(nrow(tile_design$cs_tiles))
  }

  val <- validate_conscape_tiles_idx(
    target_dir = target_tiles$asc_dir,
    source_dir = src_tiles$asc_dir,
    move_dir   = mov_tiles$asc_dir,
    delete_bad = TRUE,
    quiet      = FALSE
  )

  # Subset polygons by row index, not filename-derived IDs
  tile_design$cs_tiles <- tile_design$cs_tiles[val$ok_idx, , drop = FALSE]

  # Save bookkeeping
  tile_design$tile_idx <- val$ok_idx
  tile_design$tile_validation <- list(
    n_ok  = length(val$ok_idx),
    n_bad = length(val$bad_idx),
    bad   = val$bad_files
  )

  ## ----- write tile polygons for reference -----
  writeVector(tile_design$cs_tiles,
              file.path(asc_dir, "tiles.shp"),
              overwrite = TRUE)

  ## ----- return prep object -----
  setGDALconfig("GDAL_PAM_ENABLED", "TRUE")

  out <- list(
    cs_tiles  = tile_design$cs_tiles,
    tile_num  = tile_design$tile_num,
    asc_dir   = asc_dir,
    src       = src_tiles$asc_dir,
    target    = target_tiles$asc_dir,
    mov       = mov_tiles$asc_dir,
    tile_trim = tile_design$tile_trim,   # map units, for mosaic_conscape
    landmark  = tile_design$landmark
  )
  class(out) <- "ConScapeRtools_prep"
  out
}


