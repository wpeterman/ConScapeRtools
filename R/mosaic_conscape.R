#' Mosaic ConScape tiles into a single raster
#'
#' @description
#' After running [run_conscape()] on tiled landscapes, this function
#' reassembles the tile-level outputs into a single `SpatRaster`. Tiles
#' are optionally trimmed to remove overlapping margins, combined either
#' by merging or averaging, and then (optionally) masked to the original
#' analysis area.
#'
#' @param out_dir Directory where ConScape tile outputs (e.g., `"btwn"`
#'   or `"fcon"` subdirectories produced by [run_conscape()]) were
#'   written. All `*.asc` and `*.tif` files in this directory are
#'   treated as tiles belonging to a single surface.
#' @param mask Optional binary `SpatRaster` indicating regions of the
#'   landscape that are considered valid (e.g., potential habitat).
#'   Cells with 0 are set to 0 in the output, and cells with `NA` in the
#'   mask become `NA` in the output. If provided, it should typically be
#'   the same mask written by [conscape_prep()].
#' @param tile_trim Width (in map units) of the overlapping margins to
#'   remove from each tile before mosaicking. This should match the
#'   `tile_trim` value stored in the `"ConScapeRtools_prep"` object
#'   returned by [conscape_prep()], or an equivalent value used when
#'   creating tiles manually.
#' @param method Character string indicating how to combine overlapping
#'   tiles. `"merge"` (default) fills gaps using the first non-`NA` tile
#'   encountered, while `"mosaic"` averages values where tiles overlap.
#' @param crs Optional coordinate reference system to assign to the
#'   merged raster. Can be a proj4string, EPSG code (e.g. `"EPSG:4326"`),
#'   or a CRS taken from a `SpatRaster`. If `NULL` (default), the CRS of
#'   the first tile is retained.
#'
#' @details
#' For each tile, `tile_trim` is converted to a number of rows and
#' columns and trimmed symmetrically from all sides. To avoid trimming
#' away entire tiles, the trim width is capped at at most half the tile
#' width/height in cells. If `tile_trim` is too large relative to tile
#' dimensions, no trimming is applied to that tile.
#'
#' After trimming, tiles are combined into a single `SpatRaster` using
#' either [terra::merge()] (`method = "merge"`) or [terra::mosaic()]
#' (`method = "mosaic", fun = "mean"`). If a `mask` is supplied, the
#' function checks for compatible CRS, then crops and/or resamples the
#' mosaicked raster to match the mask's extent and resolution (using
#' nearest-neighbour resampling if needed) before applying the mask via
#' cell-wise multiplication.
#'
#' This function is typically called internally by [run_conscape()] when
#' `mosaic = TRUE`, but it can also be used directly to post-process
#' ConScape tile outputs.
#'
#' @return
#' A `SpatRaster` representing the mosaicked ConScape surface, with
#' overlapping tile margins removed, tiles combined according to
#' `method`, an optional CRS set, and (if provided) masked to the input
#' `mask`.
#'
#' @export
#' @examples
#' # example code
#'
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
#'
## Put output tiles together
#' cs_btwn <- mosaic_conscape(out_dir = cs_run$outdir_btwn,
#'                            method = 'mosaic',
#'                            crs = crs(source))
#' cs_fcon <- mosaic_conscape(out_dir = cs_run$outdir_fcon,
#'                            tile_trim = tile_trim,
#'                            crs = crs(source))
#' plot(c(cs_btwn, cs_fcon))
#' }#' @seealso [conscape_prep()], [run_conscape()], [make_tiles()]
#' @author Bill Peterman
#' @importFrom terra merge


mosaic_conscape <- function(out_dir,
                            mask = NULL,
                            tile_trim,
                            method = c("merge", "mosaic"),
                            crs = NULL) {
  method <- match.arg(method)

  # Check tile_trim
  if (!is.numeric(tile_trim) || length(tile_trim) != 1 || tile_trim < 0) {
    stop("tile_trim must be a single non-negative numeric value")
  }

  # List all raster files
  tile_files <- list.files(out_dir, pattern = "\\.asc$|\\.tif$", full.names = TRUE)
  if (length(tile_files) == 0) stop("No raster files found in out_dir")

  # Read all tiles
  r_list <- lapply(tile_files, rast)

  if (tile_trim > 0) {
    resx <- res(r_list[[1]])[1]
    resy <- res(r_list[[1]])[2]

    trim_cols <- round(tile_trim / resx)
    trim_rows <- round(tile_trim / resy)

    r_list <- lapply(r_list, function(x) {
      e   <- ext(x)
      rx  <- res(x)[1]
      ry  <- res(x)[2]
      nc  <- ncol(x)
      nr  <- nrow(x)

      # cap trim so we never remove more than half the tile
      tc <- min(trim_cols, floor((nc - 1L) / 2))
      tr <- min(trim_rows, floor((nr - 1L) / 2))

      if (tc <= 0 && tr <= 0) {
        return(x)  # nothing to trim safely
      }

      ext_trim <- ext(
        e$xmin + tc * rx,
        e$xmax - tc * rx,
        e$ymin + tr * ry,
        e$ymax - tr * ry
      )
      crop(x, ext_trim)
    })
  }

  # Combine tiles
  rc <- sprc(r_list)
  r_out <- if (method == "mosaic") {
    mosaic(rc, fun = "mean")
  } else {
    merge(rc)  # no algo / method tweaks
  }

  # Apply CRS if provided
  if (!is.null(crs)) crs(r_out) <- crs

  # Apply mask if provided
  # if (!is.null(mask)) r_out <- r_out * mask
  if (!is.null(mask)) {

    # Ensure same CRS
    if (!is.null(crs(mask)) && !is.null(crs(r_out)) &&
        crs(mask) != crs(r_out)) {
      stop("CRS of mask and r_out differ")
    }

    # If extents/resolutions differ slightly, bring r_out onto the mask grid
    same_geom <- isTRUE(all.equal(ext(r_out), ext(mask))) &&
      isTRUE(all.equal(res(r_out), res(mask)))

    if (!same_geom) {
      # try to bring r_out onto mask grid
      r_out <- crop(r_out, mask)
      same_geom <- isTRUE(all.equal(ext(r_out), ext(mask))) &&
        isTRUE(all.equal(res(r_out), res(mask)))

      if (!same_geom) {
        # last resort: resample r_out to mask grid
        r_out <- resample(r_out, mask, method = "near")
      }
    }

    # Now extents & res match, safe to multiply or use terra::mask()
    r_out <- r_out * mask
    # r_out <- mask(r_out, mask, maskvalues = 0, updatevalue = NA)
  }


  # Clean names
  names(r_out) <- sub("-.*", "", names(r_out))
  return(r_out)
}
