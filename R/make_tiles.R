#' Compute landmark-aligned tile layout for a raster
#'
#' @description
#' Internal helper to design a regular grid of interior tiles for a
#' `SpatRaster`. Tile dimensions are expressed in map units but internally
#' converted to a number of cells that is a multiple of the landmark
#' coarse-grain window (`npix`). Tile origins are aligned to the original
#' raster grid so that ConScape landmark nodes line up across tiles and
#' with a full, untiled run.
#'
#' @param r `SpatRaster`. Raster to be partitioned into interior tiles.
#'   Must be in a projected coordinate reference system where the cell
#'   size is expressed in the same units as `tile_d` and `tile_trim`
#'   (typically meters). Non-square cells are allowed but only the
#'   x-resolution is used when computing tile sizes and overlap.
#' @param tile_d Numeric scalar giving the *target* interior tile width
#'   in map units (e.g., meters). The actual tile width in cells is
#'   computed as the smallest multiple of `landmark` (in cells) that is
#'   at least `tile_d / res(r)[1]`.
#' @param tile_trim Numeric scalar giving the minimum requested overlap
#'   width (in map units) to be available on each side of a tile before
#'   trimming. This is treated as a lower bound; the effective overlap
#'   may be larger to satisfy landmark alignment.
#' @param landmark Integer giving the number of pixels in the
#'   coarse-graining window (`npix`) used by ConScape. Tiles are sized
#'   and their overlaps are constrained to be multiples of this value so
#'   that coarse-grained nodes (landmarks) fall on consistent grid
#'   locations across tiles. Default is `10L`.
#'
#' @details
#' The function works entirely in cell space but keeps all calculations
#' anchored to the original raster grid:
#'
#' * The target tile width `tile_d` is converted to a raw number of
#'   cells (`tile_cells_raw = round(tile_d / resx)`), then rounded *up*
#'   to the nearest multiple of `landmark`, with a minimum of
#'   `landmark` cells (`tile_cells`).
#' * Starting column and row indices are generated with
#'   `seq(1, ncol(r), by = tile_cells)` and
#'   `seq(1, nrow(r), by = tile_cells)`. The last tile in each row/column
#'   is truncated to the raster edge so that the union of tiles covers
#'   the full extent of `r`.
#' * For each tile, the interior extent is computed directly from
#'   `(cs, rs, ce, re)` and the raster origin, using `xmin(r)`, `ymax(r)`
#'   and the cell size. No padding or overlap is added at this step; the
#'   returned polygons describe only the interior region.
#'
#' Overlap for ConScape coarse-graining is summarised as follows:
#'
#' * `min_cells` is the minimum required overlap in cells:
#'   `floor(landmark / 2)` (half a coarse window on each side) and the
#'   user-requested `tile_trim` converted to cells are both honoured by
#'   taking their maximum.
#' * The actual overlap in cells (`overlap_cells`) is then rounded up to
#'   the nearest multiple of `landmark` so that landmark centres align
#'   across tiles.
#' * The effective overlap width in map units (`tile_trim` in the return
#'   list) is `overlap_cells * resx`. This value can differ from the
#'   user-supplied `tile_trim` but will always be at least as large.
#'
#' The interior polygons returned here are typically passed internally to
#' `tile_rast()` to construct extended tiles that include the overlap on
#' all sides. Those extended tiles are then sent to ConScape, and
#' [mosaic_conscape()] later trims and mosaics them using the same
#' landmark-aligned overlap.
#'
#' @return
#' A named list with components:
#'
#' \item{cs_tiles}{`SpatVector` of interior tile polygons covering the
#'   full extent of `r` without overlap.}
#' \item{tile_num}{Integer vector giving the tile index
#'   `seq_along(cs_tiles)`.}
#' \item{overlap_cells}{Integer giving the number of cells to use as
#'   overlap on each side of a tile (in both x and y directions). This
#'   is a multiple of `landmark`.}
#' \item{tile_trim}{Effective overlap width in map units
#'   (`overlap_cells * resx`). This is the value that should be passed to
#'   [mosaic_conscape()] and stored in `"ConScapeRtools_prep"`.}
#' \item{landmark}{Integer; the landmark window size `npix` used for
#'   alignment and coarse-graining.}
#'
#' @keywords internal
#' @noRd
#' @author Bill Peterman

make_tiles <- function(r,
                       tile_d,
                       tile_trim,
                       landmark = 10L) {

  resx <- res(r)[1]
  resy <- res(r)[2]
  if (!isTRUE(all.equal(resx, resy))) {
    warning("Non-square cells; using x-resolution.")
  }

  npix <- as.integer(landmark)
  ncols <- ncol(r)
  nrows <- nrow(r)

  # tile size in cells, multiple of npix
  tile_cells_raw <- round(tile_d / resx)
  tile_cells <- max(npix, ceiling(tile_cells_raw / npix) * npix)

  col_starts <- seq(1L, ncols, by = tile_cells)
  row_starts <- seq(1L, nrows, by = tile_cells)

  xmin_r <- xmin(r); xmax_r <- xmax(r)
  ymin_r <- ymin(r); ymax_r <- ymax(r)

  tiles_int <- vector("list", length(col_starts) * length(row_starts))
  idx <- 1L

  for (ic in seq_along(col_starts)) {
    for (ir in seq_along(row_starts)) {
      cs <- col_starts[ic]
      rs <- row_starts[ir]
      ce <- min(cs + tile_cells - 1L, ncols)
      re <- min(rs + tile_cells - 1L, nrows)

      x_min <- xmin_r + (cs - 1L) * resx
      x_max <- xmin_r + ce        * resx
      y_max <- ymax_r - (rs - 1L) * resy
      y_min <- ymax_r - re        * resy

      tiles_int[[idx]] <- ext(x_min, x_max, y_min, y_max)
      idx <- idx + 1L
    }
  }

  cs_tiles <- vect(do.call(c, lapply(tiles_int, as.polygons)))

  # npix = landmark (number of pixels in coarse-grain window)
  # tile_trim = user’s requested *minimum* overlap in map units

  min_cells <- max(
    floor(npix / 2),          # need at least half-window on each side
    ceiling(tile_trim / resx) # user minimum in cells
  )

  # force overlap_cells to be a multiple of npix
  overlap_cells <- ceiling(min_cells / npix) * npix
  overlap_dist  <- overlap_cells * resx      # effective tile_trim in map units


  list(
    cs_tiles      = cs_tiles,
    tile_num      = seq_along(cs_tiles),
    overlap_cells = overlap_cells,
    tile_trim     = overlap_dist,
    landmark      = npix
  )
}
