#' Break a raster into overlapped tiles using a predefined tiling scheme
#'
#' @description
#' This function takes a `SpatRaster` and a tiling scheme created internally by `make_tiles()`
#' and writes out overlapped raster tiles for use with ConScape. The tiling layout
#' (interior tile footprints and overlap width) is defined entirely by `make_tiles`,
#' ensuring consistent alignment across rasters and with the landmark grid used for
#' coarse graining.
#'
#' @param r `SpatRaster` to be broken into tiles.
#' @param tiles A list returned by `make_tiles()`, containing at minimum
#'   `cs_tiles` (interior tile polygons), `tile_num` (integer tile IDs),
#'   and `overlap_cells` (overlap width in cells).
#' @param out_dir Directory where `.asc` tiles will be written. If `NULL`,
#'   tiles are written to a temporary directory for the current R session.
#' @param clear_dir Logical; if `FALSE` (default) and `out_dir` is not empty,
#'   the function stops rather than overwriting existing files. If `TRUE`,
#'   existing files in `out_dir` are removed before writing tiles.
#' @param progress Logical; if `TRUE`, a text progress bar is displayed while
#'   tiles are written.
#'
#' @details
#' For each interior tile polygon in `make_tiles$cs_tiles`, the function expands
#' the extent by `make_tiles$overlap_cells` on all sides (in units of raster
#' cells), crops the extended raster (`extend()` is called once on `r` using the
#' same overlap), and writes the result as an ESRI ASCII grid with filenames of
#' the form `"r_<tile_id>.asc"`. Because the tiling scheme is generated in cell
#' space by `make_tiles()`, multiple rasters tiled with the same `make_tiles`
#' object share identical tile boundaries and overlap, which is required for
#' consistent landmark placement in ConScape coarse graining.
#' @return
#' A named list with:
#' \item{asc_dir}{Absolute path to the directory where `.asc` tiles were written.}
#' @keywords internal
#' @noRd
tile_rast <- function(r,
                      tiles,
                      out_dir,
                      clear_dir = FALSE,
                      progress = FALSE) {

  cs_tiles      <- tiles$cs_tiles
  overlap_cells <- tiles$overlap_cells

  resx <- res(r)[1]
  resy <- res(r)[2]

  # extend raster by overlap_cells cells on all sides
  r_ext <- extend(r, c(overlap_cells, overlap_cells,
                       overlap_cells, overlap_cells))

  # output directory
  if (is.null(out_dir)) {
    write_dir <- file.path(tempdir(), "asc")
  } else {
    write_dir <- out_dir
  }
  if (!dir.exists(write_dir)) dir.create(write_dir, recursive = TRUE)
  write_dir <- normalizePath(write_dir)

  if (length(list.files(write_dir)) > 0 && !clear_dir) {
    stop("Files currently exist in the specified directory. Set `clear_dir = TRUE` to remove these files.")
  } else {
    unlink(file.path(write_dir, "*"), force = TRUE)
  }

  n_tiles <- length(cs_tiles)  # number of geometries

  # build extended extents for each interior tile
  tiles_ext <- lapply(seq_len(n_tiles), function(i) {
    te <- cs_tiles[i]  # subset SpatVector
    ext(
      xmin(te) - overlap_cells * resx,
      xmax(te) + overlap_cells * resx,
      ymin(te) - overlap_cells * resy,
      ymax(te) + overlap_cells * resy
    )
  })

  # optional progress bar
  if (isTRUE(progress)) {
    pb <- txtProgressBar(min = 0, max = n_tiles, style = 3)
  }

  for (idx in seq_len(n_tiles)) {
    if (exists("pb", inherits = FALSE)) {
      setTxtProgressBar(pb, idx,
                        title = paste0("Processing ", basename(write_dir), " rasters..."))
    }

    r_tile <- crop(r_ext, tiles_ext[[idx]])
    writeRaster(
      r_tile,
      filename = file.path(write_dir, paste0("r_", idx, ".asc")),
      NAflag   = -9999,
      overwrite = TRUE
    )
  }

  if (exists("pb", inherits = FALSE)) close(pb)

  list(
    asc_dir = write_dir
  )
}
