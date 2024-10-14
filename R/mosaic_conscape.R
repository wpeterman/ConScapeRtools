#' Mosaic ConScape tiles into single contiguous raster
#'
#' @description After running [ConScapeRtools:run_conscape()], this function will reassemble tiles into a single `SpatRaster`.
#'
#' @param out_dir Directory where [ConScapeRtools:run_conscape()] results were written
#' @param tile_trim Provide `tile_trim` value used when running [ConScapeRtools:make_tiles()]

#' @return `SpatRaster`

#' @export
#' @examples examples/mosaic_conscape_example.R
#' @seealso [make_tiles()] and [run_conscape()] \cr
#' @author Bill Peterman

mosaic_conscape <- function(asc_dir,
                            tile_trim) {
  asc_files <- list.files(asc_dir,
                          pattern = "\\.asc$",
                          full.names = T)

  asc <- basename(asc_files) |>
    strsplit(split = "_", fixed = T)
  asc <- data.frame(matrix(unlist(asc), nrow=length(asc), byrow=TRUE))[,2]
  asc <- (as.numeric(gsub('.asc', '', asc)))

  crop_rast <- vector('list', length(asc_files))
  for(i in 1:length(asc_files)){
    r_ <- rast(asc_files[i])
    r_ext <- ext(r_)[1:4]

    e <- ext(r_ext[1] + tile_trim, r_ext[2] - tile_trim,
             r_ext[3] + tile_trim, r_ext[4] - tile_trim)
    t <- as.polygons(e)

    crop_rast[[i]] <- terra::mask(crop(r_, ext(t)), t)
  }

  rc <- sprc(crop_rast)
  r_mosaic <- mosaic(rc, fun = 'mean')
  return(r_mosaic)
}
