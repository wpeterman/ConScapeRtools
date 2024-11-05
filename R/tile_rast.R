#' Break raster into small tiles
#'
#' @description This uses objects created by the [ConScapeRtools:make_tiles()] function to create identical tiles from other `SpatRaster` objects.
#'
#' @param r Raster file as `SpatRaster` to be broken up into smaller tiles
#' @param make_tiles Object created from  [make_tiles()] function
#' @param out_dir Directory where .asc files of tiles will be written. This should be different directory from the `make_tiles` `asc_dir`
#' @param clear_dir Should existing files in the `out_dir` be overwritten? This function must have an empty `out_dir` to proceed
#' @return A named list with the path to the directory where .asc tiles were written

#' @export
#' @example examples/tile_rast_example.R
#' @seealso [make_tile()] for initial creation of tiles
#' @author Bill Peterman

tile_rast <- function(r,
                      make_tiles,
                      out_dir,
                      clear_dir = FALSE){
  ## Extend raster
  extnd <- ext(r) + trim
  r_e <- extend(r, extnd)
  r_e[is.na(r_e)] <- 0

  cs_tiles <- make_tiles$cs_tiles
  tile_num <- make_tiles$tile_num

  ## Break up raster
  r_list <- lapply(1:length(cs_tiles), function(x)
    terra::mask(crop(r_e, ext(cs_tiles[x])), cs_tiles[x]))

  select_rast <- tile_num

  if(is.null(out_dir)){
    write_dir <- paste0(tempdir(),"\\asc\\")
  } else {
    write_dir <- out_dir
  }

  if(!dir.exists(write_dir))
    dir.create(write_dir, recursive = T)

  write_dir <- normalizePath(write_dir)

  if(length(list.files(write_dir)) > 0 & isFALSE(clear_dir)){
    stop("Files currently exist in the specified directory. Set `clear_dir = TRUE` to remove these files.")
  } else {
    unlink(paste0(write_dir, "\\*"), force = T)
  }

  cnt <- 0
  for(i in select_rast){
    cnt <- cnt + 1
    writeRaster(r_list[[cnt]],
                filename = paste0(write_dir, '\\r_', i, '.asc'),
                NAflag = -9999,
                overwrite = TRUE)
  }

  out <- list(asc_dir = write_dir)
}
