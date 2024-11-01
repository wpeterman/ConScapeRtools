#' Wrapper function to prepare data for `run_conscape`
#'
#' @description This function will run teh `make_tiles` and `tile_rast` functions to break large raster into tiles of specified dimensions. This is a convenience wrapper to process all layers in one function call.
#'
#' @param tile_d Dimensions (in meters) of tiles
#' @param tile_trim The amount of border to be trimmed from tiles after running ConScape (meters)
#' @param asc_dir Directory where .asc files of tiles will be written to. If NULL (Default), then files will be written to the temporary directory of the R session
#' @param r_target `SpatRaster` representing target qualities. May be the same as `r_src`
#' @param r_mov `SpatRaster` representing movement probabilities
#' @param r_src `SpatRaster` representing source qualities. May be the same as `r_target`
#' @param clear_dir Logical (Default = FALSE). Should existing files in the `asc_dir` be overwritten? This function must have an empty `asc_dir` to proceed
#' @param landmark The landmark value used for 'coarse_graining' with ConScape (Default = 10L). Used to determine which landscape tiles have data to be processed with ConScape
#' @return A named list of class `ConScapeRtools_prep` containing the `SpatVector`tiles created, the numeric identifier of tiles with usable data for ConScape, the path to the directories where .asc tiles were written, and the `tile_trim` value specified.
#' @details
#' The smaller the tiles created, the faster each can be processed. The width of the `tile_trim` parameter will depend upon the movement settings of your ConScape run. If there are obvious tiling edges and artifacts in your final surfaces, then `tile_trim` and potentially `tile_d` need to be increased.
#'
#'
#' @export
#' @examples examples/conscape_prep_example.R
#' @seealso [tile_rast()]  and [make_tiles()] for individual function calls
#' @author Bill Peterman

conscape_prep <- function(tile_d,
                          tile_trim,
                          asc_dir = NULL,
                          r_target,
                          r_mov,
                          r_src,
                          clear_dir = FALSE,
                          landmark = 10L) {
  r_ext <- ext(r_target)
  extnd <- r_ext + tile_trim
  r_e <- extend(r_target, extnd)
  r_e[is.na(r_e)] <- 0
  r_ext <- ext(r_e)

  e <- ext(r_ext[1], r_ext[1] + tile_d,
           r_ext[4] - tile_d, r_ext[4])
  t <- as.polygons(e, crs = crs(r_target))

  f_shift.x <- ext(t)

  shift <- tile_d - (2*tile_trim) - (2*res(r_target)[1])

  while(f_shift.x[2] < r_ext[2]){
    f_shift.x[1:2] <- (f_shift.x[1:2] + shift)

    if(f_shift.x[2] < r_ext[2]){
      t2 <- as.polygons(f_shift.x, crs = crs(r_target))
      t <- c(t, t2)
    } else {
      f_shift.x[2] <- r_ext[2]
      f_shift.x[1] <- (f_shift.x[2] -  tile_d)
      t2 <- as.polygons(f_shift.x, crs = crs(r_target))
      t <- c(t, t2)
    }
  } ## End while

  r1_poly <- vect(t)

  ## Shift y
  all_r <- r2_poly <- r1_poly
  while(ext(r2_poly)[3] > r_ext[3]){
    y_min <- (ext(r2_poly)[3] -  shift)

    if(y_min > r_ext[3]){
      r2_poly <- shift(r2_poly, dx = 0,
                       dy = (-1 * shift))

      all_r <- c(all_r, r2_poly)
    } else {

      r2_poly <- shift(r2_poly, dx = 0,
                       dy = r_ext[3] - ext(r2_poly)[3])
      all_r <- c(all_r, r2_poly)
    }
  } ## End y while

  all_r <- vect(all_r)

  ## Get centroids
  pts <- centroids(all_r)

  ## Focal buffer
  ## Alt method
  overlap <- (3*res(r_target)[1])
  s <- vect()
  for(i in 1:length(all_r)){
    fc <- all_r[i]
    fc_ <- crop(fc, ext(ext(fc)[] + c(tile_trim - overlap, -tile_trim + overlap,
                                      tile_trim - overlap, -tile_trim + overlap)))
    s <- c(s, fc_)
    # plot(fc); plot(s[i], add=T, border = 'red')
  }

  s <- vect(s)

  ## Break up raster
  r_list <- lapply(1:length(all_r), function(x)
    terra::mask(crop(r_e, ext(all_r[x])), all_r[x]))

  r_list.crop <- lapply(1:length(all_r), function(x)
    terra::mask(crop(r_e, ext(s[x])), s[x]))

  na_rast <- lapply(1:length(r_list.crop), function(x)
    minmax(r_list.crop[[x]])[1]) |> unlist()


  agg_list <- lapply(1:length(r_list), function(x)
    terra::aggregate(r_list[[x]], fact = I(floor(landmark*2)), fun = 'mean'))

  na_rast2 <- lapply(1:length(agg_list), function(x)
    minmax(agg_list[[x]])[1]) |> unlist()

  select_rast <- which(!is.nan(na_rast) & !is.nan(na_rast2))

  if(is.null(asc_dir)){
    write_dir <- paste0(tempdir(),"\\asc\\")
  } else {
    write_dir <- asc_dir
  }

  if(!dir.exists(write_dir))
    dir.create(write_dir, recursive = T)

  write_dir <- normalizePath(write_dir)

  if(length(list.files(write_dir)) > 0 & isFALSE(clear_dir)){
    stop("Files currently exist in the specified directory. Set `clear_dir = TRUE` to remove these files.")
  } else {
    unlink(paste0(write_dir, "\\*"), force = T)
  }

  for(i in 1:length(select_rast)){
    writeRaster(r_list[[select_rast[i]]],
                filename = paste0(write_dir, '\\r_', select_rast[i], '.asc'),
                NAflag = -9999,
                overwrite = TRUE)
  }

  out <- list(cs_tiles = all_r[select_rast],
              tile_num = select_rast,
              asc_dir = write_dir,
              tile_trim = tile_trim)

  mov_tile <- tile_rast(r = r_mov,
                        make_tiles = out,
                        out_dir = file.path(out$asc_dir, 'mov'),
                        clear_dir = T)

  src_tile <- tile_rast(r = r_src,
                        make_tiles = out,
                        out_dir = file.path(out$asc_dir, 'src'),
                        clear_dir = T)

  out2 <- list(cs_tiles = all_r[select_rast],
               tile_num = select_rast,
               asc_dir = write_dir,
               src = src_tile$asc_dir,
               target = write_dir,
               mov = mov_tile$asc_dir,
               tile_trim = tile_trim)

  class(out2) <- 'ConScapeRtools_prep'
  return(out2)
} ## End function

#' @rdname conscape_prep

