#' Run ConScape
#'
#' @description This function runs ConScape (optionally in parallel) over all created landscape tiles
#'
#' @param conscape_prep Object of class `ConScapeRtools_prep` created using the `conscape_prep` function. If NULL (Default), then `hab_target`, `hab_src`, and `mov_prob` must individually be specified.
#' @param out_dir Directory where final ConScape outputs will be written
#' @param hab_target Directory with habitat target tiles
#' @param hab_src Directory with habitat source tiles
#' @param mov_prob Directory with movement probability tiles
#' @param clear_dir Should existing files in the `asc_dir` be overwritten? This function must have an empty `asc_dir` to proceed
#' @param landmark The landmark value used for 'coarse_graining' with ConScape (Default = 10L). Used to determine which landscape tiles have data to be processed with ConScape
#' @param theta Parameter to control the amount of randomness in the paths. As theta approaches 0 movement is random, whereas theta approaching infinity is optimal movement. (Default = 0.01)
#' @param exp_d Numerator of the exponential decay function used with the distance transformation of the movement grid (Default = 50)
#' @param NA_val Value to assign to NA cells to ensure ConScape can run (Default = 1e-50)
#' @param mosaic Logical. Default = TRUE. The tiles from the ConScape run will be combined into a single raster.
#' @param jl_home Path to the `bin` directory where Julia is installed
#' @param parallel Logical. If FALSE, processing will not be done in parallel.
#' @param workers If `parallel = TRUE`, provide the number of parallel workers to create. Default 0.5*Number of available cores
#' @param end_multisession Logical. (Default = FALSE)
#' @return A named list with the betweenness and functional connectivity rasters as well as the directories where output tiles were written.
#' @details
#' In most instances, it will be easiest to prepare data for analysis using the `conscape_prep` function. Provide the object created from `conscape_prep` to the `conscape_prep` parameter of `run_conscape`. Doing this eliminates the need to manually specify `landmark`, `hab_target`, `hab_src`, or `mov_prob`.

#' @importFrom JuliaConnectoR juliaImport juliaEval juliaSetupOk juliaGet juliaCall
#' @export
#' @example examples/run_conscape_example.R
#' @author Bill Peterman

run_conscape <- function(conscape_prep = NULL,
                         out_dir,
                         hab_target = NULL,
                         hab_src = NULL,
                         mov_prob = NULL,
                         landmark = 10L,
                         theta = 0.01,
                         exp_d = 150,
                         NA_val = 1e-50,
                         mosaic = TRUE,
                         jl_home,
                         parallel = FALSE,
                         workers = availableCores()/2,
                         end_multisession = FALSE){
  if(!juliaSetupOk())
    Sys.setenv(JULIA_BINDIR = jl_home)
  if(!juliaSetupOk()){
    stop("Check that the path to the Julia binary directory is correct")
  }
  suppressMessages({juliaEval("Base.redirect_stdout(devnull); Base.redirect_stderr(devnull)")})

  juliaEval('using ConScape, SparseArrays, Statistics, Plots')
  juliaEval(paste0('include("', system.file('data/conscape.jl', package = 'ConScapeRtools'), '")'))

  if(dir.exists(out_dir)){
    unlink(out_dir,
           recursive = T,
           force = T)
  }

  single_rast <- FALSE

  if(!is.null(conscape_prep) & class(conscape_prep) == 'ConScapeRtools_prep'){
    target_dir <- conscape_prep$target
    src_dir <- conscape_prep$src
    mov_dir <- conscape_prep$mov
    hab_target <- list.files(target_dir, pattern = "\\.asc$")
    hab_src <- list.files(src_dir, pattern = "\\.asc$")
    mov_prob <- list.files(mov_dir, pattern = "\\.asc$")
    landmark <- conscape_prep$landmark
    target_mask <- rast(file.path(conscape_prep$asc_dir, 'mask', "mask.asc"))
    tile_trim <- conscape_prep$tile_trim
  }

  if(class(hab_target)[[1]] == 'SpatRaster' |
     class(hab_src)[[1]] == 'SpatRaster' |
     class(mov_prob)[[1]] == 'SpatRaster'){
    if(class(hab_target)[[1]] == 'SpatRaster'){
      single_rast <- TRUE
      hab_target_dir <- basename(terra::sources(hab_target))
      if(isFALSE(dir.exists(hab_target_dir))){
        hab_target_dir <- paste0(tempdir(),'\\hab_target.asc')
        suppressWarnings(try(writeRaster(hab_target, hab_target_dir, overwrite = T, NAflag = -9999), silent = T))
        hab_target <- hab_target_dir
      }
    }

    if(class(hab_src)[[1]] == 'SpatRaster'){
      hab_src_dir <- basename(terra::sources(hab_src))
      if(isFALSE(dir.exists(hab_src_dir))){
        hab_src_dir <- paste0(tempdir(),'\\hab_src.asc')
        suppressWarnings(try(writeRaster(hab_src, hab_src_dir, overwrite = T, NAflag = -9999), silent = T))
        hab_src <- hab_src_dir
      }
    }

    if(class(mov_prob)[[1]] == 'SpatRaster'){
      mov_prob_dir <- basename(terra::sources(mov_prob))
      if(isFALSE(dir.exists(mov_prob_dir))){
        mov_prob_dir <- paste0(tempdir(),'\\mov_prob.asc')
        suppressWarnings(try(writeRaster(mov_prob, mov_prob_dir, overwrite = T, NAflag = -9999), silent = T))
        mov_prob <- mov_prob_dir
      }
    }


  } else {
    if(isFALSE(file.exists(hab_target))){
      stop("Specify either a SpatRaster object or path to directory containing *.asc files")
    }

    if(is.null(conscape_prep)){
      target_dir <- hab_target
      src_dir <- hab_src
      mov_dir <- mov_prob
      hab_target <- ifelse(grepl(".asc", basename(hab_target)), hab_target, list.files(hab_target, pattern = "\\.asc$"))
      hab_src <- ifelse(grepl(".asc", basename(hab_src)), hab_src, list.files(hab_src, pattern = "\\.asc$"))
      mov_prob <- ifelse(grepl(".asc", basename(mov_prob)), mov_prob, list.files(mov_prob, pattern = "\\.asc$"))
    }
  }
  if(isTRUE(parallel)) {

    # Parallel ----------------------------------------------------------------
    if(nbrOfWorkers() != workers){
      plan(sequential)
      options(connectionObserver = NULL)  # Disables connection tracking
      plan(multisession, workers = workers)
      options(future.globals.packages = c("stats"))

      # Ensure Julia is initialized on each worker

      future_lapply(1:workers, function(x) {
        if(!juliaSetupOk())
          Sys.setenv(JULIA_BINDIR = jl_home)
        if(!juliaSetupOk())
          stop("Check that the path to the Julia binary directory is correct")
        suppressMessages({juliaEval('using ConScape, SparseArrays, Statistics, Plots')})
        juliaEval(paste0('include("', system.file('data/conscape.jl', package = 'ConScapeRtools'), '")'))
        juliaEval("Base.redirect_stdout(devnull); Base.redirect_stderr(devnull)")
        NULL  # Just initialize, no return needed
      }
      )
    }

    suppressWarnings({
    try(cs_out <- future_lapply(1:length(hab_target),
                                FUN = cs_par.func,
                                .jl_home = jl_home,
                                .hab_target = hab_target,
                                .hab_src = hab_src,
                                .mov_prob = mov_prob,
                                .src_dir = src_dir,
                                .mov_dir = mov_dir,
                                .target_dir = target_dir,
                                .out_dir = out_dir,
                                .landmark = landmark,
                                .theta = theta,
                                .exp_d = exp_d,
                                .NA_val = NA_val,
                                future.seed = TRUE
    ),
    silent = T)
      })## End future_lapply


    ## Check for success
    btwn_files <- length(list.files(normalizePath(file.path(out_dir, "btwn")),
                                    pattern = "\\.asc$"))
    fcon_files <- length(list.files(normalizePath(file.path(out_dir, "fcon")),
                                    pattern = "\\.asc$"))
    attempt <- 1

    while(isTRUE((length(hab_target) != btwn_files) | (length(hab_target) != fcon_files))){
      attempt <- attempt + 1
      cat(paste0("\n\nParallel execution failed! Trying again...attempt #", attempt))


      # DEBUG -------------------------------------------------------------------
      browser()

      julia_parallel(workers = workers,
                     # workers = availableCores()/2,
                     jl_home = jl_home)

      try(cs_out <- future_lapply(1:length(hab_target),
                                  FUN = cs_par.func,
                                  .jl_home = jl_home,
                                  .hab_target = hab_target,
                                  .hab_src = hab_src,
                                  .mov_prob = mov_prob,
                                  .src_dir = src_dir,
                                  .mov_dir = mov_dir,
                                  .target_dir = target_dir,
                                  .out_dir = out_dir,
                                  .landmark = landmark,
                                  .theta = theta,
                                  .exp_d = exp_d,
                                  .NA_val = NA_val,
                                  future.seed = TRUE
      ),
      silent = T) ## End future_lapply

      ## Check for success
      btwn_files <- length(list.files(normalizePath(file.path(out_dir, "btwn")),
                                      pattern = "\\.asc$"))
      fcon_files <- length(list.files(normalizePath(file.path(out_dir, "fcon")),
                                      pattern = "\\.asc$"))
    } ## End while loop


    if(isTRUE(end_multisession)){
    plan(sequential)
    }

  } else {  ## End parallel

    iters <- max(length(hab_target),1)

    pb <- txtProgressBar(min = 0, max = iters,
                         initial = 0, char = "+", style = 3)

    for(i in 1:iters){

      iter <- paste0('-iter_', i)
      r_target <- hab_target[i]
      r_source <- hab_src[i]
      r_res <- mov_prob[i]

      suppressMessages({juliaCall('conscape',
                                  src_dir, mov_dir, target_dir, out_dir,
                                  r_target, r_source, r_res,
                                  landmark, theta, exp_d, NA_val, iter)})

      setTxtProgressBar(pb,i)
    } ## End for loop
    close(pb)
  } ## End ifelse

  out <- list(outdir_btwn = normalizePath(file.path(out_dir, "btwn")),
              outdir_fcon = normalizePath(file.path(out_dir, "fcon")))

  if(isTRUE(mosaic) & isFALSE(single_rast)){
    btwn <- mosaic_conscape(out_dir = out$outdir_btwn,
                            tile_trim = tile_trim,
                            method = 'mosaic')
    fcon <- mosaic_conscape(out_dir = out$outdir_fcon,
                            tile_trim = tile_trim,
                            method = 'mosaic')

    if(!is.null(conscape_prep) & class(conscape_prep) == 'ConScapeRtools_prep'){
      crs(btwn) <- crs(fcon) <- crs(target_mask)
      target_mask <- crop(target_mask, fcon)
      target_mask[target_mask == 0] <- NA
      btwn[is.na(target_mask)] <- NA
      fcon[is.na(target_mask)] <- NA
    }

    out <- list(btwn = btwn,
                fcon = fcon,
                outdir_btwn = normalizePath(file.path(out_dir, "btwn")),
                outdir_fcon = normalizePath(file.path(out_dir, "fcon")))
  } else {
    out <- list(outdir_btwn = normalizePath(file.path(out_dir, "btwn")),
                outdir_fcon = normalizePath(file.path(out_dir, "fcon")))

    if(isTRUE(single_rast)){
      out <- rast(list(btwn = rast(list.files(normalizePath(file.path(out_dir, "btwn")),
                                              full.names = T)),
                       fcon = rast(list.files(normalizePath(file.path(out_dir, "fcon")),
                                              full.names = T))))
    }
  }
  juliaEval("Base.redirect_stdout(Base.stdout); Base.redirect_stderr(Base.stderr)")

  return(out)
} ## End function

#' @importFrom parallelly availableCores
#' @importFrom future plan multisession sequential nbrOfWorkers
#' @importFrom future.apply future_lapply
#' @importFrom utils setTxtProgressBar txtProgressBar
NULL

cs_par.func <- function(i,
                        .jl_home,
                        .hab_target,
                        .hab_src,
                        .mov_prob,
                        .src_dir,
                        .mov_dir,
                        .target_dir,
                        .out_dir,
                        .landmark,
                        .theta,
                        .exp_d,
                        .NA_val){

  # if(!juliaSetupOk())
  #   Sys.setenv(JULIA_BINDIR = .jl_home)
  # if(!juliaSetupOk())
  #   stop("Check that the path to the Julia binary directory is correct")
  # suppressMessages({juliaEval('using ConScape, SparseArrays, Statistics, Plots')})
  # juliaEval(paste0('include("', system.file('data/conscape.jl', package = 'ConScapeRtools'), '")'))
  # juliaEval("Base.redirect_stdout(devnull); Base.redirect_stderr(devnull)")

  iter <- paste0('-iter_', i)

  suppressWarnings({juliaCall('conscape',
                              .src_dir, .mov_dir, .target_dir, .out_dir,
                              .hab_target[i], .hab_src[i], .mov_prob[i],
                              .landmark, .theta, .exp_d, .NA_val, iter)})
  juliaEval("GC.gc()")
  # juliaEval("Base.redirect_stdout(Base.stdout); Base.redirect_stderr(Base.stderr)")
}
