#' Run ConScape
#'
#' @description This function runs ConScape (optionally in parallel) over all created landscape tiles
#'
#' @param conscape_prep Object of class `ConScapeRtools_prep` created using the `conscape_prep` function. If NULL (Default), then `hab_target`, `hab_src`, and `mov_prob` must individually be specified.
#' @param out_dir Directory where final ConScape outputs will be written
#' @param hab_target Directory with habitat target tiles. Alternatively, a `SpatRaster`.
#' @param hab_src Directory with habitat source tiles. Alternatively, a `SpatRaster`.
#' @param mov_prob Directory with movement probability tiles. Alternatively, a `SpatRaster`.
#' @param clear_dir Should existing files in the `asc_dir` be overwritten? This function must have an empty `asc_dir` to proceed
#' @param landmark The landmark value used for 'coarse_graining' with ConScape (Default = 10L). Used to determine which landscape tiles have data to be processed with ConScape
#' @param theta Parameter to control the amount of randomness in the paths. As theta approaches 0 movement is random, whereas theta approaching infinity is optimal movement. (Default = 0.01)
#' @param exp_d Numerator of the exponential decay function used with the distance transformation of the movement grid (Default = 50)
#' @param NA_val Value to assign to NA cells to ensure ConScape can run (Default = 1e-50)
#' @param mosaic Logical. Default = TRUE. The tiles from the ConScape run will be combined into a single raster.
#' @param jl_home Path to the `bin` directory where Julia is installed
#' @param parallel Logical. If FALSE, processing will not be done in parallel.
#' @param parallel_R If `TRUE`, parallelization will be done within R.
#' @param workers If `parallel = TRUE`, provide the number of parallel workers to create. Default = 2
#' @param distributed Logical. If `parallel = TRUE`, you can choose parallelize using Julia distributed processing (`TRUE`) or threading (`FALSE`)
#' @param progress Logical. If `TRUE`, progress of processing will be reported.
#' @param parallel_R Logical. If `TRUE`, parallel processign will be done from R.
#' @return A named list with the betweenness and functional connectivity rasters as well as the directories where output tiles were written.
#' @details
#' In most instances, it will be easiest to prepare data for analysis using the `conscape_prep` function. Provide the object created from `conscape_prep` to the `conscape_prep` parameter of `run_conscape`. Doing this eliminates the need to manually specify `landmark`, `hab_target`, `hab_src`, or `mov_prob`.
#'
#' If `parallel = TRUE` and `progress = TRUE`, progress updates will be reported, but may be inconsistent.

#' @importFrom JuliaConnectoR juliaImport juliaEval juliaSetupOk juliaGet juliaCall
#' @export
#' @example examples/run_conscape_example.R
#' @author Bill Peterman

run_conscape <- function(conscape_prep = NULL,
                         out_dir,
                         hab_target = NULL,
                         hab_src = NULL,
                         mov_prob = NULL,
                         clear_dir = TRUE,
                         landmark = 10L,
                         theta = 0.01,
                         exp_d = 150,
                         NA_val = 1e-50,
                         mosaic = TRUE,
                         jl_home,
                         parallel = FALSE,
                         workers = 1,
                         distributed = FALSE,
                         progress = TRUE,
                         parallel_R = FALSE){
  stopJulia()
  if(parallel & isFALSE(parallel_R)){
    Sys.setenv(JULIA_NUM_THREADS = workers)
  }

  if(!juliaSetupOk())
    Sys.setenv(JULIA_BINDIR = jl_home)
  if(!juliaSetupOk()){
    stop("Check that the path to the Julia binary directory is correct")
  }

  if(dir.exists(out_dir)){
    if(length(list.files(out_dir)) > 0 & isFALSE(clear_dir)){
      stop("\nFiles currently exist in `out_dir`.\nEither manually delete them, or set `clear_dir = TRUE`")
    } else {
      unlink(out_dir,
             recursive = T,
             force = T)
    }
  } else {
    dir.create(out_dir, recursive = TRUE)
  }

  single_rast <- FALSE

  if(!is.null(conscape_prep) & inherits(conscape_prep, 'ConScapeRtools_prep')){
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

  if(inherits(hab_target,'SpatRaster') |
     inherits(hab_src,'SpatRaster') |
     inherits(mov_prob,'SpatRaster')){
    if(inherits(hab_target,'SpatRaster')){
      single_rast <- TRUE
      iter <- 1
      hab_target_file <- basename(terra::sources(hab_target))
      if(isFALSE(dir.exists(hab_target_file))){
        hab_target_write <- paste0(tempdir(),'\\hab_target.asc')
        suppressWarnings(try(writeRaster(hab_target, hab_target_write, overwrite = T, NAflag = -9999), silent = T))
        hab_target_dir <- dirname(hab_target_write)
      }
    }

    if(inherits(hab_src,'SpatRaster')){
      hab_src_file <- basename(terra::sources(hab_src))
      if(isFALSE(dir.exists(hab_src_file))){
        hab_src_write <- paste0(tempdir(),'\\hab_src.asc')
        suppressWarnings(try(writeRaster(hab_src, hab_src_write, overwrite = T, NAflag = -9999), silent = T))
        hab_src_dir <- dirname(hab_src_write)
      }
    }

    if(inherits(mov_prob,'SpatRaster')){
      mov_prob_file <- basename(terra::sources(mov_prob))
      if(isFALSE(dir.exists(mov_prob_file))){
        mov_prob_write <- paste0(tempdir(),'\\mov_prob.asc')
        suppressWarnings(try(writeRaster(mov_prob, mov_prob_write, overwrite = T, NAflag = -9999), silent = T))
        mov_prob_dir <- dirname(mov_prob_write)
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
  # Parallel ----------------------------------------------------------------
  if(isTRUE(parallel)) {

    # **In R ------------------------------------------------------------------

    if(parallel_R){
      plan(sequential)

      plan(multisession, workers = workers)

      suppressWarnings({
        cs_out <- future_lapply(1:length(hab_target), function(i) {
          Sys.setenv(JULIA_BINDIR = jl_home)
          conscape_file <- normalizePath(system.file('extdata', 'conscape.jl', package = 'ConScapeRtools'), winslash = "/")
          conscape_file <- paste0('include("', conscape_file, '")')

          invisible(juliaEval(conscape_file))

          cs_par.func(i,
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
                      .NA_val = NA_val)
        })
      })
      cs_out <- do.call(c, cs_out)
      plan(sequential)

    } else{
      if(!juliaSetupOk() && isFALSE(distributed)){
        Sys.setenv(JULIA_NUM_THREADS = workers)
        Sys.setenv(JULIA_BINDIR = jl_home)
        # juliaCall("addprocs", workers)
      }
      if(!juliaSetupOk())
        stop("Check that the path to the Julia binary directory is correct")


      # Threaded ----------------------------------------------------------------
      if(isFALSE(distributed)){
        juliaEval('println("Available threads: ", Threads.nthreads())')
        conscape_batch_file <- normalizePath(system.file('extdata', 'conscape_batch.jl', package = 'ConScapeRtools'), winslash = "/")
        conscape_batch_file <- paste0('include("', conscape_batch_file, '")')

        conscape_file <- normalizePath(system.file('extdata', 'conscape.jl', package = 'ConScapeRtools'), winslash = "/")
        conscape_file <- paste0('include("', conscape_file, '")')

        invisible(juliaEval(conscape_file))
        invisible(juliaEval(conscape_batch_file))

        max_retries <- 5L
        invisible(cs_out <- juliaCall("conscape_batch",
                                      src_dir, mov_dir, target_dir, out_dir,
                                      hab_target, hab_src, mov_prob,
                                      landmark, theta, exp_d, NA_val,
                                      max_retries, progress))

      } else {

        # Distributed -------------------------------------------------------------
        invisible(juliaEval("using Distributed"))
        invisible(juliaEval(paste0("addprocs(",workers,")")))

        conscape_batch_file <- normalizePath(system.file('extdata', 'conscape_batch_distributed2.jl', package = 'ConScapeRtools'), winslash = "/")
        conscape_batch_file <- paste0('include("', conscape_batch_file, '")')

        conscape_file <- normalizePath(system.file('extdata', 'conscape.jl', package = 'ConScapeRtools'), winslash = "/")
        conscape_file <- paste0('@everywhere include("', conscape_file, '")')

        invisible(juliaEval(conscape_file))
        invisible(juliaEval(conscape_batch_file))

        max_retries <- 5L
        invisible(cs_out <- juliaCall("conscape_batch_distributed",
                                      src_dir, mov_dir, target_dir, out_dir,
                                      hab_target, hab_src, mov_prob,
                                      landmark, theta, exp_d, NA_val,
                                      max_retries, progress))
      }
    }


    ## Check for success
    btwn_files <- length(list.files(normalizePath(file.path(out_dir, "btwn")),
                                    pattern = "\\.asc$"))
    fcon_files <- length(list.files(normalizePath(file.path(out_dir, "fcon")),
                                    pattern = "\\.asc$"))


    if(isTRUE((length(hab_target) != btwn_files) | (length(hab_target) != fcon_files))){
      # attempt <- attempt + 1
      # cat(paste0("\n\nParallel execution failed! Trying again...attempt #", attempt))
      warning("\n\nParallel execution failed for some or all tiles!\nInspect results carefully.\n")

      # DEBUG -------------------------------------------------------------------
      # browser()

    } ## End while loop

  } else {  ## End parallel
    # Serial ----------------------------------------------------------------

    if(!juliaSetupOk())
      Sys.setenv(JULIA_BINDIR = jl_home)
    if(!juliaSetupOk())
      stop("Check that the path to the Julia binary directory is correct")

    # Load dependencies
    juliaEval('using ConScape, SparseArrays, Statistics')

    # Include files with proper path handling
    conscape_file <- normalizePath(system.file('extdata', 'conscape.jl', package = 'ConScapeRtools'), winslash = "/")
    conscape_file <- paste0('include("', conscape_file, '")')

    invisible(juliaEval(conscape_file))

    if(single_rast){
      ## Julia order
      iter <- ""

      suppressMessages({cs_out <- juliaCall('conscape',
                                            normalizePath(hab_src_dir, winslash = "\\"),
                                            normalizePath(mov_prob_dir, winslash = "\\") ,
                                            normalizePath(hab_target_dir, winslash = "\\"),
                                            out_dir,
                                            "hab_target.asc", "hab_src.asc", "mov_prob.asc",
                                            # hab_target_file, hab_src_file, mov_prob_file,
                                            landmark, theta, exp_d, NA_val, iter)})
    } else {


      iters <- max(length(hab_target),1)

      if(isTRUE(progress)){
        pb <- txtProgressBar(min = 0, max = iters,
                             initial = 0, char = "+", style = 3)
      }

      for(i in 1:iters){

        iter <- paste0('-iter_', i)
        r_target <- hab_target[i]
        r_source <- hab_src[i]
        r_res <- mov_prob[i]

        suppressMessages({cs_out <- juliaCall('conscape',
                                              src_dir, mov_dir, target_dir, out_dir,
                                              r_target, r_source, r_res,
                                              landmark, theta, exp_d, NA_val, iter)})

        if(isTRUE(progress)){
          setTxtProgressBar(pb,i)
        }
      } ## End for loop
    } ## End single ifelse
    if(exists("pb")) close(pb)
  } ## End parallel ifelse

  out <- list(outdir_btwn = normalizePath(file.path(out_dir, "btwn")),
              outdir_fcon = normalizePath(file.path(out_dir, "fcon")))

  if(isTRUE(mosaic) & isFALSE(single_rast)){
    btwn <- mosaic_conscape(out_dir = out$outdir_btwn,
                            tile_trim = tile_trim,
                            method = 'mosaic')
    fcon <- mosaic_conscape(out_dir = out$outdir_fcon,
                            tile_trim = tile_trim,
                            method = 'mosaic')

    if(!is.null(conscape_prep) & inherits(conscape_prep, 'ConScapeRtools_prep')){
      crs(btwn) <- crs(fcon) <- crs(target_mask)
      target_mask <- crop(target_mask, fcon)
      target_mask[target_mask == 0] <- NA
      btwn[is.na(target_mask)] <- NA
      fcon[is.na(target_mask)] <- NA
    }

    out <- list(btwn = btwn,
                fcon = fcon,
                outdir_btwn = normalizePath(file.path(out_dir, "btwn")),
                outdir_fcon = normalizePath(file.path(out_dir, "fcon"))#,
                # cs_log = cs_out
                )
    class(out) <- "ConScapeResults"
  } else {
    out <- list(outdir_btwn = normalizePath(file.path(out_dir, "btwn")),
                outdir_fcon = normalizePath(file.path(out_dir, "fcon"))#,
                # cs_log = cs_out
                )
    class(out) <- "ConScapeResults"

    if(isTRUE(single_rast)){
      btwn <- rast(list.files(normalizePath(file.path(out_dir, "btwn")),
                             full.names = T))
      fcon <- rast(list.files(normalizePath(file.path(out_dir, "fcon")),
                             full.names = T))
      crs(btwn) <- crs(fcon) <- crs(mov_prob)

      out <- rast(list(btwn = btwn,
                       fcon = fcon))
    }
  }

  # juliaEval("GC.gc()")
  stopJulia()
  gc()
  return(out)
} ## End function

#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom future plan sequential multisession
#' @importFrom future.apply future_lapply
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
  iter <- paste0('-iter_', i)
  invisible({cs_out <- juliaCall('conscape',
                                 .src_dir, .mov_dir, .target_dir, .out_dir,
                                 .hab_target[i], .hab_src[i], .mov_prob[i],
                                 .landmark, .theta, .exp_d, .NA_val, iter)})
  # juliaEval("GC.gc()")
}
