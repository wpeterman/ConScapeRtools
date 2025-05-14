#' Setup Julia for parallel processing
#'
#' @description Set up parallel processing to run ConScape with parallel.
#'
#' @param workers Number of parallel workers to create
#' @param julia_home Path to the Julia 'bin' directory
#' @param attempts Number of attempts to try and process files

#' @return NULL.

#' @keywords internal
#' @author Bill Peterman

julia_parallel <- function(workers,
                           jl_home,
                           attempts = 9){
  # conscape.jl <- system.file('data/conscape.jl', package = 'ConScapeRtools')
  cat("\nSetting up Julia...\n")
  # julia <- juliaSetup(JULIA_BINDIR = .jl_home)
  # juliaEval('using ConScape, SparseArrays, Statistics, Plots')
  # juliaEval(paste0('include("', system.file('data/conscape.jl', package = 'ConScapeRtools'), '")'))
  #
  # # iter <- paste0('-iter_', i)
  # juliaEval("conscape", .src_dir, .mov_dir, .target_dir, .out_dir,
  #           .hab_target[i], .hab_src[i], .mov_prob[i],
  #           .landmark, .theta, .exp_d, .NA_val, iter)


  test_func <- function(i,
                        x,
                        jl_home){
    if(!juliaSetupOk())
      Sys.setenv(JULIA_BINDIR = jl_home)
    if(!juliaSetupOk())
      stop("Check that the path to the Julia binary directory is correct")
    juliaEval('using ConScape, SparseArrays, Statistics, Plots')
    juliaEval(paste0('include("', system.file('data/conscape.jl', package = 'ConScapeRtools'), '")'))

    iter <- paste0('-iter_', i)

    juliaSqrt <- juliaFun("sqrt")
    juliaCall("map", juliaSqrt, x[i] )
  }

  x <- 0:100
  attempt <- 0
  while(isFALSE(exists('y')) & attempt <= attempts) {
    attempt <- attempt + 1
    cat(paste0("Attempt #", attempt, " to create parallel workers.\n\n"))
    plan(multisession, workers = workers)
    suppressMessages(try(y <- future_lapply(1:length(x),
                                            test_func, x, jl_home,
                                            future.seed = 1), silent = T))
  }
  plan(sequential) #### If we set to sequential here does that carry for the rest of the script?

  if(exists('y')){
    cat(paste0("Parallel workers successfully created.\n\n"))
    # return(z)
    #^^^ artifact???
  } else {
    stop("Failed to create parallel workers!!!\n\n")
  }
}


## JuliaConnectoR
# julia_parallel <- function(workers,
#                            julia_home,
#                            attempts = 10){
#   # cl <- parallelly::makeClusterPSOCK(workers)
#   cat("\nSetting up Julia...\n")
#   ConScapeR_setup(julia_home)
#
#   test_func <- function(i,
#                         x,
#                         julia_home){
#     ConScapeR_setup(julia_home)
#     juliaCall("sqrt", x[i] )
#   }
#
#   x <- 0:10
#   attempt <- 0
#   while(isFALSE(exists('y')) & attempt <= attempts) {
#     attempt <- attempt + 1
#     cat(paste0("Attempt #", attempt, " to create parallel workers.\n\n"))
#     plan(multisession, workers = workers)
#
#     suppressMessages(try(y <- future_lapply(1:length(x),
#                                             test_func, x, julia_home,
#                                             future.seed = 1), silent = T))
#   }
#   plan(sequential) #### If we set to sequential here does that carry for the rest of the script?
#
#   if(exists('y')){
#     cat(paste0("Parallel workers successfully created.\n\n"))
#
#   } else {
#     stop("Failed to create parallel workers!!!\n\n")
#   }
# }
