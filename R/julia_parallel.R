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
                           julia_home,
                           attempts = 10){
  # conscape.jl <- system.file('data/conscape.jl', package = 'ConScapeRtools')
  cat("\nSetting up Julia...\n")
  julia_setup(julia_home)
  julia_library("ConScape")
  julia_library("SparseArrays")
  julia_library("Statistics")
  julia_library("Plots")
  julia_source(system.file('data/conscape.jl', package = 'ConScapeRtools'))

  test_func <- function(i,
                        x,
                        julia_home,
                        conscape.jl){
    julia_setup(julia_home)
    julia_library("ConScape")
    julia_library("SparseArrays")
    julia_library("Statistics")
    julia_library("Plots")
    julia_source(system.file('data/conscape.jl', package = 'ConScapeRtools'))
    julia_call("sqrt", x[i] )
  }

  x <- 0:100
  attempt <- 0
  while(isFALSE(exists('y')) & attempt <= attempts) {
    attempt <- attempt + 1
    cat(paste0("Attempt #", attempt, " to create parallel workers.\n\n"))
    plan(multisession, workers = workers)
    suppressMessages(try(y <- future_lapply(1:length(x),
                                            test_func, x, julia_home,
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
