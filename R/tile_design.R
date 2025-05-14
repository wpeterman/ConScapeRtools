#' Function to determine suitable tiling and analysis parameters
#'
#' @description Used to determine the tile size, trim width, exp_d decay value given properties of data inputs.
#'
#' @param r_mov `SpatRaster` representing movement probabilities
#' @param r_source `SpatRaster` representing landscape source quality (Default = NULL)
#' @param r_target `SpatRaster` representing landscape target quality (Default = NULL)
#' @param max_d Presumed maximum movement distance (measured in meters) through ideal habitat
#' @param theta Parameter to control the amount of randomness in the paths. As theta approaches 0 movement is random, whereas theta approaching infinity is optimal movement. (Default = 0.1)
#' @param jl_home Path to the bin directory where Julia resides
#' @return A named list with parameters for creating tiles with [conscape_prep()] or [tile_rast()]

#' @export
#' @example examples/tile_design_example.R
#' @author Bill Peterman
#' @details
#' `r_mov` SpatRaster must be in a projected coordinate reference system.
#' Relies on functions borrowed from the `ConScapeR` package (https://github.com/ConScape/ConScapeR).
#'
#' @rdname tile_design
#' @importFrom crayon %+% green red bold cyan
#' @importFrom JuliaConnectoR juliaEval juliaImport juliaSetupOk juliaFun

tile_design <- function(r_mov,
                        r_source = NULL,
                        r_target = NULL,
                        max_d,
                        theta = 0.1,
                        jl_home) {
  Sys.setenv(JULIA_BINDIR = jl_home)
  if(!juliaSetupOk())
    stop("Check that the path to the Julia binary directory is correct")

  if (!is.null(r_source) && is.null(r_target)) {
    r_target <- r_source
    message("r_target not provided. Using r_source as r_target.")
  } else if (is.null(r_source) && !is.null(r_target)) {
    r_source <- r_target
    message("r_source not provided. Using r_target as r_source.")
  }

  threshold <- 0.025
  r_res <- res(r_mov)
  dim <- ceiling((4 * max_d) / r_res)
  r_crs <- crs(r_mov)
  cntr_cell <- floor(median(1:dim[1]))
  cntr_coord <- c(cntr_cell * r_res[1],
                  cntr_cell * r_res[1])

  ## Max Values
  mx_p <- global(r_mov, 'max', na.rm = T)
  mx_src <- global(r_source, 'max', na.rm = T)
  mx_target <- global(r_target, 'max', na.rm = T)

  ## Create rast
  rmat <- matrix(NaN, dim[1], dim[1]) #mx_p
  diag(rmat) <- mx_p
  src <- target <- mov <- rast(nrows = dim[1], ncols = dim[1],
                               vals = unlist(rmat),
                               resolution = r_res, crs = r_crs,
                               xmin = 0, xmax = r_res[1] * dim[1],
                               ymin = 0, ymax = r_res[1] * dim[1])
  target[target > 0] <- mx_target[[1]]
  src[src > 0] <- mx_src[[1]]
  e_dist <- as.matrix(dist(crds(target)))[,1]

  max_cell <- which.min(abs(e_dist - max_d))

  ## Julia
  juliaEval("Base.redirect_stdout(devnull); Base.redirect_stderr(devnull)")
  suppressMessages({ConScapeR_setup(jl_home)})

  # Create ConScape Grid
  g <- Grid(affinities = mov,
            sources = src,
            targets = target,
            costs = "x -> -log(x)")

  # Create a ConScape GridRSP by providing the randomness parameter theta
  h <- GridRSP(g, theta = theta)

  dists <- expected_cost(h)

  exp_d <- optimize(exp_opt,
                    dists = dists,
                    cell = max_cell,
                    threshold = threshold,
                    interval = c(1,2500),
                    tol = 0.001,
                    maximum = T)

  ## Tile size
  prox <- exp(-dists[,1] * (1/exp_d$maximum))
  tile_d <- ceiling((e_dist[which.min(abs(prox - 0.01))] * 1.65) / r_res[1])[[1]] * r_res[1]
  tile_max <- ceiling((e_dist[which.min(abs(prox - 0.001))] * 1.7) / r_res[1])[[1]] * r_res[1]
  tile_trim <- ceiling((tile_max - tile_d) / r_res[1])[[1]] * r_res[1]

  out <- list(exp_d = as.numeric(round(exp_d$maximum,1)),
              tile_d = as.numeric(tile_max),
              tile_trim = as.numeric(tile_trim),
              theta = as.numeric(theta))

  cat(green("\n" %+% cyan("     *** Tile Design Parameters ***") %+%"\n",
            "Given a maximum dispersal of " %+% red$bold(max_d) %+%" meters,\n",
            "The `exp_d` parameter should be set to: " %+% red$bold(out$exp_d) %+% "\n\n",
            "`tile_d` should be at least: " %+% red$bold(out$tile_d) %+%",\n\n",
            "`tile_trim` should be at least: " %+% red$bold(out$tile_trim) %+%",\n\n"))
  class(out) <- 'ConScapeRtools_design'
  juliaEval("Base.redirect_stdout(Base.stdout); Base.redirect_stderr(Base.stderr)")
  return(out)
}

# Optimization function ------------------------------------------------

exp_opt <- function(x, dists, cell, threshold){
  prox <- exp(-dists[,1] * (1/x))
  y <- threshold - prox[cell]
  if(y > 0){
    y <- -9999
  }
  return(y)
}
