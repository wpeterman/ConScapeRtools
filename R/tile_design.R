#' Design tiling and decay parameters for ConScape
#'
#' @description
#' Uses the resolution and values of a movement–probability raster to
#' derive a suitable exponential decay parameter (`exp_d`) and suggested
#' tile width (`tile_d`) and overlap (`tile_trim`) for ConScape analyses,
#' given a presumed maximum dispersal distance in ideal habitat.
#' Internally, this builds a small synthetic landscape, runs ConScape's
#' randomized shortest–path model, and calibrates the decay so that
#' proximity at `max_d` matches a user–defined threshold.
#'
#' @param r_mov `SpatRaster` representing movement probabilities. Values
#'   should be in a projected coordinate reference system (map units in
#'   meters) and represent per–cell movement affinity or probability.
#' @param r_source Optional `SpatRaster` representing source quality
#'   (default `NULL`). If only `r_target` is provided, it is used for both
#'   source and target.
#' @param r_target Optional `SpatRaster` representing target quality
#'   (default `NULL`). If only `r_source` is provided, it is used for both
#'   source and target.
#' @param max_d Presumed maximum movement distance through ideal habitat,
#'   in the same map units as the raster CRS (typically meters). Used to
#'   calibrate the exponential distance–decay parameter.
#' @param theta Parameter controlling path randomness in the randomized
#'   shortest–path model (`GridRSP`). As `theta` approaches 0, movement
#'   becomes increasingly random; large values approach deterministic
#'   least–cost paths. Default is `0.1`.
#' @param jl_home Path to the `bin` directory where the Julia executable
#'   resides. This is passed to `JuliaConnectoR` and used to initialize
#'   the ConScape / ConScapeR environment.
#'
#' @details
#' The input `r_mov` must be in a projected CRS so that cell resolution and
#' `max_d` are expressed in compatible distance units. Only the cell
#' resolution and maximum values of `r_mov`, `r_source`, and `r_target`
#' are used: the function constructs a synthetic square landscape whose
#' size is chosen to span roughly four times `max_d`, sets the diagonal
#' of the movement grid to the maximum movement probability, and assigns
#' constant maximum source and target quality.
#'
#' A ConScape `Grid` and `GridRSP` object are then created, and expected
#' costs are computed from the center cell to all others. The decay
#' parameter `exp_d` is chosen by 1D optimization so that the exponential
#' proximity at the cell closest to `max_d` drops to a specified
#' threshold. This `exp_d` is then used to derive:
#'
#' * `tile_d` – a suggested minimum interior tile width for ConScape runs,
#' * `tile_trim` – a suggested minimum overlap between tiles, chosen so
#'   that contributions beyond the trimmed edge are negligible under the
#'   calibrated decay.
#'
#' These outputs are intended to be passed to [conscape_prep()] (for
#' `tile_d` and `tile_trim`) and to [run_conscape()] (for `exp_d` and
#' `theta`).
#'
#' @return
#' A named list of class `"ConScapeRtools_design"` with elements:
#'
#' * `exp_d` – calibrated exponential decay parameter for distance
#'   transformation in ConScape.
#' * `tile_d` – suggested minimum tile width (map units).
#' * `tile_trim` – suggested minimum tile overlap / trim width (map units).
#' * `theta` – the `theta` value used when calibrating `exp_d`.
#'
#' The function also prints a short, colorized summary of these design
#' parameters to the console.
#'
#' @export
#' @examples
#' \dontrun{
#' library(ConScapeRtools)
#'
#' ## Import data
#' s <- system.file("data/suitability.asc", package = "ConScapeRtools")
#' source <- terra::rast(s)
#'
#' a <- system.file("data/affinity.asc", package = "ConScapeRtools")
#' resist <- terra::rast(a)
#'
#' jl_home <- "/path/to/julia/bin"
#'
#' td <- tile_design(r_mov = resist,
#'                   r_target = source,
#'                   max_d = 7000,
#'                   theta = 0.1,
#'                   jl_home = jl_home)
#' }
#' @seealso [conscape_prep()], [tile_rast()], [run_conscape()]
#' @author Bill Peterman
#' @details
#' `r_mov` must be in a projected coordinate reference system.
#' The function relies on helper routines adapted from the `ConScapeR`
#' package (https://github.com/ConScape/ConScapeR) and uses Julia via
#' `JuliaConnectoR`.
#'
#' @rdname tile_design
#' @importFrom crayon %+% green red bold cyan
#' @importFrom JuliaConnectoR juliaEval juliaImport juliaSetupOk juliaFun stopJulia
#' @importFrom stats dist median optimize


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
  invisible({ConScapeR_setup(jl_home)})

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
  # if(method == "empirical"){
  trim_d  <- floor(exp_distance(p = 0.99, lambda = 1/exp_d$maximum)) * r_res[1]
  tile_d  <- floor(exp_distance(p = 0.9901, lambda = 1/exp_d$maximum)) * r_res[1]
  tile_max <- ceiling(exp_distance(p = 0.999, lambda = 1/exp_d$maximum)) * (r_res[1])
  # } else {
  #   prox <- exp(-dists[,1] * (1/exp_d$maximum))
  #   tile_d <- ceiling((e_dist[which.min(abs(prox - 0.01))] * 1.65) / r_res[1])[[1]] * r_res[1]
  #   tile_max <- ceiling((e_dist[which.min(abs(prox - 0.001))] * 1.7) / r_res[1])[[1]] * r_res[1]
  # }
  tile_trim <- ceiling((tile_max - trim_d) / r_res[1])[[1]] * r_res[1]

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
  invisible(juliaEval('1+1'))
  stopJulia()
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

#' Calculate the distance for a given cumulative density in a negative exponential distribution
#' @noRd

exp_distance <- function(p, lambda) {
  if (p <= 0 || p >= 1) {
    stop("Cumulative density `p` must be between 0 and 1 (exclusive).")
  }
  if (lambda <= 0) {
    stop("Rate parameter `lambda` must be positive.")
  }
  x <- -log(1 - p) / lambda
  return(x)
}

