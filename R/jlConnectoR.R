## Functions from `ConScapeR` package
#' Launch the Julia installation from R and load the `ConScape` library in Julia
#'
#' This function starts a Julia session from R and imports the `ConScape` library
#' to be used from R. It assumes that Julia is already installed in the system
#' and the path to its executables is given as the `julia_path` argument.
#' The first time the function is ran, it is best to set `install_libraries = TRUE`
#' to install the `ConScape` library in Julia.
#'
#' @param julia_path `[character]` \cr The directory for the Julia bin, e.g.
#' "C:/Programs/Julia-1.9.3/bin".
#' @param install_libraries `[logical]` \cr If `FALSE` (default), Julia will be
#' launched and the required libraries will be loaded without installing them.
#' If `TRUE`, the libraries will be (re-)installed in Julia.
#'
#' @return NULL.
#' @keywords internal
#'
ConScapeR_setup <- function(julia_path, install_libraries = FALSE) {

  Sys.setenv(JULIA_BINDIR = julia_path)

  if (install_libraries){
    Pkg <- JuliaConnectoR::juliaImport("Pkg")
    JuliaConnectoR::juliaEval("Pkg.add(\"ConScape\")")
    JuliaConnectoR::juliaEval("Pkg.add(\"SparseArrays\")")
  }
  SA <- JuliaConnectoR::juliaImport("SparseArrays")
  CS <- JuliaConnectoR::juliaImport("ConScape")
}

#' Wrapper for the Grid function of `ConScape`
#'
#' Creates a Grid from affinities, sources, targets and costs.
#' Affinities, sources and targets can be provided either as
#' `SpatRaster` from the `terra` library or as `matrix`.
#' The costs can be provided as `SpatRaster` or `matrix`,
#' but also as a string describing a transformation from the affinity matrix
#' (e.g. `"x -> -log(x)"`).
#'
#' @param affinities `[SpatRaster, matrix]` \cr The affinities, represent the likelihood
#' of movement. They can be provided as a `SpatRaster` representing the affinity,
#' permeability, or resistance of each pixel to movement (in which case it is
#' internally transformed into a matrix), or as a matrix with affinities between
#' pairs of sources-targets directly.
#' @param sources `[SpatRaster, matrix]` \cr Map or matrix presenting the habitat
#' suitability/quality of source cells.
#' @param targets `[SpatRaster, matrix]` \cr Map or matrix presenting the habitat
#' suitability/quality of target cells.
#' @param costs `[SpatRaster, matrix, character]` \cr Map or matrix presenting the
#' cost of movement through each cell. Alternatively, a string describing a transformation
#' from the affinity matrix (e.g. `"x -> -log(x)"`)
#'
#' @return A `ConScape.Grid` object within Julia.
#' @keywords internal
#'
Grid <- function(affinities, sources, targets, costs) {
  # affinities
  if (class(affinities)[1] == "SpatRaster") {
    affinities <- terra::as.matrix(affinities, wide = T)
  }
  # sources
  if (class(sources)[1] == "SpatRaster") {
    sources <- terra::as.matrix(sources, wide = T)
  }
  # targets
  if (class(targets)[1] == "SpatRaster") {
    targets <- terra::as.matrix(targets, wide = T)
  }
  # costs
  if (class(costs)[1] == "SpatRaster") {
    costs <- terra::as.matrix(costs, wide = T)
  }

  # check NaNs
  if (class(costs)[1] == "matrix") {
    nans <- is.nan(sources) | is.nan(targets) | is.nan(affinities) | is.nan(costs)
    costs[nans] <- NaN
  } else {
    nans <- is.nan(sources) | is.nan(targets) | is.nan(affinities)
  }
  sources[nans] <- NaN
  targets[nans] <- NaN
  affinities[nans] <- NaN

  # create Grid
  if (class(costs) == "character"){
    g <- JuliaConnectoR::juliaLet(
      paste0("ConScape.Grid(size(affinities)...,
              affinities=ConScape.graph_matrix_from_raster(affinities),
              source_qualities=sources,
              target_qualities=SparseArrays.sparse(targets),
              costs=ConScape.mapnz(", costs, ", ConScape.graph_matrix_from_raster(affinities)))"),
      affinities=affinities, sources=sources, targets=targets)
  } else {
    g <- JuliaConnectoR::juliaLet("ConScape.Grid(size(affinities)...,
                  affinities=ConScape.graph_matrix_from_raster(affinities),
                  source_qualities=sources,
                  target_qualities=SparseArrays.sparse(targets),
                  costs=ConScape.graph_matrix_from_raster(costs))",
                                  affinities=affinities, sources=sources, targets=targets, costs=costs)
  }
  return(g)
}


#' Convert a matrix to raster
#'
#' @param mat `[matrix]` \cr Matrix to be converted.
#' @param rast `[SpatRaster]` \cr Template raster, usually one of those used in the
#' [ConScapeR::Grid()] function.
#'
#' @return `[SpatRaster]`
#' @keywords internal
#'
mat2rast <- function(mat, rast){
  rast2 <- terra::rast(mat, extent = terra::ext(rast), crs = terra::crs(rast))
  return(rast2)
}

#' Convert a node vector to matrix
#'
#' @param vec `[matrix]` \cr Matrix to be converted
#' @param g `[Grid]` \cr Template raster, usually one of those used in the
#' [ConScapeR::Grid()] function.
#'
#' @return `[matrix]`
#' @keywords internal
#'
vec2mat <- function(vec, g){
  mat <- JuliaConnectoR::juliaLet("ConScape._vec_to_grid(g, vec)", g=g, vec=vec)
  return(mat)
}

#' Wrapper for the GridRSP function of `ConScape`
#'
#' Creates a GridRSP from a ConScape.Grid and the randomness parameter. Theta closer to zero gives a more random movement,
#' as theta increases the RSP distribution will converge to the optimal or least-cost paths.
#' Note that the computation becomes numerically unstable as theta becomes really small or large.
#' See Van Moorter et al. (2023, Methods in Ecology and Evolution) for more details.
#'
#' @param g `[Grid]` \cr The output from the [ConScapeR::Grid()] function.
#' @param theta `[numeric]` \cr The randomness parameter theta. Lower `theta` values
#' represent a more random walk, and higher values a more least-cost path behavior.
#'
#' @return A `ConScape.GridRSP` object within Julia.
#' @keywords internal
#'
GridRSP <- function(g, theta) {
  h <- JuliaConnectoR::juliaLet("ConScape.GridRSP(g, θ=theta)", g=g, theta=theta)
  return(h)
}

#' Wrapper for the `betweenness_qweighted` function of `ConScape`
#'
#' Computes the quality-weighted betweenness for a ConScape.GridRSP object.
#' See Van Moorter et al. (2023, Methods in Ecology and Evolution) for more details.
#'
#' @param h `[GridRSP]` \cr The output from the [ConScapeR::GridRSP()] function.
#'
#' @return `matrix`
#' @keywords internal
#'
betweenness_qweighted <- function(h) {
  betw <- JuliaConnectoR::juliaLet("ConScape.betweenness_qweighted(h)", h=h)
  return(betw)
}

#' Wrapper for the `betweenness_kweighted` function of `ConScape`
#'
#' Computes the quality- and proximity-weighted betweenness for a ConScape.GridRSP object.
#' See Van Moorter et al. (2023, Methods in Ecology and Evolution) for more details.
#'
#' @param h `[GridRSP]` \cr The output from the [ConScapeR::GridRSP()] function.
#' @param alpha `[numeric]` \cr The distance scaling for the exponential distance
#' to proximity transformation.
#'
#' @return `matrix`
#' @keywords internal
#'
betweenness_kweighted <- function(h, alpha) {
  betw <- JuliaConnectoR::juliaLet("ConScape.betweenness_kweighted(h, distance_transformation=x -> exp(-x*alpha))", h=h, alpha=alpha)
  return(betw)
}

#' Wrapper for the `connected_habitat` function of `ConScape`
#'
#' Computes the habitat functionality for a ConScape.GridRSP object.
#' See Van Moorter et al. (2023, Methods in Ecology and Evolution) for more details.
#'
#' @param h `[GridRSP]` \cr The output from the [ConScapeR::GridRSP()] function.
#' @param alpha `[numeric]` \cr The distance scaling for the exponential distance to
#' proximity transformation.
#'
#' @return A `matrix` with loca (pixel) measures of habitat functionality.
#' @keywords internal
#'
connected_habitat <- function(h, alpha) {
  func <- JuliaConnectoR::juliaLet("ConScape.connected_habitat(h, distance_transformation=x -> exp(-x*alpha))", h=h, alpha=alpha)
  return(func)
}

#' Wrapper for the `expected_cost` function of `ConScape`
#'
#' Computes the RSP expected cost between all sources and non-zero targets.
#' See Van Moorter et al. (2023, Methods in Ecology and Evolution) for more details.
#'
#' @param h `[GridRSP]` \cr The output from the [ConScapeR::GridRSP()] function.
#'
#' @keywords internal
#'
expected_cost <- function(h) {
  dists = JuliaConnectoR::juliaLet("ConScape.expected_cost(h)", h=h)
  return(dists)
}




