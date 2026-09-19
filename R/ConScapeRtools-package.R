#' ConScapeRtools: Tiled ConScape Workflows from R
#'
#' @description
#' ConScapeRtools prepares raster inputs, runs randomized shortest-path
#' connectivity analyses through Julia and ConScape.jl, and combines tiled
#' outputs into continuous raster surfaces. The stable workflow separates the
#' disposable tile directory created by [conscape_prep()] from the run output
#' directory used by [run_conscape()].
#'
#' Start with `vignette("ConScapeRtools_Guide", package = "ConScapeRtools")`
#' for the complete workflow. Use
#' `vignette("Tiled_ConScape_Validation_and_Performance", package =
#' "ConScapeRtools")` to choose windows and evaluate buffer convergence.
#'
#' @section Main workflow:
#' * [tile_design()] calibrates tile dimensions and distance decay.
#' * [conscape_prep()] writes matched target, source, and affinity tiles.
#' * [run_conscape()] runs ConScape and returns mosaicked results plus output
#'   directories and diagnostics.
#' * [mosaic_conscape()] reconstructs a raster manually from saved tile output.
#' * [conscape_efficiency_assessment()] compares candidate tiling designs before
#'   running Julia.
#'
#' @section External software:
#' Julia and ConScape.jl are external requirements. Normal analysis functions
#' never install or update Julia packages. Use `install_libraries = TRUE` or an
#' explicit setup helper only when you intend to modify a Julia environment.
#'
#' @name ConScapeRtools-package
#' @aliases ConScapeRtools
#' @keywords internal
"_PACKAGE"
