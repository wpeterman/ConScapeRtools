#' Run ConScape
#'
#' @description
#' Run ConScape on tiled landscapes produced by [conscape_prep()] or on
#' single, untiled rasters. Optionally uses parallel execution (in R or
#' Julia) and can mosaic tile-level outputs back into continuous rasters.
#'
#' @param conscape_prep Optional object of class `"ConScapeRtools_prep"`
#'   created by [conscape_prep()]. When supplied, this provides the tile
#'   layout, directories for source/target/movement tiles, the mask, and
#'   `landmark` / `tile_trim` values. If `NULL` (default), the user must
#'   supply `target_qualities`, `source_qualities`, and `affinities` directly.
#' @param out_dir Directory where ConScape outputs will be written. Metric
#'   subdirectories (e.g., `"btwn"`, `"fcon"`) are created inside `out_dir`.
#' @param target_qualities Either (i) a path to a directory containing
#'   target-quality tiles (`*.asc`) or (ii) a single `SpatRaster` used as
#'   ConScape's `target_qualities` layer, cells that can *receive*
#'   connectivity. The legacy name `hab_target` is accepted as a
#'   backwards-compatible alias.
#' @param source_qualities Either (i) a path to a directory containing
#'   source-quality tiles (`*.asc`) or (ii) a single `SpatRaster` used as
#'   ConScape's `source_qualities` layer, cells that *generate* connectivity.
#'   The legacy name `hab_src` is accepted as a backwards-compatible alias.
#' @param affinities Either (i) a path to a directory containing
#'   affinity/permeability tiles (`*.asc`) or (ii) a single `SpatRaster` used
#'   to build ConScape's `affinities` adjacency matrix via
#'   `graph_matrix_from_raster()`. Values should be higher where movement is
#'   easier (i.e., permeability, not resistance). The legacy name `mov_prob`
#'   is accepted as a backwards-compatible alias.
#' @param clear_dir Logical. If `TRUE` (default), any existing contents
#'   of `out_dir` are removed before writing new results. If `FALSE` and
#'   `out_dir` is not empty, the function stops with an error.
#' @param landmark Integer coarse-graining window size passed to ConScape's
#'   `coarse_graining()` (default `10L`). Larger values aggregate target
#'   qualities to fewer cells, reducing computation. When `conscape_prep` is
#'   provided, `landmark` is taken from that object. For sensitivity runs with
#'   `conscape_prep`, create the prep object with `landmark = 1L` so that target
#'   qualities are not altered by coarse graining.
#' @param theta Parameter controlling the amount of randomness in paths.
#'   As `theta` approaches 0, movement is increasingly random; as `theta`
#'   increases, paths approach deterministic least-cost paths. Default is `0.01`.
#'   Obtain a landscape-specific value from [tile_design()].
#' @param distance_scale Numerator of the exponential distance-to-proximity
#'   transformation used by ConScape (i.e., `exp(-dist / distance_scale)`).
#'   Default is `150`. Use the value returned by [tile_design()] for a
#'   landscape-calibrated decay. The legacy name `exp_d` is accepted as a
#'   backwards-compatible alias.
#' @param NA_val Backwards-compatible argument retained for older workflows.
#'   The current Julia runner treats `NA` affinities as absent graph cells and
#'   `NA` source/target qualities as zero quality.
#' @param mosaic Logical (default `TRUE`). When running on tiles
#'   (via `conscape_prep` or directory inputs), tile-level outputs are
#'   mosaicked into continuous rasters using [mosaic_conscape()]. When
#'   running on a single raster (`SpatRaster` inputs), `mosaic` is ignored.
#' @param jl_home Path to the `bin` directory where the Julia executable
#'   is installed. Used by `JuliaConnectoR` to initialize ConScape.
#' @param parallel Logical. If `FALSE` (default), tiles are processed
#'   serially. If `TRUE`, tiles are processed in parallel (either via R
#'   or Julia, depending on `parallel_R` and `distributed`).
#' @param workers Integer number of parallel workers to use when
#'   `parallel = TRUE`. Default is `1`.
#' @param blas_threads Integer number of BLAS threads to allow inside each
#'   Julia tile worker. The default, `1`, avoids oversubscribing CPU cores when
#'   ConScape is already parallelized across tiles.
#' @param distributed Logical. When `parallel = TRUE` and
#'   `parallel_R = FALSE`, controls whether Julia uses distributed
#'   workers (`TRUE`) or multithreading (`FALSE`, default).
#' @param progress Logical. If `TRUE` (default), progress information is
#'   printed for serial and R-parallel runs. Progress reporting from Julia
#'   parallel runs may be less detailed.
#' @param parallel_R Logical. If `TRUE`, tiles are processed in parallel
#'   using R (`future`/`future.apply`). If `FALSE` (default), parallelism
#'   is handled entirely within Julia.
#' @param backend Execution backend. `"classic"` uses the stable bundled Julia
#'   wrappers. `"conscape_dev"` is experimental and requires a ConScape
#'   development installation exposing `Problem`, `WindowedProblem`, and
#'   `solve`; `dev_mode = "batch"` additionally requires `BatchProblem`.
#' @param tile_trim Width of the overlapping border (in map units) to
#'   trim from each tile when mosaicking. When `conscape_prep` is supplied,
#'   this value is taken from that object and passed automatically to
#'   [mosaic_conscape()]. When running on a single `SpatRaster`, `tile_trim`
#'   controls how far the rasters are buffered on all sides before ConScape
#'   runs; the output is then cropped to the original extent so that tiled
#'   and untiled results are comparable.
#' @param cost_function Cost transformation applied to `affinities` when
#'   building the ConScape graph (`ConScape.Grid(costs = ...)`). Named
#'   options are `"minuslog"` (default), `"inverse"`, `"odds_against"`,
#'   `"odds_for"`, and `"expminus"`. A Julia lambda string such as
#'   `"x -> -log(x)"` may also be supplied, in which case costs are
#'   constructed with `ConScape.mapnz()`.
#' @param connectivity_function ConScape distance or proximity function
#'   used to build the landscape matrix for metrics. Supported values are
#'   `"expected_cost"` (default), `"least_cost_distance"`,
#'   `"free_energy_distance"`, `"survival_probability"`, and
#'   `"power_mean_proximity"`.
#' @param metrics Character vector of raster-valued ConScape metrics to
#'   compute. Supported values are `"connected_habitat"`,
#'   `"betweenness_kweighted"`, `"betweenness_qweighted"`, and
#'   `"criticality"`. The short aliases `"fcon"` and `"btwn"` are also
#'   accepted. Defaults to `c("betweenness_kweighted", "connected_habitat")`.
#' @param sensitivity Optional sensitivity specification created by
#'   [conscape_sensitivity()]. Pass `TRUE` for the default quality and linked
#'   affinity-cost elasticities; pass a character vector to specify `wrt`
#'   values directly. See [conscape_sensitivity()] for full options.
#' @param centersize Center window size in cells for the experimental
#'   `"conscape_dev"` backend.
#' @param buffer Buffer width in cells for the experimental `"conscape_dev"`
#'   backend.
#' @param window_shape Window shape for the experimental `"conscape_dev"`
#'   backend. The current ConScape dev API supports `"square"` windows.
#' @param dev_mode Experimental ConScape dev execution mode. `"windowed"` runs
#'   `WindowedProblem` directly and returns a mosaicked raster stack. `"batch"`
#'   wraps the problem in `BatchProblem`, writes intermediate batch rasters, and
#'   mosaics them with ConScape after all batches complete.
#' @param batch_grain Optional target thinning value passed to the experimental
#'   ConScape `BatchProblem`. Leave `NULL` to preserve target qualities.
#' @param batch_ext Raster extension for experimental `BatchProblem`
#'   intermediates. Defaults to `".tif"`.
#' @param dev_project Optional Julia project directory containing the ConScape
#'   development backend dependencies. When supplied with
#'   `backend = "conscape_dev"`, the project is activated via `JULIA_PROJECT`
#'   before Julia starts.
#' @param install_dev_conscape Logical. If `TRUE` and
#'   `backend = "conscape_dev"`, create or update `dev_project` by installing
#'   ConScape from `dev_conscape_url` at `dev_conscape_rev` before running.
#'   This is opt-in because it downloads a development Git dependency.
#' @param dev_conscape_rev Git branch, tag, or commit used when
#'   `install_dev_conscape = TRUE`. Defaults to `"alg_efficiency"`.
#' @param dev_conscape_url Git URL used when `install_dev_conscape = TRUE`.
#' @param mosaic_chunk_size Maximum number of output tiles to hold in memory
#'   during each mosaic batch. Smaller values reduce peak memory at the cost of
#'   more temporary files.
#' @param stop_julia Logical. If `TRUE` (default), close the Julia session when
#'   `run_conscape()` exits. Set to `FALSE` to keep Julia open for faster
#'   repeated ConScape calls with the same Julia path and process settings. Use
#'   [conscape_julia_stop()] when finished.
#' @param hab_target Deprecated. Use `target_qualities` instead.
#' @param hab_src Deprecated. Use `source_qualities` instead.
#' @param mov_prob Deprecated. Use `affinities` instead.
#' @param exp_d Deprecated. Use `distance_scale` instead.
#'
#' @details
#' In typical workflows, rasters are prepared with [conscape_prep()] and
#' the resulting object is supplied via `conscape_prep`. This ensures that
#' `target_qualities`, `source_qualities`, and `affinities` tiles share the
#' same tiling scheme, that coarse-graining landmarks are aligned across tiles,
#' and that a habitat mask and `tile_trim` value are available for mosaicking.
#'
#' When `conscape_prep` is `NULL`, the function can operate on:
#' * Directories of `*.asc` tiles passed to `target_qualities`,
#'   `source_qualities`, and `affinities`, or
#' * Three `SpatRaster` objects of identical extent and resolution passed to
#'   `target_qualities`, `source_qualities`, and `affinities`, in which case
#'   ConScape is run on the full (untiled) landscape.
#'
#' When `parallel = TRUE`, parallelization is handled either in R (with
#' `parallel_R = TRUE`, using `future`/`future.apply`) or in Julia (with
#' `parallel_R = FALSE`). For Julia-level parallelism, the `distributed`
#' argument selects between multithreading (`FALSE`, default) and distributed
#' workers (`TRUE`).
#'
#' The default run computes `connected_habitat` and `betweenness_kweighted`,
#' corresponding to the historical `fcon` and `btwn` output layers. Additional
#' metrics are requested with `metrics`. All metrics are passed to the bundled
#' Julia runner, which builds a single `GridRSP` object and then evaluates
#' each requested surface on it.
#'
#' Sensitivity outputs are requested by passing a [conscape_sensitivity()]
#' specification to `sensitivity`. These surfaces are computed after the
#' `GridRSP` object is built and written to directories named for the output,
#' for example `"elasticity_quality"` or `"elasticity_affinity_cost_linked"`.
#' Analytical sensitivity requires a ConScape installation that provides
#' `ConScape.sensitivity`; when that function is unavailable the Julia runner
#' stops with an explicit message.
#'
#' Tiled sensitivity uses the same tiling machinery as the default metric
#' outputs but with two important constraints:
#'
#' 1. The [conscape_prep()] object **must** be built with
#'    `target_mode = "full"`. The corrected `target_mode = "center"`
#'    (the new default for fcon / btwn) sets per-tile target qualities
#'    to a strict subset of source qualities, which violates the
#'    `target_equal_source = TRUE` precondition of ConScape's analytical
#'    sensitivity inside every tile. `run_conscape()` refuses
#'    `conscape_prep$target_mode == "center"` with sensitivity requested
#'    and tells you to recreate the prep with `target_mode = "full"`.
#' 2. Build the prep with `landmark = 1L` so that ConScape's
#'    `coarse_graining()` leaves target qualities equal to source
#'    qualities.
#'
#' Even with both constraints satisfied, tiled sensitivity is a
#' tile-local approximation to landscape-level sensitivity: each tile
#' computes the derivative of its own per-tile landscape summary, not
#' the global one, and the mosaic step averages overlapping per-tile
#' estimates. The approximation converges to the untiled sensitivity
#' as `tile_trim` (buffer) grows -- see
#' `tests/testthat/test-sensitivity-convergence.R`. For strict
#' correctness, run sensitivity untiled when the landscape fits in a
#' single graph.
#'
#' `run_conscape()` mosaics sensitivity surfaces with `method = "mosaic"`
#' (mean) regardless of `target_mode`, because sensitivity outputs are
#' tile-local landscape-summary derivatives rather than pairwise sums.
#' The mosaic dispatch is per output layer, so fcon / btwn can use
#' sum mosaic while sensitivity uses mean mosaic in the same run when
#' `target_mode = "full"` is in effect.
#'
#' @return
#' When tiles are used and `mosaic = TRUE`, returns an object of class
#' `"ConScapeResults"`, a named list containing:
#'
#' * `btwn` – `SpatRaster` of k-weighted betweenness (`betweenness_kweighted`),
#'   when requested.
#' * `fcon` – `SpatRaster` of connected habitat (`connected_habitat`), when
#'   requested.
#' * Additional metric or sensitivity `SpatRaster` layers named by their output
#'   identifier (e.g., `"btwn_qweighted"`, `"elasticity_quality"`).
#' * `outdirs` – named list mapping each output layer to its tile-level
#'   directory.
#' * `outdir_btwn` and `outdir_fcon` – backwards-compatible named paths for
#'   the default betweenness and connected-habitat tile directories.
#'
#' When run on a single untiled `SpatRaster`, returns a multi-layer
#' `SpatRaster` with one layer per requested metric and/or sensitivity output,
#' named by their output identifiers (e.g., `"btwn"`, `"fcon"`,
#' `"elasticity_quality"`).
#'
#' @seealso [conscape_prep()], [tile_design()], [mosaic_conscape()],
#'   [conscape_sensitivity()], [conscape_api_slots()]
#' @export
#' @examples
#' \dontrun{
#' library(ConScapeRtools)
#'
#' ## Import example data
#' s <- system.file("extdata", "suitability.asc", package = "ConScapeRtools")
#' habitat <- terra::rast(s)
#'
#' a <- system.file("extdata", "affinity.asc", package = "ConScapeRtools")
#' affinity <- terra::rast(a)
#'
#' jl_home <- "/path/to/julia/bin" # Update to your Julia installation
#'
#' ## Calibrate decay and tiling parameters
#' td <- tile_design(r_mov    = affinity,
#'                   r_target = habitat,
#'                   max_d    = 8500,
#'                   theta    = 0.15,
#'                   jl_home  = jl_home,
#'                   landmark = 5L)
#'
#' ## Prepare tiled rasters
#' prep <- conscape_prep(tile_d         = td$tile_d,
#'                       tile_trim      = td$tile_trim,
#'                       r_target       = habitat,
#'                       r_mov          = affinity,
#'                       r_src          = habitat,
#'                       clear_dir      = TRUE,
#'                       landmark       = td$landmark)
#'
#' ## Serial run (tiled)
#' cs_serial <- run_conscape(out_dir        = file.path(prep$asc_dir, "results"),
#'                           conscape_prep  = prep,
#'                           theta          = td$theta,
#'                           distance_scale = td$distance_scale,
#'                           jl_home        = jl_home)
#'
#' ## R-parallel run
#' cs_r <- run_conscape(out_dir        = file.path(prep$asc_dir, "results"),
#'                      conscape_prep  = prep,
#'                      theta          = td$theta,
#'                      distance_scale = td$distance_scale,
#'                      jl_home        = jl_home,
#'                      parallel       = TRUE,
#'                      workers        = 4,
#'                      parallel_R     = TRUE)
#'
#' ## Julia threaded parallel
#' cs_thread <- run_conscape(out_dir        = file.path(prep$asc_dir, "results"),
#'                           conscape_prep  = prep,
#'                           theta          = td$theta,
#'                           distance_scale = td$distance_scale,
#'                           jl_home        = jl_home,
#'                           parallel       = TRUE,
#'                           workers        = 4)
#'
#' ## Julia distributed parallel
#' cs_dist <- run_conscape(out_dir        = file.path(prep$asc_dir, "results"),
#'                         conscape_prep  = prep,
#'                         theta          = td$theta,
#'                         distance_scale = td$distance_scale,
#'                         jl_home        = jl_home,
#'                         parallel       = TRUE,
#'                         workers        = 4,
#'                         distributed    = TRUE)
#' plot(cs_dist)
#'
#' ## Untiled run, only suitable for small to moderate rasters
#' cs_single <- run_conscape(out_dir          = file.path(prep$asc_dir, "results"),
#'                           target_qualities = habitat,
#'                           source_qualities = habitat,
#'                           affinities       = affinity,
#'                           theta            = td$theta,
#'                           distance_scale   = td$distance_scale,
#'                           landmark         = 5L,
#'                           jl_home          = jl_home)
#'
#' plot(cs_serial$fcon - cs_single$fcon, main = "Tiled minus untiled")
#' layerCor(c(cs_single$fcon, cs_serial$fcon), fun = "cor")
#' layerCor(c(cs_single$btwn, cs_serial$btwn), fun = "cor")
#'
#' ## Experimental ConScape dev backend using package data.
#' ## This creates or reuses a dedicated Julia project with ConScape's
#' ## development backend API.
#' habitat_demo <- terra::crop(
#'   habitat,
#'   terra::ext(
#'     terra::xmin(habitat),
#'     terra::xmin(habitat) + 40 * terra::res(habitat)[1],
#'     terra::ymax(habitat) - 40 * terra::res(habitat)[2],
#'     terra::ymax(habitat)
#'   ),
#'   snap = "near"
#' )
#' affinity_demo <- terra::crop(affinity, terra::ext(habitat_demo), snap = "near")
#'
#' dev_project <- conscape_dev_backend_setup(
#'   jl_home = jl_home,
#'   rev = "alg_efficiency"
#' )
#'
#' cs_dev <- run_conscape(
#'   target_qualities = habitat_demo,
#'   source_qualities = habitat_demo,
#'   affinities = affinity_demo,
#'   out_dir = file.path(tempdir(), "conscape_dev_windowed"),
#'   jl_home = jl_home,
#'   backend = "conscape_dev",
#'   dev_mode = "windowed",
#'   centersize = 10,
#'   buffer = 10,
#'   dev_project = dev_project,
#'   landmark = 1L,
#'   theta = 0.15,
#'   distance_scale = 150,
#'   metrics = c("btwn", "fcon")
#' )
#'
#' cs_dev_batch <- run_conscape(
#'   target_qualities = habitat_demo,
#'   source_qualities = habitat_demo,
#'   affinities = affinity_demo,
#'   out_dir = file.path(tempdir(), "conscape_dev_batch"),
#'   jl_home = jl_home,
#'   backend = "conscape_dev",
#'   dev_mode = "batch",
#'   centersize = 10,
#'   buffer = 10,
#'   dev_project = dev_project,
#'   landmark = 1L,
#'   theta = 0.15,
#'   distance_scale = 150,
#'   metrics = c("btwn", "fcon")
#' )
#' }
#' @author Bill Peterman
#' @importFrom JuliaConnectoR juliaImport juliaEval juliaSetupOk juliaGet juliaCall stopJulia
run_conscape <- function(conscape_prep = NULL,
                         out_dir,
                         target_qualities = NULL,
                         source_qualities = NULL,
                         affinities = NULL,
                         clear_dir = TRUE,
                         landmark = 10L,
                         theta = 0.01,
                         distance_scale = 150,
                         NA_val = 1e-50,
                         mosaic = TRUE,
                         jl_home,
                         parallel = FALSE,
                         workers = 1,
                         distributed = FALSE,
                         progress = TRUE,
                         parallel_R = FALSE,
                         backend = c("classic", "conscape_dev"),
                         tile_trim = 0,
                         cost_function = "minuslog",
                         connectivity_function = "expected_cost",
                         metrics = c("betweenness_kweighted", "connected_habitat"),
                         sensitivity = NULL,
                         centersize = NULL,
                         buffer = NULL,
                         window_shape = c("square", "circle"),
                         dev_mode = c("windowed", "batch"),
                         batch_grain = NULL,
                         batch_ext = ".tif",
                         dev_project = NULL,
                         install_dev_conscape = FALSE,
                         dev_conscape_rev = "alg_efficiency",
                         dev_conscape_url = "https://github.com/ConScape/ConScape.jl",
                         blas_threads = 1L,
                         mosaic_chunk_size = 64L,
                         stop_julia = TRUE,
                         hab_target = NULL,
                         hab_src = NULL,
                         mov_prob = NULL,
                         exp_d = NULL){
  ## --- resolve deprecated aliases ---
  if (!is.null(hab_target)) {
    if (!is.null(target_qualities)) {
      stop("Use either target_qualities or hab_target, not both.", call. = FALSE)
    }
    target_qualities <- hab_target
  }
  if (!is.null(hab_src)) {
    if (!is.null(source_qualities)) {
      stop("Use either source_qualities or hab_src, not both.", call. = FALSE)
    }
    source_qualities <- hab_src
  }
  if (!is.null(mov_prob)) {
    if (!is.null(affinities)) {
      stop("Use either affinities or mov_prob, not both.", call. = FALSE)
    }
    affinities <- mov_prob
  }
  if (!is.null(exp_d)) {
    if (!missing(distance_scale)) {
      stop("Use either distance_scale or exp_d, not both.", call. = FALSE)
    }
    distance_scale <- exp_d   # legacy: override the default 150
  }
  ## Map internal names used throughout body
  hab_target <- target_qualities
  hab_src    <- source_qualities
  mov_prob   <- affinities
  exp_d      <- distance_scale
  if (!is.numeric(exp_d) || length(exp_d) != 1L || is.na(exp_d) || exp_d <= 0) {
    stop("distance_scale/exp_d must be a single positive numeric value.", call. = FALSE)
  }
  if (!is.numeric(theta) || length(theta) != 1L || is.na(theta) || theta <= 0) {
    stop("theta must be a single positive numeric value.", call. = FALSE)
  }
  if (!is.numeric(landmark) || length(landmark) != 1L || is.na(landmark) || landmark < 1) {
    stop("landmark must be a positive integer.", call. = FALSE)
  }
  if (!is.logical(stop_julia) || length(stop_julia) != 1L || is.na(stop_julia)) {
    stop("stop_julia must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(blas_threads) || length(blas_threads) != 1L ||
      is.na(blas_threads) || blas_threads < 1) {
    stop("blas_threads must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(mosaic_chunk_size) || length(mosaic_chunk_size) != 1L ||
      is.na(mosaic_chunk_size) || mosaic_chunk_size < 1) {
    stop("mosaic_chunk_size must be a positive integer.", call. = FALSE)
  }

  landmark <- as.integer(landmark)
  blas_threads <- as.integer(blas_threads)
  mosaic_chunk_size <- as.integer(mosaic_chunk_size)
  backend <- match.arg(backend)
  window_shape <- match.arg(window_shape)
  dev_mode <- match.arg(dev_mode)
  if (!is.null(batch_grain) &&
      (!is.numeric(batch_grain) || length(batch_grain) != 1L ||
       is.na(batch_grain) || batch_grain < 1)) {
    stop("batch_grain must be NULL or a positive integer.", call. = FALSE)
  }
  if (!is.character(batch_ext) || length(batch_ext) != 1L ||
      !nzchar(batch_ext)) {
    stop("batch_ext must be a non-empty character string.", call. = FALSE)
  }
  if (!is.null(dev_project) &&
      (!is.character(dev_project) || length(dev_project) != 1L ||
       is.na(dev_project) || !nzchar(dev_project))) {
    stop("dev_project must be NULL or a single non-empty character string.", call. = FALSE)
  }
  if (!is.logical(install_dev_conscape) || length(install_dev_conscape) != 1L ||
      is.na(install_dev_conscape)) {
    stop("install_dev_conscape must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.character(dev_conscape_rev) || length(dev_conscape_rev) != 1L ||
      is.na(dev_conscape_rev) || !nzchar(dev_conscape_rev)) {
    stop("dev_conscape_rev must be a single non-empty character string.", call. = FALSE)
  }
  if (!is.character(dev_conscape_url) || length(dev_conscape_url) != 1L ||
      is.na(dev_conscape_url) || !nzchar(dev_conscape_url)) {
    stop("dev_conscape_url must be a single non-empty character string.", call. = FALSE)
  }
  cost_function <- normalize_conscape_cost_function(cost_function)
  connectivity_function <- normalize_conscape_connectivity_function(connectivity_function)
  metrics <- normalize_conscape_metrics(metrics)
  sensitivity <- normalize_conscape_sensitivity(sensitivity)
  output_specs <- conscape_output_specs(metrics, sensitivity)
  sensitivity_wrt <- if (is.null(sensitivity)) "" else sensitivity$wrt
  sensitivity_method <- if (is.null(sensitivity)) "analytical" else sensitivity$method
  sensitivity_landscape_measure <- if (is.null(sensitivity)) "sum" else sensitivity$landscape_measure
  sensitivity_unitless <- if (is.null(sensitivity)) TRUE else sensitivity$unitless
  sensitivity_one_out_of <- if (is.null(sensitivity)) 1L else sensitivity$one_out_of
  sensitivity_diagvalue <- if (is.null(sensitivity)) NULL else sensitivity$diagvalue
  sensitivity_target_equal_source <- if (is.null(sensitivity)) TRUE else sensitivity$target_equal_source
  sensitivity_require_landmark_one <- if (is.null(sensitivity)) TRUE else sensitivity$require_landmark_one

  # Sensitivity precondition: target_equal_source = TRUE.
  #
  # ConScape's analytical sensitivity (and its simulation counterpart) is only
  # well-defined when each tile's per-cell target qualities equal its per-cell
  # source qualities. The corrected `target_mode = "center"` semantics
  # deliberately sets target = (center cells only) while source = (full
  # buffered window), which violates that precondition inside every tile and
  # also makes the sum-mosaic reduction inappropriate for per-cell
  # sensitivity surfaces (they are tile-local landscape-summary derivatives,
  # not pairwise contributions that partition disjointly over tiles).
  #
  # Until a center-mode-compatible sensitivity reduction is proven, refuse
  # the combination and tell the user to recreate the prep with
  # target_mode = "full". This must fire before we start Julia / clear
  # out_dir so the failure is recoverable.
  if (!is.null(sensitivity) &&
      !is.null(conscape_prep) &&
      inherits(conscape_prep, "ConScapeRtools_prep") &&
      identical(conscape_prep$target_mode, "center")) {
    stop(
      "ConScape sensitivity is not valid under target_mode = \"center\".\n",
      "Center-mode tiles use target = center cells and source = full buffered\n",
      "window, which violates ConScape.sensitivity's `target_equal_source = TRUE`\n",
      "precondition inside every tile, and the sum-mosaic reduction is not\n",
      "appropriate for per-cell sensitivity surfaces.\n",
      "\n",
      "Workaround: recreate the prep with target_mode = \"full\":\n",
      "  prep_sens <- conscape_prep(..., target_mode = \"full\", landmark = 1L)\n",
      "\n",
      "Note that tiled sensitivity (even under full mode) is a tile-local\n",
      "approximation to landscape-level sensitivity; for strict correctness,\n",
      "run sensitivity untiled when the landscape fits in a single graph.",
      call. = FALSE
    )
  }

  if (identical(backend, "conscape_dev")) {
    return(run_conscape_dev_backend(
      conscape_prep = conscape_prep,
      out_dir = out_dir,
      target_qualities = hab_target,
      source_qualities = hab_src,
      affinities = mov_prob,
      clear_dir = clear_dir,
      landmark = landmark,
      theta = theta,
      distance_scale = exp_d,
      jl_home = jl_home,
      parallel = parallel,
      workers = workers,
      progress = progress,
      metrics = metrics,
      connectivity_function = connectivity_function,
      cost_function = cost_function,
      sensitivity = sensitivity,
      centersize = centersize,
      buffer = buffer,
      window_shape = window_shape,
      dev_mode = dev_mode,
      batch_grain = batch_grain,
      batch_ext = batch_ext,
      dev_project = dev_project,
      install_dev_conscape = install_dev_conscape,
      dev_conscape_rev = dev_conscape_rev,
      dev_conscape_url = dev_conscape_url,
      blas_threads = blas_threads,
      stop_julia = stop_julia
    ))
  }

  if (isTRUE(stop_julia)) {
    stop_conscape_julia()
    on.exit(stop_conscape_julia(), add = TRUE)
  } else if (isTRUE(parallel) && isFALSE(parallel_R) && conscape_julia_status()) {
    current_threads <- try(suppressMessages(juliaEval("Threads.nthreads()")), silent = TRUE)
    if (inherits(current_threads, "try-error") ||
        !identical(as.integer(current_threads), as.integer(workers))) {
      warning(
        "Reusing an existing Julia session may not honor the requested workers value. ",
        "Julia thread count is fixed when Julia starts. Call conscape_julia_stop() ",
        "or set stop_julia = TRUE before changing threaded Julia settings.",
        call. = FALSE
      )
    }
  }

  if(parallel & isFALSE(parallel_R)){
    Sys.setenv(JULIA_NUM_THREADS = workers)
  }

  conscape_julia_start(jl_home, quiet = TRUE)

  if(dir.exists(out_dir)){
    if(length(list.files(out_dir)) > 0 & isFALSE(clear_dir)){
      stop("\nFiles currently exist in `out_dir`.\nEither manually delete them, or set `clear_dir = TRUE`")
    } else {
      unlink(out_dir,
             recursive = T,
             force = T)
    }
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  single_rast <- FALSE
  target_mask <- NULL

  if(!is.null(conscape_prep) & inherits(conscape_prep, 'ConScapeRtools_prep')){
    target_dir <- conscape_prep$target
    src_dir <- conscape_prep$src
    mov_dir <- conscape_prep$mov
    hab_target <- list.files(target_dir, pattern = "\\.asc$")
    hab_src <- list.files(src_dir, pattern = "\\.asc$")
    mov_prob <- list.files(mov_dir, pattern = "\\.asc$")
    landmark <- conscape_prep$landmark
    target_mask <- terra::rast(file.path(conscape_prep$asc_dir, 'mask', "mask.asc"))
    tile_trim <- conscape_prep$tile_trim
  }

  landmark <- as.integer(landmark)
  if (!is.null(sensitivity) &&
      isTRUE(sensitivity_require_landmark_one) &&
      !identical(landmark, 1L)) {
    stop(
      "ConScape sensitivity currently assumes source and target qualities are identical. ",
      "Use landmark = 1L so coarse_graining() leaves target qualities unchanged. ",
      "When using conscape_prep, create that object with landmark = 1L; ",
      "or set require_landmark_one = FALSE in conscape_sensitivity() for an experimental run.",
      call. = FALSE
    )
  }

  if (inherits(hab_target, 'SpatRaster') |
      inherits(hab_src,    'SpatRaster') |
      inherits(mov_prob,   'SpatRaster')) {

    single_rast <- TRUE
    iter <- 1

    ## optional consistency check
    if (inherits(hab_target, 'SpatRaster') &&
        inherits(hab_src,    'SpatRaster') &&
        inherits(mov_prob,   'SpatRaster')) {

      if (!isTRUE(all.equal(terra::ext(hab_target), terra::ext(hab_src))) ||
          !isTRUE(all.equal(terra::ext(hab_target), terra::ext(mov_prob)))) {
        stop("When passing SpatRasters directly, hab_target, hab_src, and mov_prob must share extent.")
      }
      if (!isTRUE(all.equal(terra::res(hab_target), terra::res(hab_src))) ||
          !isTRUE(all.equal(terra::res(hab_target), terra::res(mov_prob)))) {
        stop("When passing SpatRasters directly, hab_target, hab_src, and mov_prob must share resolution.")
      }
    }

    ## helper to extend a SpatRaster by tile_trim (map units) on all sides
    extend_if_needed <- function(r) {
      if (!inherits(r, "SpatRaster")) return(NULL)
      if (tile_trim <= 0) return(r)
      e0    <- terra::ext(r)
      e_ext <- e0 + tile_trim   # same pattern as original prep
      terra::extend(r, e_ext, fill = NA)
    }

    ## target
    if (inherits(hab_target, 'SpatRaster')) {
      r_use <- extend_if_needed(hab_target)
      if (is.null(r_use)) r_use <- hab_target

      hab_target_write <- file.path(tempdir(), "hab_target.asc")
      suppressWarnings(
        try(terra::writeRaster(r_use, hab_target_write,
                               overwrite = TRUE, NAflag = -9999),
            silent = TRUE)
      )
      hab_target_dir <- dirname(hab_target_write)
    }

    ## source
    if (inherits(hab_src, 'SpatRaster')) {
      r_use <- extend_if_needed(hab_src)
      if (is.null(r_use)) r_use <- hab_src

      hab_src_write <- file.path(tempdir(), "hab_src.asc")
      suppressWarnings(
        try(terra::writeRaster(r_use, hab_src_write,
                               overwrite = TRUE, NAflag = -9999),
            silent = TRUE)
      )
      hab_src_dir <- dirname(hab_src_write)
    }

    ## movement
    if (inherits(mov_prob, 'SpatRaster')) {
      r_use <- extend_if_needed(mov_prob)
      if (is.null(r_use)) r_use <- mov_prob

      mov_prob_write <- file.path(tempdir(), "mov_prob.asc")
      suppressWarnings(
        try(terra::writeRaster(r_use, mov_prob_write,
                               overwrite = TRUE, NAflag = -9999),
            silent = TRUE)
      )
      mov_prob_dir <- dirname(mov_prob_write)
    }

  } else {
    ## directory-based path
    if (is.null(conscape_prep) && isFALSE(file.exists(hab_target))) {
      stop("Specify either a SpatRaster object or path to directory containing *.asc files")
    }

    if (is.null(conscape_prep)) {
      target_dir <- hab_target
      src_dir    <- hab_src
      mov_dir    <- mov_prob
      hab_target <- ifelse(grepl(".asc", basename(hab_target)), hab_target,
                           list.files(hab_target, pattern = "\\.asc$"))
      hab_src    <- ifelse(grepl(".asc", basename(hab_src)), hab_src,
                           list.files(hab_src,    pattern = "\\.asc$"))
      mov_prob   <- ifelse(grepl(".asc", basename(mov_prob)), mov_prob,
                           list.files(mov_prob,   pattern = "\\.asc$"))
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
                      .NA_val = NA_val,
                      .metrics = metrics,
                      .connectivity_function = connectivity_function,
                      .cost_function = cost_function,
                      .sensitivity_wrt = sensitivity_wrt,
                      .sensitivity_method = sensitivity_method,
                      .sensitivity_landscape_measure = sensitivity_landscape_measure,
                      .sensitivity_unitless = sensitivity_unitless,
                      .sensitivity_one_out_of = sensitivity_one_out_of,
                      .sensitivity_diagvalue = sensitivity_diagvalue,
                      .sensitivity_target_equal_source = sensitivity_target_equal_source,
                      .sensitivity_require_landmark_one = sensitivity_require_landmark_one)
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
        invisible(cs_out <- juliaCall_conscape("conscape_batch",
                                                src_dir, mov_dir, target_dir, out_dir,
                                                hab_target, hab_src, mov_prob,
                                                landmark, theta, exp_d, NA_val,
                                                max_retries, progress, blas_threads,
                                                metrics, connectivity_function, cost_function,
                                                sensitivity_wrt, sensitivity_method,
                                                sensitivity_landscape_measure, sensitivity_unitless,
                                                sensitivity_one_out_of, sensitivity_diagvalue,
                                                sensitivity_target_equal_source,
                                                sensitivity_require_landmark_one))

      } else {

        # Distributed -------------------------------------------------------------
        invisible(juliaEval("using Distributed"))
        invisible(juliaEval(paste0("addprocs(",workers,")")))

        conscape_batch_file <- normalizePath(system.file('extdata', 'conscape_batch_distributed.jl', package = 'ConScapeRtools'), winslash = "/")
        conscape_batch_file <- paste0('include("', conscape_batch_file, '")')

        conscape_file <- normalizePath(system.file('extdata', 'conscape.jl', package = 'ConScapeRtools'), winslash = "/")
        conscape_file <- paste0('@everywhere include("', conscape_file, '")')

        invisible(juliaEval(conscape_file))
        invisible(juliaEval(conscape_batch_file))

        max_retries <- 5L
        invisible(cs_out <- juliaCall_conscape("conscape_batch_distributed",
                                                src_dir, mov_dir, target_dir, out_dir,
                                                hab_target, hab_src, mov_prob,
                                                landmark, theta, exp_d, NA_val,
                                                max_retries, progress, blas_threads,
                                                metrics, connectivity_function, cost_function,
                                                sensitivity_wrt, sensitivity_method,
                                                sensitivity_landscape_measure, sensitivity_unitless,
                                                sensitivity_one_out_of, sensitivity_diagvalue,
                                                sensitivity_target_equal_source,
                                                sensitivity_require_landmark_one))
      }
    }


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

      suppressMessages({cs_out <- juliaCall_conscape('conscape',
                                                     normalizePath(hab_src_dir, winslash = "\\"),
                                                     normalizePath(mov_prob_dir, winslash = "\\") ,
                                                     normalizePath(hab_target_dir, winslash = "\\"),
                                                     out_dir,
                                                     "hab_target.asc", "hab_src.asc", "mov_prob.asc",
                                                     # hab_target_file, hab_src_file, mov_prob_file,
                                                     landmark, theta, exp_d, NA_val, iter,
                                                     metrics, connectivity_function, cost_function,
                                                     sensitivity_wrt, sensitivity_method,
                                                     sensitivity_landscape_measure, sensitivity_unitless,
                                                     sensitivity_one_out_of, sensitivity_diagvalue,
                                                     sensitivity_target_equal_source,
                                                     sensitivity_require_landmark_one)})
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

        suppressMessages({cs_out <- juliaCall_conscape('conscape',
                                                       src_dir, mov_dir, target_dir, out_dir,
                                                       r_target, r_source, r_res,
                                                       landmark, theta, exp_d, NA_val, iter,
                                                       metrics, connectivity_function, cost_function,
                                                       sensitivity_wrt, sensitivity_method,
                                                       sensitivity_landscape_measure, sensitivity_unitless,
                                                       sensitivity_one_out_of, sensitivity_diagvalue,
                                                       sensitivity_target_equal_source,
                                                       sensitivity_require_landmark_one)})

        if(isTRUE(progress)){
          setTxtProgressBar(pb,i)
        }
      } ## End for loop
    } ## End single ifelse
    if(exists("pb")) close(pb)
  } ## End parallel ifelse

  expected_outputs <- if (isTRUE(single_rast)) 1L else length(hab_target)
  output_validation <- validate_conscape_outputs(
    out_dir = out_dir,
    output_specs = output_specs,
    expected_tiles = expected_outputs
  )
  if (any(!output_validation$ok)) {
    missing_layers <- paste(output_validation$layer[!output_validation$ok], collapse = ", ")
    if (isTRUE(parallel)) {
      warning(
        "\n\nParallel execution failed for some or all requested output tiles: ",
        missing_layers,
        "\nInspect results carefully.\n",
        call. = FALSE
      )
    } else {
      warning(
        "ConScape execution did not create all requested output tiles: ",
        missing_layers,
        call. = FALSE
      )
    }
  }

  batch_diagnostics <- read_conscape_batch_diagnostics(out_dir)
  run_diagnostics <- list(
    backend = backend,
    parallel = isTRUE(parallel),
    parallel_R = isTRUE(parallel_R),
    distributed = isTRUE(distributed),
    workers = as.integer(workers),
    blas_threads = blas_threads,
    output_validation = output_validation,
    batch = batch_diagnostics
  )

  output_dirs <- lapply(output_specs, function(spec) {
    normalizePath(file.path(out_dir, spec$dir), mustWork = FALSE)
  })
  names(output_dirs) <- vapply(output_specs, `[[`, character(1), "layer")

  make_result_shell <- function() {
    shell <- list(outdirs = output_dirs, diagnostics = run_diagnostics)
    if ("btwn" %in% names(output_dirs)) {
      shell$outdir_btwn <- output_dirs$btwn
    }
    if ("fcon" %in% names(output_dirs)) {
      shell$outdir_fcon <- output_dirs$fcon
    }
    class(shell) <- "ConScapeResults"
    shell
  }

  read_surface_dir <- function(spec) {
    files <- list.files(
      normalizePath(file.path(out_dir, spec$dir), mustWork = FALSE),
      pattern = "\\.asc$|\\.tif$",
      full.names = TRUE
    )
    if (length(files) == 0L) {
      stop("No raster output found for ", spec$layer, call. = FALSE)
    }
    terra::rast(files)
  }

  if(isTRUE(mosaic) & isFALSE(single_rast)){
    # Mosaic reduction is now picked per output spec via pick_mosaic_method().
    # Background:
    #   - "additive" outputs (fcon, btwn, btwn_qweighted, criticality): per-tile
    #     values are partial sums over the disjoint per-tile center targets.
    #     Under center mode they sum across tiles to reconstruct untiled
    #     ("sum"); under legacy full mode they are biased per-tile estimates
    #     that get smoothed by mean mosaic ("mosaic").
    #   - "averageable" outputs (sensitivity / elasticity_*): per-tile values
    #     are tile-local landscape-summary derivatives, NOT pairwise sums.
    #     They are always mean-mosaicked ("mosaic") regardless of target_mode.
    target_mode_for_mosaic <- if (!is.null(conscape_prep) &&
                                  inherits(conscape_prep, "ConScapeRtools_prep")) {
      conscape_prep$target_mode
    } else {
      NULL
    }

    rasters <- lapply(output_specs, function(spec) {
      spec_method <- pick_mosaic_method(spec, target_mode_for_mosaic)
      mosaic_conscape(out_dir = file.path(out_dir, spec$dir),
                      tile_trim = tile_trim,
                      method = spec_method,
                      mask = target_mask,
                      crs = terra::crs(target_mask),
                      chunk_size = mosaic_chunk_size)
    })
    names(rasters) <- vapply(output_specs, `[[`, character(1), "layer")

    if(!is.null(conscape_prep) & inherits(conscape_prep, 'ConScapeRtools_prep')){
      target_mask[target_mask == 0] <- NA
      rasters <- lapply(rasters, function(x) {
        x[is.na(target_mask)] <- NA
        x
      })
    }
    # Record the per-spec methods that were used so audits can see the
    # dispatch decision for every layer (sum vs mosaic). This replaces the
    # old single mosaic_method scalar.
    run_diagnostics$mosaic_method <- vapply(output_specs, function(spec) {
      pick_mosaic_method(spec, target_mode_for_mosaic)
    }, character(1))
    names(run_diagnostics$mosaic_method) <- vapply(output_specs, `[[`,
                                                    character(1), "layer")

    out <- c(rasters, list(outdirs = output_dirs, diagnostics = run_diagnostics))
    if ("btwn" %in% names(output_dirs)) out$outdir_btwn <- output_dirs$btwn
    if ("fcon" %in% names(output_dirs)) out$outdir_fcon <- output_dirs$fcon
    class(out) <- "ConScapeResults"
  } else {
    out <- make_result_shell()

    if (isTRUE(single_rast)) {
      rasters <- lapply(output_specs, read_surface_dir)
      names(rasters) <- vapply(output_specs, `[[`, character(1), "layer")

      ## determine original (pre-extension) extent
      orig_ext <- NULL
      if (inherits(hab_target, "SpatRaster")) {
        orig_ext <- terra::ext(hab_target)
      } else if (inherits(hab_src, "SpatRaster")) {
        orig_ext <- terra::ext(hab_src)
      } else if (inherits(mov_prob, "SpatRaster")) {
        orig_ext <- terra::ext(mov_prob)
      }

      if (!is.null(orig_ext)) {
        rasters <- lapply(rasters, terra::crop, y = orig_ext)
      }

      if (inherits(mov_prob, "SpatRaster")) {
        rasters <- lapply(rasters, function(x) {
          terra::crs(x) <- terra::crs(mov_prob)
          x
        })
      }

      out <- terra::rast(rasters)
      names(out) <- names(rasters)
      attr(out, "ConScapeRtools_diagnostics") <- run_diagnostics
    }
  }

  # juliaEval("GC.gc()")
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
                        .NA_val,
                        .metrics,
                        .connectivity_function,
                        .cost_function,
                        .sensitivity_wrt,
                        .sensitivity_method,
                        .sensitivity_landscape_measure,
                        .sensitivity_unitless,
                        .sensitivity_one_out_of,
                        .sensitivity_diagvalue,
                        .sensitivity_target_equal_source,
                        .sensitivity_require_landmark_one){
  iter <- paste0('-iter_', i)
  invisible({cs_out <- juliaCall_conscape('conscape',
                                          .src_dir, .mov_dir, .target_dir, .out_dir,
                                          .hab_target[i], .hab_src[i], .mov_prob[i],
                                          .landmark, .theta, .exp_d, .NA_val, iter,
                                          .metrics, .connectivity_function, .cost_function,
                                          .sensitivity_wrt, .sensitivity_method,
                                          .sensitivity_landscape_measure,
                                          .sensitivity_unitless, .sensitivity_one_out_of,
                                          .sensitivity_diagvalue,
                                          .sensitivity_target_equal_source,
                                          .sensitivity_require_landmark_one)})
  # juliaEval("GC.gc()")
}
