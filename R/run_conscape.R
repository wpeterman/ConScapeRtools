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
#'   ConScape's `target_qualities` layer — cells that can *receive*
#'   connectivity. The legacy name `hab_target` is accepted as a
#'   backwards-compatible alias.
#' @param source_qualities Either (i) a path to a directory containing
#'   source-quality tiles (`*.asc`) or (ii) a single `SpatRaster` used as
#'   ConScape's `source_qualities` layer — cells that *generate* connectivity.
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
#'   provided, `landmark` is taken from that object. Set `landmark = 1L` for
#'   sensitivity runs so that target qualities are not altered by coarse graining.
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
#' @param distributed Logical. When `parallel = TRUE` and
#'   `parallel_R = FALSE`, controls whether Julia uses distributed
#'   workers (`TRUE`) or multithreading (`FALSE`, default).
#' @param progress Logical. If `TRUE` (default), progress information is
#'   printed for serial and R-parallel runs. Progress reporting from Julia
#'   parallel runs may be less detailed.
#' @param parallel_R Logical. If `TRUE`, tiles are processed in parallel
#'   using R (`future`/`future.apply`). If `FALSE` (default), parallelism
#'   is handled entirely within Julia.
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
#' Tiled sensitivity uses the same tiling and mosaicking machinery as the
#' default metric outputs. It is a tile-local approximation to full-landscape
#' sensitivity, so choose a `tile_trim` large enough that paths crossing tile
#' edges have negligible influence on interior retained cells. For sensitivity
#' runs, use `landmark = 1L` so that coarse-graining does not alter target
#' qualities (the ConScape sensitivity API currently requires target quality
#' to equal source quality).
#'
#' @return
#' When tiles are used and `mosaic = TRUE`, returns an object of class
#' `"ConScapeResults"` — a named list containing:
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
#' ## Untiled run — only suitable for small to moderate rasters
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
                         tile_trim = 0,
                         cost_function = "minuslog",
                         connectivity_function = "expected_cost",
                         metrics = c("betweenness_kweighted", "connected_habitat"),
                         sensitivity = NULL,
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

  landmark <- as.integer(landmark)
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
  } else {
    dir.create(out_dir, recursive = TRUE)
  }

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
      "Use landmark = 1L so coarse_graining() leaves target qualities unchanged, ",
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
                                                max_retries, progress,
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
                                                max_retries, progress,
                                                metrics, connectivity_function, cost_function,
                                                sensitivity_wrt, sensitivity_method,
                                                sensitivity_landscape_measure, sensitivity_unitless,
                                                sensitivity_one_out_of, sensitivity_diagvalue,
                                                sensitivity_target_equal_source,
                                                sensitivity_require_landmark_one))
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

  output_dirs <- lapply(output_specs, function(spec) {
    normalizePath(file.path(out_dir, spec$dir), mustWork = FALSE)
  })
  names(output_dirs) <- vapply(output_specs, `[[`, character(1), "layer")

  make_result_shell <- function() {
    shell <- list(outdirs = output_dirs)
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
    rasters <- lapply(output_specs, function(spec) {
      mosaic_conscape(out_dir = file.path(out_dir, spec$dir),
                      tile_trim = tile_trim,
                      method = 'mosaic',
                      mask = target_mask,
                      crs = terra::crs(target_mask))
    })
    names(rasters) <- vapply(output_specs, `[[`, character(1), "layer")

    if(!is.null(conscape_prep) & inherits(conscape_prep, 'ConScapeRtools_prep')){
      target_mask[target_mask == 0] <- NA
      rasters <- lapply(rasters, function(x) {
        x[is.na(target_mask)] <- NA
        x
      })
    }

    out <- c(rasters, list(outdirs = output_dirs))
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
