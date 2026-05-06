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
#'   supply `hab_target`, `hab_src`, and `mov_prob` directly.
#' @param out_dir Directory where final ConScape outputs will be written.
#'   Subdirectories `"btwn"` and `"fcon"` are created inside `out_dir`.
#' @param hab_target,target_qualities Either (i) a path to a directory
#'   containing target-quality tiles (`*.asc`) or (ii) a single `SpatRaster`
#'   used as the ConScape `target_qualities` layer. `hab_target` is retained
#'   as a backwards-compatible alias.
#' @param hab_src,source_qualities Either (i) a path to a directory containing
#'   source-quality tiles (`*.asc`) or (ii) a single `SpatRaster` used as the
#'   ConScape `source_qualities` layer. `hab_src` is retained as a
#'   backwards-compatible alias.
#' @param mov_prob,affinities Either (i) a path to a directory containing
#'   affinity/permeability tiles (`*.asc`) or (ii) a single `SpatRaster` used
#'   to build ConScape `affinities` with `graph_matrix_from_raster()`.
#'   `mov_prob` is retained as a backwards-compatible alias.
#' @param clear_dir Logical. If `TRUE` (default), any existing contents
#'   of `out_dir` are removed before writing new results. If `FALSE` and
#'   `out_dir` is not empty, the function stops with an error.
#' @param landmark Landmark value used for ConScape's `coarse_graining`
#'   (default `10L`). When `conscape_prep` is provided, this is taken
#'   from that object and does not need to be set manually.
#' @param theta Parameter controlling the amount of randomness in paths.
#'   As `theta` approaches 0, movement is increasingly random; as `theta`
#'   becomes large, paths approach deterministic least-cost paths.
#'   Default is `0.01`.
#' @param exp_d Numerator of the exponential decay function used when
#'   transforming distance in the movement grid (i.e. `exp(-dist / exp_d)`).
#'   Default is `150`. This is typically set from [tile_design()].
#' @param distance_scale Clearer alias for `exp_d`. Use one or the other.
#' @param NA_val Backwards-compatible argument retained for older workflows.
#'   The current Julia runner treats `NA` affinities as absent graph cells and
#'   `NA` source/target qualities as zero quality.
#' @param mosaic Logical (default `TRUE`). When running on tiles
#'   (via `conscape_prep` or directory inputs), tile-level outputs are
#'   mosaicked into continuous rasters using [mosaic_conscape()]. When
#'   running on a single raster (`SpatRaster` inputs), `mosaic` is
#'   ignored.
#' @param jl_home Path to the `bin` directory where the Julia executable
#'   is installed. Used by `JuliaConnectoR` to initialize ConScape.
#' @param parallel Logical. If `FALSE` (default), tiles are processed
#'   serially. If `TRUE`, tiles are processed in parallel (either via R
#'   or Julia, depending on `parallel_R` and `distributed`).
#' @param workers Integer number of parallel workers to use when
#'   `parallel = TRUE`. Default is `1`.
#' @param distributed Logical. When `parallel = TRUE` and
#'   `parallel_R = FALSE`, controls whether Julia uses distributed
#'   workers (`TRUE`) or multithreading (`FALSE`).
#' @param progress Logical. If `TRUE`, progress information is printed
#'   for serial and R-parallel runs. Progress reporting from Julia
#'   parallel runs may be less consistent.
#' @param parallel_R Logical. If `TRUE`, tiles are processed in parallel
#'   using R (`future/future.apply`). If `FALSE` (default), parallelism
#'   is handled entirely within Julia.
#' @param tile_trim Width of the overlapping border (in map units) that
#'   will be trimmed from each tile when mosaicking. When `conscape_prep`
#'   is supplied, this parameter is taken from that object and passed to
#'   [mosaic_conscape()]. When running on a single `SpatRaster`, `tile_trim`
#'   controls how far the input rasters are extended (buffered) on all
#'   sides before ConScape is run; the output is subsequently cropped
#'   back to the original extent so that tiled and untiled runs are
#'   comparable.
#' @param cost_function Cost transformation used by `ConScape.Grid(costs = ...)`.
#'   Supported named values are `"minuslog"` (default), `"inverse"`,
#'   `"odds_against"`, `"odds_for"`, and `"expminus"`. A Julia lambda string
#'   such as `"x -> -log(x)"` can also be supplied, in which case a cost matrix
#'   is created with `ConScape.mapnz()`.
#' @param connectivity_function ConScape distance/proximity function used by
#'   metrics that require a landscape matrix. Supported values are
#'   `"expected_cost"` (default), `"least_cost_distance"`,
#'   `"free_energy_distance"`, `"survival_probability"`, and
#'   `"power_mean_proximity"`.
#' @param metrics Character vector of raster-valued ConScape metrics to write.
#'   Supported values are `"connected_habitat"`, `"betweenness_kweighted"`,
#'   `"betweenness_qweighted"`, and `"criticality"`. The aliases `"fcon"` and
#'   `"btwn"` are accepted.
#' @param sensitivity Optional sensitivity specification created by
#'   [conscape_sensitivity()]. `TRUE` uses the default quality and linked
#'   affinity-cost elasticities; a character vector is interpreted as `wrt`
#'   values passed to [conscape_sensitivity()].
#'
#' @details
#' In typical workflows, data are prepared using [conscape_prep()] and
#' the resulting object is supplied via `conscape_prep`. This ensures
#' that source, target, and movement tiles share the same tiling scheme,
#' that coarse-graining landmarks are aligned, and that a mask and
#' `tile_trim` value are available for mosaicking.
#'
#' If `conscape_prep` is not provided, the function can operate on:
#' * directories of `*.asc` tiles (`hab_target`, `hab_src`, `mov_prob`),
#'   or
#' * three `SpatRaster` objects of identical extent and resolution
#'   (`hab_target`, `hab_src`, `mov_prob`), in which case ConScape is
#'   run on the full (untiled) landscape.
#'
#' When `parallel = TRUE`, either R-level parallelization (with
#' `parallel_R = TRUE`) or Julia-level parallelization (with
#' `parallel_R = FALSE`) is used. The choice between Julia threading and
#' Julia distributed workers is controlled by `distributed`.
#'
#' The default run writes `connected_habitat` and `betweenness_kweighted`,
#' which correspond to the historical `fcon` and `btwn` outputs. Additional
#' metrics are requested with `metrics`. All metric names are passed to the
#' bundled Julia runner, which computes each requested surface on the same
#' `GridRSP` object.
#'
#' Sensitivity outputs are requested with [conscape_sensitivity()]. These
#' surfaces are computed after the `GridRSP` object is built and are written to
#' directories named for the output, for example `"elasticity_quality"` or
#' `"elasticity_affinity_cost_linked"`. Analytical sensitivity requires a
#' ConScape installation that provides `ConScape.sensitivity`; the current CRAN
#' release of ConScape may not yet include that function. When the function is
#' unavailable, the Julia runner stops with an explicit message rather than
#' silently returning an approximate substitute.
#'
#' Tiled sensitivity uses the same tiling and mosaicking machinery as the
#' default betweenness and connected-habitat outputs. It should be treated as a
#' tile-local approximation to the full-landscape sensitivity, so choose a
#' `tile_trim` large enough that paths crossing tile edges have negligible
#' influence on the interior cells retained by [mosaic_conscape()]. For
#' sensitivity runs, use `landmark = 1L` unless you have a specific reason to
#' test coarse-grained target qualities, because the ConScape sensitivity API
#' currently assumes target quality equals source quality.
#'
#' @return
#' If `mosaic = TRUE` and tiles are used (`!single_rast` internally),
#' returns an object of class `"ConScapeResults"` containing:
#'
#' * `btwn` – `SpatRaster` of k-weighted betweenness when requested,
#' * `fcon` – `SpatRaster` of connected habitat when requested,
#' * additional metric or sensitivity rasters named by output layer,
#' * `outdirs` – named list of tile-level output directories,
#' * `outdir_btwn` and `outdir_fcon` – backwards-compatible paths for the
#'   default betweenness and connected-habitat outputs when requested.
#'
#' When run on a single `SpatRaster` (no tiling), the function returns a
#' two-layer `SpatRaster` with layers named `"btwn"` and `"fcon"`.
#'
#' @seealso [conscape_prep()], [tile_design()], [mosaic_conscape()],
#'   [conscape_sensitivity()], [conscape_api_slots()]
#' @export
#' @examples
#' \dontrun{
#' library(ConScapeRtools)
#'
#' ## Import data
#' s <- system.file("extdata", "suitability.asc", package = "ConScapeRtools")
#' source <- terra::rast(s)
#'
#' a <- system.file("extdata", "affinity.asc", package = "ConScapeRtools")
#' resist <- terra::rast(a)
#'
#' jl_home <- "/path/to/julia/bin" # Update
#'
#' td <- tile_design(r_mov = resist,
#'                   r_target = source,
#'                   max_d = 8500,
#'                   theta = 0.15,
#'                   jl_home = jl_home)
#'
#' ## Tile dimension
#' tile_d <- td$tile_d
#'
#' # How much to trim tiles
#' tile_trim <- td$tile_trim
#'
#' # Makes computation more efficient
#' landmark <- 5L # Must be an integer, not numeric
#'
#' # Controls level of randomness of paths
#' theta <- td$theta
#'
#' # Controls rate of decay with distance
#' exp_d <- td$exp_d
#'
#' ## Prepare data for analysis
#' prep <- conscape_prep(tile_d = tile_d,
#'                       tile_trim = tile_trim,
#'                       r_target = source,
#'                       r_mov = resist,
#'                       r_src = source,
#'                       clear_dir = T,
#'                       landmark = landmark)
#'
#' ## Run ConScape
#' ## No parallelization
#' cs_run.serial <- run_conscape(out_dir = file.path(prep$asc_dir,"results"),
#'                               conscape_prep = prep,
#'                               theta = theta,
#'                               exp_d = exp_d,
#'                               jl_home = jl_home,
#'                               parallel = F)
#'
#' ## Parallel within R
#' cs_run.r <- run_conscape(out_dir = file.path(prep$asc_dir,"results"),
#'                          conscape_prep = prep,
#'                          theta = theta,
#'                          exp_d = exp_d,
#'                          jl_home = jl_home,
#'                          parallel = T,
#'                          workers = 4,
#'                          parallel_R = TRUE)
#' ## Threaded parallel
#' cs_run.thread <- run_conscape(out_dir = file.path(prep$asc_dir,"results"),
#'                               conscape_prep = prep,
#'                               theta = theta,
#'                               exp_d = exp_d,
#'                               jl_home = jl_home,
#'                               parallel = T,
#'                               workers = 4)
#'
#' ## Distributed parallel
#' cs_run.dist <- run_conscape(out_dir = file.path(prep$asc_dir,"results"),
#'                             conscape_prep = prep,
#'                             theta = theta,
#'                             exp_d = exp_d,
#'                             jl_home = jl_home,
#'                             parallel = T,
#'                             workers = 4,
#'                             distributed = TRUE)
#' plot(cs_run.dist)
#'
#' ## No tiling --> Only attempt with small to moderate sized rasters!
#' cs_run <- run_conscape(out_dir = file.path(prep$asc_dir,"results"),
#'                        hab_target = source,
#'                        hab_src = source,
#'                        mov_prob = resist,
#'                        theta = theta,
#'                        exp_d = exp_d,
#'                        landmark = landmark,
#'                        jl_home = jl_home)
#'
#' plot(cs_run.serial$fcon - cs_run$fcon, main = "Tiled minus untiled")
#' plot(cs_run.serial$fcon, main = "Tiled")
#' plot(cs_run$fcon, main = "Untiled")
#' layerCor(c(cs_run$fcon, cs_run.serial$fcon), fun = 'cor')
#' layerCor(c(cs_run$btwn, cs_run.serial$btwn), fun = 'cor')
#' }
#' @author Bill Peterman
#' @importFrom JuliaConnectoR juliaImport juliaEval juliaSetupOk juliaGet juliaCall stopJulia
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
                         parallel_R = FALSE,
                         tile_trim = 0,
                         target_qualities = NULL,
                         source_qualities = NULL,
                         affinities = NULL,
                         distance_scale = NULL,
                         cost_function = "minuslog",
                         connectivity_function = "expected_cost",
                         metrics = c("betweenness_kweighted", "connected_habitat"),
                         sensitivity = NULL){
  if (!is.null(target_qualities)) {
    if (!is.null(hab_target)) {
      stop("Use either hab_target or target_qualities, not both.", call. = FALSE)
    }
    hab_target <- target_qualities
  }
  if (!is.null(source_qualities)) {
    if (!is.null(hab_src)) {
      stop("Use either hab_src or source_qualities, not both.", call. = FALSE)
    }
    hab_src <- source_qualities
  }
  if (!is.null(affinities)) {
    if (!is.null(mov_prob)) {
      stop("Use either mov_prob or affinities, not both.", call. = FALSE)
    }
    mov_prob <- affinities
  }
  if (!is.null(distance_scale)) {
    if (!missing(exp_d)) {
      stop("Use either exp_d or distance_scale, not both.", call. = FALSE)
    }
    exp_d <- distance_scale
  }
  if (!is.numeric(exp_d) || length(exp_d) != 1L || is.na(exp_d) || exp_d <= 0) {
    stop("distance_scale/exp_d must be a single positive numeric value.", call. = FALSE)
  }
  if (!is.numeric(theta) || length(theta) != 1L || is.na(theta) || theta <= 0) {
    stop("theta must be a single positive numeric value.", call. = FALSE)
  }
  if (!is.numeric(landmark) || length(landmark) != 1L || is.na(landmark) || landmark < 1) {
    stop("landmark must be a positive integer.", call. = FALSE)
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

  stopJulia()
  if(parallel & isFALSE(parallel_R)){
    Sys.setenv(JULIA_NUM_THREADS = workers)
  }

  if(!juliaSetupOk()){
    Sys.setenv(JULIA_BINDIR = jl_home)
  }

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
    if (isFALSE(file.exists(hab_target))) {
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
        invisible(cs_out <- juliaCall("conscape_batch",
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
        invisible(cs_out <- juliaCall("conscape_batch_distributed",
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

      suppressMessages({cs_out <- juliaCall('conscape',
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

        suppressMessages({cs_out <- juliaCall('conscape',
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
                      method = 'merge',
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
  invisible({cs_out <- juliaCall('conscape',
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
