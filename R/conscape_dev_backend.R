#' Experimental ConScape dev backend
#' @noRd
run_conscape_dev_backend <- function(conscape_prep,
                                     out_dir,
                                     target_qualities,
                                     source_qualities,
                                     affinities,
                                     clear_dir,
                                     landmark,
                                     theta,
                                     distance_scale,
                                     jl_home,
                                     parallel,
                                     workers,
                                     progress,
                                     metrics,
                                     connectivity_function,
                                     cost_function,
                                     sensitivity,
                                     centersize,
                                     buffer,
                                     window_shape,
                                     dev_mode,
                                     batch_grain,
                                     batch_ext,
                                     dev_project,
                                     install_dev_conscape,
                                     dev_conscape_rev,
                                     dev_conscape_url,
                                     blas_threads,
                                     stop_julia) {
  if (!is.null(conscape_prep)) {
    stop(
      "backend = \"conscape_dev\" currently requires direct SpatRaster inputs. ",
      "Use target_qualities, source_qualities, and affinities rather than conscape_prep.",
      call. = FALSE
    )
  }
  if (!inherits(target_qualities, "SpatRaster") ||
      !inherits(source_qualities, "SpatRaster") ||
      !inherits(affinities, "SpatRaster")) {
    stop(
      "backend = \"conscape_dev\" currently requires SpatRaster inputs for ",
      "target_qualities, source_qualities, and affinities.",
      call. = FALSE
    )
  }
  if (is.null(centersize) || is.null(buffer)) {
    stop("backend = \"conscape_dev\" requires centersize and buffer in cells.", call. = FALSE)
  }
  if (!is.numeric(centersize) || length(centersize) != 1L ||
      is.na(centersize) || centersize < 1) {
    stop("centersize must be a positive integer cell count.", call. = FALSE)
  }
  if (!is.numeric(buffer) || length(buffer) != 1L ||
      is.na(buffer) || buffer < 0) {
    stop("buffer must be a non-negative integer cell count.", call. = FALSE)
  }
  if (!identical(dev_mode, "windowed") && !identical(dev_mode, "batch")) {
    stop("dev_mode must be either \"windowed\" or \"batch\".", call. = FALSE)
  }
  if (!is.null(batch_grain) &&
      (!is.numeric(batch_grain) || length(batch_grain) != 1L ||
       is.na(batch_grain) || batch_grain < 1)) {
    stop("batch_grain must be NULL or a positive integer.", call. = FALSE)
  }
  if (!is.character(batch_ext) || length(batch_ext) != 1L ||
      !nzchar(batch_ext)) {
    stop("batch_ext must be a non-empty character string.", call. = FALSE)
  }
  if (!startsWith(batch_ext, ".")) batch_ext <- paste0(".", batch_ext)
  if (!is.null(sensitivity)) {
    stop("backend = \"conscape_dev\" does not yet support sensitivity outputs.", call. = FALSE)
  }
  if (!identical(connectivity_function, "expected_cost")) {
    stop("backend = \"conscape_dev\" currently supports connectivity_function = \"expected_cost\".", call. = FALSE)
  }

  if (dir.exists(out_dir)) {
    if (length(list.files(out_dir)) > 0L && isFALSE(clear_dir)) {
      stop("\nFiles currently exist in `out_dir`.\nEither manually delete them, or set `clear_dir = TRUE`")
    }
    unlink(out_dir, recursive = TRUE, force = TRUE)
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  active_dev_project <- prepare_conscape_dev_project(
    jl_home = jl_home,
    dev_project = dev_project,
    install_dev_conscape = install_dev_conscape,
    dev_conscape_rev = dev_conscape_rev,
    dev_conscape_url = dev_conscape_url
  )

  if (isTRUE(stop_julia)) {
    stop_conscape_julia()
    on.exit(stop_conscape_julia(), add = TRUE)
  }
  if (isTRUE(parallel)) {
    Sys.setenv(JULIA_NUM_THREADS = workers)
  }
  conscape_julia_start(jl_home, quiet = TRUE)

  input_dir <- file.path(out_dir, "dev_inputs")
  dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
  target_file <- file.path(input_dir, "target_qualities.tif")
  source_file <- file.path(input_dir, "source_qualities.tif")
  affinities_file <- file.path(input_dir, "affinities.tif")

  terra::writeRaster(target_qualities, target_file, overwrite = TRUE, NAflag = -9999)
  terra::writeRaster(source_qualities, source_file, overwrite = TRUE, NAflag = -9999)
  terra::writeRaster(affinities, affinities_file, overwrite = TRUE, NAflag = -9999)

  helper <- normalizePath(
    system.file("extdata", "conscape_dev_windowed.jl", package = "ConScapeRtools"),
    winslash = "/"
  )
  invisible(juliaEval(paste0('include("', helper, '"); nothing')))

  invisible(juliaLet(
    paste0(
      "Base.invokelatest(",
      "Base.invokelatest(getfield, Main, :conscape_dev_run), ",
      "target_file, source_file, affinities_file, out_dir, ",
      "landmark, theta, distance_scale, centersize, buffer, window_shape, ",
      "dev_mode, batch_grain, batch_ext, threaded, blas_threads, metrics, ",
      "cost_function, progress)"
    ),
    target_file = normalizePath(target_file, winslash = "/"),
    source_file = normalizePath(source_file, winslash = "/"),
    affinities_file = normalizePath(affinities_file, winslash = "/"),
    out_dir = normalizePath(out_dir, winslash = "/"),
    landmark = as.integer(landmark),
    theta = theta,
    distance_scale = distance_scale,
    centersize = as.integer(centersize),
    buffer = as.integer(buffer),
    window_shape = window_shape,
    dev_mode = dev_mode,
    batch_grain = if (is.null(batch_grain)) NULL else as.integer(batch_grain),
    batch_ext = batch_ext,
    threaded = isTRUE(parallel),
    blas_threads = as.integer(blas_threads),
    metrics = metrics,
    cost_function = cost_function,
    progress = isTRUE(progress)
  ))

  output_dir <- file.path(out_dir, paste0("conscape_dev_", dev_mode))
  output_files <- list.files(
    output_dir,
    pattern = "\\.tif$|\\.tiff$|\\.asc$",
    recursive = FALSE,
    full.names = TRUE
  )
  output_files <- sort(output_files)
  if (length(output_files) == 0L) {
    stop("ConScape dev backend completed without writing raster outputs.", call. = FALSE)
  }

  rasters <- terra::rast(output_files)
  diagnostics <- list(
    backend = "conscape_dev",
    dev_mode = dev_mode,
    parallel = isTRUE(parallel),
    workers = as.integer(workers),
    blas_threads = as.integer(blas_threads),
    centersize = as.integer(centersize),
    buffer = as.integer(buffer),
    window_shape = window_shape,
    batch_grain = if (is.null(batch_grain)) NA_integer_ else as.integer(batch_grain),
    batch_ext = batch_ext,
    dev_project = if (is.null(active_dev_project)) NA_character_ else active_dev_project,
    install_dev_conscape = isTRUE(install_dev_conscape),
    dev_conscape_rev = dev_conscape_rev,
    dev_conscape_url = dev_conscape_url,
    output_files = normalizePath(output_files, mustWork = FALSE)
  )
  attr(rasters, "ConScapeRtools_diagnostics") <- diagnostics
  rasters
}
