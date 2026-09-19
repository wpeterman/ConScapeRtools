test_that("experimental backend writes inputs and returns documented diagnostics", {
  r <- make_test_raster(n = 4, vals = 1)
  out_dir <- file.path(tempdir(), "direct-dev-backend")
  old_env <- Sys.getenv(
    c("JULIA_BINDIR", "JULIA_NUM_THREADS", "JULIA_PROJECT"),
    unset = NA_character_
  )

  out <- testthat::with_mocked_bindings(
    getFromNamespace("run_conscape_dev_backend", "ConScapeRtools")(
      conscape_prep = NULL,
      out_dir = out_dir,
      target_qualities = r,
      source_qualities = r,
      affinities = r,
      clear_dir = TRUE,
      landmark = 1L,
      theta = 0.1,
      distance_scale = 10,
      jl_home = "C:/Julia/bin",
      parallel = TRUE,
      workers = 2L,
      progress = FALSE,
      metrics = "connected_habitat",
      connectivity_function = "expected_cost",
      cost_function = "minuslog",
      sensitivity = NULL,
      centersize = 2L,
      buffer = 1L,
      window_shape = "square",
      dev_mode = "windowed",
      batch_grain = NULL,
      batch_ext = ".tif",
      dev_project = "C:/Julia/dev-project",
      install_dev_conscape = FALSE,
      dev_conscape_rev = "9aa05cc0b0c22b9d815d3051925010a2344eada0",
      dev_conscape_url = "https://github.com/ConScape/ConScape.jl",
      blas_threads = 1L,
      stop_julia = TRUE
    ),
    prepare_conscape_dev_project = function(...) "C:/Julia/dev-project",
    conscape_julia_start = function(...) invisible(NULL),
    stop_conscape_julia = function() invisible(NULL),
    juliaEval = function(...) invisible(NULL),
    juliaLet = function(expr, ...) {
      args <- list(...)
      output_dir <- file.path(args$out_dir, "conscape_dev_windowed")
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
      output <- terra::rast(args$target_file)
      terra::writeRaster(
        output,
        file.path(output_dir, "connected_habitat.tif"),
        overwrite = TRUE
      )
      invisible(NULL)
    },
    .package = "ConScapeRtools"
  )

  expect_s4_class(out, "SpatRaster")
  diagnostics <- attr(out, "ConScapeRtools_diagnostics")
  expect_identical(diagnostics$backend, "conscape_dev")
  expect_identical(diagnostics$dev_mode, "windowed")
  expect_identical(diagnostics$workers, 2L)
  expect_identical(diagnostics$centersize, 2L)
  expect_true(all(file.exists(diagnostics$output_files)))
  expect_true(all(file.exists(file.path(
    out_dir,
    "dev_inputs",
    c("target_qualities.tif", "source_qualities.tif", "affinities.tif")
  ))))
  expect_identical(
    Sys.getenv(names(old_env), unset = NA_character_),
    old_env
  )
})

test_that("experimental backend validates inputs before writing", {
  r <- make_test_raster(n = 3, vals = 1)
  call_backend <- function(...) {
    args <- list(
      conscape_prep = NULL,
      out_dir = file.path(tempdir(), "invalid-dev-backend"),
      target_qualities = r,
      source_qualities = r,
      affinities = r,
      clear_dir = TRUE,
      landmark = 1L,
      theta = 0.1,
      distance_scale = 10,
      jl_home = "C:/Julia/bin",
      parallel = FALSE,
      workers = 1L,
      progress = FALSE,
      metrics = "connected_habitat",
      connectivity_function = "expected_cost",
      cost_function = "minuslog",
      sensitivity = NULL,
      centersize = 2L,
      buffer = 1L,
      window_shape = "square",
      dev_mode = "windowed",
      batch_grain = NULL,
      batch_ext = ".tif",
      dev_project = NULL,
      install_dev_conscape = FALSE,
      dev_conscape_rev = "9aa05cc0b0c22b9d815d3051925010a2344eada0",
      dev_conscape_url = "https://github.com/ConScape/ConScape.jl",
      blas_threads = 1L,
      stop_julia = TRUE
    )
    replacements <- list(...)
    args[names(replacements)] <- replacements
    do.call(
      getFromNamespace("run_conscape_dev_backend", "ConScapeRtools"),
      args
    )
  }

  expect_error(call_backend(centersize = 2.5), "centersize")
  expect_error(call_backend(buffer = 1.5), "buffer")
  expect_error(call_backend(batch_grain = 1.5), "batch_grain")
  expect_error(call_backend(connectivity_function = "least_cost"), "expected_cost")
})
