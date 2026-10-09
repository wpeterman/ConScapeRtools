test_that("disabled experimental backend preserves outputs and never starts Julia", {
  r <- make_test_raster(n = 4, vals = 1)
  out_dir <- tempfile("disabled-dev-backend-")
  dir.create(out_dir)
  sentinel <- file.path(out_dir, "keep.txt")
  writeLines("must survive", sentinel)
  old_env <- Sys.getenv(
    c("JULIA_BINDIR", "JULIA_NUM_THREADS", "JULIA_PROJECT"),
    unset = NA_character_
  )

  expect_error(testthat::with_mocked_bindings(
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
    prepare_conscape_dev_project = function(...) stop("must not install"),
    conscape_julia_start = function(...) stop("must not start Julia"),
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
  ), "disabled in this release")
  expect_identical(readLines(sentinel), "must survive")
  expect_identical(list.files(out_dir), "keep.txt")
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
