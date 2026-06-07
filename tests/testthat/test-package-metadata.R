test_that("installed example data paths resolve", {
  expect_true(file.exists(system.file("extdata", "affinity.asc", package = "ConScapeRtools")))
  expect_true(file.exists(system.file("extdata", "suitability.asc", package = "ConScapeRtools")))
  expect_true(file.exists(system.file("extdata", "patches.shp", package = "ConScapeRtools")))
})

test_that("Julia helper scripts referenced by run_conscape are installed", {
  expect_true(file.exists(system.file("extdata", "conscape.jl", package = "ConScapeRtools")))
  expect_true(file.exists(system.file("extdata", "conscape_batch.jl", package = "ConScapeRtools")))
  expect_true(file.exists(system.file("extdata", "conscape_batch_distributed.jl", package = "ConScapeRtools")))
  expect_true(file.exists(system.file("extdata", "conscape_dev_windowed.jl", package = "ConScapeRtools")))
})

test_that("experimental ConScape dev helper exposes windowed and batch modes", {
  dev_file <- system.file("extdata", "conscape_dev_windowed.jl", package = "ConScapeRtools")
  dev <- paste(readLines(dev_file, warn = FALSE), collapse = "\n")

  expect_match(dev, "function conscape_dev_run")
  expect_match(dev, "mode in \\(\"windowed\", \"batch\"\\)")
  expect_match(dev, "ConScape\\.WindowedProblem")
  expect_match(dev, "ConScape\\.BatchProblem")
  expect_match(dev, "ConScape\\.assess")
  expect_match(dev, "function _dev_mosaic_batch")
  expect_match(dev, "ConScape\\.batch_paths")
  expect_match(dev, "Rasters\\.mosaic")
  expect_match(dev, "ConScape\\.coarse_graining")
  expect_match(dev, "ConScape\\.Problem")
})

test_that("worked ConScape dev backend example is available", {
  example_file <- system.file(
    "examples", "conscape_dev_backend_example.R",
    package = "ConScapeRtools"
  )
  expect_true(file.exists(example_file))
  example <- paste(readLines(example_file, warn = FALSE), collapse = "\n")

  expect_match(example, "conscape_dev_backend_setup")
  expect_match(example, "window_design")
  expect_match(example, 'backend = "conscape_dev"', fixed = TRUE)
  expect_match(example, 'dev_mode = "windowed"', fixed = TRUE)
  expect_match(example, 'dev_mode = "batch"', fixed = TRUE)
  expect_match(example, "suitability.asc")
  expect_match(example, "affinity.asc")
})

test_that("Tiled ConScape vignette documents validation and dev-backend status", {
  # This test inspects the vignette source. It only runs against the source
  # checkout (devtools::test() or R CMD check on the .tar.gz), not against an
  # installed binary, because installed packages keep vignettes under
  # `doc/` (not `vignettes/`) and the original .Rmd source is not always
  # shipped. Use testthat::local_edition(3) `skip_*` patterns to bail out
  # gracefully when the source is unreachable.
  vignette_file <- testthat::test_path(
    "..", "..", "vignettes", "Tiled_ConScape_Validation_and_Performance.Rmd"
  )
  skip_if_not(
    file.exists(vignette_file),
    paste0("Vignette source not available at ", vignette_file,
           " (expected when running against an installed package).")
  )
  vignette <- paste(readLines(vignette_file, warn = FALSE), collapse = "\n")

  # Title + structure
  expect_match(vignette, "Tiled ConScape: Validation and Performance", fixed = TRUE)
  expect_match(vignette, "# Identity Validation", fixed = TRUE)
  expect_match(vignette, "Wall-clock timing", fixed = TRUE)
  expect_match(vignette, "# Decision Tree", fixed = TRUE)
  expect_match(vignette, "# Troubleshooting", fixed = TRUE)

  # Stable-backend run examples (both parameterizations)
  expect_match(vignette, "classic_prep <- conscape_prep", fixed = TRUE)
  expect_match(vignette, "center_prep <- conscape_prep", fixed = TRUE)

  # Dev backend is documented but explicitly flagged as broken upstream
  expect_match(vignette, "# Experimental: Dev Backend", fixed = TRUE)
  expect_match(vignette, "Currently Broken Upstream", fixed = TRUE)
  expect_match(vignette, "MethodError: no method matching write", fixed = TRUE)
})

test_that("threaded Julia batch avoids per-thread output redirection", {
  conscape_file <- system.file("extdata", "conscape.jl", package = "ConScapeRtools")
  batch_file <- system.file("extdata", "conscape_batch.jl", package = "ConScapeRtools")
  distributed_file <- system.file("extdata", "conscape_batch_distributed.jl", package = "ConScapeRtools")
  conscape <- paste(readLines(conscape_file, warn = FALSE), collapse = "\n")
  batch <- paste(readLines(batch_file, warn = FALSE), collapse = "\n")
  distributed <- paste(readLines(distributed_file, warn = FALSE), collapse = "\n")

  expect_match(conscape, "function _with_conscape_logging\\(f::Function, quiet::Bool\\)")
  expect_match(batch, "sensitivity_require_landmark_one,\\s*false\\)")
  expect_match(batch, 'result in \\("done", "skipped_empty_targets"\\)', fixed = FALSE)
  expect_false(grepl("redirect_stdout\\(devnull\\)", batch))
  expect_false(grepl("redirect_stdout\\(devnull\\)", distributed))
})

test_that("Julia batch helpers expose efficiency diagnostics and scheduling controls", {
  batch_file <- system.file("extdata", "conscape_batch.jl", package = "ConScapeRtools")
  distributed_file <- system.file("extdata", "conscape_batch_distributed.jl", package = "ConScapeRtools")
  batch <- paste(readLines(batch_file, warn = FALSE), collapse = "\n")
  distributed <- paste(readLines(distributed_file, warn = FALSE), collapse = "\n")

  expect_match(batch, "BLAS\\.set_num_threads")
  expect_match(batch, "work_estimates")
  expect_match(batch, "sortperm\\(work_estimates, rev = true\\)")
  expect_match(batch, "conscape_batch_diagnostics\\.csv")

  expect_match(distributed, "BLAS\\.set_num_threads")
  expect_match(distributed, "work_estimates")
  expect_match(distributed, "sortperm\\(work_estimates, rev = true\\)")
  expect_match(distributed, "conscape_batch_diagnostics\\.csv")
})

test_that("tile_design validates source and target inputs before Julia setup", {
  r <- make_test_raster()
  expect_error(
    tile_design(r_mov = r, max_d = 10, jl_home = "not-used"),
    "At least one"
  )
})

test_that("plot.ConScapeResults validates and plots result objects", {
  r <- make_test_raster(n = 3, vals = 1)
  result <- list(btwn = r, fcon = r, elasticity_quality = r)
  class(result) <- "ConScapeResults"

  png_file <- tempfile(fileext = ".png")
  grDevices::png(png_file)
  expect_invisible(plot(result))
  expect_invisible(plot(result, layers = "elasticity_quality"))
  grDevices::dev.off()

  expect_error(plot.ConScapeResults(list(btwn = r)), "ConScapeResults")
  bad <- list(outdirs = list(btwn = "path"))
  class(bad) <- "ConScapeResults"
  expect_error(plot(bad), "at least one SpatRaster")
  expect_error(plot(result, layers = "missing"), "Requested layers")
})
