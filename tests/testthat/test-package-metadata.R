test_that("installed example data paths resolve", {
  expect_true(file.exists(system.file("extdata", "affinity.asc", package = "ConScapeRtools")))
  expect_true(file.exists(system.file("extdata", "suitability.asc", package = "ConScapeRtools")))
  expect_true(file.exists(system.file("extdata", "patches.shp", package = "ConScapeRtools")))
})

test_that("Julia helper scripts referenced by run_conscape are installed", {
  expect_true(file.exists(system.file("extdata", "conscape.jl", package = "ConScapeRtools")))
  expect_true(file.exists(system.file("extdata", "conscape_batch.jl", package = "ConScapeRtools")))
  expect_true(file.exists(system.file("extdata", "conscape_batch_distributed.jl", package = "ConScapeRtools")))
})

test_that("threaded Julia batch avoids per-thread output redirection", {
  conscape_file <- system.file("extdata", "conscape.jl", package = "ConScapeRtools")
  batch_file <- system.file("extdata", "conscape_batch.jl", package = "ConScapeRtools")
  conscape <- paste(readLines(conscape_file, warn = FALSE), collapse = "\n")
  batch <- paste(readLines(batch_file, warn = FALSE), collapse = "\n")

  expect_match(conscape, "function _with_conscape_logging\\(f::Function, quiet::Bool\\)")
  expect_match(batch, "sensitivity_require_landmark_one,\\s*false\\)")
  expect_match(batch, 'result in \\("done", "skipped_empty_targets"\\)', fixed = FALSE)
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
