test_that("tile design print, summary, and plot methods are informative", {
  design <- structure(
    list(
      distance_scale = 12.5,
      exp_d = 12.5,
      tile_d = 100,
      tile_trim = 25,
      theta = 0.2,
      landmark = 10L,
      trim_threshold = 0.01,
      overlap_area_factor = 2.25,
      diagnostics = list(
        expected_tile_count = 4L,
        proximity_at_effective_trim = 0.01
      )
    ),
    class = "ConScapeRtools_design"
  )

  printed <- capture.output(print(design))
  expect_match(printed[1], "ConScape tile design")
  expect_true(any(grepl("distance_scale", printed)))
  expect_true(any(grepl("expected tiles", printed)))

  s <- summary(design)
  expect_s3_class(s, "summary_ConScapeRtools_design")
  expect_true(all(c(
    "distance_scale", "theta", "tile_d", "tile_trim",
    "expected_tile_count", "overlap_area_factor", "proximity_at_trim"
  ) %in% s$parameter))

  png_file <- tempfile(fileext = ".png")
  grDevices::png(png_file)
  expect_invisible(plot(design))
  grDevices::dev.off()
})

test_that("prep print, summary, and plot methods describe tile bookkeeping", {
  r <- make_test_raster(n = 10, vals = 1)
  prep <- conscape_prep(
    tile_d = 5,
    tile_trim = 2,
    asc_dir = file.path(tempdir(), "s3-prep"),
    r_target = r,
    r_mov = r,
    r_src = r,
    clear_dir = TRUE,
    landmark = 5L,
    progress = FALSE
  )

  printed <- capture.output(print(prep))
  expect_match(printed[1], "prepared tile set")
  expect_true(any(grepl("tiles", printed)))

  s <- summary(prep)
  expect_s3_class(s, "summary_ConScapeRtools_prep")
  expect_true("bad_tiles_removed" %in% s$field)
  expect_equal(as.integer(s$value[s$field == "landmark"]), 5L)

  png_file <- tempfile(fileext = ".png")
  grDevices::png(png_file)
  expect_invisible(plot(prep))
  grDevices::dev.off()
})

test_that("ConScapeResults print and summary methods report raster layers", {
  r <- make_test_raster(n = 3, vals = 1)
  result <- list(
    btwn = r,
    fcon = r * 2,
    outdirs = list(btwn = "btwn-dir", fcon = "fcon-dir")
  )
  class(result) <- "ConScapeResults"

  printed <- capture.output(print(result))
  expect_match(printed[1], "ConScape results")
  expect_true(any(grepl("btwn, fcon", printed)))

  s <- summary(result)
  expect_s3_class(s, "summary_ConScapeResults")
  expect_equal(s$layer, c("btwn", "fcon"))
  expect_equal(s$nrow, c(3L, 3L))
  expect_equal(s$min, c(1, 2))
})

test_that("sensitivity specifications have compact print and summary methods", {
  spec <- conscape_sensitivity(wrt = c("Q", "A&C=f(A)"), method = "simulation")

  printed <- capture.output(print(spec))
  expect_match(printed[1], "sensitivity specification")
  expect_true(any(grepl("A&C=f\\(A\\)", printed)))

  s <- summary(spec)
  expect_s3_class(s, "summary_ConScapeSensitivitySpec")
  expect_true(all(c("wrt", "method", "unitless") %in% s$field))
})
