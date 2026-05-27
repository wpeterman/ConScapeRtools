integration_jl_home <- function() {
  candidates <- c(
    Sys.getenv("CONSCAPERTOOLS_JL_HOME", ""),
    Sys.getenv("JULIA_BINDIR", "")
  )
  candidates <- candidates[nzchar(candidates)]
  if (length(candidates)) candidates[[1]] else ""
}

crop_cells_from_top_left <- function(r, n = 40L) {
  x0 <- terra::xmin(r)
  y1 <- terra::ymax(r)
  terra::crop(
    r,
    terra::ext(
      x0,
      x0 + n * terra::res(r)[1],
      y1 - n * terra::res(r)[2],
      y1
    ),
    snap = "near"
  )
}

surface_difference_summary <- function(tiled, untiled) {
  diff <- terra::values(tiled - untiled, mat = FALSE)
  diff <- diff[is.finite(diff)]
  c(
    max_abs = max(abs(diff)),
    mean_abs = mean(abs(diff)),
    rmse = sqrt(mean(diff^2))
  )
}

test_that("tiled example-data surfaces match untiled surfaces after mosaicing", {
  skip_on_cran()
  skip_if_not(
    identical(Sys.getenv("RUN_CONSCAPERTOOLS_INTEGRATION", ""), "true"),
    "Set RUN_CONSCAPERTOOLS_INTEGRATION=true to run Julia-backed ConScape checks."
  )

  jl_home <- integration_jl_home()
  skip_if_not(
    nzchar(jl_home),
    "Set CONSCAPERTOOLS_JL_HOME or JULIA_BINDIR to the Julia bin directory."
  )

  habitat_file <- system.file("extdata", "suitability.asc", package = "ConScapeRtools")
  affinity_file <- system.file("extdata", "affinity.asc", package = "ConScapeRtools")
  habitat <- crop_cells_from_top_left(terra::rast(habitat_file), n = 40L)
  affinity <- crop_cells_from_top_left(terra::rast(affinity_file), n = 40L)

  expect_equal(dim(habitat)[1:2], c(40L, 40L))
  expect_equal(dim(affinity)[1:2], c(40L, 40L))

  cell_size <- terra::res(habitat)[1]
  tile_width <- 20 * cell_size
  tile_trim <- 20 * cell_size

  prep <- conscape_prep(
    tile_d = tile_width,
    tile_trim = tile_trim,
    asc_dir = file.path(tempdir(), "agreement-prep"),
    r_target = habitat,
    r_mov = affinity,
    r_src = habitat,
    clear_dir = TRUE,
    landmark = 5L,
    progress = FALSE
  )
  expect_equal(length(prep$tile_num), 4L)

  tiled <- run_conscape(
    conscape_prep = prep,
    out_dir = file.path(tempdir(), "agreement-tiled"),
    theta = 0.15,
    distance_scale = 150,
    jl_home = jl_home,
    progress = FALSE
  )

  untiled <- run_conscape(
    target_qualities = habitat,
    source_qualities = habitat,
    affinities = affinity,
    out_dir = file.path(tempdir(), "agreement-untiled"),
    theta = 0.15,
    distance_scale = 150,
    landmark = 5L,
    tile_trim = prep$tile_trim,
    jl_home = jl_home,
    progress = FALSE
  )

  summaries <- rbind(
    btwn = surface_difference_summary(tiled$btwn, untiled[["btwn"]]),
    fcon = surface_difference_summary(tiled$fcon, untiled[["fcon"]])
  )

  expect_equal(unname(summaries[, "max_abs"]), c(0, 0), tolerance = 1e-8)
  expect_equal(unname(summaries[, "rmse"]), c(0, 0), tolerance = 1e-10)
})
