# Worked example for the experimental ConScape dev backend.
#
# This example uses the package data shipped with ConScapeRtools. It is not run
# automatically because it installs and calls a development branch of ConScape.jl.

library(ConScapeRtools)
library(terra)

# Update this path for your Julia installation.
jl_home <- "C:/Users/peterman.73/AppData/Local/Programs/Julia-1.12.6/bin/"

habitat_file <- system.file("extdata", "suitability.asc", package = "ConScapeRtools")
affinity_file <- system.file("extdata", "affinity.asc", package = "ConScapeRtools")

habitat <- rast(habitat_file)
affinity <- rast(affinity_file)

# Keep the demonstration small while exercising real package data.
demo_extent <- ext(
  xmin(habitat),
  xmin(habitat) + 40 * res(habitat)[1],
  ymax(habitat) - 40 * res(habitat)[2],
  ymax(habitat)
)
habitat_demo <- crop(habitat, demo_extent, snap = "near")
affinity_demo <- crop(affinity, demo_extent, snap = "near")

# Translate map-unit design values into ConScape dev window values.
# In a real analysis, use tile_design() first to derive these distances from
# max_d, theta, and trim_threshold.
wd <- window_design(
  r = habitat_demo,
  tile_d = 2700,
  tile_trim = 2700,
  landmark = 1L
)

# Create a dedicated Julia project with the ConScape development backend.
# Reuse the returned path across sessions to avoid reinstalling each time.
dev_project <- conscape_dev_backend_setup(
  jl_home = jl_home,
  rev = "alg_efficiency",
  quiet = FALSE
)

# WindowedProblem mode is the first target for correctness and API stability.
windowed <- run_conscape(
  target_qualities = habitat_demo,
  source_qualities = habitat_demo,
  affinities = affinity_demo,
  out_dir = file.path(tempdir(), "conscape_dev_windowed"),
  jl_home = jl_home,
  backend = "conscape_dev",
  dev_mode = "windowed",
  centersize = wd$centersize,
  buffer = wd$buffer,
  dev_project = dev_project,
  landmark = 1L,
  theta = 0.15,
  distance_scale = 150,
  metrics = c("btwn", "fcon"),
  progress = TRUE
)

plot(windowed)

# BatchProblem mode writes intermediate batches and mosaics them after all
# batches complete. It is the throughput-oriented path for larger extents.
batched <- run_conscape(
  target_qualities = habitat_demo,
  source_qualities = habitat_demo,
  affinities = affinity_demo,
  out_dir = file.path(tempdir(), "conscape_dev_batch"),
  jl_home = jl_home,
  backend = "conscape_dev",
  dev_mode = "batch",
  centersize = wd$centersize,
  buffer = wd$buffer,
  dev_project = dev_project,
  landmark = 1L,
  theta = 0.15,
  distance_scale = 150,
  metrics = c("btwn", "fcon"),
  progress = TRUE
)

plot(batched)

# WindowedProblem and BatchProblem should agree for the same parameters.
terra::global(abs(windowed - batched), c("max", "mean"), na.rm = TRUE)
