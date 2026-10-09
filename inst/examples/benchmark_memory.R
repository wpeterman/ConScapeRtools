# Run each checkout in a fresh process with the same R and terra versions:
# Rscript benchmark_memory.R /path/to/checkout prep
# Rscript benchmark_memory.R /path/to/checkout mosaic
# R allocations exclude terra/GDAL/Julia allocations; use OS peak working set
# as well when comparing total memory. Temporary benchmark data are disposable.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L, args[2] %in% c("prep", "mosaic"))
pkgload::load_all(args[1], quiet = TRUE)
root <- tempfile("conscape-benchmark-")
dir.create(root)
set.seed(2718)
terra::terraOptions(memmax = 0.25, memmin = 0)
r <- terra::rast(nrows = 1200, ncols = 1200, xmin = 0, xmax = 1200,
                 ymin = 0, ymax = 1200, crs = "EPSG:3857")
terra::values(r) <- runif(terra::ncell(r), 0.1, 0.9)
r <- terra::writeRaster(r, file.path(root, "input.tif"), datatype = "FLT8S")
if (args[2] == "prep") {
  operation <- function() conscape_prep(
    centersize = 300, buffer = 75, landmark = 1L, r_target = r, r_mov = r,
    r_src = r, asc_dir = file.path(root, "tiles"), clear_dir = TRUE, progress = FALSE
  )
} else {
  tiles <- file.path(root, "tiles")
  dir.create(tiles)
  for (i in 1:16) {
    row <- (i - 1L) %% 4L
    col <- (i - 1L) %/% 4L
    tile <- terra::crop(r, terra::ext(max(0, col * 300 - 75), min(1200, (col + 1) * 300 + 75),
                                    max(0, row * 300 - 75), min(1200, (row + 1) * 300 + 75)))
    terra::writeRaster(tile, file.path(tiles, sprintf("tile-%02d.tif", i)), datatype = "FLT8S")
  }
  operation <- function() mosaic_conscape(tiles, tile_trim = 75, method = "sum", chunk_size = 4L)
}
# Warm dispatch and disk cache. Three measured repeats exclude input generation.
invisible(operation())
measurements <- lapply(1:3, function(i) {
  invisible(gc())
  profile <- file.path(root, "allocation.log")
  utils::Rprofmem(profile)
  elapsed <- system.time(out <- operation())[["elapsed"]]
  utils::Rprofmem(NULL)
  lines <- readLines(profile, warn = FALSE)
  sizes <- suppressWarnings(as.numeric(sub(" .*", "", lines)))
  data.frame(case = args[2], replicate = i, seconds = elapsed,
             r_allocated_mb = sum(sizes, na.rm = TRUE) / 1e6,
             largest_r_allocation_mb = max(sizes, na.rm = TRUE) / 1e6)
})
print(do.call(rbind, measurements), row.names = FALSE)
cat(R.version.string, "terra", as.character(utils::packageVersion("terra")), "\n")
unlink(root, recursive = TRUE)
