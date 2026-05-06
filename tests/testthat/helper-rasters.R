make_test_raster <- function(n = 20, vals = 1, crs = "") {
  terra::rast(
    nrows = n, ncols = n,
    xmin = 0, xmax = n,
    ymin = 0, ymax = n,
    vals = vals,
    crs = crs
  )
}

write_test_asc <- function(r, path) {
  terra::writeRaster(r, path, overwrite = TRUE, NAflag = -9999)
}
