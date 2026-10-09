mock_julia <- function(.local_envir = parent.frame()) {
  testthat::local_mocked_bindings(
    stopJulia = function() invisible(TRUE),
    juliaSetupOk = function() TRUE,
    juliaEval = function(...) invisible(TRUE),
    ConScapeR_setup = function(...) invisible(TRUE),
    .package = "ConScapeRtools", .env = .local_envir
  )
}

write_fake_conscape_outputs <- function(out_dir, source_file, iter = "") {
  r <- terra::rast(source_file)
  btwn <- r
  fcon <- r
  terra::values(btwn) <- seq_len(terra::ncell(btwn))
  terra::values(fcon) <- seq_len(terra::ncell(fcon)) * 2
  dir.create(file.path(out_dir, "btwn"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out_dir, "fcon"), recursive = TRUE, showWarnings = FALSE)
  terra::writeRaster(btwn, file.path(out_dir, "btwn", paste0("betweenness", iter, ".asc")),
                     overwrite = TRUE, NAflag = -9999)
  terra::writeRaster(fcon, file.path(out_dir, "fcon", paste0("fcon", iter, ".asc")),
                     overwrite = TRUE, NAflag = -9999)
  "done"
}
