#' Validate requested ConScape output directories
#' @noRd
validate_conscape_outputs <- function(out_dir,
                                      output_specs,
                                      expected_tiles,
                                      pattern = "\\.asc$|\\.tif$") {
  rows <- lapply(output_specs, function(spec) {
    output_dir <- file.path(out_dir, spec$dir)
    files <- if (dir.exists(output_dir)) {
      list.files(output_dir, pattern = pattern, full.names = FALSE)
    } else {
      character()
    }
    data.frame(
      layer = spec$layer,
      dir = spec$dir,
      expected = expected_tiles,
      found = length(files),
      missing = max(expected_tiles - length(files), 0L),
      ok = length(files) == expected_tiles,
      path = normalizePath(output_dir, mustWork = FALSE),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

#' Read Julia batch diagnostics when present
#' @noRd
read_conscape_batch_diagnostics <- function(out_dir) {
  diag_file <- file.path(out_dir, "conscape_batch_diagnostics.csv")
  if (!file.exists(diag_file)) {
    return(NULL)
  }
  tryCatch(
    utils::read.csv(diag_file, stringsAsFactors = FALSE),
    error = function(e) NULL
  )
}
