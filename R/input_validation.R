# Shared preflight checks run before Julia starts or output files are removed.
validate_conscape_rasters <- function(rasters) {
  if (!all(vapply(rasters, inherits, logical(1), "SpatRaster")) ||
      any(vapply(rasters, terra::nlyr, numeric(1)) != 1L)) {
    stop("Inputs must all be single-layer SpatRaster objects.", call. = FALSE)
  }
  first <- rasters[[1L]]
  for (r in rasters[-1L]) {
    if (!isTRUE(all.equal(terra::ext(first), terra::ext(r)))) {
      stop("Raster inputs must share extent.", call. = FALSE)
    }
    if (!isTRUE(all.equal(terra::res(first), terra::res(r)))) {
      stop("Raster inputs must share resolution.", call. = FALSE)
    }
    if (!terra::same.crs(first, r)) {
      stop("Raster inputs must share CRS.", call. = FALSE)
    }
  }
  if (!isTRUE(all.equal(terra::res(first)[1], terra::res(first)[2]))) {
    stop("ConScape ASCII inputs require square cells; resample explicitly first.", call. = FALSE)
  }
  invisible(NULL)
}

conscape_path <- function(path) {
  # Resolve existing ancestors too, so .. and symlink aliases cannot evade checks.
  if (file.exists(path) || dir.exists(path)) {
    out <- normalizePath(path, winslash = "/", mustWork = TRUE)
  } else {
    parent <- dirname(path)
    out <- if (identical(parent, path)) path else file.path(conscape_path(parent), basename(path))
    out <- normalizePath(out, winslash = "/", mustWork = FALSE)
  }
  out <- sub("/+$", "", out)
  if (.Platform$OS.type == "windows") out <- tolower(out)
  out
}

validate_conscape_output_dir <- function(path, clear_dir, protected = character()) {
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("Output directory must be a single non-empty path.", call. = FALSE)
  }
  resolved <- conscape_path(path)
  broad <- vapply(c(getwd(), path.expand("~"), tempdir(), dirname(resolved)),
                  conscape_path, character(1))
  if (resolved %in% broad ||
      any(startsWith(broad[1:3], paste0(resolved, "/")))) {
    stop("Use a dedicated working directory, not a workspace, home, temporary root, or filesystem root.", call. = FALSE)
  }
  protected <- protected[!is.na(protected) & nzchar(protected)]
  for (p in protected) {
    p <- conscape_path(p)
    if (identical(resolved, p) || startsWith(p, paste0(resolved, "/")) ||
        (dir.exists(p) && startsWith(resolved, paste0(p, "/")))) {
      stop("Output directory would overwrite input files. Use separate input and output directories.", call. = FALSE)
    }
  }
  if (file.exists(path) && !dir.exists(path)) {
    stop("Output directory points to an existing file.", call. = FALSE)
  }
  if (dir.exists(path) && !clear_dir &&
      length(list.files(path, all.files = TRUE, no.. = TRUE)) > 0L) {
    stop("Files currently exist in the output directory. Set clear_dir = TRUE to remove them.", call. = FALSE)
  }
  invisible(NULL)
}

conscape_input_sources <- function(rasters) {
  unlist(lapply(rasters, function(x) {
    if (inherits(x, "SpatRaster")) terra::sources(x) else character()
  }), use.names = FALSE)
}

conscape_tile_manifest <- function(target_dir, source_dir, move_dir) {
  dirs <- c(target_dir, source_dir, move_dir)
  if (length(dirs) != 3L || anyNA(dirs) || any(!dir.exists(dirs))) {
    stop("Supply three existing tile directories or three SpatRaster objects.", call. = FALSE)
  }
  files <- lapply(dirs, list.files, pattern = "\\.asc$")
  if (!length(files[[1L]]) || !identical(files[[1L]], files[[2L]]) ||
      !identical(files[[1L]], files[[3L]])) {
    stop("Tile directories must contain non-empty, identical sets of .asc filenames.", call. = FALSE)
  }
  for (f in files[[1L]]) {
    validate_conscape_rasters(lapply(file.path(dirs, f), terra::rast))
  }
  files[[1L]]
}
