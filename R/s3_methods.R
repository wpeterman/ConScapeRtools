conscape_raster_layers <- function(x) {
  is_raster <- vapply(x, inherits, logical(1), what = "SpatRaster")
  names(x)[is_raster]
}

conscape_format_value <- function(x) {
  if (is.null(x)) {
    return("NULL")
  }
  if (is.numeric(x) && length(x) == 1L) {
    return(format(x, digits = 6, trim = TRUE))
  }
  if (is.logical(x) && length(x) == 1L) {
    return(as.character(x))
  }
  if (is.character(x)) {
    return(paste(x, collapse = ", "))
  }
  paste(utils::capture.output(utils::str(x, give.attr = FALSE)), collapse = " ")
}

#' S3 display methods for ConScapeRtools objects
#'
#' @description
#' Print and summarize the main object types returned by ConScapeRtools:
#' design objects from [tile_design()], prep objects from [conscape_prep()],
#' window design objects from [window_design()], result objects from
#' [run_conscape()], and sensitivity specifications from [conscape_sensitivity()].
#'
#' @param x Object to print or plot.
#' @param object Object to summarize.
#' @param ... Additional arguments passed to methods.
#'
#' @return
#' Print and plot methods return their input invisibly. Summary methods return
#' a compact data frame or list describing the object.
#'
#' @name conscape_s3_methods
NULL

#' @rdname conscape_s3_methods
#' @method print ConScapeRtools_design
#' @export
print.ConScapeRtools_design <- function(x, ...) {
  cat("ConScape tile design\n")
  cat("  distance_scale: ", conscape_format_value(x$distance_scale), "\n", sep = "")
  cat("  theta:          ", conscape_format_value(x$theta), "\n", sep = "")
  cat("  tile_d:         ", conscape_format_value(x$tile_d), "\n", sep = "")
  cat("  tile_trim:      ", conscape_format_value(x$tile_trim), "\n", sep = "")
  if (!is.null(x$centersize) && !is.null(x$buffer)) {
    cat("  centersize:     ", conscape_format_value(x$centersize), "\n", sep = "")
    cat("  buffer:         ", conscape_format_value(x$buffer), "\n", sep = "")
  }
  if (!is.null(x$landmark)) {
    cat("  landmark:       ", conscape_format_value(x$landmark), "\n", sep = "")
  }
  if (is.list(x$diagnostics) && !is.null(x$diagnostics$expected_tile_count)) {
    cat("  expected tiles: ", conscape_format_value(x$diagnostics$expected_tile_count), "\n", sep = "")
    cat("  area factor:    ", conscape_format_value(round(x$overlap_area_factor, 3)), "\n", sep = "")
  }
  invisible(x)
}

#' @rdname conscape_s3_methods
#' @method summary ConScapeRtools_design
#' @export
summary.ConScapeRtools_design <- function(object, ...) {
  diagnostics <- object$diagnostics
  diagnostic_value <- function(name) {
    if (is.list(diagnostics) && !is.null(diagnostics[[name]])) {
      return(diagnostics[[name]])
    }
    NA_real_
  }
  out <- data.frame(
    parameter = c(
      "distance_scale",
      "theta",
      "tile_d",
      "tile_trim",
      "centersize",
      "buffer",
      "landmark",
      "trim_threshold",
      "expected_tile_count",
      "overlap_area_factor",
      "proximity_at_trim"
    ),
    value = c(
      object$distance_scale,
      object$theta,
      object$tile_d,
      object$tile_trim,
      if (is.null(object$centersize)) NA_real_ else object$centersize,
      if (is.null(object$buffer)) NA_real_ else object$buffer,
      if (is.null(object$landmark)) NA_real_ else object$landmark,
      if (is.null(object$trim_threshold)) NA_real_ else object$trim_threshold,
      diagnostic_value("expected_tile_count"),
      if (is.null(object$overlap_area_factor)) NA_real_ else object$overlap_area_factor,
      diagnostic_value("proximity_at_effective_trim")
    ),
    description = c(
      "Exponential distance-decay numerator for run_conscape()",
      "Randomized shortest-path theta used during calibration",
      "Suggested minimum interior tile width",
      "Suggested minimum tile overlap and mosaic trim width",
      "Equivalent center window size for WindowedProblem",
      "Equivalent buffer width for WindowedProblem",
      "Landmark value used to round tile width and trim",
      "Requested maximum proximity at the trim distance",
      "Expected number of interior tiles for r_mov",
      "Extended tile area divided by retained interior tile area",
      "Expected proximity at the effective trim distance"
    ),
    stringsAsFactors = FALSE
  )
  class(out) <- c("summary_ConScapeRtools_design", class(out))
  out
}

#' @rdname conscape_s3_methods
#' @method print ConScapeRtools_window_design
#' @export
print.ConScapeRtools_window_design <- function(x, ...) {
  cat("ConScape window design\n")
  cat("  centersize:       ", conscape_format_value(x$centersize), " cells\n", sep = "")
  cat("  buffer:           ", conscape_format_value(x$buffer), " cells\n", sep = "")
  cat("  tile_d:           ", conscape_format_value(x$tile_d), "\n", sep = "")
  cat("  tile_trim:        ", conscape_format_value(x$tile_trim), "\n", sep = "")
  cat("  source window:    ",
      conscape_format_value(x$diagnostics$total_window_width_cells),
      " x ",
      conscape_format_value(x$diagnostics$total_window_width_cells),
      " cells\n", sep = "")
  cat("  expected windows: ", conscape_format_value(x$diagnostics$expected_windows), "\n", sep = "")
  cat("  area factor:      ",
      conscape_format_value(round(x$diagnostics$overlap_area_factor, 3)),
      "\n", sep = "")
  invisible(x)
}

#' @rdname conscape_s3_methods
#' @method summary ConScapeRtools_window_design
#' @export
summary.ConScapeRtools_window_design <- function(object, ...) {
  diagnostics <- object$diagnostics
  out <- data.frame(
    parameter = c(
      "centersize",
      "buffer",
      "tile_d",
      "tile_trim",
      "source_window_cells",
      "center_cells",
      "expected_windows",
      "overlap_area_factor",
      "workload_proxy"
    ),
    value = c(
      object$centersize,
      object$buffer,
      object$tile_d,
      object$tile_trim,
      diagnostics$source_window_cells,
      diagnostics$center_cells,
      diagnostics$expected_windows,
      diagnostics$overlap_area_factor,
      diagnostics$workload_proxy
    ),
    description = c(
      "Center window size in cells",
      "Context buffer width in cells",
      "Equivalent center width in map units",
      "Equivalent buffer width in map units",
      "Cells in the full source/context window",
      "Cells retained as center targets",
      "Approximate number of windows over r",
      "Context window area divided by center area",
      "Approximate source cells times target cells times windows"
    ),
    stringsAsFactors = FALSE
  )
  class(out) <- c("summary_ConScapeRtools_window_design", class(out))
  out
}

#' @rdname conscape_s3_methods
#' @method print summary_ConScapeRtools_window_design
#' @export
print.summary_ConScapeRtools_window_design <- function(x, ...) {
  cat("ConScape window design summary\n")
  print.data.frame(x, row.names = FALSE, ...)
  invisible(x)
}

#' @rdname conscape_s3_methods
#' @method print summary_ConScapeRtools_design
#' @export
print.summary_ConScapeRtools_design <- function(x, ...) {
  cat("ConScape tile design summary\n")
  print.data.frame(x, row.names = FALSE, ...)
  invisible(x)
}

#' @rdname conscape_s3_methods
#' @method plot ConScapeRtools_design
#' @export
plot.ConScapeRtools_design <- function(x, ...) {
  vals <- c(
    distance_scale = x$distance_scale,
    tile_d = x$tile_d,
    tile_trim = x$tile_trim
  )
  graphics::barplot(
    vals,
    ylab = "Map units",
    main = "ConScape Tile Design",
    ...
  )
  invisible(x)
}

#' @rdname conscape_s3_methods
#' @method print ConScapeRtools_prep
#' @export
print.ConScapeRtools_prep <- function(x, ...) {
  cat("ConScape prepared tile set\n")
  cat("  tiles:            ", length(x$tile_num), "\n", sep = "")
  cat("  landmark:         ", conscape_format_value(x$landmark), "\n", sep = "")
  cat("  tile_trim:        ", conscape_format_value(x$tile_trim), "\n", sep = "")
  cat("  target_threshold: ", conscape_format_value(x$target_threshold), "\n", sep = "")
  cat("  asc_dir:          ", conscape_format_value(x$asc_dir), "\n", sep = "")
  invisible(x)
}

#' @rdname conscape_s3_methods
#' @method summary ConScapeRtools_prep
#' @export
summary.ConScapeRtools_prep <- function(object, ...) {
  n_tiles <- length(object$tile_num)
  n_polygons <- if (inherits(object$cs_tiles, "SpatVector")) nrow(object$cs_tiles) else NA_integer_
  validation <- object$tile_validation
  n_bad <- if (is.list(validation) && !is.null(validation$n_bad)) validation$n_bad else NA_integer_

  out <- data.frame(
    field = c(
      "tiles",
      "tile_polygons",
      "bad_tiles_removed",
      "landmark",
      "tile_trim",
      "target_threshold",
      "asc_dir",
      "target_dir",
      "source_dir",
      "movement_dir"
    ),
    value = c(
      n_tiles,
      n_polygons,
      n_bad,
      object$landmark,
      object$tile_trim,
      object$target_threshold,
      object$asc_dir,
      object$target,
      object$src,
      object$mov
    ),
    stringsAsFactors = FALSE
  )
  class(out) <- c("summary_ConScapeRtools_prep", class(out))
  out
}

#' @rdname conscape_s3_methods
#' @method print summary_ConScapeRtools_prep
#' @export
print.summary_ConScapeRtools_prep <- function(x, ...) {
  cat("ConScape prepared tile set summary\n")
  print.data.frame(x, row.names = FALSE, ...)
  invisible(x)
}

#' @rdname conscape_s3_methods
#' @method plot ConScapeRtools_prep
#' @export
plot.ConScapeRtools_prep <- function(x, ...) {
  if (!inherits(x$cs_tiles, "SpatVector")) {
    stop("x$cs_tiles must be a SpatVector to plot a ConScapeRtools_prep object.", call. = FALSE)
  }
  terra::plot(x$cs_tiles, main = "ConScape Tile Layout", ...)
  if ("tile_id" %in% names(x$cs_tiles)) {
    centers <- terra::centroids(x$cs_tiles)
    coords <- terra::crds(centers)
    graphics::text(coords[, 1], coords[, 2], labels = x$cs_tiles$tile_id, cex = 0.7)
  }
  invisible(x)
}

#' @rdname conscape_s3_methods
#' @method print ConScapeResults
#' @export
print.ConScapeResults <- function(x, ...) {
  raster_layers <- conscape_raster_layers(x)
  outdirs <- if (is.list(x$outdirs)) names(x$outdirs) else character()

  cat("ConScape results\n")
  cat("  raster layers: ", if (length(raster_layers)) paste(raster_layers, collapse = ", ") else "none", "\n", sep = "")
  cat("  output dirs:   ", if (length(outdirs)) paste(outdirs, collapse = ", ") else "none", "\n", sep = "")
  invisible(x)
}

#' @rdname conscape_s3_methods
#' @method summary ConScapeResults
#' @export
summary.ConScapeResults <- function(object, ...) {
  raster_layers <- conscape_raster_layers(object)
  if (length(raster_layers) == 0L) {
    out <- data.frame(layer = character(), nrow = integer(), ncol = integer(),
                      nlyr = integer(), min = numeric(), max = numeric())
  } else {
    out <- do.call(rbind, lapply(raster_layers, function(layer) {
      r <- object[[layer]]
      rng <- tryCatch(
        terra::global(r, c("min", "max"), na.rm = TRUE),
        error = function(e) data.frame(min = NA_real_, max = NA_real_)
      )
      data.frame(
        layer = layer,
        nrow = terra::nrow(r),
        ncol = terra::ncol(r),
        nlyr = terra::nlyr(r),
        min = rng[1, "min"],
        max = rng[1, "max"],
        stringsAsFactors = FALSE
      )
    }))
  }
  class(out) <- c("summary_ConScapeResults", class(out))
  out
}

#' @rdname conscape_s3_methods
#' @method print summary_ConScapeResults
#' @export
print.summary_ConScapeResults <- function(x, ...) {
  cat("ConScape results summary\n")
  print.data.frame(x, row.names = FALSE, ...)
  invisible(x)
}

#' @rdname conscape_s3_methods
#' @method print ConScapeSensitivitySpec
#' @export
print.ConScapeSensitivitySpec <- function(x, ...) {
  cat("ConScape sensitivity specification\n")
  cat("  wrt:                  ", conscape_format_value(x$wrt), "\n", sep = "")
  cat("  method:               ", conscape_format_value(x$method), "\n", sep = "")
  cat("  landscape_measure:    ", conscape_format_value(x$landscape_measure), "\n", sep = "")
  cat("  unitless:             ", conscape_format_value(x$unitless), "\n", sep = "")
  cat("  require_landmark_one: ", conscape_format_value(x$require_landmark_one), "\n", sep = "")
  invisible(x)
}

#' @rdname conscape_s3_methods
#' @method summary ConScapeSensitivitySpec
#' @export
summary.ConScapeSensitivitySpec <- function(object, ...) {
  out <- data.frame(
    field = names(object),
    value = vapply(object, conscape_format_value, character(1)),
    stringsAsFactors = FALSE
  )
  class(out) <- c("summary_ConScapeSensitivitySpec", class(out))
  out
}

#' @rdname conscape_s3_methods
#' @method print summary_ConScapeSensitivitySpec
#' @export
print.summary_ConScapeSensitivitySpec <- function(x, ...) {
  cat("ConScape sensitivity specification summary\n")
  print.data.frame(x, row.names = FALSE, ...)
  invisible(x)
}
