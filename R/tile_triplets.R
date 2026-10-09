#' Write matched ConScape tile triplets
#' @noRd
tile_rast_triplets <- function(r_target,
                               r_mov,
                               r_src,
                               tiles,
                               root_dir,
                               clear_dir = FALSE,
                               progress = FALSE,
                               target_mode = c("full", "center")) {
  target_mode <- match.arg(target_mode)

  dirs <- list(
    target = file.path(root_dir, "target"),
    mov = file.path(root_dir, "mov"),
    src = file.path(root_dir, "src")
  )

  for (write_dir in dirs) {
    if (!dir.exists(write_dir)) dir.create(write_dir, recursive = TRUE)
    if (length(list.files(write_dir)) > 0L && !clear_dir) {
      stop("Files currently exist in the specified directory. Set `clear_dir = TRUE` to remove these files.")
    }
    unlink(file.path(write_dir, "*"), force = TRUE)
  }

  cs_tiles <- tiles$cs_tiles
  overlap_cells <- tiles$overlap_cells
  resx <- terra::res(r_target)[1]
  resy <- terra::res(r_target)[2]
  n_tiles <- length(cs_tiles)

  tile_extents <- lapply(seq_len(n_tiles), function(i) {
    te <- cs_tiles[i]
    terra::ext(
      terra::xmin(te) - overlap_cells * resx,
      terra::xmax(te) + overlap_cells * resx,
      terra::ymin(te) - overlap_cells * resy,
      terra::ymax(te) + overlap_cells * resy
    )
  })

  tile_counts <- function(r) {
    finite <- is.finite(r)
    c(finite = terra::global(finite, "sum", na.rm = TRUE)[1, 1],
      positive = terra::global(finite & r > 0, "sum", na.rm = TRUE)[1, 1])
  }

  if (isTRUE(progress)) {
    pb <- txtProgressBar(min = 0, max = n_tiles, style = 3)
    on.exit(close(pb), add = TRUE)
  }

  diagnostics <- vector("list", n_tiles)
  ok <- rep(FALSE, n_tiles)

  for (idx in seq_len(n_tiles)) {
    if (exists("pb", inherits = FALSE)) {
      setTxtProgressBar(pb, idx,
                        title = "Writing matched ConScape tile triplets...")
    }

    # Crop from the original backing file; pad only the current tile.
    target_tile <- terra::crop(r_target, tile_extents[[idx]], extend = TRUE)
    mov_tile <- terra::crop(r_mov, tile_extents[[idx]], extend = TRUE)
    src_tile <- terra::crop(r_src, tile_extents[[idx]], extend = TRUE)

    # Target qualities are part of the RSP problem, not merely an output mask:
    # they define the absorbing-state set in ConScape's GridRSP, which controls
    # the columns of Z = (I - W) \ I_target.
    #
    # `target_mode = "full"` keeps the full buffered target raster, so each tile
    # solves a Z whose columns include every target cell in the buffered window.
    # Tiled outputs are then averaged across overlapping tiles by
    # `mosaic_conscape(method = "mosaic")`. This is a tile-local approximation:
    # omitted targets and changed graph boundaries both affect the result,
    # so the direction of its error is not guaranteed.
    #
    # `target_mode = "center"` restricts the absorbing target set to the
    # center cells of the tile, while keeping the full buffered source qualities
    # and affinities. Adjacent tiles' centers are disjoint, so every target
    # appears in exactly one tile. Tiled outputs are then combined with
    # `mosaic_conscape(method = "sum")` and the result equals the untiled
    # surface up to truncation error from paths exiting the buffer. This is the
    # same per-window math that ConScape's dev WindowedProblem uses, executed
    # against the stable bundled Julia API.
    output_target_tile <- if (identical(target_mode, "center")) {
      terra::extend(terra::crop(target_tile, terra::ext(cs_tiles[idx])),
                    terra::ext(target_tile))
    } else {
      target_tile
    }

    tile_name <- paste0("r_", idx, ".asc")
    target_file <- file.path(dirs$target, tile_name)
    mov_file <- file.path(dirs$mov, tile_name)
    src_file <- file.path(dirs$src, tile_name)

    target_counts <- tile_counts(target_tile)
    source_counts <- tile_counts(src_tile)
    move_counts <- tile_counts(mov_tile)
    output_counts <- if (target_mode == "center") tile_counts(output_target_tile) else target_counts
    n_finite_target <- unname(target_counts["finite"])
    n_finite_source <- unname(source_counts["finite"])
    n_finite_move <- unname(move_counts["finite"])
    n_positive_target <- unname(target_counts["positive"])
    n_positive_source <- unname(source_counts["positive"])
    n_positive_move <- unname(move_counts["positive"])
    n_finite_output_target <- unname(output_counts["finite"])
    n_positive_output_target <- unname(output_counts["positive"])

    # In center mode we require at least one center target to keep the tile;
    # otherwise the tile would solve over an empty absorbing set.
    ok[idx] <- n_finite_source > 0L && n_finite_move > 0L &&
      (if (identical(target_mode, "center")) n_positive_output_target > 0L
       else n_finite_target > 0L)

    if (ok[idx]) {
      # Write the target raster that ConScape will actually ingest: center-only
      # in center mode, full buffered in full mode.
      tile_target_to_write <- if (identical(target_mode, "center")) {
        output_target_tile
      } else {
        target_tile
      }
      terra::writeRaster(tile_target_to_write, filename = target_file,
                         NAflag = -9999, overwrite = TRUE)
      terra::writeRaster(mov_tile, filename = mov_file,
                         NAflag = -9999, overwrite = TRUE)
      terra::writeRaster(src_tile, filename = src_file,
                         NAflag = -9999, overwrite = TRUE)
    }

    diagnostics[[idx]] <- data.frame(
      idx = idx,
      file = tile_name,
      nrow = nrow(target_tile),
      ncol = ncol(target_tile),
      ncell = terra::ncell(target_tile),
      target_mode = target_mode,
      n_finite_target = n_finite_target,
      n_finite_source = n_finite_source,
      n_finite_move = n_finite_move,
      n_positive_target = n_positive_target,
      n_positive_source = n_positive_source,
      n_positive_move = n_positive_move,
      n_finite_output_target = n_finite_output_target,
      n_positive_output_target = n_positive_output_target,
      written = ok[idx],
      stringsAsFactors = FALSE
    )
  }

  diagnostics <- do.call(rbind, diagnostics)
  ok_idx <- unname(which(ok))
  bad_idx <- unname(which(!ok))
  bad_files <- if (length(bad_idx) == 0L) {
    data.frame(
      idx = integer(),
      target = character(),
      source = character(),
      move = character(),
      n_finite_target = numeric(),
      n_finite_source = numeric(),
      n_finite_move = numeric(),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      idx = bad_idx,
      target = file.path(dirs$target, paste0("r_", bad_idx, ".asc")),
      source = file.path(dirs$src, paste0("r_", bad_idx, ".asc")),
      move = file.path(dirs$mov, paste0("r_", bad_idx, ".asc")),
      n_finite_target = diagnostics$n_finite_target[bad_idx],
      n_finite_source = diagnostics$n_finite_source[bad_idx],
      n_finite_move = diagnostics$n_finite_move[bad_idx],
      stringsAsFactors = FALSE
    )
  }

  list(
    target_dir = normalizePath(dirs$target),
    mov_dir = normalizePath(dirs$mov),
    src_dir = normalizePath(dirs$src),
    ok_idx = ok_idx,
    bad_idx = bad_idx,
    bad_files = bad_files,
    diagnostics = diagnostics
  )
}
