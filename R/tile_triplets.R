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

  r_target_ext <- terra::extend(r_target, c(overlap_cells, overlap_cells,
                                            overlap_cells, overlap_cells))
  r_mov_ext <- terra::extend(r_mov, c(overlap_cells, overlap_cells,
                                      overlap_cells, overlap_cells))
  r_src_ext <- terra::extend(r_src, c(overlap_cells, overlap_cells,
                                      overlap_cells, overlap_cells))

  tile_extents <- lapply(seq_len(n_tiles), function(i) {
    te <- cs_tiles[i]
    terra::ext(
      terra::xmin(te) - overlap_cells * resx,
      terra::xmax(te) + overlap_cells * resx,
      terra::ymin(te) - overlap_cells * resy,
      terra::ymax(te) + overlap_cells * resy
    )
  })

  count_finite <- function(r) {
    vals <- terra::values(r, mat = FALSE)
    sum(is.finite(vals))
  }

  count_positive <- function(r) {
    vals <- terra::values(r, mat = FALSE)
    sum(is.finite(vals) & vals > 0)
  }

  mask_to_center <- function(r, core_ext) {
    xy <- terra::xyFromCell(r, seq_len(terra::ncell(r)))
    vals <- terra::values(r, mat = FALSE)
    inside <- xy[, 1] >= core_ext$xmin & xy[, 1] < core_ext$xmax &
      xy[, 2] > core_ext$ymin & xy[, 2] <= core_ext$ymax
    vals[!inside] <- NA_real_
    terra::values(r) <- vals
    r
  }

  if (isTRUE(progress)) {
    pb <- txtProgressBar(min = 0, max = n_tiles, style = 3)
  }

  diagnostics <- vector("list", n_tiles)
  ok <- rep(FALSE, n_tiles)

  for (idx in seq_len(n_tiles)) {
    if (exists("pb", inherits = FALSE)) {
      setTxtProgressBar(pb, idx,
                        title = "Writing matched ConScape tile triplets...")
    }

    target_tile <- terra::crop(r_target_ext, tile_extents[[idx]])
    mov_tile <- terra::crop(r_mov_ext, tile_extents[[idx]])
    src_tile <- terra::crop(r_src_ext, tile_extents[[idx]])

    # Target qualities are part of the RSP problem, not merely an output mask:
    # they define the absorbing-state set in ConScape's GridRSP, which controls
    # the columns of Z = (I - W) \ I_target.
    #
    # `target_mode = "full"` keeps the full buffered target raster, so each tile
    # solves a Z whose columns include every target cell in the buffered window.
    # Tiled outputs are then averaged across overlapping tiles by
    # `mosaic_conscape(method = "mean")`. This is a smoothing heuristic; each
    # tile's per-cell sum is biased low (it cannot see targets in other tiles),
    # so the mean is also biased.
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
      mask_to_center(target_tile, terra::ext(cs_tiles[idx]))
    } else {
      target_tile
    }

    tile_name <- paste0("r_", idx, ".asc")
    target_file <- file.path(dirs$target, tile_name)
    mov_file <- file.path(dirs$mov, tile_name)
    src_file <- file.path(dirs$src, tile_name)

    n_finite_target <- count_finite(target_tile)
    n_finite_source <- count_finite(src_tile)
    n_finite_move <- count_finite(mov_tile)
    n_positive_target <- count_positive(target_tile)
    n_positive_source <- count_positive(src_tile)
    n_positive_move <- count_positive(mov_tile)
    n_finite_output_target <- count_finite(output_target_tile)
    n_positive_output_target <- count_positive(output_target_tile)

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

  if (exists("pb", inherits = FALSE)) close(pb)

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
