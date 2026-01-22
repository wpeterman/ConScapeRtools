#' Validate ConScape tile triplets and return indices to keep
#' @noRd
validate_conscape_tiles_idx <- function(target_dir,
                                        source_dir,
                                        move_dir,
                                        pattern_target = "\\.asc$",
                                        pattern_source = "\\.asc$",
                                        pattern_move   = "\\.asc$",
                                        delete_bad = TRUE,
                                        quiet = TRUE) {
  stopifnot(dir.exists(target_dir), dir.exists(source_dir), dir.exists(move_dir))

  .finite_n <- function(f) {
    r <- terra::rast(f)
    as.numeric(terra::global(r, fun = "notNA", na.rm = TRUE)[1, 1])
  }

  t_files <- sort(list.files(target_dir, pattern_target, full.names = TRUE))
  s_files <- sort(list.files(source_dir, pattern_source, full.names = TRUE))
  m_files <- sort(list.files(move_dir,   pattern_move,   full.names = TRUE))

  if (length(t_files) == 0L) stop("No target tiles found in: ", target_dir)
  if (length(s_files) == 0L) stop("No source tiles found in: ", source_dir)
  if (length(m_files) == 0L) stop("No movement tiles found in: ", move_dir)

  # Require 1:1:1 pairing by order
  n <- length(t_files)
  if (length(s_files) != n || length(m_files) != n) {
    stop("Unequal number of tiles across dirs: ",
         "target=", length(t_files), ", source=", length(s_files), ", move=", length(m_files), ". ",
         "This validator assumes files are paired by sort order.")
  }

  n_fin_t <- vapply(t_files, .finite_n, numeric(1))
  n_fin_s <- vapply(s_files, .finite_n, numeric(1))
  n_fin_m <- vapply(m_files, .finite_n, numeric(1))

  ok <- (n_fin_t > 0) & (n_fin_s > 0) & (n_fin_m > 0)
  ok_idx  <- which(ok)
  bad_idx <- which(!ok)

  bad_df <- data.frame(
    idx = bad_idx,
    target = t_files[bad_idx],
    source = s_files[bad_idx],
    move   = m_files[bad_idx],
    n_finite_target = n_fin_t[bad_idx],
    n_finite_source = n_fin_s[bad_idx],
    n_finite_move   = n_fin_m[bad_idx],
    stringsAsFactors = FALSE
  )

  if (delete_bad && nrow(bad_df) > 0L) {
    file.remove(bad_df$target, bad_df$source, bad_df$move)
  }

  if (!quiet) {
    message(sprintf("Tile validation: %d ok, %d bad (deleted=%s)",
                    length(ok_idx), length(bad_idx), delete_bad))
  }

  list(
    ok_idx = ok_idx,
    bad_idx = bad_idx,
    bad_files = bad_df
  )
}


#' Ensure cs_tiles has a tile id column
#' @noRd
ensure_tile_id <- function(cs_tiles, id_col = NULL) {
  n <- nrow(cs_tiles)
  if (n == 0L) stop("cs_tiles has 0 rows.")

  # If user supplied an id_col, trust but verify
  if (!is.null(id_col)) {
    if (!id_col %in% names(cs_tiles)) stop("id_col not found in cs_tiles: ", id_col)
    return(list(cs_tiles = cs_tiles, id_col = id_col))
  }

  # Heuristics: common names
  candidates <- c("tile_num", "tile_id", "tile", "id", "ID", "TileID", "Tile_Id", "name")
  candidates <- candidates[candidates %in% names(cs_tiles)]

  if (length(candidates)) {
    # pick the first that looks usable (unique-ish)
    for (cc in candidates) {
      v <- cs_tiles[[cc]]
      # accept numeric/integer, or character that can be coerced
      if (is.numeric(v) || is.integer(v)) {
        if (length(unique(v)) == n) return(list(cs_tiles = cs_tiles, id_col = cc))
      } else if (is.character(v)) {
        vv <- suppressWarnings(as.integer(v))
        if (!anyNA(vv) && length(unique(vv)) == n) {
          cs_tiles[[cc]] <- vv
          return(list(cs_tiles = cs_tiles, id_col = cc))
        }
      }
    }
  }

  # Otherwise: create sequential tile_id (stable within this prep run)
  cs_tiles$tile_id <- seq_len(n)
  list(cs_tiles = cs_tiles, id_col = "tile_id")
}

