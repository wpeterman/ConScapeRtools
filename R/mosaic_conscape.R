#' Mosaic ConScape tiles into a single raster
#'
#' @description
#' After running [run_conscape()] on tiled landscapes, this function
#' reassembles the tile-level outputs into a single `SpatRaster`. Tiles
#' are optionally trimmed to remove overlapping margins, combined either
#' by merging or averaging, and then (optionally) masked to the original
#' analysis area.
#'
#' @param out_dir Directory where ConScape tile outputs (e.g., `"btwn"`
#'   or `"fcon"` subdirectories produced by [run_conscape()]) were
#'   written. All `*.asc` and `*.tif` files in this directory are
#'   treated as tiles belonging to a single surface.
#' @param mask Optional binary `SpatRaster` indicating regions of the
#'   landscape that are considered valid (e.g., potential habitat).
#'   Cells with 0 are set to 0 in the output, and cells with `NA` in the
#'   mask become `NA` in the output. If provided, it should typically be
#'   the same mask written by [conscape_prep()].
#' @param tile_trim Width (in map units) of the overlapping margins to
#'   remove from each tile before mosaicking. This should match the
#'   `tile_trim` value stored in the `"ConScapeRtools_prep"` object
#'   returned by [conscape_prep()], or an equivalent value used when
#'   creating tiles manually.
#' @param method Character string indicating how to combine overlapping
#'   tiles. `"sum"` adds tile values where they overlap, treating `NA` as 0;
#'   this is the mathematically correct reduction for ConScape outputs when
#'   each cell appears as a *target* in exactly one tile (i.e., when
#'   `conscape_prep()` was called with `target_mode = "center"`). `"mosaic"`
#'   averages values where tiles overlap; this is a smoothing heuristic that
#'   should only be used with `target_mode = "full"` tiles. `"merge"` fills
#'   gaps using the first non-`NA` tile encountered. The default is `"sum"`
#'   because the current `conscape_prep()` default is `target_mode = "center"`.
#' @param crs Optional coordinate reference system to assign to the
#'   merged raster. Can be a proj4string, EPSG code (e.g. `"EPSG:4326"`),
#'   or a CRS taken from a `SpatRaster`. If `NULL` (default), the CRS of
#'   the first tile is retained.
#' @param chunk_size Maximum number of tile rasters to read and combine at
#'   once. If there are more tiles than `chunk_size`, temporary intermediate
#'   mosaics are written and then combined in a final pass.
#'
#' @details
#' For each tile, `tile_trim` is converted to a number of rows and
#' columns and trimmed symmetrically from all sides. To avoid trimming
#' away entire tiles, the trim width is capped at at most half the tile
#' width/height in cells. If `tile_trim` is too large relative to tile
#' dimensions, no trimming is applied to that tile.
#'
#' After trimming, tiles are combined into a single `SpatRaster` according to
#' `method`. Choice of `method` reflects how the per-tile ConScape problem was
#' posed:
#'
#' * `"sum"` adds tile values cell-wise, treating `NA` as 0. This is the
#'   correct reduction when each cell appears as a *target* in exactly one
#'   tile (the `target_mode = "center"` case in [conscape_prep()]). Each
#'   tile contributes `q_src(i) * Σ_{j in tile center} q_tgt(j) * proximity(i,j)`
#'   and similar partial sums for betweenness; summing over tiles recovers the
#'   full source-to-target sum up to paths that exit the tile buffer.
#' * `"mosaic"` averages tile values where they overlap (`terra::mosaic` with
#'   `fun = "mean"`). This is a smoothing heuristic appropriate only for
#'   `target_mode = "full"` tiles. Because each tile undercounts targets that
#'   live in other tiles, the average is biased; use this only when reproducing
#'   the legacy ConScapeRtools workflow.
#' * `"merge"` (`terra::merge`) fills cells with the first non-`NA` tile value
#'   encountered. Useful for diagnostics and for non-overlapping tiles.
#'
#' If a `mask` is supplied, the function checks for compatible CRS, then crops
#' and/or resamples the mosaicked raster to match the mask's extent and
#' resolution (using nearest-neighbour resampling if needed) before applying
#' the mask via cell-wise multiplication.
#'
#' This function is typically called internally by [run_conscape()] when
#' `mosaic = TRUE`, but it can also be used directly to post-process
#' ConScape tile outputs.
#'
#' @return
#' A `SpatRaster` representing the mosaicked ConScape surface, with
#' overlapping tile margins removed, tiles combined according to
#' `method`, an optional CRS set, and (if provided) masked to the input
#' `mask`.
#'
#' @export
#' @examples
#' \dontrun{
#' library(ConScapeRtools)
#'
#' ## Import example data
#' s <- system.file("extdata", "suitability.asc", package = "ConScapeRtools")
#' habitat <- terra::rast(s)
#'
#' a <- system.file("extdata", "affinity.asc", package = "ConScapeRtools")
#' affinity <- terra::rast(a)
#'
#' jl_home <- "/path/to/julia/bin"
#'
#' td <- tile_design(r_mov    = affinity,
#'                   r_target = habitat,
#'                   max_d    = 7000,
#'                   theta    = 0.1,
#'                   jl_home  = jl_home,
#'                   landmark = 5L)
#'
#' prep <- conscape_prep(tile_d    = td$tile_d,
#'                       tile_trim = td$tile_trim,
#'                       r_target  = habitat,
#'                       r_mov     = affinity,
#'                       r_src     = habitat,
#'                       clear_dir = TRUE,
#'                       landmark  = td$landmark)
#'
#' cs_res <- run_conscape(conscape_prep  = prep,
#'                        out_dir        = "conscape_out",
#'                        theta          = td$theta,
#'                        distance_scale = td$distance_scale,
#'                        jl_home        = jl_home)
#'
#' ## Manually reassemble tile outputs (done automatically when mosaic = TRUE)
#' mask <- terra::rast(file.path(prep$asc_dir, "mask", "mask.asc"))
#' cs_btwn <- mosaic_conscape(out_dir   = cs_res$outdir_btwn,
#'                            mask      = mask,
#'                            tile_trim = prep$tile_trim,
#'                            method    = "mosaic",
#'                            crs       = terra::crs(habitat))
#' cs_fcon <- mosaic_conscape(out_dir   = cs_res$outdir_fcon,
#'                            mask      = mask,
#'                            tile_trim = prep$tile_trim,
#'                            crs       = terra::crs(habitat))
#' plot(c(cs_btwn, cs_fcon))
#' }
#' @seealso [conscape_prep()], [run_conscape()]
#' @author Bill Peterman
#' @importFrom terra merge


mosaic_conscape <- function(out_dir,
                            mask = NULL,
                            tile_trim,
                            method = c("sum", "mosaic", "merge"),
                            crs = NULL,
                            chunk_size = 64L) {
  method <- match.arg(method)

  # Check tile_trim
  if (!is.numeric(tile_trim) || length(tile_trim) != 1 || tile_trim < 0) {
    stop("tile_trim must be a single non-negative numeric value")
  }
  if (!is.numeric(chunk_size) || length(chunk_size) != 1 ||
      is.na(chunk_size) || chunk_size < 1) {
    stop("chunk_size must be a single positive integer")
  }
  chunk_size <- as.integer(chunk_size)

  # List all raster files
  tile_files <- list.files(out_dir, pattern = "\\.asc$|\\.tif$", full.names = TRUE)
  if (length(tile_files) == 0) stop("No raster files found in out_dir")

  first_tile <- terra::rast(tile_files[1])

  # Trimming semantics depend on the reduction method.
  #
  # `mosaic` (mean) and `merge`: each tile has full-buffered output and
  # overlapping tiles produce redundant estimates of the same per-cell value.
  # Trimming the buffer keeps only the interior cells the tile is most
  # confident about and avoids the seam where buffered neighbors disagree.
  #
  # `sum`: per-tile output values in the buffer region are NOT redundant
  # estimates of the same quantity. With `target_mode = "center"`, each
  # buffer cell's fcon/btwn value represents that cell's contribution to its
  # tile's *center* targets via the RSP. Different tiles compute disjoint
  # contributions (because centers tile the landscape), so to recover the
  # untiled per-cell sum we must sum *every* per-tile contribution,
  # including buffer-cell contributions from neighboring tiles. Trimming
  # would discard those contributions, leaving each cell with only its own
  # tile's contribution. We therefore set the effective trim to zero for
  # the sum reduction.
  effective_tile_trim <- if (identical(method, "sum")) 0 else tile_trim
  if (effective_tile_trim > 0) {
    resx <- terra::res(first_tile)[1]
    resy <- terra::res(first_tile)[2]
    trim_cols <- round(effective_tile_trim / resx)
    trim_rows <- round(effective_tile_trim / resy)
  } else {
    trim_cols <- 0L
    trim_rows <- 0L
  }

  read_tile <- function(file, trim = TRUE) {
    x <- terra::rast(file)
    if (isTRUE(trim) && (trim_cols > 0L || trim_rows > 0L)) {
      e   <- terra::ext(x)
      rx  <- terra::res(x)[1]
      ry  <- terra::res(x)[2]
      nc  <- ncol(x)
      nr  <- nrow(x)

      # cap trim so we never remove more than half the tile
      tc <- min(trim_cols, floor((nc - 1L) / 2))
      tr <- min(trim_rows, floor((nr - 1L) / 2))

      if (tc <= 0 && tr <= 0) {
        return(x)
      }

      ext_trim <- terra::ext(
        e$xmin + tc * rx,
        e$xmax - tc * rx,
        e$ymin + tr * ry,
        e$ymax - tr * ry
      )
      x <- terra::crop(x, ext_trim)
    }
    x
  }

  # NA-as-zero helper used by the "sum" reduction. ConScape writes NA into
  # cells outside the tile's analysis area (e.g., cells masked out by
  # target_threshold or outside the buffered window). For a sum mosaic those
  # absent contributions must be 0, not NA, or terra::mosaic would propagate
  # NA across the whole landscape.
  zero_out_na <- function(r) {
    vals <- terra::values(r, mat = FALSE)
    vals[!is.finite(vals)] <- 0
    out <- r
    terra::values(out) <- vals
    out
  }

  combine_files <- function(files, trim = TRUE) {
    r_list <- lapply(files, read_tile, trim = trim)
    rc <- terra::sprc(r_list)
    if (method == "sum") {
      rc_zero <- terra::sprc(lapply(r_list, zero_out_na))
      terra::mosaic(rc_zero, fun = "sum")
    } else if (method == "mosaic") {
      terra::mosaic(rc, fun = "mean")
    } else {
      terra::merge(rc)
    }
  }

  combine_sum_count <- function(files, trim = TRUE) {
    r_list <- lapply(files, read_tile, trim = trim)
    sum_list <- lapply(r_list, function(x) {
      y <- x
      vals <- terra::values(y, mat = FALSE)
      vals[is.na(vals)] <- 0
      terra::values(y) <- vals
      y
    })
    count_list <- lapply(r_list, function(x) {
      y <- x
      vals <- terra::values(y, mat = FALSE)
      vals <- ifelse(is.na(vals), 0, 1)
      terra::values(y) <- vals
      y
    })
    list(
      sum = terra::mosaic(terra::sprc(sum_list), fun = "sum"),
      count = terra::mosaic(terra::sprc(count_list), fun = "sum")
    )
  }

  if (length(tile_files) <= chunk_size) {
    r_out <- combine_files(tile_files, trim = TRUE)
  } else if (method %in% c("merge", "sum")) {
    # `merge` and `sum` are both associative: combining chunks of partial
    # results is equivalent to one-shot combine. For `sum`, NAs are already
    # treated as zeros inside `combine_files`, so the same recursion works.
    chunks <- split(tile_files, ceiling(seq_along(tile_files) / chunk_size))
    tmp_dir <- tempfile("conscape_mosaic_chunks_")
    dir.create(tmp_dir, recursive = TRUE)
    on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)

    tmp_files <- character(length(chunks))
    for (i in seq_along(chunks)) {
      chunk_rast <- combine_files(chunks[[i]], trim = TRUE)
      tmp_files[i] <- file.path(tmp_dir, paste0("chunk_", i, ".tif"))
      terra::writeRaster(chunk_rast, tmp_files[i], overwrite = TRUE)
    }
    r_out <- combine_files(tmp_files, trim = FALSE)
  } else {
    # method == "mosaic" (mean). Chunked mean = chunked sum / chunked count.
    chunks <- split(tile_files, ceiling(seq_along(tile_files) / chunk_size))
    tmp_dir <- tempfile("conscape_mosaic_chunks_")
    dir.create(tmp_dir, recursive = TRUE)
    on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)

    sum_files <- character(length(chunks))
    count_files <- character(length(chunks))
    for (i in seq_along(chunks)) {
      chunk <- combine_sum_count(chunks[[i]], trim = TRUE)
      sum_files[i] <- file.path(tmp_dir, paste0("sum_", i, ".tif"))
      count_files[i] <- file.path(tmp_dir, paste0("count_", i, ".tif"))
      terra::writeRaster(chunk$sum, sum_files[i], overwrite = TRUE)
      terra::writeRaster(chunk$count, count_files[i], overwrite = TRUE)
    }
    total <- combine_sum_count(sum_files, trim = FALSE)$sum
    counts <- combine_sum_count(count_files, trim = FALSE)$sum
    r_out <- total / counts
    r_out[counts == 0] <- NA
  }

  # Apply CRS if provided
  if (!is.null(crs)) terra::crs(r_out) <- crs

  # Apply mask if provided
  # if (!is.null(mask)) r_out <- r_out * mask
  if (!is.null(mask)) {

    # Ensure same CRS when both rasters carry a non-empty CRS.
    mask_crs <- terra::crs(mask)
    out_crs <- terra::crs(r_out)
    has_mask_crs <- !is.null(mask_crs) && length(mask_crs) == 1L && !is.na(mask_crs) && nzchar(mask_crs)
    has_out_crs <- !is.null(out_crs) && length(out_crs) == 1L && !is.na(out_crs) && nzchar(out_crs)
    if (has_mask_crs && has_out_crs && mask_crs != out_crs) {
      stop("CRS of mask and r_out differ")
    }
    if (has_mask_crs && !has_out_crs) terra::crs(r_out) <- mask_crs
    if (!has_mask_crs && has_out_crs) terra::crs(mask) <- out_crs
    if (!has_mask_crs && !has_out_crs) {
      terra::crs(r_out) <- ""
      terra::crs(mask) <- ""
    }

    # If extents/resolutions differ slightly, bring r_out onto the mask grid
    same_geom <- isTRUE(all.equal(terra::ext(r_out), terra::ext(mask))) &&
      isTRUE(all.equal(terra::res(r_out), terra::res(mask)))

    if (!same_geom) {
      # try to bring r_out onto mask grid
      r_out <- terra::crop(r_out, mask)
      same_geom <- isTRUE(all.equal(terra::ext(r_out), terra::ext(mask))) &&
        isTRUE(all.equal(terra::res(r_out), terra::res(mask)))

      if (!same_geom) {
        # last resort: resample r_out to mask grid
        r_out <- terra::resample(r_out, mask, method = "near")
      }
    }

    # Now extents & res match, safe to multiply or use terra::mask()
    r_out <- r_out * mask
    # r_out <- mask(r_out, mask, maskvalues = 0, updatevalue = NA)
  }


  # Clean names
  names(r_out) <- sub("-.*", "", names(r_out))
  return(r_out)
}
