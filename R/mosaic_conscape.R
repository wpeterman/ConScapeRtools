#' Mosaic ConScape tiles into a single raster
#'
#' @description
#' After running [run_conscape()] on tiled landscapes, this function
#' reassembles one directory of tile-level outputs into a single `SpatRaster`.
#' Depending on the reduction method, tiles are retained in full or trimmed to
#' remove overlapping margins, then combined and optionally masked to the
#' original analysis area.
#'
#' @param out_dir Directory containing one ConScape output layer. For results
#'   returned by [run_conscape()], use a path from the canonical `outdirs` map,
#'   such as `result$outdirs$btwn` or `result$outdirs$fcon`. All `*.asc` and
#'   `*.tif` files directly in this directory are treated as tiles belonging to
#'   that surface.
#' @param mask Optional binary `SpatRaster` indicating regions of the
#'   landscape that are considered valid (e.g., potential habitat).
#'   Cells with 0 are set to 0 in the output, and cells with `NA` in the
#'   mask become `NA` in the output. If provided, it should typically be
#'   the same mask written by [conscape_prep()].
#' @param tile_trim Width (in map units) of the overlapping margins. This
#'   should match the `tile_trim` value stored in the
#'   `"ConScapeRtools_prep"` object returned by [conscape_prep()], or an
#'   equivalent value used when creating tiles manually. The margins are
#'   removed for `method = "mosaic"` and `method = "merge"`. They are retained
#'   for `method = "sum"` because overlapping buffer cells contain distinct
#'   contributions that must be added.
#' @param method Character string indicating how to combine overlapping
#'   tiles. `"sum"` adds tile values where they overlap, treating `NA` as 0;
#'   this is the mathematically correct reduction for additive connectivity
#'   metrics when [conscape_prep()] used `target_mode = "center"`. `"mosaic"`
#'   averages values where tiles overlap. [run_conscape()] uses this mean for
#'   additive metrics under `target_mode = "full"` and for all sensitivity or
#'   elasticity surfaces. `"merge"` fills gaps using the first non-`NA` tile
#'   encountered and is mainly useful for diagnostics or non-overlapping
#'   tiles. The default is `"sum"` because the current [conscape_prep()] default
#'   is `target_mode = "center"`.
#' @param crs Optional coordinate reference system to assign to the
#'   merged raster. Can be a proj4string, EPSG code (e.g. `"EPSG:4326"`),
#'   or a CRS taken from a `SpatRaster`. If `NULL` (default), the CRS of
#'   the first tile is retained.
#' @param chunk_size Maximum number of raster sources per reduction batch.
#'   Intermediate batches are reduced recursively, using disk-backed,
#'   double-precision rasters. A value of 1 uses pairwise reduction (2 sources).
#'   Smaller batches limit memory pressure but increase temporary-file I/O.
#'
#' @details
#' For `method = "mosaic"` and `method = "merge"`, `tile_trim` is converted to
#' rows and columns and applied symmetrically to all sides. To avoid removing
#' an entire tile, the trim width is capped below half the tile width and
#' height. For `method = "sum"`, the effective trim is zero so contributions
#' from neighboring tile buffers are preserved.
#'
#' After trimming, tiles are combined into a single `SpatRaster` according to
#' `method`. Choice of `method` reflects how the per-tile ConScape problem was
#' posed:
#'
#' * `"sum"` adds tile values cell-wise, treating `NA` as 0. This is the
#'   correct reduction when each cell appears as a *target* in exactly one
#'   tile (the `target_mode = "center"` case in [conscape_prep()]). Each
#'   tile contributes the source quality times the sum of target quality
#'   and proximity over targets in that tile's center
#'   and similar partial sums for betweenness; summing over tiles recovers the
#'   full source-to-target sum up to paths that exit the tile buffer.
#' * `"mosaic"` averages tile values where they overlap (`terra::mosaic` with
#'   `fun = "mean"`). For additive connectivity metrics, this is the legacy
#'   reduction for `target_mode = "full"`. Omitted targets and changed graph
#'   boundaries affect this approximation; its error need not be negative.
#'   Sensitivity and elasticity surfaces
#'   are tile-local estimates rather than partial sums, so [run_conscape()]
#'   also uses this mean reduction for those layers and for tile-local
#'   criticality. Criticality cannot use center-target sum reduction.
#' * `"merge"` (`terra::merge`) fills cells with the first non-`NA` tile value
#'   encountered. Useful for diagnostics and for non-overlapping tiles.
#'
#' If a `mask` is supplied, the function checks for compatible CRS, then crops
#' and/or resamples the mosaicked raster to match the mask's extent and
#' resolution (using nearest-neighbor resampling if needed) before applying
#' the mask via cell-wise multiplication.
#'
#' This function is normally called internally by [run_conscape()] when
#' `mosaic = TRUE`. For manual post-processing, call [run_conscape()] with
#' `mosaic = FALSE`, select a layer directory from `result$outdirs`, and use
#' the corresponding method recorded in `result$diagnostics$mosaic_method`.
#'
#' @return
#' A `SpatRaster` representing the combined ConScape surface. Overlapping
#' margins are retained for `method = "sum"` and trimmed for the other methods.
#' The optional CRS and mask are applied before the raster is returned.
#' The result is backed by a temporary GeoTIFF for the current R session.
#' Use [terra::writeRaster()] to save a durable result outside the tile and
#' run working directories.
#'
#' @export
#' @examples
#' tile_dir <- tempfile("conscape-output-")
#' dir.create(tile_dir)
#' tile <- terra::rast(nrows = 3, ncols = 3)
#' terra::values(tile) <- seq_len(terra::ncell(tile))
#' terra::writeRaster(tile, file.path(tile_dir, "tile.asc"), overwrite = TRUE)
#' combined <- mosaic_conscape(tile_dir, tile_trim = 0, method = "merge")
#' combined
#'
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
#' analysis_dir <- file.path(tempdir(), "conscape_analysis")
#' tile_dir <- file.path(analysis_dir, "tiles")
#' run_dir <- file.path(analysis_dir, "run")
#'
#' prep <- conscape_prep(tile_d    = td$tile_d,
#'                       tile_trim = td$tile_trim,
#'                       asc_dir   = tile_dir,
#'                       r_target  = habitat,
#'                       r_mov     = affinity,
#'                       r_src     = habitat,
#'                       clear_dir = TRUE,
#'                       landmark  = td$landmark)
#'
#' cs_res <- run_conscape(conscape_prep  = prep,
#'                        out_dir        = run_dir,
#'                        theta          = td$theta,
#'                        distance_scale = td$distance_scale,
#'                        jl_home        = jl_home,
#'                        mosaic         = FALSE)
#'
#' ## Manually reassemble one layer using run_conscape()'s recorded method.
#' mask <- terra::rast(file.path(prep$asc_dir, "mask", "mask.asc"))
#' cs_btwn <- mosaic_conscape(out_dir   = cs_res$outdirs$btwn,
#'                            mask      = mask,
#'                            tile_trim = prep$tile_trim,
#'                            method    = cs_res$diagnostics$mosaic_method[["btwn"]],
#'                            crs       = terra::crs(habitat))
#' cs_fcon <- mosaic_conscape(out_dir   = cs_res$outdirs$fcon,
#'                            mask      = mask,
#'                            tile_trim = prep$tile_trim,
#'                            method    = cs_res$diagnostics$mosaic_method[["fcon"]],
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

  # Validate the layer directory and reduction controls.
  if (!is.character(out_dir) || length(out_dir) != 1L || is.na(out_dir) ||
      !nzchar(out_dir) || !dir.exists(out_dir)) {
    stop("out_dir must be an existing directory containing one output layer", call. = FALSE)
  }
  if (!is.numeric(tile_trim) || length(tile_trim) != 1L ||
      is.na(tile_trim) || !is.finite(tile_trim) || tile_trim < 0) {
    stop("tile_trim must be a single non-negative numeric value", call. = FALSE)
  }
  if (!is.numeric(chunk_size) || length(chunk_size) != 1L ||
      is.na(chunk_size) || !is.finite(chunk_size) || chunk_size < 1 ||
      chunk_size != floor(chunk_size) || chunk_size > .Machine$integer.max) {
    stop("chunk_size must be a single positive integer", call. = FALSE)
  }
  chunk_size <- as.integer(chunk_size)

  # List all raster files
  tile_files <- list.files(out_dir, pattern = "\\.asc$|\\.tif$", full.names = TRUE)
  if (length(tile_files) == 0) {
    stop("No raster files found in out_dir", call. = FALSE)
  }

  first_tile <- terra::rast(tile_files[1])
  if (!is.null(mask) && (!inherits(mask, "SpatRaster") || terra::nlyr(mask) != 1L)) {
    stop("mask must be a single-layer SpatRaster", call. = FALSE)
  }

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
    if (terra::nlyr(x) != 1L || !terra::same.crs(first_tile, x) ||
        !isTRUE(all.equal(terra::res(first_tile), terra::res(x))) ||
        !isTRUE(all.equal(terra::origin(first_tile), terra::origin(x)))) {
      stop("Tiles must be single-layer rasters with matching CRS, resolution, and grid origin.",
           call. = FALSE)
    }
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

  # Keep only bounded batches of file-backed rasters. Intermediates use doubles
  # so chunking does not silently round contributions to 32-bit floats.
  tmp_dir <- tempfile("conscape-mosaic-")
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
  temp_raster <- function() tempfile(tmpdir = tmp_dir, fileext = ".tif")
  wopt <- list(datatype = "FLT8S")
  fan_in <- max(2L, chunk_size)

  combine <- function(files, fun, trim = FALSE) {
    layers <- lapply(files, read_tile, trim = trim)
    if (identical(fun, "merge")) {
      terra::merge(terra::sprc(layers), filename = temp_raster(), wopt = wopt)
    } else {
      terra::mosaic(terra::sprc(layers), fun = fun,
                    filename = temp_raster(), wopt = wopt)
    }
  }
  reduce_files <- function(files, fun) {
    while (length(files) > 1L) {
      groups <- split(files, ceiling(seq_along(files) / fan_in))
      files <- vapply(groups, function(group) {
        if (length(group) == 1L) return(group)
        terra::sources(combine(group, fun))[1L]
      }, character(1))
    }
    terra::rast(files[[1L]])
  }

  groups <- split(tile_files, ceiling(seq_along(tile_files) / fan_in))
  sums <- counts <- character(length(groups))
  for (i in seq_along(groups)) {
    if (method == "merge") {
      part <- combine(groups[[i]], "merge", trim = TRUE)
    } else {
      layers <- lapply(groups[[i]], read_tile)
      # Raster arithmetic streams through terra; no full values() or xy matrix.
      clean <- lapply(layers, function(x) {
        valid <- if (method == "sum") is.finite(x) else !is.na(x)
        terra::ifel(valid, x, 0, filename = temp_raster(), wopt = wopt)
      })
      part <- terra::mosaic(terra::sprc(clean), fun = "sum",
                            filename = temp_raster(), wopt = wopt)
      if (method == "mosaic") {
        present <- lapply(layers, function(x) {
          terra::ifel(is.na(x), 0, 1, filename = temp_raster(), wopt = wopt)
        })
        count <- terra::mosaic(terra::sprc(present), fun = "sum",
                               filename = temp_raster(), wopt = wopt)
        counts[i] <- terra::sources(count)[1L]
      }
    }
    sums[i] <- terra::sources(part)[1L]
  }
  r_out <- reduce_files(sums, if (method == "merge") "merge" else "sum")
  if (method == "mosaic") {
    count <- reduce_files(counts, "sum")
    r_out <- terra::ifel(count > 0, r_out / count, NA,
                         filename = temp_raster(), wopt = wopt)
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
    if (has_mask_crs && has_out_crs && !terra::same.crs(mask, r_out)) {
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
  # The returned raster must not reference files removed by on.exit(). Let terra
  # manage this final backing file for the remainder of the R session.
  terra::writeRaster(r_out, tempfile("conscape-mosaic-result-", fileext = ".tif"),
                      wopt = wopt)
}
