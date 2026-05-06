#' @title Plot Method for ConScapeResults Objects
#' @description Plot raster outputs stored in a `"ConScapeResults"` object.
#'
#' @param x A `"ConScapeResults"` object returned by [run_conscape()] when
#'   tiled outputs are mosaicked.
#' @param layers Optional character vector naming raster outputs to plot. If
#'   `NULL` (default), all `SpatRaster` elements in `x` are plotted. Common
#'   layer names include `"btwn"`, `"fcon"`, `"btwn_qweighted"`,
#'   `"criticality"`, and sensitivity names such as `"elasticity_quality"`.
#' @param ... Additional arguments passed to `terra::plot`
#'
#' @return Invisibly returns `x` after plotting the selected rasters.
#'
#' @examples
#' \dontrun{
#' # Assuming 'cs_results' is a ConScapeResults object
#' plot(cs_results)
#' plot(cs_results, layers = c("fcon", "elasticity_quality"))
#' }
#'
#' @importFrom terra plot
#' @importFrom graphics par mtext title
#' @method plot ConScapeResults
#' @export
plot.ConScapeResults <- function(x, layers = NULL, ...) {
  # Input validation
  if (!inherits(x, "ConScapeResults"))
    stop("Object must be of class 'ConScapeResults'")

  raster_layers <- vapply(x, inherits, logical(1), what = "SpatRaster")
  raster_names <- names(x)[raster_layers]
  if (length(raster_names) == 0L) {
    stop("Object must contain at least one SpatRaster output")
  }

  if (!is.null(layers)) {
    if (!is.character(layers) || length(layers) == 0L) {
      stop("layers must be NULL or a non-empty character vector")
    }
    missing_layers <- setdiff(layers, raster_names)
    if (length(missing_layers) > 0L) {
      stop(
        "Requested layers not found: ",
        paste(missing_layers, collapse = ", ")
      )
    }
    raster_names <- layers
  }

  n_layers <- length(raster_names)
  n_col <- min(2L, n_layers)
  n_row <- ceiling(n_layers / n_col)

  op <- par(mfrow = c(n_row, n_col), mar = c(3, 3, 3, 5))
  on.exit(par(op))

  display_names <- c(
    btwn = "Betweenness",
    fcon = "Connected Habitat",
    btwn_qweighted = "Quality-weighted Betweenness",
    criticality = "Criticality",
    elasticity_quality = "Quality Elasticity",
    elasticity_affinity = "Affinity Elasticity",
    elasticity_cost = "Cost Elasticity",
    elasticity_affinity_cost_linked = "Linked Affinity-Cost Elasticity",
    elasticity_cost_affinity_linked = "Linked Cost-Affinity Elasticity",
    sensitivity_quality = "Quality Sensitivity",
    sensitivity_affinity = "Affinity Sensitivity",
    sensitivity_cost = "Cost Sensitivity",
    sensitivity_affinity_cost_linked = "Linked Affinity-Cost Sensitivity",
    sensitivity_cost_affinity_linked = "Linked Cost-Affinity Sensitivity"
  )

  for (nm in raster_names) {
    main <- if (nm %in% names(display_names)) display_names[[nm]] else nm
    plot(x[[nm]], main = main, ...)
  }

  invisible(x)
}

