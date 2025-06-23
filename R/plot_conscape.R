#' @title Plot Method for ConScapeResults Objects
#' @description Creates a two-panel plot showing betweenness and functional connectivity rasters.
#'
#' @param x A ConScapeResults object containing `$btwn` and `$fcon` SpatRaster objects
#' @param ... Additional arguments passed to `terra::plot`
#'
#' @return A two-panel plot (side-by-side) of the spatial rasters
#'
#' @examples
#' \dontrun{
#' # Assuming 'cs_results' is a ConScapeResults object
#' plot(cs_results)
#' plot(cs_results, main = "My Connectivity Results", col = terrain.colors(50))
#' }
#'
#' @importFrom terra plot
#' @importFrom graphics par mtext title
#' @export
plot.ConScapeResults <- function(x, ...) {
  # Input validation
  if (!inherits(x, "ConScapeResults"))
    stop("Object must be of class 'ConScapeResults'")
  if (!all(c("btwn", "fcon") %in% names(x)))
    stop("Object must contain 'btwn' and 'fcon' elements")

  # Set up plot layout
  op <- par(mfrow = c(1, 2), mar = c(3, 3, 3, 5))
  on.exit(par(op))

  # Plot each raster
  plot(x$btwn, main = "Betweenness", ...)
  plot(x$fcon, main = "Functional Connectivity", ...)

  invisible(x)
}

