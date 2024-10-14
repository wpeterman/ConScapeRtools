#------------------------------------------------------------------------------
# Sample datasets

#' Sample map with affinities
#'
#' Map representing the affinities (i.e. the likelihood of movement) across
#' each pixel in an area.
#'
#' @format A ASC file. Original CRS UTM zone 17N.
#'
#' @examples
#' (f <- system.file("data/affinity", package = "ConScapeRtools"))
#' r <- terra::rast(f)
#' plot(r)
#'
#' @name affinity.asc
#' @seealso [ConScapeRtools::suitability.asc] \cr
#'

NULL

#' Sample map with suitability
#'
#' Map representing the habitat quality in each pixel in an area.
#'
#' @format A ASC file. Original CRS: UTM zone 17N.
#'
#' @examples
#' (f <- system.file("data/suitability.asc", package = "ConScapeRtools"))
#' r <- terra::rast(f)
#' plot(r)
#'
#' @name suitability.asc
#' @seealso [ConScapeRtools::affinity.asc] \cr
#'
NULL

#' Sample shapefile of target patches
#'
#' A shapefile delineating potential habitat patches.
#'
#' @format A SHP file. Original CRS: UTM zone 17N.
#'
#' @examples
#' (f <- system.file("data/patches.shp", package = "ConScapeRtools"))
#' p <- terra::vect(f)
#' plot(p)
#'
#' @name patches.shp
#'
NULL
