#' Build a spatial grid for use with TwilightFreeGrid
#'
#' Creates a \code{SpatRaster} over a lon/lat bounding box at a specified
#' resolution, optionally masked to retain only sea or land cells.
#'
#' @param lon Numeric vector of length 2: \code{c(min_lon, max_lon)} in decimal degrees.
#' @param lat Numeric vector of length 2: \code{c(min_lat, max_lat)} in decimal degrees.
#' @param cell.size Cell size in decimal degrees (both longitude and latitude).
#' @param mask Character string: \code{"none"} (default, no masking), \code{"sea"} (retain
#'   only ocean cells), or \code{"land"} (retain only land cells).
#'
#' @return A \code{SpatRaster} (WGS84) with value 1 in retained cells and \code{NA} in
#'   masked cells. Pass directly to the \code{grid} argument of \link{TwilightFreeGrid}.
#'
#' @details
#' Land/sea masking is performed by rasterising
#' \code{rnaturalearth::ne_countries(scale = "medium")} onto the grid.  Cells
#' whose centre falls inside a country polygon are treated as land; all others
#' are treated as sea.  Small islands may be missed at coarse resolutions.
#'
#' @importFrom terra rast values rasterize vect
#' @export
#' @examples
#' \dontrun{
#'   # 2-degree grid over the Southern Ocean, sea cells only
#'   grid <- makeGrid(lon = c(-180, 180), lat = c(-70, -30), cell.size = 2, mask = "sea")
#'   TwilightFreeGrid(date_time, light, grid = grid, ...)
#' }
makeGrid <- function(lon, lat, cell.size, mask = c("none", "sea", "land")) {
  mask <- match.arg(mask)
  if (length(lon) != 2 || length(lat) != 2) stop("lon and lat must each be length-2 vectors: c(min, max)")
  if (lon[1] >= lon[2]) stop("lon[1] must be less than lon[2]")
  if (lat[1] >= lat[2]) stop("lat[1] must be less than lat[2]")
  if (cell.size <= 0) stop("cell.size must be positive")

  r <- terra::rast(
    xmin = lon[1], xmax = lon[2],
    ymin = lat[1], ymax = lat[2],
    resolution = cell.size,
    crs = "EPSG:4326"
  )
  terra::values(r) <- 1

  if (mask != "none") {
    land <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
    land_r <- terra::rasterize(terra::vect(land), r, field = 1, background = 0)
    if (mask == "sea") {
      r[terra::values(land_r) == 1] <- NA
    } else {
      r[terra::values(land_r) == 0] <- NA
    }
  }

  r
}
