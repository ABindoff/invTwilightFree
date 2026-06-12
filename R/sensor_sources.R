#' Bathymetry source: seafloor depth from ETOPO 2022 (NOAA, via marmap)
#'
#' Builds a static location source that returns seafloor depth (positive-down,
#' metres) at candidate cells, for use as a dive-depth constraint with
#' [floor_rule()]. On first evaluation it downloads the ETOPO 2022 15 arc-second
#' relief grid for the bounding box of the supplied coordinates (via
#' \code{marmap::getNOAA.bathy}), caches it, and reuses it for every subsequent
#' knot (bathymetry does not change with date).
#'
#' Land cells (positive elevation) are set to depth `0`, so a cell on land
#' cannot satisfy a non-trivial dive under [floor_rule()]. For a clean land/sea
#' boundary independent of dive depth, pair with a sea mask as well.
#'
#' `marmap` is an optional dependency; this function errors with an install hint
#' if it is not available. Supply your own raster via `rast` (for example a
#' GEBCO 2024 tile) to avoid the download entirely or to work offline.
#'
#' @param resolution Grid resolution in arc-minutes passed to
#'   \code{marmap::getNOAA.bathy}. Default `4` (about 7 km) is ample for a depth
#'   feasibility constraint and downloads quickly; lower it (down to `0.25`) for
#'   finer coastlines.
#' @param cache_dir Directory for the cached grid (marmap `keep = TRUE` writes a
#'   CSV here). Defaults to a per-user cache via
#'   \code{tools::R_user_dir("invTwilightFree", "cache")}.
#' @param antimeridian If `TRUE`, fetch across the 180 degree line; your grid
#'   longitudes must then use the same convention (the seal example spanning
#'   130 to 185 needs this). Default `FALSE`.
#' @param pad Degrees of padding added around the coordinate bounding box when
#'   fetching, to avoid edge `NA`s. Default `1`.
#' @param rast Optional user-supplied \code{RasterLayer} of seafloor depth
#'   (positive-down, land as `0` or `NA`). If given, no download occurs and this
#'   raster is used directly. Mainly for offline use, custom bathymetry, and
#'   testing.
#'
#' @return A source closure of class `tf_source`: `function(lon, lat, date)`
#'   returning a numeric vector of per-cell seafloor depths (positive-down,
#'   metres); `date` is ignored. Pair with [floor_rule()] inside
#'   [location_term()].
#'
#' @seealso [floor_rule()]
#' @examples
#' \dontrun{
#' # Download ETOPO 2022 for the candidate region on first use.
#' bs <- bathy_source(resolution = 4)
#' depth <- bs(lon = c(150, 155), lat = c(-45, -40), date = NULL)
#' }
#'
#' # Offline: supply a synthetic depth raster (positive-down).
#' r <- raster::raster(nrows = 2, ncols = 2,
#'                     xmn = 149, xmx = 151, ymn = -46, ymx = -44)
#' raster::values(r) <- c(50, 1500, 0, 3000)   # metres, land = 0
#' bs <- bathy_source(rast = r)
#' bs(lon = c(149.5, 150.5), lat = c(-45.5, -44.5))
#' @export
bathy_source <- function(resolution = 4, cache_dir = NULL,
                         antimeridian = FALSE, pad = 1, rast = NULL) {
  if (!is.null(rast) && !inherits(rast, "RasterLayer")) {
    stop("`rast` must be a RasterLayer of seafloor depth (positive-down)")
  }
  if (is.null(cache_dir)) {
    cache_dir <- tools::R_user_dir("invTwilightFree", "cache")
  }

  state <- new.env(parent = emptyenv())
  state$rast <- rast

  fetch <- function(lon, lat) {
    if (!requireNamespace("marmap", quietly = TRUE)) {
      stop("bathy_source() needs the 'marmap' package. ",
           "Install it with install.packages('marmap'), ",
           "or pass your own raster via rast=.")
    }
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    b <- marmap::getNOAA.bathy(
      lon1 = min(lon) - pad, lon2 = max(lon) + pad,
      lat1 = min(lat) - pad, lat2 = max(lat) + pad,
      resolution = resolution, keep = TRUE, path = cache_dir,
      antimeridian = antimeridian
    )
    r <- marmap::as.raster(b)        # elevation: negative = depth, positive = land
    depth <- -1 * r                  # seafloor depth, positive-down
    depth[depth < 0] <- 0            # land (elevation > 0) -> depth 0
    state$rast <- depth
    invisible(depth)
  }

  src <- function(lon, lat, date = NULL) {
    if (is.null(state$rast)) fetch(lon, lat)
    as.numeric(raster::extract(state$rast, cbind(lon, lat)))
  }
  class(src) <- c("tf_source", "function")
  src
}


#' Sea surface temperature source: NOAA OISST v2.1 (daily, via rerddap/ERDDAP)
#'
#' Builds a time-indexed location source that returns sea surface temperature
#' (degrees C) at candidate cells for each knot's date, for matching against
#' tag-recorded temperature with [student_rule()]. SST layers are fetched lazily
#' per unique date from the NOAA OISST v2.1 daily product on an ERDDAP server
#' (via \code{rerddap::griddap}) and cached for the session; OISST is a
#' gap-filled (interpolated) product, so cells are rarely missing except over
#' land.
#'
#' `rerddap` is an optional dependency; this function errors with an install
#' hint if it is missing. Supply your own data via `rast` (a single
#' \code{RasterLayer} used for all dates, e.g. a climatology, or a named list of
#' \code{RasterLayer}s keyed by \code{"YYYY-MM-DD"}) to work offline or use a
#' different product.
#'
#' Longitude convention: the default dataset uses -180..180. For a grid that
#' crosses the 180 degree line (e.g. the seal example spanning 130 to 185), use
#' the 0..360 dataset (`dataset = "ncdcOisst21Agg"`) and a matching grid.
#'
#' Note for the calling framework: SST need only be evaluated for knots that
#' have a surface temperature reading. `build_aux_matrix()` should skip the
#' source (and the fetch) for knots whose reduced tag value is `NA`, since
#' [student_rule()] returns no information there anyway.
#'
#' @param dataset ERDDAP dataset id. Default `"ncdcOisst21Agg_LonPM180"`
#'   (NOAA OISST v2.1, daily, 0.25 degree, longitudes -180..180).
#' @param url ERDDAP server base URL. Default NOAA CoastWatch West Coast node.
#' @param field SST variable name in the dataset. Default `"sst"`.
#' @param pad Degrees of padding around the coordinate bounding box when
#'   fetching, to avoid edge `NA`s. Default `0.5`.
#' @param cache_dir Reserved for on-disk caching of fetched layers; `rerddap`
#'   also keeps its own HTTP cache. Default `NULL`.
#' @param rast Optional offline data: a single \code{RasterLayer} (used for all
#'   dates) or a named list of \code{RasterLayer}s keyed by `"YYYY-MM-DD"`. If
#'   given, no download occurs. Mainly for offline use and testing.
#'
#' @return A source closure of class `tf_source`: `function(lon, lat, date)`
#'   returning a numeric vector of per-cell SST (degrees C) for the knot date.
#'   Pair with [student_rule()] inside [location_term()].
#'
#' @seealso [student_rule()]
#' @examples
#' \dontrun{
#' ss <- sst_source()                       # NOAA OISST via ERDDAP
#' sst <- ss(lon = c(150, 152), lat = c(-45, -44), date = as.Date("2024-01-15"))
#' }
#'
#' # Offline: a single static layer used for every date.
#' r <- raster::raster(nrows = 2, ncols = 2,
#'                     xmn = 149, xmx = 151, ymn = -46, ymx = -44)
#' raster::values(r) <- c(11.2, 11.8, 12.1, 12.6)
#' ss <- sst_source(rast = r)
#' ss(lon = c(149.5, 150.5), lat = c(-44.5, -45.5), date = Sys.Date())
#' @export
sst_source <- function(dataset = "ncdcOisst21Agg_LonPM180",
                       url = "https://coastwatch.pfeg.noaa.gov/erddap/",
                       field = "sst", pad = 0.5,
                       cache_dir = NULL, rast = NULL) {
  if (!is.null(rast) &&
      !inherits(rast, "RasterLayer") &&
      !(is.list(rast) && length(rast) > 0L)) {
    stop("`rast` must be a RasterLayer or a non-empty named list of RasterLayers")
  }

  state <- new.env(parent = emptyenv())
  state$static <- if (inherits(rast, "RasterLayer")) rast else NULL
  state$layers <- if (is.list(rast) && !inherits(rast, "RasterLayer")) rast else list()
  state$info   <- NULL

  fetch_day <- function(day, lon, lat) {
    if (!requireNamespace("rerddap", quietly = TRUE)) {
      stop("sst_source() needs the 'rerddap' package. ",
           "Install it with install.packages('rerddap'), or pass your own data via rast=.")
    }
    if (is.null(state$info)) state$info <- rerddap::info(dataset, url = url)
    g <- rerddap::griddap(
      state$info,
      time      = c(day, day),
      latitude  = c(min(lat) - pad, max(lat) + pad),
      longitude = c(min(lon) - pad, max(lon) + pad),
      fields    = field,
      fmt       = "csv"
    )
    df <- g$data[, c("longitude", "latitude", field)]
    df <- df[stats::complete.cases(df), ]
    raster::rasterFromXYZ(df, crs = "+proj=longlat +datum=WGS84")
  }

  src <- function(lon, lat, date) {
    if (!is.null(state$static)) {
      return(as.numeric(raster::extract(state$static, cbind(lon, lat))))
    }
    key <- as.character(as.Date(date))
    if (is.null(state$layers[[key]])) {
      state$layers[[key]] <- fetch_day(key, lon, lat)
    }
    as.numeric(raster::extract(state$layers[[key]], cbind(lon, lat)))
  }
  class(src) <- c("tf_source", "function")
  src
}


#' Raster source: wrap a user-supplied raster or named list of rasters
#'
#' A thin adapter that turns a \code{RasterLayer} (static, date-independent) or
#' a named list of \code{RasterLayer}s (keyed by \code{"YYYY-MM-DD"}) into a
#' `tf_source` closure for use in [location_term()]. This is the generic
#' "bring your own raster" entry point; [sst_source()] and [bathy_source()] are
#' specialised versions that also handle fetching.
#'
#' @param r A \code{RasterLayer} (used for every date) or a non-empty named list
#'   of \code{RasterLayer}s keyed by \code{"YYYY-MM-DD"} strings.
#'
#' @return A source closure of class `tf_source`: `function(lon, lat, date)`
#'   returning a numeric vector of field values at the candidate cells.
#'
#' @seealso [function_source()], [prior_raster()], [sst_source()]
#' @examples
#' r <- raster::raster(nrows = 1, ncols = 2,
#'                     xmn = 0, xmx = 2, ymn = 0, ymx = 1)
#' raster::values(r) <- c(5.0, 10.0)
#' rs <- raster_source(r)
#' rs(lon = c(0.5, 1.5), lat = c(0.5, 0.5), date = Sys.Date())
#' @export
raster_source <- function(r) {
  is_layer <- inherits(r, "RasterLayer")
  is_list  <- is.list(r) && length(r) > 0L && !is_layer
  if (!is_layer && !is_list) {
    stop("`r` must be a RasterLayer or a non-empty named list of RasterLayers")
  }
  src <- function(lon, lat, date = NULL) {
    if (is_layer) {
      return(as.numeric(raster::extract(r, cbind(lon, lat))))
    }
    key <- as.character(as.Date(date))
    layer <- r[[key]]
    if (is.null(layer)) stop("no raster layer for date ", key)
    as.numeric(raster::extract(layer, cbind(lon, lat)))
  }
  class(src) <- c("tf_source", "function")
  src
}


#' Function source: wrap a closure as a location source
#'
#' Promotes any function `fn(lon, lat, date)` to a `tf_source` for use in
#' [location_term()]. This is the escape hatch for arbitrary environmental
#' fields (geomagnetic models, tidal models, custom climatologies) that you can
#' evaluate programmatically at any coordinate and date.
#'
#' @param fn A function accepting `(lon, lat, date)` — equal-length numeric
#'   vectors of coordinates and a single date (or anything the function
#'   accepts) — and returning a numeric vector of field values the same length
#'   as `lon`.
#'
#' @return A source closure of class `tf_source`.
#'
#' @seealso [raster_source()], [custom_rule()]
#' @examples
#' # A synthetic "field" returning latitude as the value (useful for testing).
#' fs <- function_source(function(lon, lat, date) lat)
#' fs(lon = c(10, 20), lat = c(-30, 30), date = Sys.Date())
#' @export
function_source <- function(fn) {
  if (!is.function(fn)) stop("`fn` must be a function(lon, lat, date)")
  src <- function(lon, lat, date) fn(lon, lat, date)
  class(src) <- c("tf_source", "function")
  src
}


#' Prior raster source: a spatial log-prior (or density) surface
#'
#' A convenience wrapper around [raster_source()] for user-supplied prior
#' surfaces, documenting the expected convention. Use with [identity_rule()].
#' If `is_log = FALSE` the values are treated as non-negative densities and
#' logged before being passed to the rule (zeros become `-Inf`).
#'
#' @param r A \code{RasterLayer} whose values are log-prior densities (if
#'   `is_log = TRUE`, the default) or non-negative densities.
#' @param is_log If `TRUE` (default) values are already on the log scale. If
#'   `FALSE` they are logged internally.
#'
#' @return A source closure of class `tf_source`. Pair with [identity_rule()].
#'
#' @seealso [identity_rule()], [hemisphere_prior()]
#' @examples
#' r <- raster::raster(nrows = 1, ncols = 2,
#'                     xmn = 0, xmx = 2, ymn = 0, ymx = 1)
#' raster::values(r) <- c(log(0.8), log(0.2))   # log-prior
#' ps <- prior_raster(r)
#' ps(lon = c(0.5, 1.5), lat = c(0.5, 0.5), date = NULL)
#' @export
prior_raster <- function(r, is_log = TRUE) {
  if (!inherits(r, "RasterLayer")) stop("`r` must be a RasterLayer")
  src <- function(lon, lat, date = NULL) {
    vals <- as.numeric(raster::extract(r, cbind(lon, lat)))
    if (!isTRUE(is_log)) {
      vals <- ifelse(is.na(vals) | vals <= 0, -Inf, log(vals))
    }
    vals
  }
  class(src) <- c("tf_source", "function")
  src
}


#' Sea mask source: land exclusion from Natural Earth polygons
#'
#' Rasterises `rnaturalearth::ne_countries()` onto the grid and returns a static
#' source that is `1` over sea cells and `0` over land cells. Pair with
#' [mask_rule()] to hard-exclude land candidate locations in [TwilightFreeGrid()].
#'
#' The rasterisation is performed once on the first call (using the supplied
#' coordinate bounding box) and cached for the life of the source object.
#' At coarse resolutions, small islands may not appear. For fine resolution
#' grids, the initial rasterisation may take a few seconds.
#'
#' @return A source closure of class `tf_source`: `function(lon, lat, date)`
#'   returning a numeric vector (`1` = sea, `0` = land). `date` is ignored.
#'   Pair with [mask_rule()].
#'
#' @seealso [mask_rule()]
#' @examples
#' \dontrun{
#' sm <- sea_mask_source()
#' val <- sm(lon = c(151, 148), lat = c(-34, -33))  # Sydney offshore vs land
#' }
#' @export
sea_mask_source <- function() {
  state <- new.env(parent = emptyenv())
  state$rast <- NULL

  build <- function(lon, lat) {
    r <- raster::raster(
      xmn = min(lon) - 1, xmx = max(lon) + 1,
      ymn = min(lat) - 1, ymx = max(lat) + 1,
      resolution = max(diff(range(lon)), diff(range(lat))) / 200,
      crs = "+proj=longlat +datum=WGS84"
    )
    land <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sp")
    land_r <- raster::rasterize(land, r, field = 1, background = 0)
    sea_r <- land_r
    sea_r[] <- ifelse(land_r[] == 1, 0, 1)
    state$rast <- sea_r
    invisible(sea_r)
  }

  src <- function(lon, lat, date = NULL) {
    if (is.null(state$rast)) build(lon, lat)
    as.numeric(raster::extract(state$rast, cbind(lon, lat)))
  }
  class(src) <- c("tf_source", "function")
  src
}
