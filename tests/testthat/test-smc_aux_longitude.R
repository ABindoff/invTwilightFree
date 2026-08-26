# Regression: sensor-fusion terms must keep their longitude dependence in the
# SMC engine whatever longitude convention the auxiliary grid uses.
#
# The defect this pins: particle longitudes are held in -180..180, while the aux
# raster carries the caller's convention. An analysis working in 0..360 supplies
# an extent like 150..250, so `(p_lon - xmin) / cell_w` was negative for every
# particle and `negative f64 as usize` saturates to 0 in Rust. Every particle
# read column 0 and the terms silently lost all longitude dependence -- with no
# warning, and with latitude still working, so the feature looked alive.
#
# Found while setting up the elephant-seal reporting run, where two fusion arms
# came back bit-identical to their light-only counterparts.

make_track <- function(n_days = 12, lon0 = 200, lat0 = 40) {
  times <- seq(as.POSIXct("2021-06-01", tz = "UTC"), by = "30 min",
               length.out = n_days * 48)
  z <- solar_zenith(as.numeric(times), rep(lon0, length(times)), rep(lat0, length(times)))
  list(time = times, light = pmax(0, pmin(64, 558.5 - 5.818 * z)),
       lon = lon0, lat = lat0)
}

fit_with <- function(trk, terms, mask) {
  TwilightFreeSMC(date_time = trk$time, light = trk$light,
                  calibration = c(558.5, 5.818),
                  likelihood_params = c(1 / 32, 64, 0.1),
                  n_particles = 400L,
                  start_lon = trk$lon, start_lat = trk$lat,
                  end_lon = trk$lon, end_lat = trk$lat,
                  method = "ffbs", step_hours = 12, diffusion = 110,
                  spatial_mask = mask, terms = terms, seed = 42)
}

flat_grid <- function(xmin, xmax, ymin = 20, ymax = 70) {
  m <- terra::rast(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                   resolution = 1, crs = "EPSG:4326")
  terra::values(m) <- 1
  m
}

# A term that depends on longitude ALONE, pulling hard toward one meridian.
lon_pull <- function(target) {
  location_term(
    name = "lon_pull",
    source = function_source(function(lon, lat, date = NULL) lon),
    rule = custom_rule(function(obs, expected) {
      -0.5 * ((as.numeric(expected) - target) / 3)^2
    }))
}

test_that("a longitude-only term moves the SMC fit in the 0-360 convention", {
  skip_if_not_installed("terra")
  trk  <- make_track()
  mask <- flat_grid(150, 250)                     # 0-360, as a Pacific study uses
  base <- fit_with(trk, list(), mask)
  pull <- fit_with(trk, list(lon_pull(215)), mask)
  # The engine reports longitude in -180..180; compare on a common convention.
  lon360 <- function(x) (x %% 360 + 360) %% 360
  expect_false(identical(base$lon, pull$lon))
  expect_gt(max(abs(lon360(pull$lon) - lon360(base$lon))), 0.5)
  # and it must move TOWARD the target, not merely differ
  expect_lt(abs(mean(lon360(pull$lon)) - 215), abs(mean(lon360(base$lon)) - 215))
})

test_that("the same term still works in the -180..180 convention", {
  skip_if_not_installed("terra")
  trk  <- make_track(lon0 = -160)
  mask <- flat_grid(-180, -100)
  base <- fit_with(trk, list(), mask)
  pull <- fit_with(trk, list(lon_pull(-145)), mask)
  expect_false(identical(base$lon, pull$lon))
  expect_lt(abs(mean(pull$lon) - (-145)), abs(mean(base$lon) - (-145)))
})

test_that("a latitude-only term still works (it always did; guard the fix)", {
  skip_if_not_installed("terra")
  trk  <- make_track()
  mask <- flat_grid(150, 250)
  lat_pull <- location_term(
    name = "lat_pull",
    source = function_source(function(lon, lat, date = NULL) lat),
    rule = custom_rule(function(obs, expected) -0.5 * ((as.numeric(expected) - 52) / 3)^2))
  base <- fit_with(trk, list(), mask)
  pull <- fit_with(trk, list(lat_pull), mask)
  expect_false(identical(base$lat, pull$lat))
  expect_lt(abs(mean(pull$lat) - 52), abs(mean(base$lat) - 52))
})

test_that("no terms is bit-identical before and after the wrap", {
  skip_if_not_installed("terra")
  trk  <- make_track()
  mask <- flat_grid(150, 250)
  a <- fit_with(trk, list(), mask)
  b <- fit_with(trk, list(), mask)
  expect_identical(a$lon, b$lon)
  expect_identical(a$lat, b$lat)
})
