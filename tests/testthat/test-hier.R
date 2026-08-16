test_that("TwilightFreeHier runs end-to-end and returns a well-formed object", {
  skip_on_cran()

  mk_tag <- function(days, lat0, lon0, seed) {
    set.seed(seed)
    times <- seq(as.POSIXct("2024-05-01", tz = "UTC"), by = "10 min", length.out = days * 144)
    lat <- lat0 + cumsum(stats::rnorm(length(times), 0, 0.02))
    lon <- lon0 + cumsum(stats::rnorm(length(times), 0, 0.02))
    z <- solar_zenith(as.numeric(times), lon, lat)
    light <- pmin(pmax(558.5 - 5.818 * z, 0), 64)
    list(df = data.frame(Date = times, Light = light), lat = lat, lon = lon)
  }
  a <- mk_tag(5, -48, 150, 1)
  b <- mk_tag(7, -46, 148, 2)          # different length
  data <- list(A = a$df, B = b$df)
  locations <- data.frame(
    id = c("A", "B"),
    deploy_lon = c(a$lon[1], b$lon[1]), deploy_lat = c(a$lat[1], b$lat[1]),
    retrieve_lon = c(a$lon[length(a$lon)], NA),   # A fixed, B free
    retrieve_lat = c(a$lat[length(a$lat)], NA))

  fit <- TwilightFreeHier(data, locations, step_hours = 12,
                          mesh_pad = 8, coarse_res = 2, surrogate_diffusion = 60,
                          sweeps = 200L, burn = 50L, thin = 2L, seed = 1L)

  expect_s3_class(fit, "TwilightFreeHier")
  expect_equal(nrow(fit$movement), 2L)
  expect_equal(fit$movement$id, c("A", "B"))
  expect_true(all(fit$movement$sigma_kmday > 0))
  expect_named(fit$tracks, c("A", "B"))
  expect_true(all(is.finite(fit$tracks$A$lon)))

  # deployment is fixed: first knot has zero posterior sd and equals the input
  expect_equal(fit$tracks$A$lon_sd[1], 0)
  expect_equal(fit$tracks$A$lon[1], locations$deploy_lon[1], tolerance = 1e-6)
  # A has a fixed retrieval (last knot pinned, ~zero sd); B is free (positive sd)
  expect_lt(fit$tracks$A$lat_sd[nrow(fit$tracks$A)], 1e-6)
  expect_gt(fit$tracks$B$lat_sd[nrow(fit$tracks$B)], 0)
})

test_that("TwilightFreeHier accepts a long data frame with an id column", {
  skip_on_cran()
  set.seed(3)
  one <- function(id, lat0, lon0) {
    times <- seq(as.POSIXct("2024-05-01", tz = "UTC"), by = "10 min", length.out = 5 * 144)
    lat <- lat0 + cumsum(stats::rnorm(length(times), 0, 0.02))
    lon <- lon0 + cumsum(stats::rnorm(length(times), 0, 0.02))
    z <- solar_zenith(as.numeric(times), lon, lat)
    data.frame(id = id, Date = times, Light = pmin(pmax(558.5 - 5.818 * z, 0), 64),
               dlon = lon[1], dlat = lat[1])
  }
  long <- rbind(one("X", -48, 150), one("Y", -45, 149))
  locs <- data.frame(id = c("X", "Y"),
                     deploy_lon = c(long$dlon[1], long$dlon[long$id == "Y"][1]),
                     deploy_lat = c(long$dlat[1], long$dlat[long$id == "Y"][1]))
  fit <- TwilightFreeHier(long[, c("id", "Date", "Light")], locs, id = "id",
                          mesh_pad = 8, coarse_res = 2, sweeps = 150L, burn = 50L, thin = 2L, seed = 2L)
  expect_s3_class(fit, "TwilightFreeHier")
  expect_equal(fit$movement$id, c("X", "Y"))
})

test_that("a hemisphere_prior term shifts the latitude posterior in its direction", {
  skip_on_cran()
  mk <- function(seed, lon0) {
    set.seed(seed)
    times <- seq(as.POSIXct("2024-09-20", tz = "UTC"), by = "10 min", length.out = 8 * 144)
    lat <- cumsum(stats::rnorm(length(times), 0, 0.03))
    lon <- lon0 + cumsum(stats::rnorm(length(times), 0, 0.03))
    z <- solar_zenith(as.numeric(times), lon, lat)
    data.frame(Date = times, Light = pmin(pmax(558.5 - 5.818 * z, 0), 64))
  }
  data <- list(T1 = mk(1, -30))
  locs <- data.frame(id = "T1", deploy_lon = -30, deploy_lat = 0,
                     retrieve_lon = NA, retrieve_lat = NA)
  hemi <- function(side) location_term("hemi", source = hemisphere_prior(function(d) side),
                                       rule = identity_rule())
  common <- list(data = data, locations = locs, step_hours = 12, mesh_pad = 25, coarse_res = 3,
                 surrogate_diffusion = 80, sweeps = 1500L, burn = 500L, thin = 3L, seed = 11L)
  mlat <- function(fit) mean(fit$tracks$T1$lat)
  lat_N <- mlat(do.call(TwilightFreeHier, c(common, list(terms = hemi("N")))))
  lat_S <- mlat(do.call(TwilightFreeHier, c(common, list(terms = hemi("S")))))
  expect_gt(lat_N, lat_S)     # 'N' prior pulls north of the 'S' prior
})

test_that("the spherical metric carries the area element and points equatorward", {
  # The block sampler's target is a density w.r.t. d(lon) d(lat), so it owes the
  # sphere's cos(lat) area element. Omitting it is the grid engine's uncorrected
  # kernel in continuous form: mass accumulates as 1/cos(lat) and tracks lean
  # poleward. The factor does not cancel -- the proposal is Gaussian in degrees.
  #
  # No stored pre-fix reference exists, so this uses a contrast that needs none.
  # Mute the emission (constant light carries no position) and free the retrieval
  # end, so the movement prior alone sets latitude, then compare metrics:
  #
  #   flat       symmetric RW in degrees under a FIXED reference metric. It owes
  #              no Jacobian (the degrees->km map is affine, so its constant
  #              Jacobian cancels) and has no latitude direction. The control.
  #   spherical  the same walk plus the great-circle upgrade and the area element,
  #              which is the only term here with a latitude DIRECTION.
  #
  # So spherical must sit EQUATORWARD of flat, and the sign must FLIP across the
  # equator. Both hemispheres are run because a same-signed shift at both would be
  # some generic pull rather than the measure. The sign is the assertion; the
  # magnitude is not (an uninformative posterior spreads to the mesh edge).
  skip_on_cran()
  mk_flat <- function() {
    tm <- seq(as.POSIXct("2024-03-18", tz = "UTC"), by = "30 min",
              length.out = 30 * 48)
    data.frame(time = tm, light = rep(20, length(tm)))
  }
  drift <- function(dep_lat, metric) {
    loc <- data.frame(id = "a", deploy_lon = 150, deploy_lat = dep_lat,
                      retrieve_lon = NA, retrieve_lat = NA)
    f <- TwilightFreeHier(data = list(a = mk_flat()), locations = loc,
                          step_hours = 12, calibration = c(64, 64 / 90),
                          likelihood_params = c(0.5, 64, 0.05), metric = metric,
                          surrogate_diffusion = 200, mesh_pad = c(30, 30),
                          sweeps = 1200L, burn = 400L, thin = 4L, seed = 7L)
    mean(f$tracks[[1]]$lat) - dep_lat
  }
  expect_lt(drift( 45, "spherical"), drift( 45, "flat"))   # north: pulled south
  expect_gt(drift(-45, "spherical"), drift(-45, "flat"))   # south: pulled north
})
