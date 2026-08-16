make_grid <- function() {
  g <- terra::rast(xmin = 140, xmax = 160, ymin = -50, ymax = -35,
                   resolution = 1, crs = "EPSG:4326")
  terra::values(g) <- 1
  g
}

sim_light <- function(lon = 150, lat = -45, days = 12) {
  t <- seq(as.POSIXct("2021-06-01", tz = "UTC"), by = "20 min",
           length.out = days * 72)
  z <- solar_zenith(as.numeric(t), rep(lon, length(t)), rep(lat, length(t)))
  list(time = t, light = pmin(pmax(558.5 - 5.818 * z, 0), 64))
}

test_that("the default is bit-identical to drift_correction = FALSE", {
  # `drift_correction` defaults to FALSE precisely so that every number in the
  # elephant-seal campaign stays reproducible. Checked on log_z as well as the
  # track: the shift was added to the forward and backward passes separately, so
  # an inconsistency between them would surface in the marginal likelihood first.
  skip_if_not_installed("terra")
  d <- sim_light(); g <- make_grid()
  a <- TwilightFreeGrid(d$time, d$light, grid = g, start_lon = 150, start_lat = -45,
                        step_hours = 12, diffusion = 80)
  b <- TwilightFreeGrid(d$time, d$light, grid = g, start_lon = 150, start_lat = -45,
                        step_hours = 12, diffusion = 80, drift_correction = FALSE)
  expect_identical(a$fit$lon, b$fit$lon)
  expect_identical(a$fit$lat, b$fit$lat)
  expect_identical(a$log_z, b$log_z)
})

test_that("drift_correction = TRUE reaches the engine but barely moves the fit", {
  # This guards the NEGATIVE RESULT, in both directions.
  #
  # It must reach the engine: `log_z` has to move, otherwise the argument is not
  # plumbed through and a future reader would think the question had been settled
  # when it had not. An earlier version of this work was nearly derailed by
  # exactly that -- a coarse check showed identical MODES and looked inert.
  #
  # And it must barely move the fit: the correction nulls the one-step kernel
  # drift exactly, yet on a 476-knot record the posterior mean latitude shifted
  # 0.00009 deg. The reason is that a forward-backward smoother cancels a drift
  # entering from both time directions. If a future change makes this argument
  # suddenly matter, the cancellation argument in
  # notes/latitude_bias_investigation.md section 4a needs revisiting.
  skip_if_not_installed("terra")
  d <- sim_light(); g <- make_grid()
  off <- TwilightFreeGrid(d$time, d$light, grid = g, start_lon = 150, start_lat = -45,
                          step_hours = 12, diffusion = 80, drift_correction = FALSE)
  on  <- TwilightFreeGrid(d$time, d$light, grid = g, start_lon = 150, start_lat = -45,
                          step_hours = 12, diffusion = 80, drift_correction = TRUE)
  expect_false(identical(off$log_z, on$log_z))
  expect_lt(abs(mean(on$fit$lat) - mean(off$fit$lat)), 0.05)
})
