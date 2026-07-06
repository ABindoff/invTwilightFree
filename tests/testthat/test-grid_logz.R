# Tests that run_grid_hmm returns a usable log marginal likelihood (log_z) and
# that it behaves as model evidence for hemisphere class selection.

test_that("TwilightFreeGrid returns a finite scalar log_z", {
  skip_if_not_installed("raster")
  r <- makeGrid(lon = c(140, 160), lat = c(-60, -40), cell.size = 2)
  times <- seq(as.POSIXct("2024-06-01", tz = "UTC"),
               as.POSIXct("2024-06-08", tz = "UTC"), by = "20 min")
  t_unix <- as.numeric(times)
  cal <- c(64, 64 / 90); lp <- c(0.5, 64, 0.05)
  zen <- solar_zenith(t_unix, rep(150, length(t_unix)), rep(-50, length(t_unix)))
  obs <- pmax(0, pmin(64, cal[1] - cal[2] * zen))

  fit <- TwilightFreeGrid(times, obs, r, calibration = cal, likelihood_params = lp,
                          step_hours = 12, diffusion = 150)
  expect_true(is.numeric(fit$log_z))
  expect_length(fit$log_z, 1L)
  expect_true(is.finite(fit$log_z))
})

test_that("log_z prefers the hemisphere that generated the data (solstice)", {
  skip_if_not_installed("raster")
  # Symmetric grid spanning both hemispheres; data generated from the south.
  r <- makeGrid(lon = c(140, 160), lat = c(-60, 60), cell.size = 2)
  times <- seq(as.POSIXct("2024-06-15", tz = "UTC"),       # near solstice: strong signal
               as.POSIXct("2024-06-29", tz = "UTC"), by = "20 min")
  t_unix <- as.numeric(times)
  cal <- c(64, 64 / 90); lp <- c(0.5, 64, 0.05)
  zen <- solar_zenith(t_unix, rep(150, length(t_unix)), rep(-50, length(t_unix)))
  obs <- pmax(0, pmin(64, cal[1] - cal[2] * zen))

  south <- TwilightFreeGrid(times, obs, r, calibration = cal, likelihood_params = lp,
    step_hours = 12, diffusion = 150,
    terms = list(location_term("S", source = hemisphere_prior(function(d) "S", softness = 0),
                               rule = identity_rule())))
  north <- TwilightFreeGrid(times, obs, r, calibration = cal, likelihood_params = lp,
    step_hours = 12, diffusion = 150,
    terms = list(location_term("N", source = hemisphere_prior(function(d) "N", softness = 0),
                               rule = identity_rule())))

  expect_gt(south$log_z, north$log_z)            # evidence favours the true branch
})
