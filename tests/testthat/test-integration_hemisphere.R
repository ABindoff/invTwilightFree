# End-to-end integration test: hemisphere prior fixes equinox hemisphere confusion.
#
# Physical setup:
#   A southern-hemisphere species has a GLS light record generated from near the
#   September equinox.  The solar zenith at +25 N on Sept 22-24 is nearly
#   identical to the zenith at -25 S (sun declination ~0; max diff ~0.35 deg),
#   so the light likelihood is bimodal with peaks at approximately +-25 deg lat.
#   The light data here were generated at the NORTHERN mirror (+25 N), so the
#   light-only MAP is in the wrong hemisphere.  A hemisphere_prior("S") adds
#   ~-6.9 nats to every northern cell and flips the MAP to the correct southern
#   hemisphere.
#
# Regression value: if the bimodal peak shifts sides across package versions,
# this test will catch it, prompting investigation rather than silent drift.

test_that("hemisphere_prior flips equinox MAP from northern to southern hemisphere", {
  skip_if_not_installed("raster")

  # --- Synthetic light data generated at the NORTHERN mirror (+25 N, -30 E) ---
  # With set.seed, the obs vector is identical every run; TwilightFreeGrid is
  # deterministic (Grid HMM), so the full test is reproducible with no seed arg.
  times <- seq(
    as.POSIXct("2024-09-22", tz = "UTC"),
    as.POSIXct("2024-09-24", tz = "UTC"),
    by = "10 min"
  )
  t_unix <- as.numeric(times)
  cal    <- c(64, 64 / 90)          # intercept, slope; 64 units full light
  lp     <- c(0.5, 64, 0.05)        # lambda, max_light, prob_slab

  set.seed(42)
  zenith_n <- solar_zenith(t_unix, rep(-30, length(t_unix)), rep(25, length(t_unix)))
  obs <- pmax(0, pmin(64, cal[1] - cal[2] * zenith_n + rnorm(length(t_unix), 0, 1)))

  # Symmetric grid: equal extent in N and S hemispheres
  grid <- makeGrid(lon = c(-35, -25), lat = c(-35, 35), cell.size = 2.5)

  # ---- Light only (no prior, no fixed endpoints) ----
  fit_lo <- TwilightFreeGrid(
    date_time         = times,
    light             = obs,
    grid              = grid,
    calibration       = cal,
    likelihood_params = lp,
    step_hours        = 12,
    diffusion         = 200
  )
  lat_lo <- fit_lo$fit$lat

  # ---- With southern-hemisphere prior ----
  term_s <- location_term(
    "hemisphere",
    source = hemisphere_prior(function(d) "S", softness = 1e-3),
    rule   = identity_rule()
  )
  fit_pr <- TwilightFreeGrid(
    date_time         = times,
    light             = obs,
    grid              = grid,
    calibration       = cal,
    likelihood_params = lp,
    step_hours        = 12,
    diffusion         = 200,
    terms             = list(term_s)
  )
  lat_pr <- fit_pr$fit$lat

  # Light-only MAP is in the wrong (northern) hemisphere
  expect_true(all(lat_lo > 0),
    label = paste("light-only lat should all be >0; got", paste(round(lat_lo, 1), collapse = " ")))

  # Prior flips every knot to the correct (southern) hemisphere
  expect_true(all(lat_pr <= 0),
    label = paste("prior lat should all be <=0; got", paste(round(lat_pr, 1), collapse = " ")))

  # The shift should be substantial (at least 30 degrees, not just numerical noise)
  expect_gt(mean(lat_lo) - mean(lat_pr), 30)
})


test_that("terms = list() is byte-identical to a no-terms call", {
  skip_if_not_installed("raster")

  times <- seq(
    as.POSIXct("2024-09-22", tz = "UTC"),
    as.POSIXct("2024-09-24", tz = "UTC"),
    by = "10 min"
  )
  t_unix <- as.numeric(times)
  cal <- c(64, 64 / 90); lp <- c(0.5, 64, 0.05)

  set.seed(42)
  zenith_n <- solar_zenith(t_unix, rep(-30, length(t_unix)), rep(25, length(t_unix)))
  obs <- pmax(0, pmin(64, cal[1] - cal[2] * zenith_n + rnorm(length(t_unix), 0, 1)))

  grid <- makeGrid(lon = c(-35, -25), lat = c(-35, 35), cell.size = 2.5)

  fit_base <- TwilightFreeGrid(
    date_time = times, light = obs, grid = grid,
    calibration = cal, likelihood_params = lp, step_hours = 12, diffusion = 200
  )
  fit_empty <- TwilightFreeGrid(
    date_time = times, light = obs, grid = grid,
    calibration = cal, likelihood_params = lp, step_hours = 12, diffusion = 200,
    terms = list()
  )
  expect_identical(fit_base$fit$lat, fit_empty$fit$lat)
  expect_identical(fit_base$fit$lon, fit_empty$fit$lon)
})
