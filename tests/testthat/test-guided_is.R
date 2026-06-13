# Verifies that the guided particle filter IS correction produces estimates
# consistent with the unguided forward filter on the same light data.
#
# Setup: noiseless observations from a known stationary location (lon=147, lat=-43).
# We run (a) unguided forward SMC and (b) guided SMC with the TRUE endpoint supplied.
# After the IS correction, the guided MAP should not be pulled off the true location;
# both methods should agree within ~3 degrees.

make_guided_smc <- function(method, seed = 42) {
  times <- seq(
    as.POSIXct("2020-06-01", tz = "UTC"),
    as.POSIXct("2020-06-04", tz = "UTC"),
    by = "10 min"
  )
  t_unix <- as.numeric(times)
  true_lon <- 147.0; true_lat <- -43.0
  cal <- c(64, 64 / 90)

  zenith <- solar_zenith(t_unix, rep(true_lon, length(t_unix)), rep(true_lat, length(t_unix)))
  obs    <- pmax(0, pmin(64, cal[1] - cal[2] * zenith))

  TwilightFreeSMC(
    date_time    = times,
    light        = obs,
    start_lat    = true_lat,
    start_lon    = true_lon,
    end_lat      = if (method == "guided") true_lat else NA_real_,
    end_lon      = if (method == "guided") true_lon else NA_real_,
    n_particles  = 500,
    step_hours   = 12,
    diffusion    = 100,
    calibration  = cal,
    method       = method,
    seed         = seed
  )
}

test_that("guided IS correction: MAP is consistent with unguided forward filter", {
  fwd  <- make_guided_smc("forward")
  guid <- make_guided_smc("guided")

  # Both should recover approximately the true location
  true_lon <- 147.0; true_lat <- -43.0
  fwd_lon_err  <- abs(mean(fwd$lon)  - true_lon)
  fwd_lat_err  <- abs(mean(fwd$lat)  - true_lat)
  guid_lon_err <- abs(mean(guid$lon) - true_lon)
  guid_lat_err <- abs(mean(guid$lat) - true_lat)

  expect_lt(fwd_lon_err,  5, label = "forward SMC lon error < 5 deg")
  expect_lt(fwd_lat_err,  5, label = "forward SMC lat error < 5 deg")
  expect_lt(guid_lon_err, 5, label = "guided SMC lon error < 5 deg")
  expect_lt(guid_lat_err, 5, label = "guided SMC lat error < 5 deg")

  # After IS correction the guided MAP must not be pulled significantly further
  # from the truth than the unguided MAP (was biased toward endpoint before fix;
  # endpoint IS the true location here, so both should be similar).
  expect_lt(
    abs(guid_lat_err - fwd_lat_err), 3,
    label = "guided and forward lat MAPs agree within 3 deg after IS correction"
  )
})

test_that("guided IS correction: last-step fallback does not panic", {
  # A 2-step track (start + 1 knot + end) triggers the degenerate last-step
  # case where sigma_g -> 0. The do_guided guard should fall back to unguided
  # without a panic.
  times <- seq(
    as.POSIXct("2020-06-01", tz = "UTC"),
    as.POSIXct("2020-06-02", tz = "UTC"),
    by = "10 min"
  )
  t_unix <- as.numeric(times)
  cal    <- c(64, 64 / 90)
  zenith <- solar_zenith(t_unix, rep(147, length(t_unix)), rep(-43, length(t_unix)))
  obs    <- pmax(0, pmin(64, cal[1] - cal[2] * zenith))

  expect_no_error(
    TwilightFreeSMC(
      date_time   = times,
      light       = obs,
      start_lat   = -43, start_lon = 147,
      end_lat     = -43, end_lon   = 147,
      n_particles = 100,
      step_hours  = 12,
      diffusion   = 100,
      calibration = cal,
      method      = "guided",
      seed        = 1
    )
  )
})
