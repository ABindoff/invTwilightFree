# Tests that TwilightFreeGrid respects NA cells produced by makeGrid(mask=).
# Uses a manually constructed raster so no rnaturalearthdata download is needed.

test_that("TwilightFreeGrid excludes NA (masked) grid cells from HMM candidate set", {
  skip_if_not_installed("raster")

  # 1-row, 4-column raster: cells at lon = -15, -5, +5, +15 (10-degree spacing).
  # Middle two cells (lon = -5, +5) are valid (value = 1); outer two are masked (NA).
  # After the fix, lon_vec/lat_vec passed to Rust should contain only {-5, +5}.
  r <- raster::raster(
    matrix(c(NA, 1, 1, NA), nrow = 1),
    xmn = -20, xmx = 20, ymn = 27, ymx = 33
  )

  times  <- seq(as.POSIXct("2024-06-01", tz = "UTC"),
                as.POSIXct("2024-06-02", tz = "UTC"),
                by = "10 min")
  t_unix <- as.numeric(times)
  cal    <- c(64, 64 / 90)
  lp     <- c(0.5, 64, 0.05)

  # Noiseless obs from lon=+5 (inside the valid zone)
  zenith <- solar_zenith(t_unix, rep(5, length(t_unix)), rep(30, length(t_unix)))
  obs    <- pmax(0, pmin(64, cal[1] - cal[2] * zenith))

  fit <- TwilightFreeGrid(
    date_time         = times,
    light             = obs,
    grid              = r,
    calibration       = cal,
    likelihood_params = lp,
    step_hours        = 12,
    diffusion         = 200
  )

  # Valid cell centres: lon = -5 and lon = +5.  Masked cell centres: -15 and +15.
  # Every MAP longitude must fall within the valid zone, not at a masked cell.
  expect_true(
    all(fit$fit$lon > -10 & fit$fit$lon < 10),
    label = paste("MAP lon should be in valid zone (-10, 10); got:",
                  paste(round(fit$fit$lon, 1), collapse = " "))
  )
})

test_that("unmasked grid (all cells valid) is unaffected by the NA filter", {
  skip_if_not_installed("raster")

  # All four cells have value 1 — no NAs — so behaviour must be identical to
  # the pre-fix code path.
  r_full <- raster::raster(
    matrix(c(1, 1, 1, 1), nrow = 1),
    xmn = -20, xmx = 20, ymn = 27, ymx = 33
  )
  r_na <- raster::raster(
    matrix(c(NA, 1, 1, NA), nrow = 1),
    xmn = -20, xmx = 20, ymn = 27, ymx = 33
  )

  times  <- seq(as.POSIXct("2024-06-01", tz = "UTC"),
                as.POSIXct("2024-06-02", tz = "UTC"),
                by = "10 min")
  t_unix <- as.numeric(times)
  cal    <- c(64, 64 / 90); lp <- c(0.5, 64, 0.05)
  zenith <- solar_zenith(t_unix, rep(5, length(t_unix)), rep(30, length(t_unix)))
  obs    <- pmax(0, pmin(64, cal[1] - cal[2] * zenith))

  fit_full <- TwilightFreeGrid(date_time = times, light = obs, grid = r_full,
    calibration = cal, likelihood_params = lp, step_hours = 12, diffusion = 200)
  fit_na   <- TwilightFreeGrid(date_time = times, light = obs, grid = r_na,
    calibration = cal, likelihood_params = lp, step_hours = 12, diffusion = 200)

  # With all 4 cells valid, obs from lon=+5 picks the +5 cell in both cases.
  # With the NA mask, the HMM can only pick from {-5, +5} so the result
  # should still be +5. Both fits should agree on the MAP longitude.
  expect_equal(fit_full$fit$lon, fit_na$fit$lon)
})
