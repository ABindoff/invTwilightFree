# Tests for the dual (time-conditioned) primitives. Use the compiled solar model,
# so they run after devtools::load_all().

cal <- c(64, 64 / 90)
lp  <- c(0.5, 64, 0.05)

# noiseless light generated from a known location/time
gen_light <- function(lon, lat, times) {
  z <- solar_zenith(as.numeric(times), rep(lon, length(times)), rep(lat, length(times)))
  pmin(pmax(cal[1] - cal[2] * z, 0), lp[2])
}

test_that("eval_logpt_loc peaks at zero offset for the true location and clock", {
  times <- seq(as.POSIXct("2024-05-01", tz = "UTC"),
               as.POSIXct("2024-05-04", tz = "UTC"), by = "10 min")
  light <- gen_light(150, -50, times)
  d <- eval_logpt_loc(150, -50, times, light, calibration = cal, likelihood_params = lp)
  peak <- d$offset[which.max(d$loglik)]
  expect_lt(abs(peak), 240)                          # within a couple of steps of 0
})

test_that("eval_logpt_loc recovers an injected clock offset", {
  times <- seq(as.POSIXct("2024-05-01", tz = "UTC"),
               as.POSIXct("2024-05-04", tz = "UTC"), by = "10 min")
  light <- gen_light(150, -50, times)
  delta <- 25 * 60                                   # tag clock 25 min fast
  wrong_times <- times + delta                       # timestamps carry the error
  d <- eval_logpt_loc(150, -50, wrong_times, light, calibration = cal, likelihood_params = lp)
  peak <- d$offset[which.max(d$loglik)]
  expect_lt(abs(peak - (-delta)), 180)               # correction ~ -delta
})

# NOTE: an earlier test here asserted the dual is "sharper in time than the
# spatial likelihood is in latitude". It was removed deliberately: the
# conditioning analysis (notes/topology/conditioning_finding.md) showed that, for
# the continuous-amplitude model, latitude is well determined by midday amplitude
# and is NOT reliably the softer direction, so the comparison is not a true
# invariant. (It was also degenerate on noiseless data, where both widths
# collapse to one grid step.) The dual's sharpness in time is instead checked
# directly by the "peaks at zero offset" and "recovers an injected offset" tests.

test_that("calibrate_clock recovers a linear drift and zeroes the interior", {
  times <- seq(as.POSIXct("2024-05-01", tz = "UTC"),
               as.POSIXct("2024-08-29", tz = "UTC"), by = "20 min")   # ~120 days
  light <- gen_light(150, -50, times)                # stationary known site for the test
  t0 <- as.numeric(times[1]); t1 <- as.numeric(times[length(times)])
  # injected drift: 0 at start, +30 min at end
  drift_true <- function(t) (30 * 60) * (as.numeric(t) - t0) / (t1 - t0)
  wrong_times <- times + drift_true(times)

  cc <- calibrate_clock(
    deploy   = list(lon = 150, lat = -50, time = times[1]),
    retrieve = list(lon = 150, lat = -50, time = times[length(times)]),
    times = wrong_times, light = light,
    calibration = cal, likelihood_params = lp, window_days = 4)

  corrected <- cc$correct(wrong_times)
  resid <- as.numeric(corrected) - as.numeric(times)  # should be ~0 across interior
  expect_lt(max(abs(resid)), 300)                     # within 5 minutes everywhere
  expect_lt(abs(cc$offset_retrieve - cc$offset_deploy - (-30 * 60)), 300)
})

test_that("calibrate_clock_from_endpoints derives a rough calibration and recovers an offset", {
  times <- seq(as.POSIXct("2024-05-01", tz = "UTC"),
               as.POSIXct("2024-05-10", tz = "UTC"), by = "20 min")
  light <- gen_light(150, -50, times)
  delta <- 20 * 60
  wrong <- times + delta
  cc <- calibrate_clock_from_endpoints(wrong, light,
          start_lon = 150, start_lat = -50, end_lon = 150, end_lat = -50)  # no calibration supplied
  # both endpoints share the location here, so both offsets ~ -delta
  expect_lt(abs(cc$offset_deploy - (-delta)), 240)
})

test_that("calibrate_clock_from_endpoints returns NULL when an endpoint is unknown", {
  times <- seq(as.POSIXct("2024-05-01", tz = "UTC"),
               as.POSIXct("2024-05-03", tz = "UTC"), by = "30 min")
  light <- gen_light(150, -50, times)
  expect_warning(
    res <- calibrate_clock_from_endpoints(times, light, 150, -50, NA_real_, NA_real_),
    "known start AND end")
  expect_null(res)
})

test_that("TwilightFreeGrid calibrate = TRUE corrects an injected clock offset", {
  skip_if_not_installed("raster")
  grid <- makeGrid(lon = c(140, 160), lat = c(-60, -40), cell.size = 2)
  times <- seq(as.POSIXct("2024-06-10", tz = "UTC"),
               as.POSIXct("2024-06-24", tz = "UTC"), by = "20 min")
  light <- gen_light(150, -50, times)
  wrong <- times + 20 * 60                              # 20-min clock error
  fit <- TwilightFreeGrid(wrong, light, grid,
    start_lat = -50, start_lon = 150, end_lat = -50, end_lon = 150,
    calibration = cal, likelihood_params = lp, step_hours = 12, diffusion = 150,
    calibrate = TRUE)
  expect_false(is.null(fit$clock))
  expect_lt(abs(fit$clock$offset_deploy - (-20 * 60)), 240)
  out <- paste(capture.output(print(fit)), collapse = "\n")
  expect_match(out, "Clock drift")                    # the calibrated fit reports its drift
})

test_that("format_clock summarises a fitted clock and is empty for NULL", {
  clock <- list(offset_deploy = -20 * 60, offset_retrieve = -50 * 60,
                drift = c(intercept = -20 * 60, slope_per_sec = (-30 * 60) / (100 * 86400)))
  s <- format_clock(clock)
  expect_length(s, 3)
  expect_match(s[1], "deploy")
  expect_match(paste(s, collapse = " "), "min/day")
  expect_equal(format_clock(NULL), character(0))
})
