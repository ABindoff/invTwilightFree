mk <- function(lon = 200, lat = 45, days = 14) {
  t <- seq(as.POSIXct("2021-06-01", tz = "UTC"), by = "30 min", length.out = days * 48)
  z <- solar_zenith(as.numeric(t), rep(lon, length(t)), rep(lat, length(t)))
  list(time = t, light = pmin(pmax(160 - 1.6 * z, 0), 100))
}
CAL <- c(160, 1.6); LP <- c(1/50, 100, 0.10)
gr <- function() {
  g <- terra::rast(xmin = 180, xmax = 220, ymin = 30, ymax = 60,
                   resolution = 2, crs = "EPSG:4326")
  terra::values(g) <- 1; g
}
fit1 <- function(...) TwilightFreeGrid(mk()$time, mk()$light, grid = gr(),
  start_lon = 200, start_lat = 45, step_hours = 12, diffusion = 80,
  calibration = CAL, likelihood_params = LP, ...)

test_that("lambda_scale = NULL and a constant 1 are bit-identical", {
  # The default must not move any existing result. A constant schedule also has to
  # be exactly inert, because that is what makes `slope = 0` a genuine control
  # rather than an approximate one.
  skip_if_not_installed("terra")
  a <- fit1()
  b <- fit1(lambda_scale = 1)
  k <- length(a$fit$lat)
  expect_identical(a$fit$lat, b$fit$lat)
  expect_identical(a$fit$lon, b$fit$lon)
  expect_identical(a$log_z, b$log_z)
  expect_identical(fit1(lambda_scale = rep(1, k))$fit$lat, a$fit$lat)
})

test_that("a constant scale c is the same as scaling lambda by c", {
  # Fixes the MEANING of the multiplier: it multiplies lambda, nothing else.
  # Without this the argument could silently be doing something adjacent, such as
  # scaling the whole log-likelihood, and the schedule would be uninterpretable.
  skip_if_not_installed("terra")
  scaled <- TwilightFreeGrid(mk()$time, mk()$light, grid = gr(), start_lon = 200,
    start_lat = 45, step_hours = 12, diffusion = 80, calibration = CAL,
    likelihood_params = c(LP[1] * 3, LP[2], LP[3]))
  viascale <- fit1(lambda_scale = 3)
  expect_equal(viascale$fit$lat, scaled$fit$lat, tolerance = 1e-10)
  expect_equal(viascale$log_z, scaled$log_z, tolerance = 1e-8)
})

test_that("a non-constant schedule changes the fit", {
  # Guards against the argument being accepted and ignored -- a failure mode this
  # package has had before, where `mask_matrix` in the SMC path is documented and
  # inert.
  #
  # Assert on log_z and the posterior, NOT on the MAP track: the MAP is snapped to
  # cell centres, so on a well-identified synthetic it can be bit-identical while
  # the likelihood underneath has changed a great deal. An earlier version of this
  # test checked the track and passed a working implementation as broken.
  skip_if_not_installed("terra")
  a <- fit1()
  k <- length(a$fit$lat)
  s <- seq(0.4, 2.5, length.out = k)
  b <- fit1(lambda_scale = s / exp(mean(log(s))))
  expect_false(isTRUE(all.equal(a$log_z, b$log_z)))
  pa <- grid_posterior(a)$P; pb <- grid_posterior(b)$P
  expect_gt(max(abs(pa - pb)), 1e-8)
})

test_that("lambda_scale rejects bad input", {
  skip_if_not_installed("terra")
  expect_error(fit1(lambda_scale = c(1, 2, 3)), "one value per knot")
  expect_error(fit1(lambda_scale = -1), "positive")
  expect_error(fit1(lambda_scale = 0), "positive")
  expect_error(fit1(lambda_scale = NA_real_), "positive|finite")
})

test_that("declination_lambda_scale is tighter at the equinox than the solstice", {
  # The whole point of the schedule. Latitude is inferred from day length, which
  # varies least with latitude near an equinox, so that is where the likelihood must
  # work hardest.
  k <- seq(as.POSIXct("2021-01-01", tz = "UTC"), by = "12 hours", length.out = 730)
  s <- declination_lambda_scale(k)
  dec <- abs(solar_declination(k))
  expect_gt(mean(s[dec < 3]), mean(s[dec > 20]))     # equinox tighter than solstice
  expect_equal(exp(mean(log(s))), 1, tolerance = 1e-10)   # level unchanged
  expect_lte(max(s) / min(s), 6 + 1e-8)                    # cap respected
})

test_that("slope = 0 gives an exactly constant schedule", {
  # so that `slope = 0` is a clean control for the simulation study
  k <- seq(as.POSIXct("2021-01-01", tz = "UTC"), by = "12 hours", length.out = 200)
  s <- declination_lambda_scale(k, slope = 0)
  expect_equal(diff(range(s)), 0, tolerance = 1e-12)
  expect_equal(unname(s[1]), 1, tolerance = 1e-12)
})

test_that("max_ratio caps an aggressive slope", {
  k <- seq(as.POSIXct("2021-01-01", tz = "UTC"), by = "12 hours", length.out = 730)
  s <- declination_lambda_scale(k, slope = -0.5, max_ratio = 2)
  expect_lte(max(s) / min(s), 2 + 1e-8)
  expect_equal(exp(mean(log(s))), 1, tolerance = 1e-10)
})

test_that("solar_declination has the right sign and magnitude", {
  d <- solar_declination(as.POSIXct(c("2021-06-21", "2021-12-21", "2021-03-20"),
                                    tz = "UTC"))
  expect_gt(d[1], 22); expect_lt(d[1], 24)     # northern solstice
  expect_lt(d[2], -22); expect_gt(d[2], -24)   # southern solstice
  expect_lt(abs(d[3]), 2)                       # equinox
})
