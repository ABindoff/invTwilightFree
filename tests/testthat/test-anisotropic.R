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

test_that("omitting diffusion_lon is bit-identical to the isotropic engine", {
  # The whole point of the default being NULL: existing analyses must not move.
  # Checked on log_z as well as the track, because the two passes of the smoother
  # were edited separately and an inconsistency between them would show up in the
  # marginal likelihood before it showed up in the MAP path.
  skip_if_not_installed("terra")
  d <- sim_light(); g <- make_grid()
  a <- TwilightFreeGrid(d$time, d$light, grid = g, start_lon = 150, start_lat = -45,
                        step_hours = 12, diffusion = 80)
  b <- TwilightFreeGrid(d$time, d$light, grid = g, start_lon = 150, start_lat = -45,
                        step_hours = 12, diffusion = 80, diffusion_lon = NULL)
  expect_identical(a$fit$lon, b$fit$lon)
  expect_identical(a$fit$lat, b$fit$lat)
  expect_identical(a$log_z, b$log_z)
})

test_that("equal anisotropic scales reproduce the isotropic fit closely", {
  # Not bit-identical: isotropic uses great-circle distance, anisotropic resolves
  # the step into components on a local flat-earth approximation. Over one step
  # the two agree to well under a kilometre, so the tracks should coincide on a
  # 1 degree grid and log_z should be close.
  skip_if_not_installed("terra")
  d <- sim_light(); g <- make_grid()
  iso <- TwilightFreeGrid(d$time, d$light, grid = g, start_lon = 150, start_lat = -45,
                          step_hours = 12, diffusion = 80)
  ani <- TwilightFreeGrid(d$time, d$light, grid = g, start_lon = 150, start_lat = -45,
                          step_hours = 12, diffusion = 80, diffusion_lon = 80)
  expect_equal(ani$fit$lon, iso$fit$lon, tolerance = 1e-6)
  expect_equal(ani$fit$lat, iso$fit$lat, tolerance = 1e-6)
  expect_equal(ani$log_z, iso$log_z, tolerance = 0.05)
})

test_that("a loose east-west scale widens longitude and not latitude", {
  # The mechanism the argument exists for: each axis takes its own prior width.
  # A prior loosened only in longitude must widen the longitude posterior while
  # leaving latitude where it was.
  skip_if_not_installed("terra")
  d <- sim_light(); g <- make_grid()
  post_sd <- function(f) {
    gp <- grid_posterior(f)
    sl <- sp <- numeric(nrow(gp$P))
    for (i in seq_len(nrow(gp$P))) {
      w <- gp$P[i, ]; s <- sum(w)
      if (!is.finite(s) || s <= 0) { sl[i] <- sp[i] <- NA; next }
      w <- w / s
      sl[i] <- sqrt(sum(w * (gp$lon - sum(w * gp$lon))^2))
      sp[i] <- sqrt(sum(w * (gp$lat - sum(w * gp$lat))^2))
    }
    c(lon = stats::median(sl, na.rm = TRUE), lat = stats::median(sp, na.rm = TRUE))
  }
  base <- post_sd(TwilightFreeGrid(d$time, d$light, grid = g, start_lon = 150,
                                   start_lat = -45, step_hours = 12, diffusion = 60,
                                   diffusion_lon = 60))
  wide <- post_sd(TwilightFreeGrid(d$time, d$light, grid = g, start_lon = 150,
                                   start_lat = -45, step_hours = 12, diffusion = 60,
                                   diffusion_lon = 240))
  expect_gt(wide[["lon"]], base[["lon"]])
  expect_lt(abs(wide[["lat"]] - base[["lat"]]), 0.5 * base[["lat"]])
})

test_that("diffusion_lon must match diffusion in length", {
  skip_if_not_installed("terra")
  d <- sim_light(); g <- make_grid()
  expect_error(
    TwilightFreeGrid(d$time, d$light, grid = g, start_lon = 150, start_lat = -45,
                     step_hours = 12, diffusion = c(40, 200),
                     trans_prob = c(0.9, 0.1, 0.1, 0.9), diffusion_lon = 80),
    "diffusion_lon")
})

test_that("area_correction weights cells by cos(latitude) and is optional", {
  # The movement kernel is a density on the sphere evaluated on a grid uniform in
  # degrees, so a cell's probability needs its area. Without the correction every
  # cell counts equally and the smaller high-latitude cells collect weight they
  # have not earned, which biases tracks poleward.
  #
  # The effect is only visible where the PRIOR has room to speak. Given light
  # strong enough to pin the position, the two agree, and that is a property
  # worth asserting rather than a nuisance: a correction to the prior must not
  # disturb a well-determined fit. So this checks both regimes, using `temper`
  # to weaken the emission in the second.
  skip_if_not_installed("terra")
  g <- terra::rast(xmin = 190, xmax = 210, ymin = 20, ymax = 70,
                   resolution = 2, crs = "EPSG:4326")
  terra::values(g) <- 1
  t <- seq(as.POSIXct("2021-06-01", tz = "UTC"), by = "30 min", length.out = 20 * 48)
  z <- solar_zenith(as.numeric(t), rep(200, length(t)), rep(45, length(t)))
  lig <- pmin(pmax(160 - 1.6 * z, 0), 100)
  # Both endpoints anchored and a 78 km step on a 222 km grid leaves the track no
  # room to move, so the two settings agree trivially. Free the far end and widen
  # the step, or this tests nothing.
  run <- function(area, temper) {
    TwilightFreeGrid(t, lig, grid = g, start_lon = 200, start_lat = 45,
                     step_hours = 12, diffusion = 500, calibration = c(160, 1.6),
                     likelihood_params = c(1/50, 100, 0.10, 0.35, 2, 0.10, temper),
                     area_correction = area)
  }
  # strong light: the correction must not move a fit the data already determine
  expect_equal(mean(run(TRUE, 1)$fit$lat), mean(run(FALSE, 1)$fit$lat),
               tolerance = 0.05)

  # weak light: now the prior decides, and the uncorrected kernel drifts north
  off <- run(FALSE, 0.02); on <- run(TRUE, 0.02)
  expect_gt(mean(off$fit$lat), 45.5)          # uncorrected leaves the truth
  expect_lt(mean(on$fit$lat), mean(off$fit$lat))
  expect_lt(abs(mean(on$fit$lat) - 45), abs(mean(off$fit$lat) - 45))

  # NOTE: this test constrains the SIGN of the correction and nothing else. A
  # coefficient sweep through the engine is monotone in the weight with no plateau
  # at 1, so every assertion above also passes for a half-Jacobian or a double one.
  # The coefficient is pinned by the equivalence test below, not here.
})

test_that("area_correction is exactly the deprecated area_prior() term", {
  # This is what pins the COEFFICIENT. `area_correction` and an `area_prior()`
  # term routed through `identity_rule()` are the same operator -- the same
  # log(cos(lat)) with the same 1e-6 clamp, applied to the same destination cell
  # in the initial, forward and backward passes -- so they must agree bit for bit.
  #
  # It does double duty. Any change to the weight, the clamp, or which cell index
  # the term lands on breaks the identity, and because `area_prior()` is plain R
  # the reference side cannot drift silently with the engine. It is also the proof
  # that `area_correction = FALSE` is the pre-correction engine: FALSE plus the
  # term reproduces TRUE, so FALSE is TRUE with the term removed.
  skip_if_not_installed("terra")
  g <- terra::rast(xmin = 190, xmax = 210, ymin = 20, ymax = 70,
                   resolution = 2, crs = "EPSG:4326")
  terra::values(g) <- 1
  t <- seq(as.POSIXct("2021-06-01", tz = "UTC"), by = "30 min", length.out = 20 * 48)
  z <- solar_zenith(as.numeric(t), rep(200, length(t)), rep(45, length(t)))
  lig <- pmin(pmax(160 - 1.6 * z, 0), 100)
  args <- list(t, lig, grid = g, start_lon = 200, start_lat = 45,
               step_hours = 12, diffusion = 500, calibration = c(160, 1.6),
               likelihood_params = c(1/50, 100, 0.10, 0.35, 2, 0.10, 0.02))

  flag <- do.call(TwilightFreeGrid, c(args, list(area_correction = TRUE)))
  term <- do.call(TwilightFreeGrid, c(args, list(
    area_correction = FALSE,
    terms = list(location_term(name = "area",
                               source = suppressWarnings(area_prior()),
                               rule = identity_rule())))))

  expect_identical(flag$fit$lat, term$fit$lat)
  expect_identical(flag$fit$lon, term$fit$lon)
  expect_equal(flag$log_z, term$log_z, tolerance = 1e-10)
})

test_that("supplying area_prior() alongside area_correction is an error", {
  # The double-count is silent and moves the fit a long way equatorward, so the
  # engine refuses instead of warning.
  skip_if_not_installed("terra")
  d <- sim_light(); g <- make_grid()
  tm <- location_term(name = "area", source = suppressWarnings(area_prior()),
                      rule = identity_rule())
  expect_error(
    TwilightFreeGrid(d$time, d$light, grid = g, start_lon = 150, start_lat = -45,
                     step_hours = 12, diffusion = 80, terms = list(tm),
                     area_correction = TRUE),
    "same operator")
  # ...and is allowed when the engine is not also applying it
  expect_no_error(
    TwilightFreeGrid(d$time, d$light, grid = g, start_lon = 150, start_lat = -45,
                     step_hours = 12, diffusion = 80, terms = list(tm),
                     area_correction = FALSE))
})

test_that("area_prior() is deprecated and warns", {
  expect_warning(area_prior(), "deprecated")
})

test_that("the 1/cos(lat) lean follows from the kernel arithmetic", {
  # HONESTY NOTE: this is a derivation check, not a test of the shipped engine.
  # It replicates the kernel in R and would pass if run_grid_hmm did not exist.
  # It is kept because it is what licenses the coefficient -- the empirical fits
  # cannot distinguish a weight of 1 from 0.6 -- but the engine's own behaviour is
  # constrained by the equivalence test above, not by this.
  #
  # The defect, isolated from any fitting: incoming kernel mass for a
  # destination cell, from a uniform source over a lon-lat grid.
  lon <- seq(150.5, 249.5, by = 1); lat <- seq(20.5, 69.5, by = 1)
  gg <- expand.grid(lon = lon, lat = lat)
  gcd <- function(lo1, la1, lo2, la2) {
    p1 <- la1 * pi / 180; p2 <- la2 * pi / 180
    a <- sin((p2 - p1) / 2)^2 + cos(p1) * cos(p2) * sin((lo2 - lo1) * pi / 360)^2
    2 * 6371 * asin(pmin(1, sqrt(a)))
  }
  incoming <- function(la) {
    d <- gcd(200, la, gg$lon, gg$lat)
    k <- d <= 550
    sum(exp(-(d[k]^2) / (2 * 110^2)))
  }
  las <- seq(26, 62, by = 4)
  w <- vapply(las, incoming, 0)
  expect_gt(cor(w / w[1], (1 / cos(las * pi / 180)) / (1 / cos(las[1] * pi / 180))), 0.999)
  # and the cos(lat) weighting removes it
  expect_lt(diff(range(w * cos(las * pi / 180))) / mean(w * cos(las * pi / 180)), 0.05)
})
