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
