mk_track <- function(lon = 200, lat = 45, days = 14) {
  t <- seq(as.POSIXct("2021-06-01", tz = "UTC"), by = "30 min", length.out = days * 48)
  z <- solar_zenith(as.numeric(t), rep(lon, length(t)), rep(lat, length(t)))
  list(time = t, light = pmin(pmax(160 - 1.6 * z, 0), 100))
}
CAL <- c(160, 1.6); LP <- c(1/50, 100, 0.10)

test_that("posterior_extent brackets the posterior and shrinks as epsilon grows", {
  skip_if_not_installed("terra")
  d <- mk_track()
  g <- makeGrid(lon = c(170, 230), lat = c(25, 65), cell.size = 2)
  f <- TwilightFreeGrid(d$time, d$light, grid = g, start_lon = 200, start_lat = 45,
                        step_hours = 12, diffusion = 60, calibration = CAL,
                        likelihood_params = LP)
  tight <- posterior_extent(f, epsilon = 1e-9)
  loose <- posterior_extent(f, epsilon = 1e-2)

  # the truth must be inside, or the box is not fit for refitting
  expect_lt(tight[["lon_min"]], 200); expect_gt(tight[["lon_max"]], 200)
  expect_lt(tight[["lat_min"]], 45);  expect_gt(tight[["lat_max"]], 45)
  # a laxer epsilon discards more mass, so the box cannot grow
  expect_lte(loose[["lon_max"]] - loose[["lon_min"]],
             tight[["lon_max"]] - tight[["lon_min"]] + 1e-8)
  expect_true(attr(tight, "retained") <= 1)
})

test_that("the margin is applied and scales with the movement scale", {
  # Not cosmetic: the transition kernel reaches 5 sigma and is NOT renormalised per
  # source cell, so a box drawn tight to the posterior would drop outgoing mass at
  # every knot near the new edge -- a position-dependent tilt, not a rounding error.
  skip_if_not_installed("terra")
  d <- mk_track()
  g <- makeGrid(lon = c(170, 230), lat = c(25, 65), cell.size = 2)
  run <- function(diff) TwilightFreeGrid(d$time, d$light, grid = g, start_lon = 200,
      start_lat = 45, step_hours = 12, diffusion = diff, calibration = CAL,
      likelihood_params = LP)
  slow <- posterior_extent(run(40), epsilon = 1e-6)
  fast <- posterior_extent(run(160), epsilon = 1e-6)
  expect_gt(attr(fast, "pad_km"), attr(slow, "pad_km"))
  expect_gte(fast[["lat_max"]] - fast[["lat_min"]], slow[["lat_max"]] - slow[["lat_min"]])
})

test_that("refining reproduces the full-domain fit", {
  # THE ACCEPTANCE TEST. Refinement is only worth having if it is an optimisation
  # rather than an approximation. It cannot be bit-identical -- a different extent
  # moves the cell centres, so the refined grid is not a subset of the original -- so
  # the standard is agreement to within a cell, at a conservative epsilon.
  skip_if_not_installed("terra")
  d <- mk_track()
  full <- TwilightFreeGrid(d$time, d$light,
    grid = makeGrid(lon = c(170, 230), lat = c(25, 65), cell.size = 1),
    start_lon = 200, start_lat = 45, step_hours = 12, diffusion = 60,
    calibration = CAL, likelihood_params = LP)
  ref <- TwilightFreeGridRefine(d$time, d$light, start_lon = 200, start_lat = 45,
    step_hours = 12, diffusion = 60, calibration = CAL, likelihood_params = LP,
    lon = c(170, 230), lat = c(25, 65), cell.size = 1, coarse.size = 4,
    epsilon = 1e-9, verbose = FALSE)

  expect_s3_class(ref, "TwilightFreeGrid")
  expect_equal(length(ref$fit$lat), length(full$fit$lat))
  expect_lt(max(abs(ref$fit$lat - full$fit$lat)), 1.01)   # within one 1-degree cell
  expect_lt(max(abs(ref$fit$lon - full$fit$lon)), 1.01)
  expect_lt(abs(mean(ref$fit$lat) - mean(full$fit$lat)), 0.5)
  # and it must actually have shrunk the domain, or it saved nothing
  expect_lt(length(ref$cell_lon), length(full$cell_lon))
  expect_true(is.numeric(attr(ref, "extent")))
  expect_s3_class(attr(ref, "coarse"), "TwilightFreeGrid")
})

test_that("refine rejects incoherent arguments", {
  skip_if_not_installed("terra")
  d <- mk_track()
  expect_error(
    TwilightFreeGridRefine(d$time, d$light, grid = 1, lon = c(170, 230),
                           lat = c(25, 65), cell.size = 1, verbose = FALSE),
    "grid")
  expect_error(
    TwilightFreeGridRefine(d$time, d$light, lon = c(170, 230), lat = c(25, 65),
                           cell.size = 4, coarse.size = 1, verbose = FALSE),
    "coarse")
  g <- makeGrid(lon = c(170, 230), lat = c(25, 65), cell.size = 4)
  f <- TwilightFreeGrid(d$time, d$light, grid = g, start_lon = 200, start_lat = 45,
                        step_hours = 12, diffusion = 60, calibration = CAL,
                        likelihood_params = LP)
  expect_error(posterior_extent(f, epsilon = 0), "epsilon")
  expect_error(posterior_extent(f, epsilon = 1), "epsilon")
})

test_that("the fit records the movement settings the margin needs", {
  # posterior_extent() sizes its margin in km, so it needs diffusion and step_hours.
  # Before this they were not returned and the margin could not be computed from a
  # fit alone.
  skip_if_not_installed("terra")
  d <- mk_track()
  f <- TwilightFreeGrid(d$time, d$light,
    grid = makeGrid(lon = c(180, 220), lat = c(35, 55), cell.size = 4),
    start_lon = 200, start_lat = 45, step_hours = 12, diffusion = 77,
    calibration = CAL, likelihood_params = LP)
  expect_equal(f$diffusion, 77)
  expect_equal(f$step_hours, 12)
  expect_true(isTRUE(f$area_correction))
})
