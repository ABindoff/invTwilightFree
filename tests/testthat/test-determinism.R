make_smc_call <- function(seed) {
  # Minimal synthetic light series: 48 hours at 10-minute intervals, Hobart
  times <- seq(
    as.POSIXct("2020-06-01 00:00:00", tz = "UTC"),
    as.POSIXct("2020-06-03 00:00:00", tz = "UTC"),
    by = "10 min"
  )
  # Simple sinusoidal light proxy (not calibrated, but enough to exercise the engine)
  hour  <- as.numeric(format(times, "%H")) + as.numeric(format(times, "%M")) / 60
  light <- pmax(0, 10 * sin((hour - 6) * pi / 12))

  TwilightFreeSMC(
    date_time    = times,
    light        = light,
    start_lat    = -42.88,
    start_lon    = 147.33,
    n_particles  = 200,
    step_hours   = 12,
    diffusion    = 50,
    method       = "forward",
    seed         = seed
  )
}

test_that("TwilightFreeSMC is reproducible with a fixed seed", {
  r1 <- make_smc_call(seed = 42)
  r2 <- make_smc_call(seed = 42)
  expect_identical(r1$lat, r2$lat)
  expect_identical(r1$lon, r2$lon)
})

test_that("TwilightFreeSMC gives different results with different seeds", {
  r1 <- make_smc_call(seed = 1)
  r2 <- make_smc_call(seed = 2)
  # Different seeds should (almost certainly) produce different tracks
  expect_false(identical(r1$lat, r2$lat))
})

test_that("TwilightFreeGrid is deterministic without a seed", {
  skip_if_not_installed("rnaturalearth")
  times <- seq(
    as.POSIXct("2020-06-01 00:00:00", tz = "UTC"),
    as.POSIXct("2020-06-03 00:00:00", tz = "UTC"),
    by = "10 min"
  )
  hour  <- as.numeric(format(times, "%H")) + as.numeric(format(times, "%M")) / 60
  light <- pmax(0, 10 * sin((hour - 6) * pi / 12))

  grid <- makeGrid(lon = c(140, 155), lat = c(-48, -37), cell.size = 1)

  run_grid <- function() {
    TwilightFreeGrid(
      date_time  = times,
      light      = light,
      grid       = grid,
      start_lat  = -42.88,
      start_lon  = 147.33,
      step_hours = 12,
      diffusion  = 50
    )
  }

  g1 <- run_grid()
  g2 <- run_grid()
  expect_identical(g1$fit$lat, g2$fit$lat)
  expect_identical(g1$fit$lon, g2$fit$lon)
})
