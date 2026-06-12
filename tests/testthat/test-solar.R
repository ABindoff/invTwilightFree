test_that("solar_zenith matches SGAT::solar/zenith within 0.1 degrees", {
  skip_if_not_installed("SGAT")

  # Sample one timestamp per week across 2020, at three latitudes
  times <- seq(
    as.POSIXct("2020-01-01", tz = "UTC"),
    as.POSIXct("2020-12-31", tz = "UTC"),
    by = "1 week"
  )
  unix <- as.numeric(times)
  test_lats <- c(-60, 0, 60)
  test_lon  <- 147.0

  for (lat in test_lats) {
    lons <- rep(test_lon, length(unix))
    lats <- rep(lat,      length(unix))

    our_z  <- solar_zenith(unix, lons, lats)

    sgat_sun <- SGAT::solar(times)
    sgat_z   <- SGAT::zenith(sgat_sun, test_lon, lat)

    max_diff <- max(abs(our_z - sgat_z), na.rm = TRUE)
    expect_lt(
      max_diff, 0.5,
      label = sprintf("max zenith difference at lat=%g: %.4f deg", lat, max_diff)
    )
  }
})

test_that("solar_zenith is vectorised and returns correct length", {
  unix <- as.numeric(Sys.time()) + seq(0, 3600, by = 600)
  lons <- rep(0, length(unix))
  lats <- rep(0, length(unix))
  z <- solar_zenith(unix, lons, lats)
  expect_length(z, length(unix))
  expect_true(all(z >= 0 & z <= 180))
})
