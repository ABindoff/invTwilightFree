test_that("TwilightFreeHier runs end-to-end and returns a well-formed object", {
  skip_on_cran()

  mk_tag <- function(days, lat0, lon0, seed) {
    set.seed(seed)
    times <- seq(as.POSIXct("2024-05-01", tz = "UTC"), by = "10 min", length.out = days * 144)
    lat <- lat0 + cumsum(stats::rnorm(length(times), 0, 0.02))
    lon <- lon0 + cumsum(stats::rnorm(length(times), 0, 0.02))
    z <- solar_zenith(as.numeric(times), lon, lat)
    light <- pmin(pmax(558.5 - 5.818 * z, 0), 64)
    list(df = data.frame(Date = times, Light = light), lat = lat, lon = lon)
  }
  a <- mk_tag(5, -48, 150, 1)
  b <- mk_tag(7, -46, 148, 2)          # different length
  data <- list(A = a$df, B = b$df)
  locations <- data.frame(
    id = c("A", "B"),
    deploy_lon = c(a$lon[1], b$lon[1]), deploy_lat = c(a$lat[1], b$lat[1]),
    retrieve_lon = c(a$lon[length(a$lon)], NA),   # A fixed, B free
    retrieve_lat = c(a$lat[length(a$lat)], NA))

  fit <- TwilightFreeHier(data, locations, step_hours = 12,
                          mesh_pad = 8, coarse_res = 2, surrogate_diffusion = 60,
                          sweeps = 200L, burn = 50L, thin = 2L, seed = 1L)

  expect_s3_class(fit, "TwilightFreeHier")
  expect_equal(nrow(fit$movement), 2L)
  expect_equal(fit$movement$id, c("A", "B"))
  expect_true(all(fit$movement$sigma_kmday > 0))
  expect_named(fit$tracks, c("A", "B"))
  expect_true(all(is.finite(fit$tracks$A$lon)))

  # deployment is fixed: first knot has zero posterior sd and equals the input
  expect_equal(fit$tracks$A$lon_sd[1], 0)
  expect_equal(fit$tracks$A$lon[1], locations$deploy_lon[1], tolerance = 1e-6)
  # A has a fixed retrieval (last knot pinned, ~zero sd); B is free (positive sd)
  expect_lt(fit$tracks$A$lat_sd[nrow(fit$tracks$A)], 1e-6)
  expect_gt(fit$tracks$B$lat_sd[nrow(fit$tracks$B)], 0)
})

test_that("TwilightFreeHier accepts a long data frame with an id column", {
  skip_on_cran()
  set.seed(3)
  one <- function(id, lat0, lon0) {
    times <- seq(as.POSIXct("2024-05-01", tz = "UTC"), by = "10 min", length.out = 5 * 144)
    lat <- lat0 + cumsum(stats::rnorm(length(times), 0, 0.02))
    lon <- lon0 + cumsum(stats::rnorm(length(times), 0, 0.02))
    z <- solar_zenith(as.numeric(times), lon, lat)
    data.frame(id = id, Date = times, Light = pmin(pmax(558.5 - 5.818 * z, 0), 64),
               dlon = lon[1], dlat = lat[1])
  }
  long <- rbind(one("X", -48, 150), one("Y", -45, 149))
  locs <- data.frame(id = c("X", "Y"),
                     deploy_lon = c(long$dlon[1], long$dlon[long$id == "Y"][1]),
                     deploy_lat = c(long$dlat[1], long$dlat[long$id == "Y"][1]))
  fit <- TwilightFreeHier(long[, c("id", "Date", "Light")], locs, id = "id",
                          mesh_pad = 8, coarse_res = 2, sweeps = 150L, burn = 50L, thin = 2L, seed = 2L)
  expect_s3_class(fit, "TwilightFreeHier")
  expect_equal(fit$movement$id, c("X", "Y"))
})
