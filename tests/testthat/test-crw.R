test_that("movement='brownian' is unchanged by the CRW addition", {
  skip_on_cran()
  mk <- function(days, lat0, lon0, seed) {
    set.seed(seed)
    t <- seq(as.POSIXct("2024-05-01", tz = "UTC"), by = "10 min", length.out = days * 144)
    lat <- lat0 + cumsum(stats::rnorm(length(t), 0, 0.02))
    lon <- lon0 + cumsum(stats::rnorm(length(t), 0, 0.02))
    z <- solar_zenith(as.numeric(t), lon, lat)
    list(df = data.frame(Date = t, Light = pmin(pmax(558.5 - 5.818 * z, 0), 64)),
         lon = lon, lat = lat)
  }
  a <- mk(5, -48, 150, 1); b <- mk(6, -46, 148, 2)
  data <- list(A = a$df, B = b$df)
  locs <- data.frame(id = c("A", "B"),
                     deploy_lon = c(a$lon[1], b$lon[1]),
                     deploy_lat = c(a$lat[1], b$lat[1]),
                     retrieve_lon = c(a$lon[length(a$lon)], NA),
                     retrieve_lat = c(a$lat[length(a$lat)], NA))
  cfg <- list(data = data, locations = locs, step_hours = 12, mesh_pad = 8,
              coarse_res = 2, surrogate_diffusion = 60,
              sweeps = 300L, burn = 100L, thin = 2L, seed = 7L)

  # The default must be identical to an explicit brownian request, and the CRW
  # code path must not perturb it: rho = 0 leaves every acceptance ratio alone.
  f_default  <- do.call(TwilightFreeHier, cfg)
  f_brownian <- do.call(TwilightFreeHier, c(cfg, list(movement = "brownian")))
  expect_equal(f_default$tracks$A$lon, f_brownian$tracks$A$lon)
  expect_equal(f_default$movement$sigma_kmday, f_brownian$movement$sigma_kmday)
  expect_null(f_default$persistence)

  f_crw <- do.call(TwilightFreeHier, c(cfg, list(movement = "crw")))
  expect_s3_class(f_crw, "TwilightFreeHier")
  expect_false(is.null(f_crw$persistence))
  expect_equal(nrow(f_crw$persistence), 2L)
  expect_true(all(abs(f_crw$persistence$rho) <= 0.99))
  expect_identical(f_crw$model, "crw")
})

test_that("the correlated random walk recovers directional persistence", {
  skip_on_cran()
  # A track built as an AR(1) on increments, so the truth is known. Persistence
  # is strong, which is the regime a memoryless walk cannot represent.
  mk_crw <- function(days, lat0, lon0, rho, seed, s = 0.02) {
    set.seed(seed)
    n <- days * 144
    dlon <- dlat <- numeric(n)
    for (i in 2:n) {
      dlon[i] <- rho * dlon[i - 1] + stats::rnorm(1, 0, s)
      dlat[i] <- rho * dlat[i - 1] + stats::rnorm(1, 0, s)
    }
    lon <- lon0 + cumsum(dlon); lat <- lat0 + cumsum(dlat)
    t <- seq(as.POSIXct("2024-05-01", tz = "UTC"), by = "10 min", length.out = n)
    z <- solar_zenith(as.numeric(t), lon, lat)
    list(df = data.frame(Date = t, Light = pmin(pmax(558.5 - 5.818 * z, 0), 64)),
         lon = lon, lat = lat)
  }
  a <- mk_crw(8, -48, 150, 0.9, 11)
  b <- mk_crw(8, -46, 148, 0.9, 12)
  data <- list(A = a$df, B = b$df)
  locs <- data.frame(id = c("A", "B"),
                     deploy_lon = c(a$lon[1], b$lon[1]),
                     deploy_lat = c(a$lat[1], b$lat[1]),
                     retrieve_lon = c(a$lon[length(a$lon)], b$lon[length(b$lon)]),
                     retrieve_lat = c(a$lat[length(a$lat)], b$lat[length(b$lat)]))
  f <- TwilightFreeHier(data, locs, step_hours = 12, mesh_pad = 8, coarse_res = 2,
                        surrogate_diffusion = 60, sweeps = 1500L, burn = 500L,
                        thin = 3L, movement = "crw", seed = 3L)
  # The knots are 12 h apart while the simulation steps every 10 min, so the
  # knot-scale persistence is weaker than the generating rho. The test is that
  # it is detected as clearly positive, not that it equals 0.9.
  expect_true(all(f$persistence$rho > 0.05))
  expect_true(all(f$persistence$lower > -0.5))
})

test_that("persistence is near zero for a memoryless track", {
  skip_on_cran()
  set.seed(21)
  n <- 8 * 144
  t <- seq(as.POSIXct("2024-05-01", tz = "UTC"), by = "10 min", length.out = n)
  lat <- -48 + cumsum(stats::rnorm(n, 0, 0.02))
  lon <- 150 + cumsum(stats::rnorm(n, 0, 0.02))
  z <- solar_zenith(as.numeric(t), lon, lat)
  df <- data.frame(Date = t, Light = pmin(pmax(558.5 - 5.818 * z, 0), 64))
  locs <- data.frame(id = "A", deploy_lon = lon[1], deploy_lat = lat[1],
                     retrieve_lon = lon[n], retrieve_lat = lat[n])
  f <- TwilightFreeHier(list(A = df), locs, step_hours = 12, mesh_pad = 8,
                        coarse_res = 2, surrogate_diffusion = 60,
                        sweeps = 1500L, burn = 500L, thin = 3L,
                        movement = "crw", seed = 5L)
  # A Brownian track has no directional memory, so the interval should admit 0.
  expect_lt(f$persistence$lower[1], 0.2)
})
