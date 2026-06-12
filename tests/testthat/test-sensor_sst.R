# Tests for student_rule() and sst_source().
# student_rule() is pure and fully tested. sst_source() is tested through
# injected rasters (no network); the live ERDDAP fetch is not exercised here.

test_that("student_rule matches a scaled, logged Student-t density", {
  sr <- student_rule(sd = 2, df = 4)
  obs <- 10
  expected <- c(10, 12, 6)
  resid <- (obs - expected) / 2
  expect_equal(sr(obs, expected),
               stats::dt(resid, df = 4, log = TRUE) - log(2))
})

test_that("student_rule peaks at a perfect match and decays with distance", {
  sr <- student_rule(sd = 1.5, df = 4)
  out <- sr(obs = 12, expected = c(12, 13, 18))
  expect_equal(which.max(out), 1L)                 # exact match is most likely
  expect_true(all(diff(out) < 0))                  # monotonically worse as SST departs
})

test_that("student_rule approaches Gaussian as df grows", {
  obs <- 10; expected <- c(9, 11, 14); sd <- 2
  big <- student_rule(sd = sd, df = 1e6)(obs, expected)
  gauss <- stats::dnorm((obs - expected) / sd, log = TRUE) - log(sd)
  expect_equal(big, gauss, tolerance = 1e-3)
})

test_that("student_rule treats missing obs or missing field as uninformative", {
  sr <- student_rule(sd = 1)
  expect_equal(sr(obs = NA_real_, expected = c(10, 11)), c(0, 0))
  expect_equal(sr(obs = numeric(0), expected = c(10, 11)), c(0, 0))
  expect_equal(sr(obs = 10, expected = c(10, NA, 12))[2], 0)
})

test_that("student_rule validates sd and df", {
  expect_error(student_rule(sd = 0), "positive")
  expect_error(student_rule(sd = 1, df = -2), "positive")
})

test_that("student_rule returns a tf_rule", {
  expect_s3_class(student_rule(sd = 1), "tf_rule")
})

test_that("sst_source with a static raster returns SST and ignores date", {
  skip_if_not_installed("raster")
  r <- raster::raster(nrows = 2, ncols = 2,
                      xmn = 149, xmx = 151, ymn = -46, ymx = -44)
  raster::values(r) <- c(11.2, 11.8, 12.1, 12.6)   # top row first
  ss <- sst_source(rast = r)
  a <- ss(lon = c(149.5, 150.5), lat = c(-44.5, -45.5), date = as.Date("2024-01-15"))
  b <- ss(lon = c(149.5, 150.5), lat = c(-44.5, -45.5), date = as.Date("2024-07-15"))
  expect_equal(a, c(11.2, 12.6))
  expect_equal(a, b)                                # static layer: date-independent
  expect_s3_class(ss, "tf_source")
})

test_that("sst_source with a date-keyed list selects the right layer", {
  skip_if_not_installed("raster")
  mk <- function(v) {
    r <- raster::raster(nrows = 1, ncols = 1, xmn = 149, xmx = 151, ymn = -46, ymx = -44)
    raster::values(r) <- v
    r
  }
  layers <- list("2024-01-15" = mk(18.0), "2024-07-15" = mk(9.0))
  ss <- sst_source(rast = layers)
  expect_equal(ss(150, -45, as.Date("2024-01-15")), 18.0)
  expect_equal(ss(150, -45, as.Date("2024-07-15")), 9.0)
})

test_that("sst_source composes with student_rule into an SST match surface", {
  skip_if_not_installed("raster")
  r <- raster::raster(nrows = 1, ncols = 3, xmn = 0, xmx = 3, ymn = 0, ymx = 1)
  raster::values(r) <- c(10, 14, 18)               # SST per cell
  ss <- sst_source(rast = r)
  sr <- student_rule(sd = 1.5, df = 4)
  sst <- ss(lon = c(0.5, 1.5, 2.5), lat = c(0.5, 0.5, 0.5), date = Sys.Date())
  ll <- sr(obs = 14, expected = sst)               # tag read 14 C
  expect_equal(which.max(ll), 2L)                  # the 14 C cell wins
})

test_that("sst_source errors helpfully when rerddap is absent and no raster given", {
  skip_if(requireNamespace("rerddap", quietly = TRUE),
          "rerddap is installed; the missing-dependency path cannot be exercised")
  ss <- sst_source()
  expect_error(ss(lon = 150, lat = -45, date = Sys.Date()), "rerddap")
})
