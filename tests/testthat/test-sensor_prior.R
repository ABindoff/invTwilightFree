# Tests for hemisphere_prior() and identity_rule().
# These two functions are pure (no I/O), so they are fully testable now,
# ahead of the rest of the sensor-fusion framework (see HANDOFF_2_SENSOR_FUSION.md).

test_that("hemisphere_prior favours the requested southern hemisphere", {
  src <- hemisphere_prior(function(d) "S", softness = 1e-3)
  lat <- c(-30, -1, 0, 1, 30)
  out <- src(rep(140, length(lat)), lat, as.Date("2024-06-21"))
  expect_equal(out[lat <= 0], rep(0, sum(lat <= 0)))
  expect_equal(out[lat > 0], rep(log(1e-3), sum(lat > 0)))
})

test_that("hemisphere_prior favours the requested northern hemisphere", {
  src <- hemisphere_prior(function(d) "N", softness = 1e-2)
  lat <- c(-10, 0, 10)
  out <- src(c(1, 1, 1), lat, as.Date("2024-01-01"))
  expect_equal(out, c(log(1e-2), 0, 0))   # boundary (0) counts as allowed for N
})

test_that("'both' applies no constraint", {
  src <- hemisphere_prior(function(d) "both")
  expect_equal(src(c(1, 2), c(-50, 50), Sys.Date()), c(0, 0))
})

test_that("softness = 0 is a hard cut (-Inf)", {
  src <- hemisphere_prior(function(d) "N", softness = 0)
  out <- src(c(1, 1), c(-10, 10), Sys.Date())
  expect_identical(out, c(-Inf, 0))
})

test_that("boundary slack shifts the dividing line", {
  src <- hemisphere_prior(function(d) "S", softness = 1e-3, boundary = 10)
  out <- src(c(1, 1), c(5, 15), Sys.Date())   # 5 <= 10 allowed; 15 not
  expect_equal(out, c(0, log(1e-3)))
})

test_that("date is actually passed through to the season function", {
  # June -> S, December -> both
  season <- function(d) if (as.integer(format(d, "%m")) == 12) "both" else "S"
  src <- hemisphere_prior(season, softness = 1e-3)
  jun <- src(c(1, 1), c(-20, 20), as.Date("2024-06-15"))
  dec <- src(c(1, 1), c(-20, 20), as.Date("2024-12-15"))
  expect_equal(jun, c(0, log(1e-3)))
  expect_equal(dec, c(0, 0))
})

test_that("invalid softness errors", {
  expect_error(hemisphere_prior(function(d) "S", softness = 2), "in \\[0, 1\\]")
  expect_error(hemisphere_prior(function(d) "S", softness = -1), "in \\[0, 1\\]")
})

test_that("invalid hemisphere_by_date return value errors at evaluation", {
  src <- hemisphere_prior(function(d) "X")
  expect_error(src(1, 1, Sys.Date()), "must return one of")
})

test_that("identity_rule passes log-densities and neutralises NA", {
  r <- identity_rule()
  expect_equal(r(obs = NULL, expected = c(0, log(1e-3), NA)), c(0, log(1e-3), 0))
})

test_that("identity_rule preserves -Inf hard cuts", {
  r <- identity_rule()
  expect_identical(r(NULL, c(-Inf, 0)), c(-Inf, 0))
})

test_that("identity_rule logs densities when is_log = FALSE", {
  r <- identity_rule(is_log = FALSE)
  expect_equal(r(obs = NULL, expected = c(1, exp(1), 0)), c(0, 1, -Inf))
})

test_that("returned objects carry the framework classes", {
  expect_s3_class(hemisphere_prior(function(d) "both"), "tf_source")
  expect_s3_class(identity_rule(), "tf_rule")
})
