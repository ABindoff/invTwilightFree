# Tests for floor_rule() and bathy_source().
# floor_rule() is pure and fully tested here. bathy_source() is tested through
# an injected raster (no network); the live marmap download is exercised only in
# an internet-gated test.

test_that("floor_rule hard constraint excludes cells shallower than the dive", {
  fr <- floor_rule()
  expect_identical(fr(obs = 300, expected = c(50, 299, 300, 1200)),
                   c(-Inf, -Inf, 0, 0))            # 300 m seabed just clears a 300 m dive
})

test_that("floor_rule soft constraint is zero when deep enough, quadratic when not", {
  fs <- floor_rule(sd = 200)
  out <- fs(obs = 300, expected = c(300, 1000, 100))
  expect_equal(out[1:2], c(0, 0))                  # deep enough -> no penalty
  expect_equal(out[3], -0.5 * ((300 - 100) / 200)^2)
})

test_that("floor_rule treats a missing observation as uninformative", {
  fr <- floor_rule()
  expect_equal(fr(obs = NA_real_, expected = c(10, 20, 30)), c(0, 0, 0))
  expect_equal(fr(obs = numeric(0), expected = c(10, 20)), c(0, 0))
})

test_that("floor_rule treats missing bathymetry as uninformative, not impossible", {
  fr <- floor_rule()                               # hard rule
  expect_equal(fr(obs = 500, expected = c(NA, 100, 1000)), c(0, -Inf, 0))
})

test_that("floor_rule covers altimetry under a positive-down convention", {
  # Tag at sea level (depth 0); terrain elevations expressed as depth-below-sea:
  # summit -2000, plain -100, shore 0, open sea 1200. Terrain above the animal
  # (negative depth, i.e. land poking up) is ruled out; shore and sea are fine.
  fr <- floor_rule()
  expect_identical(fr(obs = 0, expected = c(-2000, -100, 0, 1200)),
                   c(-Inf, -Inf, 0, 0))
  # A bird at +50 m (depth -50) can now clear the 100 m plain... no: plain depth
  # -100 is still above -50, so still ruled out; only the shore/sea remain.
  expect_identical(fr(obs = -50, expected = c(-2000, -100, 0, 1200)),
                   c(-Inf, -Inf, 0, 0))
})

test_that("floor_rule validates sd", {
  expect_error(floor_rule(sd = 0), "positive")
  expect_error(floor_rule(sd = -5), "positive")
  expect_error(floor_rule(sd = c(1, 2)), "single positive")
})

test_that("floor_rule returns a tf_rule", {
  expect_s3_class(floor_rule(), "tf_rule")
})

test_that("bathy_source returns positive-down depths from an injected raster", {
  skip_if_not_installed("raster")
  r <- raster::raster(nrows = 2, ncols = 2,
                      xmn = 149, xmx = 151, ymn = -46, ymx = -44)
  # raster fills by row from the top-left: top row (lat -44.5), then bottom (-45.5)
  raster::values(r) <- c(50, 1500, 0, 3000)
  bs <- bathy_source(rast = r)
  # point in the top-left cell (lon 149.5, lat -44.5) -> 50 m
  d <- bs(lon = c(149.5, 150.5, 149.5), lat = c(-44.5, -44.5, -45.5))
  expect_equal(d, c(50, 1500, 0))
  expect_s3_class(bs, "tf_source")
})

test_that("bathy_source + floor_rule compose into a depth feasibility surface", {
  skip_if_not_installed("raster")
  r <- raster::raster(nrows = 1, ncols = 3,
                      xmn = 0, xmx = 3, ymn = 0, ymx = 1)
  raster::values(r) <- c(100, 800, 2000)           # seabed depth per cell
  bs <- bathy_source(rast = r)
  fr <- floor_rule()
  depth <- bs(lon = c(0.5, 1.5, 2.5), lat = c(0.5, 0.5, 0.5))
  expect_equal(fr(obs = 500, expected = depth), c(-Inf, 0, 0))
})

test_that("bathy_source errors helpfully when marmap is absent and no raster given", {
  skip_if(requireNamespace("marmap", quietly = TRUE),
          "marmap is installed; the missing-dependency path cannot be exercised")
  bs <- bathy_source()
  expect_error(bs(lon = 150, lat = -45), "marmap")
})
