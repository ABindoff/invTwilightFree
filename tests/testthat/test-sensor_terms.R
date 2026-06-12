# Tests for location_term() and build_aux_matrix(), the sensor-fusion plumbing.
# The flatten-order test is the important one: build_aux_matrix returns K x n,
# TwilightFreeGrid flattens with t() to row-major (k*n + i), and run_grid_hmm
# reads aux_logl[k*n + i]. A mismatch here is silent and catastrophic.

make_knots <- function(start = "2024-06-20", K = 4, step_h = 12) {
  t0 <- as.numeric(as.POSIXct(start, tz = "UTC"))
  t0 + (0:(K - 1)) * step_h * 3600
}

test_that("location_term validates its inputs", {
  expect_error(location_term("x", source = 1, rule = identity_rule()), "tf_source")
  expect_error(location_term("x", source = hemisphere_prior(function(d) "S"), rule = 1),
               "tf_rule")
  bad_tag <- data.frame(t = 1, v = 2)
  expect_error(
    location_term("x", tag = bad_tag, source = bathy_source(rast = NULL),
                  rule = floor_rule()),
    "columns 'time' and 'value'"
  )
})

test_that("empty terms give an all-zero K x n matrix", {
  kt <- make_knots(K = 3)
  m <- build_aux_matrix(list(), lon = c(0, 1), lat = c(0, 1), knot_times = kt)
  expect_equal(dim(m), c(3, 2))
  expect_true(all(m == 0))
})

test_that("a hemisphere prior fills each knot row from the source", {
  kt <- make_knots(start = "2024-06-21", K = 2)        # austral winter
  term <- location_term("hemi",
                        source = hemisphere_prior(function(d) "S", softness = 1e-3),
                        rule = identity_rule())
  lon <- c(140, 140, 140)
  lat <- c(-30, 0, 30)
  m <- build_aux_matrix(list(term), lon, lat, kt)
  # southern and equatorial cells allowed (0), northern downweighted
  expect_equal(m[1, ], c(0, 0, log(1e-3)))
  expect_equal(m[2, ], c(0, 0, log(1e-3)))
})

test_that("flatten order matches the Rust k*n + i convention", {
  # Build a term whose value is deterministic per (knot, cell) so we can check
  # the exact mapping. Use a custom source returning lat as the 'field' and an
  # identity rule, so aux[k, i] = lat[i] for every k.
  src <- structure(function(lon, lat, date) lat, class = c("tf_source", "function"))
  term <- location_term("probe", source = src, rule = identity_rule())
  lon <- c(10, 20, 30)
  lat <- c(-1, -2, -3)
  kt <- make_knots(K = 4)
  aux <- build_aux_matrix(list(term), lon, lat, kt)

  flat <- as.numeric(t(aux))                 # what TwilightFreeGrid sends to Rust
  n <- length(lon); K <- length(kt)
  for (k in seq_len(K)) {
    for (i in seq_len(n)) {
      # Rust reads aux_logl[(k-1)*n + (i-1)] (0-based) == aux[k, i]
      expect_equal(flat[(k - 1) * n + i], aux[k, i])
      expect_equal(flat[(k - 1) * n + i], lat[i])
    }
  }
})

test_that("a missing tag reading contributes the term's `missing` value", {
  kt <- make_knots(start = "2024-06-20", K = 3, step_h = 12)
  # Tag reading only in the first interval (t0-12h, t0]; later knots are missing.
  # The first knot is at midnight June 20, so the first interval is noon-to-midnight
  # June 19. Observation at 23:00 June 19 falls inside it.
  tag <- data.frame(time = as.POSIXct("2024-06-19 23:00", tz = "UTC"), value = 12)
  r <- raster::raster(nrows = 1, ncols = 2, xmn = 0, xmx = 2, ymn = 0, ymx = 1)
  raster::values(r) <- c(12, 20)
  term <- location_term("sst", tag = tag, reduce = "median",
                        source = sst_source(rast = r),
                        rule = student_rule(sd = 1.5),
                        missing = 0)
  m <- build_aux_matrix(list(term), lon = c(0.5, 1.5), lat = c(0.5, 0.5), knot_times = kt)
  # First knot has the reading (non-zero ll); the rest are exactly the missing value.
  expect_false(all(m[1, ] == 0))
  expect_true(all(m[2, ] == 0))
  expect_true(all(m[3, ] == 0))
})

test_that("weight scales a term's contribution linearly", {
  kt <- make_knots(K = 2)
  mk <- function(w) {
    term <- location_term("hemi", weight = w,
                          source = hemisphere_prior(function(d) "S", softness = 1e-2),
                          rule = identity_rule())
    build_aux_matrix(list(term), lon = c(0, 0), lat = c(-10, 10), kt)
  }
  expect_equal(mk(2), mk(1) * 2)
})

test_that("multiple terms sum", {
  kt <- make_knots(K = 2)
  t1 <- location_term("a", source = hemisphere_prior(function(d) "S", softness = 1e-2),
                      rule = identity_rule())
  t2 <- location_term("b", source = hemisphere_prior(function(d) "S", softness = 1e-2),
                      rule = identity_rule())
  one  <- build_aux_matrix(list(t1), c(0, 0), c(-10, 10), kt)
  both <- build_aux_matrix(list(t1, t2), c(0, 0), c(-10, 10), kt)
  expect_equal(both, one * 2)
})
