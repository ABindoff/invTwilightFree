skip_if_not_installed("terra")

mk_grid <- function(xmin = 140, xmax = 180, ymin = -70, ymax = -40, res = 2) {
  g <- terra::rast(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                   resolution = res, crs = "EPSG:4326")
  terra::values(g) <- 1
  g
}

# Tracks are placed on CELL CENTRES. A point on a cell corner is equidistant
# from four cells, so a symmetric kernel ties and which.max() picks arbitrarily
# among them -- which looks like a centring bug and is not one.
mk_track <- function(lon, lat, n = 40, sd = NULL,
                     t0 = as.POSIXct("2021-06-01", tz = "UTC")) {
  d <- data.frame(time = seq(t0, by = "12 hours", length.out = n),
                  lon = rep(lon, length.out = n), lat = rep(lat, length.out = n))
  if (!is.null(sd)) { d$lon_sd <- sd; d$lat_sd <- sd }
  d
}

test_that("a point map is a normalised utilisation distribution", {
  g <- mk_grid()
  m <- occupancy_map(mk_track(151, -61), grid = g, method = "point")
  expect_s4_class(m, "SpatRaster")
  expect_equal(terra::nlyr(m), 1L)
  expect_equal(sum(terra::values(m), na.rm = TRUE), 1, tolerance = 1e-10)
  # a stationary track puts everything in one cell
  expect_equal(max(terra::values(m), na.rm = TRUE), 1, tolerance = 1e-10)
  expect_equal(sum(terra::values(m) > 0, na.rm = TRUE), 1L)
})

test_that("`by` produces one normalised layer per period", {
  g <- mk_grid()
  tr <- mk_track(151, -61, n = 200)          # 100 days, so several months
  tr$lon <- c(rep(151, 100), rep(171, 100))
  m <- occupancy_map(tr, grid = g, method = "point",
                     by = function(t) format(t, "%Y-%m"))
  expect_gt(terra::nlyr(m), 1L)
  s <- terra::global(m, "sum", na.rm = TRUE)[, 1]
  expect_true(all(abs(s - 1) < 1e-10))
})

test_that("weight = 'tag' stops a long track outvoting a short one", {
  g <- mk_grid()
  short <- mk_track(151, -61, n = 4)
  long  <- mk_track(171, -51, n = 400)
  by_tag  <- occupancy_map(list(a = short, b = long), grid = g,
                           method = "point", weight = "tag")
  by_knot <- occupancy_map(list(a = short, b = long), grid = g,
                           method = "point", weight = "knot")
  cell_short <- terra::cellFromXY(g, cbind(151, -61))
  expect_equal(terra::values(by_tag)[cell_short], 0.5, tolerance = 1e-8)
  # under knot weighting the short track is swamped
  expect_lt(terra::values(by_knot)[cell_short], 0.05)
})

test_that("the posterior form is more diffuse than the point form", {
  g <- mk_grid()
  tr <- mk_track(151, -61, n = 30, sd = 3)
  pt <- occupancy_map(tr, grid = g, method = "point")
  po <- occupancy_map(tr, grid = g, method = "posterior")
  expect_equal(sum(terra::values(po), na.rm = TRUE), 1, tolerance = 1e-8)
  expect_gt(sum(terra::values(po) > 0, na.rm = TRUE),
            sum(terra::values(pt) > 0, na.rm = TRUE))
  # both are centred on the same cell
  expect_equal(which.max(terra::values(po)), which.max(terra::values(pt)))
})

test_that("posterior is refused when a track carries no uncertainty", {
  g <- mk_grid()
  expect_error(occupancy_map(mk_track(151, -61), grid = g, method = "posterior"),
               "standard deviations")
})

test_that("occupancy_overlap is 1 for a map against itself", {
  g <- mk_grid()
  m <- occupancy_map(mk_track(151, -61, n = 30, sd = 3), grid = g,
                     method = "posterior")
  o <- occupancy_overlap(m, m)
  expect_equal(o$overlap, 1, tolerance = 1e-8)
  expect_equal(o$bhattacharyya, 1, tolerance = 1e-8)
  expect_equal(o$correlation, 1, tolerance = 1e-8)
})

test_that("occupancy_overlap is 0 for disjoint maps and ordered in between", {
  g <- mk_grid()
  a <- occupancy_map(mk_track(151, -61), grid = g, method = "point")
  b <- occupancy_map(mk_track(177, -45), grid = g, method = "point")
  near <- occupancy_map(mk_track(155, -61), grid = g, method = "point")
  expect_equal(occupancy_overlap(a, b)$overlap, 0, tolerance = 1e-12)
  # a diffuse version of `a` overlaps `a` more than a displaced point does
  a_diffuse <- occupancy_map(mk_track(151, -61, sd = 4), grid = g,
                             method = "posterior")
  expect_gt(occupancy_overlap(a, a_diffuse)$overlap,
            occupancy_overlap(a, near)$overlap)
})

test_that("bhattacharyya forgives diffuseness more than overlap does", {
  g <- mk_grid()
  truth   <- occupancy_map(mk_track(151, -61), grid = g, method = "point")
  diffuse <- occupancy_map(mk_track(151, -61, sd = 5), grid = g, method = "posterior")
  o <- occupancy_overlap(truth, diffuse)
  # a correctly centred but uncertain estimate keeps more Bhattacharyya than
  # shared mass -- that ordering is the reason both are reported
  expect_gt(o$bhattacharyya, o$overlap)
})

test_that("longitude convention is wrapped, not silently mis-indexed", {
  # the same track expressed in 0-360 and in -180..180 must give the same map
  g360 <- mk_grid(xmin = 140, xmax = 220)
  a <- occupancy_map(mk_track(201, -61), grid = g360, method = "point")
  b <- occupancy_map(mk_track(-159, -61), grid = g360, method = "point")
  expect_equal(terra::values(a), terra::values(b))
  expect_equal(which.max(terra::values(a)),
               terra::cellFromXY(g360, cbind(201, -61)))
})

test_that("masked cells are dropped and the rest still normalise", {
  g <- mk_grid()
  v <- terra::values(g); v[1:100] <- NA; terra::values(g) <- v
  m <- occupancy_map(mk_track(151, -61, n = 30, sd = 3), grid = g,
                     method = "posterior")
  expect_true(all(is.na(terra::values(m)[1:100])))
  expect_equal(sum(terra::values(m), na.rm = TRUE), 1, tolerance = 1e-8)
})

test_that("per_area rescales without changing where the mass is", {
  g <- mk_grid()
  tr <- mk_track(151, -61, n = 30, sd = 3)
  m  <- occupancy_map(tr, grid = g, method = "posterior")
  ma <- occupancy_map(tr, grid = g, method = "posterior", per_area = TRUE)
  a  <- terra::values(terra::cellSize(g, unit = "km"))[, 1]
  # per_area is exactly mass / area, cell by cell
  expect_equal(terra::values(ma)[, 1], terra::values(m)[, 1] / a,
               tolerance = 1e-10)
  expect_false(isTRUE(all.equal(sum(terra::values(ma), na.rm = TRUE), 1)))
})
