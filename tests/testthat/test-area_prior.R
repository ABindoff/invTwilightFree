test_that("area_prior returns log(cos(lat)) and down-weights poleward cells", {
  # deprecated in favour of TwilightFreeGrid(area_correction = TRUE); the warning
  # is asserted in test-anisotropic.R, so silence it here rather than let an
  # expected warning hide an unexpected one
  src <- suppressWarnings(area_prior())
  expect_s3_class(src, "tf_source")
  lat <- c(0, -60, 60, -30)
  expect_equal(src(lon = rep(0, 4), lat = lat), log(cos(lat * pi / 180)))
  v <- src(lon = c(0, 0), lat = c(-30, -60))
  expect_gt(v[1], v[2])                 # the equatorward (larger-area) cell is favoured
})

test_that("area_prior is longitude-independent and date-independent", {
  # deprecated in favour of TwilightFreeGrid(area_correction = TRUE); the warning
  # is asserted in test-anisotropic.R, so silence it here rather than let an
  # expected warning hide an unexpected one
  src <- suppressWarnings(area_prior())
  a <- src(lon = c(10, 200), lat = c(-45, -45))
  expect_equal(a[1], a[2])              # constant along a parallel
  expect_equal(src(0, -45, as.Date("2024-01-01")), src(0, -45, as.Date("2024-07-01")))
})
