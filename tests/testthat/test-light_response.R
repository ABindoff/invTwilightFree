test_that("a logger linear in irradiance is identified as narrow", {
  # The clamped-linear response this package started from: light falls from
  # saturation to zero across civil twilight, an 11 degree transition.
  t <- seq(as.POSIXct("2021-06-01", tz = "UTC"), by = "5 min", length.out = 10 * 288)
  z <- solar_zenith(as.numeric(t), rep(150, length(t)), rep(-45, length(t)))
  lig <- pmin(pmax(558.5 - 5.818 * z, 0), 64)

  r <- fit_light_response(t, lig, 150, -45)
  expect_s3_class(r, "tf_light_response")
  expect_lt(r$width_deg, 12)                        # a log channel gives tens
  expect_equal(r$zero_at, 558.5 / 5.818, tolerance = 0.02)

  # The returned slope is steeper than the 5.818 the data were built with, by
  # about the 1.374 that fitting a logistic to a clamped line costs: the
  # logistic bends at both ends, so its gradient at half amplitude exceeds the
  # line's. Asserted rather than tolerated, so that a change to how the tangent
  # is taken shows up here instead of silently moving every fitted track.
  expect_equal(r$slope / 5.818, 1.374, tolerance = 0.05)
  expect_equal(4 * r$slope * r$scale, r$max_light, tolerance = 1e-8)
})

test_that("a log-scaled channel gives a wide transition", {
  # The Mk9 case: a smooth decline over tens of degrees onto a non-zero dark
  # reading. Fitting a narrow clamped line to this is the error that motivated
  # the function, so the fit has to be able to tell the two apart.
  t <- seq(as.POSIXct("2021-06-01", tz = "UTC"), by = "5 min", length.out = 10 * 288)
  z <- solar_zenith(as.numeric(t), rep(150, length(t)), rep(-45, length(t)))
  lig <- 20 + 140 / (1 + exp((z - 92) / 6))

  r <- fit_light_response(t, lig, 150, -45)
  expect_gt(r$width_deg, 20)
  expect_equal(r$z50, 92, tolerance = 0.05)
  expect_equal(r$scale, 6, tolerance = 0.1)
  expect_equal(r$floor, 20, tolerance = 1)
})

test_that("the envelope ignores shading, which only ever subtracts", {
  # Cloud, shading and depth can only reduce the reading, so the clear-sky curve
  # is the top of the point cloud. Randomly attenuating three quarters of the
  # record must not move the fitted geometry.
  t <- seq(as.POSIXct("2021-06-01", tz = "UTC"), by = "5 min", length.out = 12 * 288)
  z <- solar_zenith(as.numeric(t), rep(150, length(t)), rep(-45, length(t)))
  clear <- 20 + 140 / (1 + exp((z - 92) / 6))
  set.seed(42)
  shaded <- clear
  hit <- sample(length(clear), floor(0.75 * length(clear)))
  shaded[hit] <- 20 + (clear[hit] - 20) * runif(length(hit), 0, 1)

  r <- fit_light_response(t, shaded, 150, -45)
  expect_equal(r$z50, 92, tolerance = 1)
  expect_equal(r$scale, 6, tolerance = 1.5)
})

test_that("scale_light sets the intensity scale without touching the geometry", {
  # The two halves of the fit want different windows: geometry needs known
  # positions, intensity needs a record whose shading matches the one being
  # tracked. Supplying `scale_light` must move the second and leave the first
  # exactly where it was.
  t <- seq(as.POSIXct("2021-06-01", tz = "UTC"), by = "5 min", length.out = 10 * 288)
  z <- solar_zenith(as.numeric(t), rep(150, length(t)), rep(-45, length(t)))
  lig <- 20 + 140 / (1 + exp((z - 92) / 6))

  base <- fit_light_response(t, lig, 150, -45)
  # a record with deep dives in it: darker at the bottom, so a wider range
  dived <- c(lig, rep(0, length(lig) %/% 4))
  alt <- fit_light_response(t, lig, 150, -45, scale_light = dived)

  expect_identical(alt$z50, base$z50)
  expect_identical(alt$scale, base$scale)
  expect_lt(alt$baseline, base$baseline)
  expect_gt(alt$max_light, base$max_light)
  # the tangent must still reach the top of the range the light is expressed on
  expect_equal(4 * alt$slope * alt$scale, alt$max_light, tolerance = 1e-8)

  expect_null(fit_light_response(t, lig, 150, -45, scale_light = numeric(0)))
})

test_that("too little record, or a flat one, returns NULL rather than nonsense", {
  t <- seq(as.POSIXct("2021-06-01", tz = "UTC"), by = "5 min", length.out = 100)
  expect_null(fit_light_response(t, rep(50, 100), 150, -45))

  t2 <- seq(as.POSIXct("2021-06-01", tz = "UTC"), by = "5 min", length.out = 2000)
  expect_null(fit_light_response(t2, rep(50, 2000), 150, -45))
})

test_that("pooling replaces geometry with the median and keeps each scale", {
  t <- seq(as.POSIXct("2021-06-01", tz = "UTC"), by = "5 min", length.out = 10 * 288)
  z <- solar_zenith(as.numeric(t), rep(150, length(t)), rep(-45, length(t)))
  mk <- function(z50, s, gain) {
    fit_light_response(t, gain * (20 + 140 / (1 + exp((z - z50) / s))), 150, -45)
  }
  rs <- list(a = mk(90, 5, 1), b = mk(92, 6, 1), c = mk(94, 7, 2))

  p <- pool_light_responses(rs)
  g <- attr(p, "pooled")
  expect_equal(unname(g[["z50"]]), 92, tolerance = 0.2)
  expect_equal(unname(g[["n"]]), 3)
  expect_identical(names(p), names(rs))

  # every tag ends up on the same geometry ...
  expect_equal(unname(vapply(p, function(r) r$z50, 0)), rep(unname(g[["z50"]]), 3))
  # ... and keeps its own intensity scale, so the doubled tag stays doubled
  expect_equal(vapply(p, function(r) r$max_light, 0),
               vapply(rs, function(r) r$max_light, 0))
  expect_gt(p$c$max_light, 1.5 * p$a$max_light)
  # derived fields must follow the pooled geometry, not the tag's own
  expect_equal(p$a$width_deg, 4.394 * g[["scale"]])
  expect_equal(4 * p$a$slope * p$a$scale, p$a$max_light, tolerance = 1e-8)
})

test_that("pooling tolerates gaps and refuses an empty geometry set", {
  t <- seq(as.POSIXct("2021-06-01", tz = "UTC"), by = "5 min", length.out = 10 * 288)
  z <- solar_zenith(as.numeric(t), rep(150, length(t)), rep(-45, length(t)))
  r <- fit_light_response(t, 20 + 140 / (1 + exp((z - 92) / 6)), 150, -45)

  # An animal that departs too fast to be calibrated contributes no geometry but
  # still needs a response: it takes the pool's, and keeps its own scale.
  p <- pool_light_responses(list(a = r, b = NULL), scale = list(a = r, b = r))
  expect_equal(p$b$z50, r$z50, tolerance = 1e-8)
  expect_equal(attr(p, "pooled")[["n"]], 1)

  # A tag with no scale of its own cannot be given a response at all.
  p2 <- pool_light_responses(list(a = r), scale = list(a = r, b = NULL))
  expect_null(p2$b)

  expect_error(pool_light_responses(list(NULL, NULL)), "no fitted responses")
  expect_error(pool_light_responses(list(a = 1)), "tf_light_response")
})
