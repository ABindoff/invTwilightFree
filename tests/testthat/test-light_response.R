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

test_that("calibration_logistic is on the engine's scale, not the tag's", {
  # The engine is handed `pmax(0, light - baseline)` and evaluates
  # `floor + amp / (1 + exp((z - z50) / scale))`. So the floor in
  # `calibration_logistic` has to be measured from `baseline`, not from zero.
  # Return the raw dark reading instead and the expected curve sits
  # `dark_level - baseline` too high at every zenith, which on real Mk9 records
  # was 11% of the response range and made the four-parameter form look far
  # worse than the tangent it was being compared against.
  t <- seq(as.POSIXct("2021-06-01", tz = "UTC"), by = "5 min", length.out = 10 * 288)
  z <- solar_zenith(as.numeric(t), rep(150, length(t)), rep(-45, length(t)))
  lig <- 20 + 140 / (1 + exp((z - 92) / 6))

  r <- fit_light_response(t, lig, 150, -45)
  expect_equal(r$calibration_logistic[1], r$dark_level - r$baseline,
               tolerance = 1e-8)
  expect_equal(r$calibration_logistic[2:4], c(r$amp, r$z50, r$scale),
               tolerance = 1e-8)

  # The curve and the baselined light must agree at both ends of the range.
  mu <- function(zz) r$calibration_logistic[1] +
    r$calibration_logistic[2] / (1 + exp((zz - r$calibration_logistic[3]) /
                                           r$calibration_logistic[4]))
  obs <- pmax(0, lig - r$baseline)
  expect_equal(mu(140), min(obs), tolerance = 1)          # night
  expect_equal(mu(20), max(obs), tolerance = 2)           # full day

  # Never negative: the observations are clamped at zero and the spike
  # normaliser is only defined over [0, max_light].
  bright <- 500 + 140 / (1 + exp((z - 92) / 6))
  expect_gte(fit_light_response(t, bright, 150, -45)$calibration_logistic[1], 0)

  # Pooling must not reintroduce the raw scale.
  p <- pool_light_responses(list(a = r, b = r))
  expect_equal(p$a$calibration_logistic[1], r$dark_level - r$baseline,
               tolerance = 1e-8)
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

# The R-side twin of the Rust `expected_light()` is internal, so name it once
# here rather than reaching into the namespace at every call site.
expected_light_at <- invTwilightFree:::.tf_expected_light

test_that("the lookup-table response evaluates as measured", {
  # A channel that neither parametric form can express: a daytime plateau, a
  # narrow collapse, and a hard non-zero floor. The clamped line has no floor
  # parameter and the logistic cannot have a narrow transition and a wide
  # shoulder at once, so this is the case the table exists for.
  t <- seq(as.POSIXct("2021-06-01", tz = "UTC"), by = "5 min", length.out = 20 * 288)
  z <- solar_zenith(as.numeric(t), rep(150, length(t)), rep(-45, length(t)))
  truth <- function(zz) 40 + 120 / (1 + exp((zz - 92) / 4))
  # A diving animal, so the 5th percentile sits near zero while the clear-sky
  # reading at night is still 40. That gap is the whole point: on the baselined
  # scale the response has a floor a quarter of the way up the range, and it is
  # what the tangent cannot express. Without dives the baseline IS the night
  # reading and the floor is zero, which tests nothing.
  set.seed(7)
  lig <- truth(z)
  hit <- sample(length(lig), floor(0.3 * length(lig)))
  lig[hit] <- lig[hit] * runif(length(hit), 0, 0.05)
  r <- fit_light_response(t, lig, 150, -45)
  expect_lt(r$baseline, 10)

  tab <- light_response_table(r, from = 30, to = 140, by = 1)
  expect_lt(tab[1], 0)                      # the sentinel the engines dispatch on
  expect_equal(tab[2:3], c(30, 1))
  expect_length(tab, 3 + length(seq(30, 140, 1)))
  expect_equal(attr(tab, "n_tags"), 1L)

  # It must reproduce the measured curve, on the baselined scale. The envelope
  # is a 95th percentile in one-degree bins, so it sits a little under the truth
  # through the steep part; the tolerance is that binning, not slack.
  zz <- seq(40, 130, by = 2)
  got <- expected_light_at(zz, tab, max_light = r$max_light)
  want <- pmin(pmax(truth(zz) - r$baseline, 0), r$max_light)
  expect_lt(max(abs(got - want)), 8)
  expect_lt(sqrt(mean((got - want)^2)), 3)

  # ... including the floor, which the tangent forces to zero.
  expect_gt(expected_light_at(130, tab, r$max_light), 20)
  expect_equal(expected_light_at(130, r$calibration, r$max_light), 0)

  # Held constant outside the grid, never extrapolated off the end.
  expect_equal(expected_light_at(0, tab, r$max_light),
               expected_light_at(30, tab, r$max_light))
  expect_equal(expected_light_at(179, tab, r$max_light),
               expected_light_at(140, tab, r$max_light))

  # The R twin and the Rust engine must agree, or the clock estimator profiles a
  # different likelihood from the one the engines evaluate. `eval_logpk_grid()`
  # is the exported route into the Rust `expected_light()`: hold position and
  # light fixed and the log-likelihood is a function of the response alone, so
  # two calibrations agree there iff they agree in the engine.
  lp <- c(1 / (r$max_light * 0.5), r$max_light, 0.05)
  fine <- seq(30, 140, by = 0.5)
  tab_fine <- c(-1, 30, 0.5,
                expected_light_at(fine, tab, r$max_light))
  expect_equal(eval_logpk_grid(150, -45, as.numeric(t[1:500]), lig[1:500],
                               tab, lp, 2),
               eval_logpk_grid(150, -45, as.numeric(t[1:500]), lig[1:500],
                               tab_fine, lp, 2),
               tolerance = 1e-6)
})

test_that("pooling a table across tags is on each tag's own scale", {
  t <- seq(as.POSIXct("2021-06-01", tz = "UTC"), by = "5 min", length.out = 20 * 288)
  z <- solar_zenith(as.numeric(t), rep(150, length(t)), rep(-45, length(t)))
  a <- fit_light_response(t, 40 + 120 / (1 + exp((z - 92) / 2.5)), 150, -45)
  # same shape, twice as bright: pooling must not let it dominate the median
  b <- fit_light_response(t, 80 + 240 / (1 + exp((z - 92) / 2.5)), 150, -45)

  ta <- light_response_table(a); tb <- light_response_table(b)
  tp <- light_response_table(list(a, b))
  expect_equal(attr(tp, "n_tags"), 2L)
  mid <- function(v) v[3 + which(seq(30, 140, 1) == 92)]
  expect_true(mid(tp) > mid(ta) && mid(tp) < mid(tb))
  expect_equal(length(tp), length(ta))
})

test_that("flat_outside holds the table constant beyond the transition", {
  # The clamped-linear response owes its accuracy to being FLAT outside the
  # transition, not to describing the channel well. `flat_outside` gives the
  # table the same property while keeping the measured shape where the signal
  # is, and reuses the window the tangent already implies, so it adds no tuning
  # constant of its own.
  t <- seq(as.POSIXct("2021-06-01", tz = "UTC"), by = "5 min", length.out = 20 * 288)
  z <- solar_zenith(as.numeric(t), rep(150, length(t)), rep(-45, length(t)))
  set.seed(11)
  lig <- 40 + 120 / (1 + exp((z - 92) / 4))
  hit <- sample(length(lig), floor(0.3 * length(lig)))
  lig[hit] <- lig[hit] * runif(length(hit), 0, 0.05)
  r <- fit_light_response(t, lig, 150, -45)

  win <- c(r$saturate_at, r$zero_at)
  plain <- light_response_table(r, from = 30, to = 140, by = 1)
  flat  <- light_response_table(r, from = 30, to = 140, by = 1, flat_outside = win)
  expect_equal(attr(flat, "flat_outside"), win)

  grid <- seq(30, 140, 1)
  yp <- plain[-(1:3)]; yf <- flat[-(1:3)]
  # constant outside, and strictly so
  expect_equal(length(unique(round(yf[grid <= win[1]], 9))), 1L)
  expect_equal(length(unique(round(yf[grid >= win[2]], 9))), 1L)
  # the plain table is NOT constant there: that is the whole point
  expect_gt(diff(range(yp[grid <= win[1]])), 1)
  # untouched inside the window
  inside <- grid > win[1] & grid < win[2]
  expect_equal(yf[inside], yp[inside])
  # and the floor survives: flattening must not send night back to zero
  expect_gt(yf[length(yf)], 20)

  expect_error(light_response_table(r, flat_outside = 90), "lo < hi")
  expect_error(light_response_table(r, flat_outside = c(110, 80)), "lo < hi")
})

test_that("pooling a table rescales to each tag's own range", {
  # Tags of the same model differ in gain. A shared SHAPE must not impose a
  # shared gain: pooled in absolute units, a dim tag gets a curve that tops out
  # above its range and saturates, and a bright one gets a curve that never
  # reaches its brightest clear-sky observations.
  t <- seq(as.POSIXct("2021-06-01", tz = "UTC"), by = "5 min", length.out = 20 * 288)
  z <- solar_zenith(as.numeric(t), rep(150, length(t)), rep(-45, length(t)))
  set.seed(3)
  mk <- function(gain) {
    v <- gain * (40 + 120 / (1 + exp((z - 92) / 4)))
    hit <- sample(length(v), floor(0.3 * length(v)))
    v[hit] <- v[hit] * runif(length(hit), 0, 0.05)
    fit_light_response(t, v, 150, -45)
  }
  dim_tag <- mk(1); bright <- mk(2)
  expect_gt(bright$max_light, 1.5 * dim_tag$max_light)

  fits <- list(dim_tag, bright)
  td <- light_response_table(fits, max_light = dim_tag$max_light)
  tb <- light_response_table(fits, max_light = bright$max_light)
  expect_equal(attr(td, "max_light"), dim_tag$max_light)

  peak <- function(v) max(v[-(1:3)])
  # each table lands on its own tag's range ...
  expect_equal(peak(td) / dim_tag$max_light, peak(tb) / bright$max_light,
               tolerance = 1e-8)
  expect_gt(peak(tb), 1.5 * peak(td))
  # ... and neither overshoots the range the engine will clamp at
  expect_lt(peak(td), 1.05 * dim_tag$max_light)
  expect_lt(peak(tb), 1.05 * bright$max_light)

  # the shape is identical once scaled out
  expect_equal(td[-(1:3)] / dim_tag$max_light, tb[-(1:3)] / bright$max_light,
               tolerance = 1e-8)

  # pooling in absolute units is what this replaces: one curve for both, so the
  # dim tag's table overshoots its own range
  abs_tab <- light_response_table(fits)
  expect_gt(peak(abs_tab), 1.15 * dim_tag$max_light)

  expect_error(light_response_table(fits, max_light = -1), "positive")
  expect_error(light_response_table(fits, max_light = c(1, 2)), "single")
})
