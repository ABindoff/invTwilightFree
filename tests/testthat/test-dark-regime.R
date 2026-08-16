# The Kuo-Mallick darkness regime: a second slab weight and a sharper upper arm
# where the model expects darkness. `likelihood_params` of length 6 enables it;
# lengths 3 and 4 must behave exactly as before.

lp3 <- function(ml = 100) c(1 / (ml * 0.5), ml, 0.10)
lp6 <- function(ml = 100, frac = 0.35, ratio = 8, ps = 0.25)
  c(1 / (ml * 0.5), ml, 0.10, frac, ratio, ps)

test_that("length-3 params are untouched by the darkness code", {
  set.seed(1)
  t <- seq(as.POSIXct("2021-06-01", tz = "UTC"), by = "30 min", length.out = 800)
  z <- solar_zenith(as.numeric(t), rep(150, length(t)), rep(-45, length(t)))
  lig <- pmin(pmax(160 - 1.6 * z, 0), 100) * runif(length(z), 0.3, 1)
  cal <- c(160, 1.6)
  a <- eval_logpk_grid(150, -45, as.numeric(t), lig, cal, lp3(), 2)
  # a six-long vector whose dark weight equals the day weight and whose dark
  # ratio equals shade_ratio is the same model, so it must agree to the bit
  b <- eval_logpk_grid(150, -45, as.numeric(t), lig, cal,
                       c(lp3(), 0.35, 2, 0.10), 2)
  expect_equal(a, b, tolerance = 1e-12)
})

test_that("the darkness regime only touches observations where dark is expected", {
  set.seed(2)
  t <- seq(as.POSIXct("2021-06-01", tz = "UTC"), by = "30 min", length.out = 1200)
  z <- solar_zenith(as.numeric(t), rep(150, length(t)), rep(-45, length(t)))
  cal <- c(160, 1.6)                       # zero above zenith 100
  mu <- pmin(pmax(160 - 1.6 * z, 0), 100)
  lig <- mu * runif(length(mu), 0.3, 1)

  # daylight only: the two models must agree, because none of these observations
  # sits in the darkness regime
  day <- mu > 0.35 * 100
  expect_gt(sum(day), 100)
  a <- eval_logpk_grid(150, -45, as.numeric(t)[day], lig[day], cal, lp3(), 2)
  b <- eval_logpk_grid(150, -45, as.numeric(t)[day], lig[day], cal, lp6(), 2)
  expect_equal(a, b, tolerance = 1e-12)

  # night only: they must NOT agree
  night <- mu <= 0.05 * 100
  expect_gt(sum(night), 100)
  c3 <- eval_logpk_grid(150, -45, as.numeric(t)[night], lig[night], cal, lp3(), 2)
  c6 <- eval_logpk_grid(150, -45, as.numeric(t)[night], lig[night], cal, lp6(), 2)
  expect_false(isTRUE(all.equal(c3, c6)))
})

test_that("a bright reading at night has a bounded cost", {
  # What the indicator guarantees is a BOUND, not a discount. However bright a
  # false reading is, the slab can always explain it, so its marginal density
  # cannot fall below pi/max_light and its cost cannot exceed -log of that. The
  # darkness regime raises the floor where artificial light and moonlight
  # actually occur, without paying for it in daylight.
  #
  # Note what this does NOT say. One anomaly inside an otherwise quiet night
  # stands out MORE under the darkness regime, not less, because the sharper arm
  # also makes the quiet readings far more informative. That is the intended
  # behaviour: the anomaly is capped while the honest observations do the work.
  t <- as.POSIXct("2021-06-01 08:00", tz = "UTC") + (0:19) * 1800
  cal <- c(160, 1.6)
  mu <- pmin(pmax(160 - 1.6 * solar_zenith(as.numeric(t), rep(150, 20),
                                           rep(-45, 20)), 0), 100)
  expect_equal(max(mu), 0)

  # one observation at a time, so the density is read directly
  den1 <- function(lp, y)
    exp(eval_logpk_grid(150, -45, as.numeric(t[10]), y, cal, lp, 2))

  for (y in c(20, 60, 95, 100)) {
    expect_gte(den1(lp6(), y), 0.25 / 100)          # floored by the dark slab
    expect_gte(den1(lp3(), y), 0.10 / 100)          # and by the day slab
  }
  # There is a crossover, and it is the point of the design. Below it the
  # darkness regime is STRICTER, because a moderate excursion above the dark
  # envelope is informative about position and should be paid for. Above it the
  # regime is more forgiving, because an extreme reading is contamination and
  # nothing is learned by fighting it.
  expect_lt(den1(lp6(), 60), den1(lp3(), 60))       # stricter where it informs
  expect_gt(den1(lp6(), 100), den1(lp3(), 100))     # forgiving where it does not
  cross <- vapply(seq(60, 100, by = 1),
                  function(y) den1(lp6(), y) - den1(lp3(), y), 0)
  expect_equal(sum(diff(sign(cross)) != 0), 1)      # exactly one crossing
  # the worst case really is the slab, so cost is capped
  expect_lt(-log(den1(lp6(), 100)), -log(0.25 / 100) + 1e-6)

  # a larger dark slab absorbs more; a sharper dark arm does not change the cap
  expect_gt(den1(lp6(ps = 0.60), 95), den1(lp6(ps = 0.10), 95))
  expect_equal(den1(lp6(ratio = 40), 95), den1(lp6(ratio = 8), 95),
               tolerance = 1e-3)

  # and the sharper arm is what buys the signal: a reading consistent with
  # darkness becomes more probable, which is where day length comes from
  expect_gt(den1(lp6(ratio = 40), 1), den1(lp6(ratio = 2), 1))
})

test_that("the sharper arm still discriminates dark from not-dark", {
  # Absorbing anomalies must not cost the day-length signal. Readings consistent
  # with darkness should prefer a position where the model expects darkness over
  # one where it expects twilight. Same times, same longitude, opposite
  # hemisphere: June is night at 45 S and twilight at 45 N.
  t <- as.POSIXct("2021-06-01 08:00", tz = "UTC") + (0:19) * 1800
  cal <- c(160, 1.6)
  mu_s <- pmin(pmax(160 - 1.6 * solar_zenith(as.numeric(t), rep(150, 20),
                                             rep(-45, 20)), 0), 100)
  mu_n <- pmin(pmax(160 - 1.6 * solar_zenith(as.numeric(t), rep(150, 20),
                                             rep(45, 20)), 0), 100)
  expect_equal(max(mu_s), 0)
  expect_gt(max(mu_n), 30)

  lig <- rep(2, length(t))                 # consistent with darkness
  expect_gt(eval_logpk_grid(150, -45, as.numeric(t), lig, cal, lp6(), 2),
            eval_logpk_grid(150,  45, as.numeric(t), lig, cal, lp6(), 2))

  # a sharper dark arm must not blunt that preference
  gap <- function(lp) eval_logpk_grid(150, -45, as.numeric(t), lig, cal, lp, 2) -
                      eval_logpk_grid(150,  45, as.numeric(t), lig, cal, lp, 2)
  expect_gt(gap(lp6(ratio = 40)), 0)
  expect_gt(gap(lp6(ratio = 8)), 0)
})

test_that("tempering scales the light likelihood and nothing else", {
  # The seventh parameter multiplies each window's light log-likelihood, standing
  # in for an effective sample size below the observation count. It must be
  # exactly linear in the log, and inert by default.
  set.seed(9)
  t <- seq(as.POSIXct("2021-06-01", tz = "UTC"), by = "30 min", length.out = 600)
  z <- solar_zenith(as.numeric(t), rep(150, length(t)), rep(-45, length(t)))
  cal <- c(160, 1.6)
  lig <- pmin(pmax(160 - 1.6 * z, 0), 100) * runif(length(z), 0.3, 1)
  f <- function(lp) eval_logpk_grid(150, -45, as.numeric(t), lig, cal, lp, 2)

  base <- f(lp6())
  expect_equal(f(c(lp6(), 1.0)), base, tolerance = 1e-12)   # default is inert
  expect_equal(f(c(lp6(), 0.5)), 0.5 * base, tolerance = 1e-9)
  expect_equal(f(c(lp6(), 0.25)), 0.25 * base, tolerance = 1e-9)

  # a non-positive value is ignored rather than silently zeroing the likelihood
  expect_equal(f(c(lp6(), 0)), base, tolerance = 1e-12)
  expect_equal(f(c(lp6(), -1)), base, tolerance = 1e-12)

  # and it must reach the length-3 path too
  expect_equal(f(c(lp3(), 0.35, 2, 0.10, 0.5)), 0.5 * f(lp3()), tolerance = 1e-9)
})

test_that("tempering widens the posterior without moving its centre", {
  # The point is calibration, not relocation: flattening every window equally
  # must leave the emission's preferred latitude alone while making it less sure.
  t <- seq(as.POSIXct("2021-09-15", tz = "UTC"), by = "30 min", length.out = 48)
  cal <- c(160, 1.6)
  lats <- seq(20, 70, by = 0.5)
  z <- solar_zenith(as.numeric(t), rep(200, length(t)), rep(45, length(t)))
  lig <- pmin(pmax(160 - 1.6 * z, 0), 100)

  prof <- function(tm) {
    ll <- eval_logpk_grid(rep(200, length(lats)), lats, as.numeric(t), lig, cal,
                          c(lp6(), tm), 2)
    w <- exp(ll - max(ll)); w <- w / sum(w)
    c(mean = sum(w * lats), sd = sqrt(sum(w * (lats - sum(w * lats))^2)))
  }
  a <- prof(1); b <- prof(0.3)
  expect_equal(unname(a["mean"]), unname(b["mean"]), tolerance = 0.5)
  expect_gt(b["sd"], a["sd"])
})
