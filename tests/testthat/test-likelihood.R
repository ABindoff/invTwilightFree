test_that("light_log_likelihood matches hand-computed value", {
  # Spike-and-slab with asymmetric exponential spike:
  #   spike(obs, exp) = lambda * exp(-lambda * (exp - obs))      if obs <= exp
  #                   = lambda * exp(-lambda * 2 * (obs - exp))  if obs > exp
  #   mix = (1 - p) * spike + p * (1/max_light)
  #   loglik = sum(log(mix))

  obs      <- c(5.0, 3.0, 8.0)  # 8.0 > expected, triggers penalty branch
  expected <- c(6.0, 6.0, 6.0)
  lambda   <- 0.5
  max_light <- 10.0
  prob_slab  <- 0.1

  slab_d <- 1.0 / max_light

  spike <- ifelse(
    obs <= expected,
    lambda * exp(-lambda * (expected - obs)),
    lambda * exp(-lambda * 2.0 * (obs - expected))
  )
  mix <- (1 - prob_slab) * spike + prob_slab * slab_d
  expected_ll <- sum(log(mix))

  result <- light_log_likelihood(obs, expected, lambda, max_light, prob_slab)

  expect_equal(result, expected_ll, tolerance = 1e-10)
})

test_that("light_log_likelihood uses same factor-2 penalty as the engine (P5 unification)", {
  # After the P5 spike_density() refactor, the exported function and the engine
  # both use factor 2 on the obs > expected branch. Verify consistency.
  obs       <- c(8.0)
  expected  <- c(6.0)
  lambda    <- 0.5
  max_light <- 10.0
  prob_slab <- 0.0  # pure spike for clarity

  spike_factor2 <- lambda * exp(-lambda * 2.0 * (obs - expected))
  ll_factor2    <- log(spike_factor2)

  result <- light_log_likelihood(obs, expected, lambda, max_light, prob_slab)
  expect_equal(result, ll_factor2, tolerance = 1e-10,
               label = "exported function uses factor-2 over-bright penalty")
})
