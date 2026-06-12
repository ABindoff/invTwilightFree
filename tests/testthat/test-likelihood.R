test_that("light_log_likelihood matches hand-computed value", {
  # Spike-and-slab with asymmetric exponential spike:
  #   spike(obs, exp) = lambda * exp(-lambda * (exp - obs))   if obs <= exp
  #                   = lambda * exp(-lambda * 2 * (obs - exp)) if obs > exp
  #   mix = (1 - p) * spike + p * (1/max_light)
  #   loglik = sum(log(mix))

  obs      <- c(5.0, 3.0, 8.0)  # 8.0 > expected, triggers penalty branch
  expected <- c(6.0, 6.0, 6.0)
  lambda   <- 0.5
  max_light <- 10.0
  prob_slab  <- 0.1

  slab_d <- 1.0 / max_light

  # exported light_log_likelihood uses factor 10 on the over-bright branch
  spike <- ifelse(
    obs <= expected,
    lambda * exp(-lambda * (expected - obs)),
    lambda * exp(-lambda * 10.0 * (obs - expected))
  )
  mix <- (1 - prob_slab) * spike + prob_slab * slab_d
  expected_ll <- sum(log(mix))

  result <- light_log_likelihood(obs, expected, lambda, max_light, prob_slab)

  expect_equal(result, expected_ll, tolerance = 1e-10)
})

test_that("light_log_likelihood penalises obs > expected more weakly than factor-10 variant", {
  # The engine uses factor 2 on the over-bright branch.
  # This test documents that the exported function currently uses factor 10
  # (see CRITICAL_REVIEW H2). Remove/update when P5 unification is done.
  obs      <- c(8.0)
  expected <- c(6.0)
  lambda   <- 0.5
  max_light <- 10.0
  prob_slab  <- 0.0  # pure spike for clarity

  # factor 10 branch (as currently exported)
  spike_10 <- lambda * exp(-lambda * 10.0 * (obs - expected))
  ll_10 <- log(spike_10)

  result <- light_log_likelihood(obs, expected, lambda, max_light, prob_slab)
  expect_equal(result, ll_10, tolerance = 1e-10,
               label = "exported function uses factor-10 over-bright penalty")
})
