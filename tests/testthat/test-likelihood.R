# Hand-computed checks of the spike-and-slab light likelihood. The spike is
# normalised over [0, max_light] by N(mu); without this the model is improper and
# biases the latitude posterior (caught by SBC). N(mu) is included below.

# Normalising constant of the asymmetric-exponential spike over [0, max_light].
spike_normaliser <- function(mu, lambda, max_light) {
  lo <- 1 - exp(-lambda * mu)
  hi <- 0.5 * (1 - exp(-2 * lambda * (max_light - mu)))
  pmax(lo + hi, 1e-12)
}

test_that("light_log_likelihood matches hand-computed value", {
  #   spike(obs, mu) = lambda * exp(-lambda * (mu - obs))      if obs <= mu
  #                  = lambda * exp(-lambda * 2 * (obs - mu))  if obs > mu
  #   spike_norm = spike / N(mu)
  #   mix = (1 - p) * spike_norm + p * (1/max_light)
  #   loglik = sum(log(mix))
  obs       <- c(5.0, 3.0, 8.0)  # 8.0 > expected, triggers penalty branch
  expected  <- c(6.0, 6.0, 6.0)
  lambda    <- 0.5
  max_light <- 10.0
  prob_slab <- 0.1

  slab_d <- 1.0 / max_light
  spike <- ifelse(
    obs <= expected,
    lambda * exp(-lambda * (expected - obs)),
    lambda * exp(-lambda * 2.0 * (obs - expected))
  )
  spike_norm <- spike / spike_normaliser(expected, lambda, max_light)
  mix <- (1 - prob_slab) * spike_norm + prob_slab * slab_d
  expected_ll <- sum(log(mix))

  result <- light_log_likelihood(obs, expected, lambda, max_light, prob_slab, 2)
  expect_equal(result, expected_ll, tolerance = 1e-10)
})

test_that("light_log_likelihood uses the factor-2 over-bright penalty, normalised", {
  # The obs > expected branch decays at twice the rate; the spike is normalised
  # by N(mu). With prob_slab = 0 the density is the normalised spike alone.
  obs       <- c(8.0)
  expected  <- c(6.0)
  lambda    <- 0.5
  max_light <- 10.0
  prob_slab <- 0.0  # pure spike for clarity

  spike_factor2 <- lambda * exp(-lambda * 2.0 * (obs - expected))
  ll <- log(spike_factor2 / spike_normaliser(expected, lambda, max_light))

  result <- light_log_likelihood(obs, expected, lambda, max_light, prob_slab, 2)
  expect_equal(result, ll, tolerance = 1e-10,
               label = "exported function uses the normalised factor-2 spike")
})
