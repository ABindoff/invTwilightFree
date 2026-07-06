# The spike-and-slab observation model must be a proper density: integrating
# over the tag's light range [0, max_light] gives 1. This guards the spike
# normalisation that SBC showed was missing (an un-normalised spike biased the
# latitude posterior). Runs after the Rust core is rebuilt.

test_that("the normalised spike-and-slab density integrates to 1 across mu", {
  lambda <- 0.5; maxl <- 64; pslab <- 0.05
  Lg <- seq(0, maxl, length.out = 4001)
  dens_integral <- function(mu) {
    d <- exp(vapply(Lg, function(L)
      light_log_likelihood(L, mu, lambda, maxl, pslab), numeric(1)))
    sum((d[-1] + d[-length(d)]) / 2 * diff(Lg))         # trapezoid
  }
  for (mu in c(5, 20, 40, 60)) {                         # interior and near both clamps
    expect_equal(dens_integral(mu), 1, tolerance = 0.01,
                 label = sprintf("integral at mu=%g", mu))
  }
})

test_that("the density is a proper mixture for the extreme expected values too", {
  lambda <- 0.5; maxl <- 64; pslab <- 0.1
  Lg <- seq(0, maxl, length.out = 4001)
  for (mu in c(0, maxl)) {                               # fully dark / fully lit cells
    d <- exp(vapply(Lg, function(L)
      light_log_likelihood(L, mu, lambda, maxl, pslab), numeric(1)))
    integral <- sum((d[-1] + d[-length(d)]) / 2 * diff(Lg))
    expect_equal(integral, 1, tolerance = 0.02)
  }
})
