# Class-aware hemisphere selection against the package.
#
# Now runnable: run_grid_hmm returns the forward-pass log marginal likelihood
# (model evidence) as `log_z`, surfaced on the fit as `fit$log_z`. The hemisphere
# Bayes factor is the difference in log_z between fits under opposite hard priors.

library(invTwilightFree)

#' Select the hemisphere branch by model evidence.
#'
#' Fits the grid HMM twice under hard hemisphere priors (south vs north) and
#' compares the log marginal likelihood. Returns the chosen class, the log Bayes
#' factor (south minus north; > 0 favours south), and the winning fit.
class_select <- function(date_time, light, grid,
                         start_lon, abs_start_lat, ...) {
  fit_S <- TwilightFreeGrid(
    date_time, light, grid,
    start_lat = -abs(abs_start_lat), start_lon = start_lon,
    terms = list(location_term("hemi-S",
      source = hemisphere_prior(function(d) "S", softness = 0),  # hard cut at equator
      rule = identity_rule())), ...)

  fit_N <- TwilightFreeGrid(
    date_time, light, grid,
    start_lat = abs(abs_start_lat), start_lon = start_lon,
    terms = list(location_term("hemi-N",
      source = hemisphere_prior(function(d) "N", softness = 0),
      rule = identity_rule())), ...)

  lbf <- fit_S$log_z - fit_N$log_z
  list(
    class = if (lbf >= 0) "S" else "N",
    log_bayes_factor = lbf,
    log_z = c(S = fit_S$log_z, N = fit_N$log_z),
    fit = if (lbf >= 0) fit_S else fit_N
  )
}

# Usage:
# grid <- makeGrid(lon = c(130, 170), lat = c(-65, -35), cell.size = 1, mask = "sea")
# res  <- class_select(track$time, track$light, grid,
#                      start_lon = 158, abs_start_lat = 54, diffusion = 80)
# res$class            # "S"
# res$log_bayes_factor # large positive = decisive
#
# Decision-curve experiment (the useful design result): truncate the deployment
# to D days around an equinox and watch res$log_bayes_factor grow with the amount
# of non-equinox season included.
