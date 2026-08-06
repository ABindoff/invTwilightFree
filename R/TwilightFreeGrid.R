#' Reconstruct animal tracks using SGAT::essie and continuous likelihood
#' 
#' @param date_time A vector of POSIXct dates
#' @param light A vector of light observations
#' @param grid A \code{SpatRaster} defining the grid (from \link{makeGrid})
#' @param start_lat Initial latitude (Decimal Degrees)
#' @param start_lon Initial longitude (Decimal Degrees)
#' @param end_lat Final latitude (optional, defaults to NA)
#' @param end_lon Final longitude (optional, defaults to NA)
#' @param fixed Optional data frame of fixed locations with columns `time` (POSIXct), `lat`, and `lon`.
#' @param step_hours Time step for the HMM knots in hours (default 12.0)
#' @param diffusion Movement diffusion in km/sqrt(day)
#' @param trans_prob Optional transition probability matrix (flattened, row-major) for behavioral states. Defaults to 0.9 diagonal if multiple diffusions are provided.
#' @param calibration Calibration parameters c(intercept, slope)
#' @param likelihood_params Likelihood parameters c(lambda, max_light, prob_slab)
#' @param overcount Effective over-counting of light observations within a knot
#'   (default 1). Light half an hour apart shares cloud, sea state and
#'   behaviour, so n observations carry roughly n/overcount observations' worth
#'   of information; the per-knot log-likelihood is divided by this, widening
#'   the posterior by about sqrt(overcount) without moving the point estimate.
#'   Use it when reported intervals are too narrow against known positions.
#' @param shade_ratio Ratio of the spike's upper-arm decay rate to its shading
#'   (lower-arm) rate. The default `2` reproduces the historical fixed ratio.
#'   The shading arm decides how cheaply the model can explain light far below
#'   the clear-sky expectation, so a value below 1 suits a continuously diving
#'   animal, whose record is mostly attenuated; a larger value makes shading
#'   more surprising.
#' @param terms Optional list of `location_term()` objects (priors, SST,
#'   bathymetry, masks). Each contributes an additive log-likelihood over grid
#'   cells that is combined with the light likelihood. Defaults to `list()`
#'   (light only), which reproduces the original behaviour exactly.
#' @param calibrate If `TRUE` and both endpoint locations are known, the tag clock
#'   is calibrated from the known deployment and retrieval fixes (via
#'   [calibrate_clock_from_endpoints()]) and the observation times are corrected
#'   before fitting. The fitted clock model is returned as `$clock`. Default
#'   `FALSE` (no correction).
#' @importFrom stats quantile lm coef
#' @export
TwilightFreeGrid <- function(date_time, light, grid,
                             start_lat = NA_real_, start_lon = NA_real_,
                             end_lat = NA_real_, end_lon = NA_real_,
                             fixed = NULL,
                             step_hours = 12.0,
                             diffusion = 50,
                             trans_prob = NULL,
                             calibration = NULL,
                             likelihood_params = NULL,
                             shade_ratio = 2,
                             overcount = 1,
                             terms = list(),
                             calibrate = FALSE) {

  if(!inherits(date_time, "POSIXct")) stop("date_time must be POSIXct")

  # Ensure sorted
  ord <- order(date_time)
  date_time <- date_time[ord]
  light <- light[ord]

  # Optional clock calibration from the known endpoints (engine-agnostic): correct
  # the timestamps before auto-calibration and knot construction.
  clock <- NULL
  if (isTRUE(calibrate)) {
    clock <- calibrate_clock_from_endpoints(date_time, light, start_lon, start_lat,
                                            end_lon, end_lat, calibration, likelihood_params)
    if (!is.null(clock)) date_time <- clock$correct(date_time)
  }
  
  if (length(diffusion) > 1 && is.null(trans_prob)) {
    # Default to a "sticky" behavior: 90% chance to stay in current state
    k <- length(diffusion)
    m <- matrix((1.0 - 0.9) / (k - 1), nrow = k, ncol = k)
    diag(m) <- 0.9
    trans_prob <- as.numeric(t(m)) # Flatten to row-major
  } else if (is.null(trans_prob)) {
    trans_prob <- c(1.0)
  }

  # Auto-Calibration
  if (is.null(calibration) || is.null(likelihood_params)) {
    min_l <- quantile(light, 0.05, na.rm = TRUE)
    max_l <- quantile(light, 0.95, na.rm = TRUE)
    
    light_shifted <- pmax(0, light - min_l)
    max_shifted <- max_l - min_l
    
    if (is.null(calibration)) {
      # Data-driven calibration using first 3 days at start location if available
      if (!is.na(start_lat) && !is.na(start_lon)) {
        cal_idx <- date_time < (date_time[1] + 3*24*3600)
        cal_times <- date_time[cal_idx]
        cal_light <- light_shifted[cal_idx]
        
        # Calculate true solar zeniths
        cal_zenith <- solar_zenith(as.numeric(cal_times), 
                                   rep(start_lon, length(cal_times)), 
                                   rep(start_lat, length(cal_times)))
        
        # Fit linear model for transitions (zenith between 85 and 100, light > 0 and < max)
        trans_idx <- which(cal_light > 0 & cal_light < max_shifted * 0.95 & cal_zenith > 85 & cal_zenith < 100)
        if (length(trans_idx) > 10) {
          fit_cal <- lm(cal_light[trans_idx] ~ cal_zenith[trans_idx])
          intercept_est <- coef(fit_cal)[1]
          slope_est <- -coef(fit_cal)[2] 
          if (is.na(slope_est) || slope_est <= 0) slope_est <- max_shifted / (96 - 85)
        } else {
          slope_est <- max_shifted / (96 - 85)
          intercept_est <- slope_est * 96
        }
        calibration <- c(intercept_est, slope_est)
      } else {
        # Fallback to arbitrary if no start location
        calibration <- c(max_shifted * 1.1, max_shifted / (96 - 85))
      }
    }
    
    if (is.null(likelihood_params)) {
      likelihood_params <- c(1.0 / (max_shifted * 0.5), max_shifted, 0.10)
    }
    process_light <- light_shifted
  } else {
    process_light <- light
  }
  
  unix_times <- as.numeric(date_time)
  
  # Define Knots
  t_start <- unix_times[1]
  t_end <- unix_times[length(unix_times)]
  k_steps <- ceiling((t_end - t_start) / (step_hours * 3600)) + 1
  t_step <- if (k_steps > 1) (t_end - t_start) / (k_steps - 1) else 0
  
  knot_times <- t_start + (0:(k_steps-1)) * t_step
  time <- as.POSIXct(knot_times, origin="1970-01-01", tz="UTC")
  
  # Initialize x0 and fixed vectors
  x0 <- matrix(0.0, nrow = k_steps, ncol = 2)
  fixed_vec <- rep(FALSE, k_steps)
  
  # Handle start/end convenience params
  if (!is.na(start_lat) && !is.na(start_lon)) {
    x0[1, ] <- c(start_lon, start_lat)
    fixed_vec[1] <- TRUE
  }
  if (!is.na(end_lat) && !is.na(end_lon)) {
    x0[k_steps, ] <- c(end_lon, end_lat)
    fixed_vec[k_steps] <- TRUE
  }
  
  # Handle any other fixed locations
  if (!is.null(fixed)) {
    if (!inherits(fixed, "data.frame") || !all(c("time", "lat", "lon") %in% names(fixed))) {
      stop("fixed must be a data frame with columns time, lat, lon")
    }
    for (i in 1:nrow(fixed)) {
      f_time <- as.numeric(fixed$time[i])
      # Find nearest knot
      k_idx <- which.min(abs(knot_times - f_time))
      x0[k_idx, ] <- c(fixed$lon[i], fixed$lat[i])
      fixed_vec[k_idx] <- TRUE
    }
  }
  
  # terra::crds() excludes NA cells by default (na.rm=TRUE), so lon/lat already
  # contain only the valid (non-masked) grid cell centres.
  lon_vec <- terra::crds(grid)[, 1]
  lat_vec <- terra::crds(grid)[, 2]
  
  # Ensure valid data
  valid_obs <- !is.na(process_light)
  obs_light_clean <- process_light[valid_obs]
  obs_times_clean <- unix_times[valid_obs]

  # Auxiliary location terms (priors, SST, bathymetry, masks) -> additive
  # k_steps x n log-likelihood, flattened row-major (k*n + i) to match Rust.
  if (length(terms) > 0) {
    aux_mat <- build_aux_matrix(terms, lon_vec, lat_vec, knot_times)
    aux_flat <- as.numeric(t(aux_mat))
  } else {
    aux_flat <- numeric(0)
  }

  cat("Running custom Grid HMM solver in Rust...\n")
  fit <- run_grid_hmm(
    lon = as.numeric(lon_vec),
    lat = as.numeric(lat_vec),
    knot_times = as.numeric(knot_times),
    obs_times = as.numeric(obs_times_clean),
    obs_light = as.numeric(obs_light_clean),
    fixed_idx = as.integer(which(fixed_vec) - 1),
    fixed_lon = as.numeric(x0[fixed_vec, 1]),
    fixed_lat = as.numeric(x0[fixed_vec, 2]),
    diffusion = as.numeric(diffusion),
    trans_prob = as.numeric(trans_prob),
    calibration = as.numeric(calibration),
    likelihood_params = as.numeric(likelihood_params),
    shade_ratio = as.numeric(shade_ratio),
    overcount = as.numeric(overcount),
    aux_logl = aux_flat
  )
  
  # Return combined object. `log_z` is the grid HMM's log marginal likelihood
  # (model evidence) from the forward pass; use it for model comparison, e.g. a
  # hemisphere Bayes factor by differencing log_z across opposite priors.
  res <- list(
    fit = fit,
    grid = grid,
    obs_light = light,
    obs_times = date_time,
    calibration = calibration,
    likelihood_params = likelihood_params,
    log_z = fit$log_z,
    clock = clock,                # NULL unless calibrate = TRUE; the fitted clock model
    cell_lon = lon_vec,           # candidate cell coordinates (after any NA-mask filter)
    cell_lat = lat_vec            # columns of fit$posterior correspond to these cells
  )
  class(res) <- "TwilightFreeGrid"
  return(res)
}

#' Per-knot posterior over grid cells from a TwilightFreeGrid fit
#'
#' Returns the exact marginal posterior of location at each knot (from the grid
#' HMM forward-backward pass), reshaped to a `K` by `n` matrix aligned with the
#' candidate cell coordinates. Use for honest uncertainty (credible regions,
#' coverage) and for simulation-based calibration.
#'
#' @param fit A `TwilightFreeGrid` object.
#' @return A list with `lon`, `lat` (length-`n` candidate cell coordinates) and
#'   `P` (a `K` by `n` matrix; row `k` is the posterior over cells at knot `k`,
#'   summing to 1 where the knot is identified).
#' @export
grid_posterior <- function(fit) {
  stopifnot(inherits(fit, "TwilightFreeGrid"))
  n <- length(fit$cell_lon)
  K <- length(fit$fit$time)
  P <- matrix(fit$fit$posterior, nrow = K, ncol = n, byrow = TRUE)
  list(lon = fit$cell_lon, lat = fit$cell_lat, P = P)
}

#' @method print TwilightFreeGrid
#' @export
print.TwilightFreeGrid <- function(x, ...) {
  n_knots <- length(x$fit$time)
  start_t <- as.POSIXct(min(x$fit$time), origin = "1970-01-01", tz = "UTC")
  end_t   <- as.POSIXct(max(x$fit$time), origin = "1970-01-01", tz = "UTC")

  cat("\nTwilightFree Grid HMM Track\n")
  cat("=============================================\n")
  cat(sprintf("Light Observations: %d\n", length(x$obs_times)))
  cat(sprintf("Track Knots:        %d\n", n_knots))
  cat(sprintf("Duration:           %.2f days\n", diff(range(x$fit$time)) / 86400))
  cat(sprintf("Start:              %s\n", format(start_t)))
  cat(sprintf("End:                %s\n", format(end_t)))
  cat("---------------------------------------------\n")
  cat(sprintf("Mean Lat:           %.2f\n", mean(x$fit$lat)))
  cat(sprintf("Mean Lon:           %.2f\n", mean(x$fit$lon)))
  if (!is.null(x$log_z)) cat(sprintf("Log evidence (logZ): %.1f\n", x$log_z))
  if (!is.null(x$clock)) {
    cat("---------------------------------------------\n")
    for (line in format_clock(x$clock)) cat(line, "\n", sep = "")
  }
  cat("=============================================\n")
  invisible(x)
}
