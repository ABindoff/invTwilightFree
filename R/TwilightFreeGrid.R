#' Reconstruct animal tracks using SGAT::essie and continuous likelihood
#' 
#' @param date_time A vector of POSIXct dates
#' @param light A vector of light observations
#' @param grid A RasterLayer or SpatialPixels object defining the grid (from TwilightFree::makeGrid)
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
#' @export
TwilightFreeGrid <- function(date_time, light, grid, 
                             start_lat = NA_real_, start_lon = NA_real_,
                             end_lat = NA_real_, end_lon = NA_real_,
                             fixed = NULL,
                             step_hours = 12.0, 
                             diffusion = 50, 
                             trans_prob = NULL,
                             calibration = NULL, 
                             likelihood_params = NULL) {
  
  if(!inherits(date_time, "POSIXct")) stop("date_time must be POSIXct")
  
  # Ensure sorted
  ord <- order(date_time)
  date_time <- date_time[ord]
  light <- light[ord]
  
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
  
  # Call Rust grid HMM solver
  lon_vec <- raster::coordinates(grid)[, 1]
  lat_vec <- raster::coordinates(grid)[, 2]
  
  # Ensure valid data
  valid_obs <- !is.na(process_light)
  obs_light_clean <- process_light[valid_obs]
  obs_times_clean <- unix_times[valid_obs]
  
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
    likelihood_params = as.numeric(likelihood_params)
  )
  
  # Return combined object
  res <- list(
    fit = fit,
    grid = grid,
    obs_light = light,
    obs_times = date_time,
    calibration = calibration,
    likelihood_params = likelihood_params
  )
  class(res) <- "TwilightFreeGrid"
  return(res)
}
