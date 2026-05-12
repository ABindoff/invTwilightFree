#' Reconstruct animal tracks using a Sequential Monte Carlo (Particle Filter)
#' 
#' This reimagined approach uses all light data and a spike-and-slab likelihood 
#' to jointly estimate locations and handle non-solar light noise.
#' 
#' @param date_time A vector of POSIXct dates
#' @param light A vector of light observations
#' @param start_time Optional start date-time (POSIXct). Data before this time will be ignored.
#' @param end_time Optional end date-time (POSIXct). Data after this time will be ignored.
#' @param n_particles Number of particles (default 1000)
#' @param start_lat Initial latitude (Decimal Degrees)
#' @param start_lon Initial longitude (Decimal Degrees)
#' @param end_lat Final latitude (optional, defaults to NA)
#' @param end_lon Final longitude (optional, defaults to NA)
#' @param method Track reconstruction method: "guided" (fast Brownian Bridge), "ffbs" (rigorous Forward-Filtering Backward-Smoothing), or "forward" (standard SMC).
#' @param step_hours Time step for the particle filter knots in hours (default 12.0). Movement is proposed every `step_hours`, while likelihood is evaluated continuously.
#' @param diffusion Movement diffusion in km/sqrt(day). Typical values range from 10 to 100.
#' @param calibration Calibration parameters c(intercept, slope) mapping zenith to light. If NULL, auto-calibrates from data.
#' @param likelihood_params Likelihood parameters c(lambda, max_light, prob_slab). If NULL, auto-calibrates from data.
#' @export
#' @return A `TwilightFreeTrack` object
TwilightFreeSMC <- function(date_time, light, 
                           start_time = NULL, end_time = NULL,
                           n_particles = 1000, 
                           start_lat, start_lon, 
                           end_lat = NA_real_, end_lon = NA_real_,
                           method = c("guided", "ffbs", "forward"),
                           step_hours = 12.0,
                           diffusion = 50, 
                           calibration = NULL,
                           likelihood_params = NULL) {
  
  if(!inherits(date_time, "POSIXct")) {
    stop("date_time must be POSIXct")
  }
  
  # Filter by time if provided
  if (!is.null(start_time)) {
    keep <- date_time >= start_time
    date_time <- date_time[keep]
    light <- light[keep]
  }
  if (!is.null(end_time)) {
    keep <- date_time <= end_time
    date_time <- date_time[keep]
    light <- light[keep]
  }

  if (length(date_time) == 0) stop("No data remains after time filtering!")
  
  method <- match.arg(method)
  
  # Ensure sorted
  ord <- order(date_time)
  date_time <- date_time[ord]
  light <- light[ord]
  
  # Auto-Calibration
  if (is.null(calibration) || is.null(likelihood_params)) {
    min_l <- quantile(light, 0.05, na.rm = TRUE)
    max_l <- quantile(light, 0.95, na.rm = TRUE)
    
    # Shift light so baseline is 0
    light_shifted <- pmax(0, light - min_l)
    max_shifted <- max_l - min_l
    
    if (is.null(calibration)) {
      # Data-driven calibration using first 3 days at start location
      cal_idx <- date_time < (date_time[1] + 3*24*3600)
      cal_times <- date_time[cal_idx]
      cal_light <- light_shifted[cal_idx]
      
      # Calculate true solar zeniths
      d <- (as.numeric(cal_times) - 946728000) / 86400
      g <- (357.529 + 0.98560028 * d) * pi / 180
      q <- 280.459 + 0.98564736 * d
      l_sun <- (q + 1.915 * sin(g) + 0.020 * sin(2 * g)) * pi / 180
      e <- (23.439 - 0.00000036 * d) * pi / 180
      delta <- asin(sin(e) * sin(l_sun))
      ra <- atan2(cos(e) * sin(l_sun), cos(l_sun))
      gmst <- (18.697374558 + 24.06570982441908 * d) * 15 * pi / 180
      lmst <- gmst + start_lon * pi / 180
      h <- lmst - ra
      phi <- start_lat * pi / 180
      cos_z <- sin(phi) * sin(delta) + cos(phi) * cos(delta) * cos(h)
      cal_zenith <- acos(pmax(-1, pmin(1, cos_z))) * 180 / pi
      
      # Fit linear model for transitions (zenith between 85 and 96, light > 0 and < max)
      trans_idx <- cal_light > 0 & cal_light < max_shifted * 0.95 & cal_zenith > 85 & cal_zenith < 100
      if (sum(trans_idx) > 10) {
        fit_cal <- lm(cal_light[trans_idx] ~ cal_zenith[trans_idx])
        intercept_est <- coef(fit_cal)[1]
        slope_est <- -coef(fit_cal)[2] # slope is negative in the model, but we want positive slope for calibration
        if (is.na(slope_est) || slope_est <= 0) slope_est <- max_shifted / (96 - 85)
      } else {
        slope_est <- max_shifted / (96 - 85)
        intercept_est <- slope_est * 96
      }
      
      calibration <- c(intercept_est, slope_est)
    }
    
    if (is.null(likelihood_params)) {
      # lambda controls the exponential decay for shading. 
      # We use the max light as the scale factor.
      lambda_est <- 1.0 / (max_shifted * 0.5)
      # We also increase prob_slab slightly to 0.1 to account for genuine outliers.
      likelihood_params <- c(lambda_est, max_shifted, 0.10)
    }
    
    # We must use the shifted light for the particle filter
    process_light <- light_shifted
  } else {
    process_light <- light
  }
  
  unix_times <- as.numeric(date_time)
  
  res <- run_particle_filter(
    unix_times = unix_times,
    obs_light = process_light,
    n_particles = as.integer(n_particles),
    start_lat = as.numeric(start_lat),
    start_lon = as.numeric(start_lon),
    end_lat = as.numeric(end_lat),
    end_lon = as.numeric(end_lon),
    method = as.character(method),
    step_hours = as.numeric(step_hours),
    diffusion = as.numeric(diffusion),
    calibration = as.numeric(calibration),
    likelihood_params = as.numeric(likelihood_params)
  )
  
  res$obs_light <- light
  class(res) <- "TwilightFreeTrack"
  return(res)
}

#' @method print TwilightFreeTrack
#' @export
print.TwilightFreeTrack <- function(x, ...) {
  n_knots <- length(x$knot_times)
  n_obs <- length(x$obs_times)
  start_t <- as.POSIXct(min(x$knot_times), origin = "1970-01-01", tz = "UTC")
  end_t <- as.POSIXct(max(x$knot_times), origin = "1970-01-01", tz = "UTC")
  
  cat("\nTwilightFree Continuous Particle Filter Track\n")
  cat("=============================================\n")
  cat(sprintf("Light Observations: %d\n", n_obs))
  cat(sprintf("Track Knots:        %d\n", n_knots))
  cat(sprintf("Duration:           %.2f days\n", diff(range(x$knot_times))/86400))
  cat(sprintf("Start:              %s\n", format(start_t)))
  cat(sprintf("End:                %s\n", format(end_t)))
  cat("---------------------------------------------\n")
  cat(sprintf("Mean Lat:           %.2f (Avg SD: %.2f)\n", mean(x$lat), mean(x$lat_sd)))
  cat(sprintf("Mean Lon:           %.2f (Avg SD: %.2f)\n", mean(x$lon), mean(x$lon_sd)))
  cat(sprintf("Avg Anomaly Prob:   %.1f%%\n", mean(x$prob_false) * 100))
  cat("=============================================\n")
  invisible(x)
}

#' @import ggplot2
#' @importFrom rnaturalearth ne_countries
#' @importFrom sf st_as_sf st_coordinates
#' @method plot TwilightFreeTrack
#' @export
plot.TwilightFreeTrack <- function(x, type = c("track", "uncertainty", "diagnostics", "image"), ...) {
  type <- match.arg(type)
  d_time <- as.POSIXct(x$knot_times, origin = "1970-01-01", tz = "UTC")
  obs_time <- as.POSIXct(x$obs_times, origin = "1970-01-01", tz = "UTC")
  
  if (type == "track") {
    # Prepare data
    df <- data.frame(lon = x$lon, lat = x$lat)
    
    # Get world map
    world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
    
    # Bounding box with padding
    lon_range <- range(df$lon, na.rm = TRUE)
    lat_range <- range(df$lat, na.rm = TRUE)
    lon_pad <- max(5, diff(lon_range) * 0.2)
    lat_pad <- max(5, diff(lat_range) * 0.2)
    
    p <- ggplot() +
      geom_sf(data = world, fill = "#e0e0e0", color = "#f8f9fa", size = 0.2) +
      geom_path(data = df, aes(x = lon, y = lat), color = "#2c3e50", size = 1.2, alpha = 0.8) +
      # Points with alpha to show density
      geom_point(data = df, aes(x = lon, y = lat), color = "#34495e", size = 1.5, alpha = 0.5) +
      # Start and end markers (sized so overlap is visible as a bullseye)
      geom_point(data = df[1, , drop=FALSE], aes(x = lon, y = lat), 
                 shape = 21, fill = "#2ecc71", color = "#c0392b", size = 5, stroke = 1.5) +
      geom_point(data = df[nrow(df), , drop=FALSE], aes(x = lon, y = lat), 
                 shape = 21, fill = "#e74c3c", color = "#c0392b", size = 2.5, stroke = 1.5) +
      coord_sf(xlim = c(lon_range[1] - lon_pad, lon_range[2] + lon_pad),
               ylim = c(lat_range[1] - lat_pad, lat_range[2] + lat_pad), 
               expand = FALSE) +
      theme_minimal(base_family = "sans") +
      theme(
        panel.background = element_rect(fill = "#ebf5fb", color = NA),
        panel.grid.major = element_line(color = "white", size = 0.5),
        plot.title = element_text(face = "bold", size = 14, color = "#2c3e50"),
        axis.title = element_text(face = "bold", color = "#34495e")
      ) +
      labs(title = "Reconstructed Track (Continuous SMC)", x = "Longitude", y = "Latitude")
    
    return(p)
    
  } else if (type == "uncertainty") {
    df <- data.frame(time = d_time, lat = x$lat, lon = x$lon, lat_sd = x$lat_sd, lon_sd = x$lon_sd)
    
    p1 <- ggplot(df, aes(x = time)) +
      geom_ribbon(aes(ymin = lat - 2*lat_sd, ymax = lat + 2*lat_sd), fill = "#3498db", alpha = 0.3) +
      geom_line(aes(y = lat), color = "#2980b9", size = 1) +
      theme_minimal() + labs(y = "Latitude", x = "") +
      ggtitle("Posterior Latitude & Longitude Uncertainty")
      
    p2 <- ggplot(df, aes(x = time)) +
      geom_ribbon(aes(ymin = lon - 2*lon_sd, ymax = lon + 2*lon_sd), fill = "#e67e22", alpha = 0.3) +
      geom_line(aes(y = lon), color = "#d35400", size = 1) +
      theme_minimal() + labs(y = "Longitude", x = "Time")
      
    # If patchwork or gridExtra is needed, but we can just use cowplot or print in viewport.
    # Using simple layout via basic graphics is impossible if returning ggplot objects.
    # We will return a list of plots and print them.
    # To keep dependencies low, we just print them on a grid.
    gridExtra::grid.arrange(p1, p2, ncol = 1)
    invisible(list(p1 = p1, p2 = p2))
    
  } else if (type == "diagnostics") {
    df_obs <- data.frame(time = obs_time, zenith = x$obs_zenith, light = x$obs_light, prob = x$prob_false)
    
    # Calculate appropriate bar width based on data sampling interval
    obs_diff <- as.numeric(difftime(df_obs$time[2], df_obs$time[1], units="secs"))
    bar_width <- if(!is.na(obs_diff) && obs_diff > 0) obs_diff else 600
    
    p1 <- ggplot(df_obs, aes(x = time, y = light)) +
      geom_col(aes(fill = prob > 0.5), width = bar_width) +
      scale_fill_manual(values = c("FALSE" = "#bdc3c7", "TRUE" = "#c0392b")) +
      theme_minimal() + labs(x = "", y = "Observed Light") +
      theme(legend.position = "none") +
      ggtitle("Observed Light (Anomalies Highlighted in Red)")
      
    p2 <- ggplot(df_obs, aes(x = time, y = zenith)) +
      geom_line(color = "#2c3e50", size = 0.5) +
      theme_minimal() + labs(x = "", y = "Est. Zenith") +
      ggtitle("Estimated Solar Zenith Angle")
      
    p3 <- ggplot(df_obs, aes(x = time, y = prob)) +
      geom_col(fill = "#e67e22", width = bar_width) +
      geom_hline(yintercept = 0.5, linetype = "dashed", color = "#7f8c8d") +
      theme_minimal() + labs(x = "Time", y = "P(Anomaly)") +
      coord_cartesian(ylim = c(0, 1))
      
    gridExtra::grid.arrange(p1, p2, p3, ncol = 1)
    invisible(list(p1 = p1, p2 = p2, p3 = p3))
  } else if (type == "image") {
    # LightImage style plot
    df <- data.frame(
      time = obs_time,
      light = x$obs_light,
      date = as.Date(obs_time),
      hour = as.numeric(format(obs_time, "%H")) + as.numeric(format(obs_time, "%M"))/60
    )
    
    p <- ggplot(df, aes(x = date, y = hour, fill = light)) +
      geom_tile() +
      scale_fill_viridis_c(option = "magma", name = "Light") +
      scale_y_continuous(breaks = seq(0, 24, by = 4), limits = c(0, 24)) +
      theme_minimal() +
      labs(title = "Light Intensity Image", x = "Date", y = "Hour of Day (UTC)") +
      theme(panel.grid = element_blank())
    
    return(p)
  }
}
