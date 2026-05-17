
devtools::load_all(".")
library(ggplot2)
library(patchwork)
library(raster)
library(dplyr)
library(tidyr)

# Helper for Distance
get_dist <- function(lon1, lat1, lon2, lat2) {
  rad <- pi/180
  a1 <- lat1 * rad; a2 <- lat2 * rad
  b1 <- lon1 * rad; b2 <- lon2 * rad
  dlon <- b2 - b1
  acos(pmin(1, sin(a1)*sin(a2) + cos(a1)*cos(a2)*cos(dlon))) * 6371.0
}

# Load the Scenarios and the Ground Truth
track_data <- readRDS("scratch/simulated_seal_light_scenarios.rds")

get_paths <- function(scenario = "cloudy") {
  cat("  -> Loading scenario:", scenario, "\n")
  prefix <- paste0("_seal_", scenario)
  light_col <- paste0("light_", scenario)
  if (scenario == "cloudy") light_col <- "light_ideal"
  
  # Generate a spatial mask
  cat("    - Generating spatial mask...\n")
  grid <- TwilightFree::makeGrid(lon=c(130, 185), lat=c(-85, -40), cell.size=2, mask="sea")
  
  # Fit Original TwilightFree
  cat("    - Fitting Original TwilightFree...\n")
  df_tf <- data.frame(Date = track_data$time, Light = track_data[[light_col]])
  # Ensure POSIXct
  df_tf$Date <- as.POSIXct(df_tf$Date, tz="UTC")
  
  df_orig_tf <- NULL
  tryCatch({
    fit_tf_orig <- TwilightFree::TwilightFree(df_tf, deployed.at=c(track_data$true_lon[1], track_data$true_lat[1]), retrieved.at=c(track_data$true_lon[nrow(track_data)], track_data$true_lat[nrow(track_data)]))
    ess <- SGAT::essie(fit_tf_orig, grid, epsilon1=0.01)
    trip_mat <- TwilightFree::trip(ess)
    df_orig_tf <- data.frame(
      Time = as.POSIXct(trip_mat[, "Date"], origin="1970-01-01", tz="UTC"), 
      Lat = trip_mat[, "Lat"], 
      Lon = trip_mat[, "Lon"], 
      Method = "TwilightFree (Orig)"
    )
  }, error = function(e) {
    cat("    - Original TwilightFree failed:", e$message, "\n")
  })
  
  # Fit invTF (Guided)
  cat("    - Fitting invTF (Guided)...\n")
  fit_inv_guided <- TwilightFreeSMC(
    date_time = track_data$time, light = track_data[[light_col]],
    start_lat = track_data$true_lat[1], start_lon = track_data$true_lon[1],
    end_lat = track_data$true_lat[nrow(track_data)], end_lon = track_data$true_lon[nrow(track_data)],
    method = "guided", n_particles = 1000,
    calibration = true_cal, likelihood_params = true_lik,
    diffusion = 200, step_hours = 2.0, spatial_mask = grid
  )
  
  # Fit invTF (FFBS)
  cat("    - Fitting invTF (FFBS)...\n")
  fit_inv_ffbs <- TwilightFreeSMC(
    date_time = track_data$time, light = track_data[[light_col]],
    start_lat = track_data$true_lat[1], start_lon = track_data$true_lon[1],
    end_lat = track_data$true_lat[nrow(track_data)], end_lon = track_data$true_lon[nrow(track_data)],
    method = "ffbs", n_particles = 1000,
    calibration = true_cal, likelihood_params = true_lik,
    diffusion = 200, step_hours = 2.0, spatial_mask = grid
  )
  
  # Fit invTF (Grid)
  cat("    - Fitting invTF (Grid)...\n")
  df_inv_grid <- NULL
  tryCatch({
    res_grid <- TwilightFreeGrid(
      date_time = as.POSIXct(track_data$time, tz="UTC"), 
      light = track_data[[light_col]], 
      grid = grid,
      step_hours = 4.0,
      diffusion = 100, 
      calibration = true_cal, 
      likelihood_params = true_lik
    )
    trip_mat_grid <- TwilightFree::trip(res_grid$fit)
    df_inv_grid <- data.frame(
      Time = as.POSIXct(trip_mat_grid[, "Date"], origin="1970-01-01", tz="UTC"), 
      Lat = trip_mat_grid[, "Lat"], 
      Lon = trip_mat_grid[, "Lon"], 
      Method = "invTF (Grid)"
    )
  }, error = function(e) {
    cat("    - invTF (Grid) failed:", e$message, "\n")
  })
  
  # Load Caches
  cat("    - Loading FLightR cache...\n")
  f_flightr <- paste0("vignettes/cache_flightr", prefix, ".rds")
  cat("    - Loading SGAT cache (this may take a while)...\n")
  f_sgat <- paste0("vignettes/cache_sgat", prefix, ".rds")
  
  if (!file.exists(f_flightr)) {
    cat("Missing cache for", scenario, "\n")
    return(NULL)
  }
  
  c_flightr <- readRDS(f_flightr)
  obj_flightr <- if ("fit" %in% names(c_flightr$fit)) c_flightr$fit$fit else c_flightr$fit
  
  df_flightr <- data.frame(Time = as.POSIXct(obj_flightr$Results$Quantiles$time, tz="UTC"), Lat = obj_flightr$Results$Quantiles$Medianlat, Lon = obj_flightr$Results$Quantiles$Medianlon, Method = "FLightR")
  
  df_inv_guided <- data.frame(Time = as.POSIXct(fit_inv_guided$knot_times, origin="1970-01-01", tz="UTC"), Lat = fit_inv_guided$lat, Lon = fit_inv_guided$lon)
  df_inv_ffbs <- data.frame(Time = as.POSIXct(fit_inv_ffbs$knot_times, origin="1970-01-01", tz="UTC"), Lat = fit_inv_ffbs$lat, Lon = fit_inv_ffbs$lon)
  
  # Fair evaluation: invTF produces arbitrary knots (e.g. noon/midnight). 
  # Twilight geolocation only has information at twilights. 
  # We interpolate invTF to the twilight times to evaluate it fairly.
  twl_times <- as.numeric(df_flightr$Time)
  
  df_inv_guided_eval <- data.frame(
    Time = df_flightr$Time,
    Lat = approx(as.numeric(df_inv_guided$Time), df_inv_guided$Lat, xout = twl_times, rule=2)$y,
    Lon = approx(as.numeric(df_inv_guided$Time), df_inv_guided$Lon, xout = twl_times, rule=2)$y,
    Method = "invTF (Guided)"
  )
  
  df_inv_ffbs_eval <- data.frame(
    Time = df_flightr$Time,
    Lat = approx(as.numeric(df_inv_ffbs$Time), df_inv_ffbs$Lat, xout = twl_times, rule=2)$y,
    Lon = approx(as.numeric(df_inv_ffbs$Time), df_inv_ffbs$Lon, xout = twl_times, rule=2)$y,
    Method = "invTF (FFBS)"
  )
  
  df_inv_grid_eval <- NULL
  if (!is.null(df_inv_grid)) {
    df_inv_grid_eval <- data.frame(
      Time = df_flightr$Time,
      Lat = approx(as.numeric(df_inv_grid$Time), df_inv_grid$Lat, xout = twl_times, rule=2)$y,
      Lon = approx(as.numeric(df_inv_grid$Time), df_inv_grid$Lon, xout = twl_times, rule=2)$y,
      Method = "invTF (Grid)"
    )
  }
  
  all_dfs <- list(df_inv_guided_eval, df_inv_ffbs_eval, df_flightr)
  if (!is.null(df_inv_grid_eval)) {
    all_dfs <- append(all_dfs, list(df_inv_grid_eval))
  }
  if (!is.null(df_orig_tf)) {
    all_dfs <- append(all_dfs, list(df_orig_tf))
  }
  
  if (file.exists(f_sgat)) {
    c_sgat <- readRDS(f_sgat)
    obj_sgat <- if ("fit" %in% names(c_sgat$fit)) c_sgat$fit$fit else c_sgat$fit
    sgat_time <- if(!is.null(obj_sgat$model$z)) obj_sgat$model$z else obj_sgat$time
    sgat_lat <- if(!is.null(obj_sgat$x)) apply(obj_sgat$x[[1]][,2,], 1, mean) else obj_sgat$lat
    sgat_lon <- if(!is.null(obj_sgat$x)) apply(obj_sgat$x[[1]][,1,], 1, mean) else obj_sgat$lon
    df_sgat <- data.frame(Time = as.POSIXct(sgat_time, tz="UTC"), Lat = sgat_lat, Lon = sgat_lon, Method = "SGAT")
    all_dfs <- append(all_dfs, list(df_sgat))
  }
  
  # Combine and calculate errors
  all_df <- bind_rows(all_dfs)
  
  all_df <- all_df %>%
    mutate(
      true_lat = approx(as.numeric(track_data$time), track_data$true_lat, xout = as.numeric(Time))$y,
      true_lon = approx(as.numeric(track_data$time), track_data$true_lon, xout = as.numeric(Time))$y,
      lat_err = Lat - true_lat,
      lon_err = Lon - true_lon,
      dist_err = get_dist(Lon, Lat, true_lon, true_lat)
    ) %>%
    mutate(
      lon_err = ifelse(lon_err > 180, lon_err - 360, ifelse(lon_err < -180, lon_err + 360, lon_err))
    )
  
  return(list(df = all_df, scenario = scenario))
}

plot_errors <- function(res) {
  df <- res$df
  scenario <- res$scenario
  
  p1 <- ggplot(df, aes(x = Time, y = lat_err, color = Method)) +
    geom_line(alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_minimal() +
    labs(title = paste("Latitude Error -", scenario), y = "Error (deg)")
    
  p2 <- ggplot(df, aes(x = Time, y = lon_err, color = Method)) +
    geom_line(alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_minimal() +
    labs(title = paste("Longitude Error -", scenario), y = "Error (deg)")
    
  p3 <- ggplot(df, aes(x = Time, y = dist_err, color = Method)) +
    geom_line(alpha = 0.7) +
    theme_minimal() +
    labs(title = paste("Distance Error -", scenario), y = "Error (km)")
    
  return(p1 / p2 / p3 + plot_layout(guides = "collect"))
}

# Analyze all scenarios
cat("Analyzing scenarios...\n")

# Provide true calibration from generate_seal_light.R
# light = 558.5 - 5.818 * z
true_cal <- c(558.5, 5.818)
# likelihood params: lambda (shading), max_light, prob_slab
true_lik <- c(1.0 / (64 * 0.5), 64, 0.1)

res_cloudy <- get_paths("cloudy")
if(!is.null(res_cloudy)) {
  cat("  -> Saving Cloudy plot...\n")
  ggsave("scratch/diagnostic_cloudy.png", plot_errors(res_cloudy), width = 12, height = 12)
}

res_shaded <- get_paths("shaded")
# res_shaded <- NULL
if(!is.null(res_shaded)) {
  cat("  -> Saving Shaded plot...\n")
  ggsave("scratch/diagnostic_shaded.png", plot_errors(res_shaded), width = 12, height = 12)
}

res_alan <- get_paths("alan")
# res_alan <- NULL
if(!is.null(res_alan)) {
  cat("  -> Saving ALAN plot...\n")
  ggsave("scratch/diagnostic_alan.png", plot_errors(res_alan), width = 12, height = 12)
}

# Print summary RMSE
cat("\nSummary RMSE (km):\n")
if(!is.null(res_cloudy)) {
  cat("Cloudy:\n")
  res_cloudy$df %>% group_by(Method) %>% summarize(
    RMSE_Dist = sqrt(mean(dist_err^2, na.rm=TRUE)), 
    RMSE_Lat_deg = sqrt(mean(lat_err^2, na.rm=TRUE)),
    RMSE_Lon_deg = sqrt(mean(lon_err^2, na.rm=TRUE)),
    MedianLonErr = median(lon_err, na.rm=TRUE)
  ) %>% print()
}
if(!is.null(res_shaded)) {
  cat("Shaded:\n")
  res_shaded$df %>% group_by(Method) %>% summarize(
    RMSE_Dist = sqrt(mean(dist_err^2, na.rm=TRUE)), 
    RMSE_Lat_deg = sqrt(mean(lat_err^2, na.rm=TRUE)),
    RMSE_Lon_deg = sqrt(mean(lon_err^2, na.rm=TRUE)),
    MedianLonErr = median(lon_err, na.rm=TRUE)
  ) %>% print()
}
if(!is.null(res_alan)) {
  cat("ALAN:\n")
  res_alan$df %>% group_by(Method) %>% summarize(
    RMSE_Dist = sqrt(mean(dist_err^2, na.rm=TRUE)), 
    RMSE_Lat_deg = sqrt(mean(lat_err^2, na.rm=TRUE)),
    RMSE_Lon_deg = sqrt(mean(lon_err^2, na.rm=TRUE)),
    MedianLonErr = median(lon_err, na.rm=TRUE)
  ) %>% print()
}

