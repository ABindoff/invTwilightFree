## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 12,
  fig.height = 12,
  out.width = "100%"
)
# Ensure the local version of the package is loaded

library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(geosphere)

# Helper for distance calculation
get_dist <- function(lon1, lat1, lon2, lat2) {
  distHaversine(cbind(lon1, lat1), cbind(lon2, lat2)) / 1000
}

# Helper to unwrap caches
find_fit <- function(obj, target) {
  curr <- obj
  for(i in 1:10) {
    if (target %in% names(curr)) return(curr)
    if ("fit" %in% names(curr)) curr <- curr$fit else break
  }
  return(NULL)
}

# Core benchmarking function
get_bench_data <- function(scenario = "ideal", track_data) {
  # Prefix for caches
  f_prefix <- if(scenario == "ideal") "_ideal" else paste0("_seal_", scenario)
  light_col <- if(scenario == "ideal") "light_ideal" else if(scenario == "cloudy") "light_ideal" else paste0("light_", scenario)
  
  # Fit invTF (Guided)
  t_inv_guided <- system.time({
    fit_inv_guided <- TwilightFreeSMC(
      date_time = track_data$time, light = track_data[[light_col]],
      start_lat = track_data$true_lat[1], start_lon = track_data$true_lon[1],
      end_lat = track_data$true_lat[nrow(track_data)], end_lon = track_data$true_lon[nrow(track_data)],
      method = "guided", n_particles = 5000,
      diffusion = 60,
      calibration = c(311.57, 3.19),
      likelihood_params = c(1.0, 64, 0.05)
    )
  })
  
  # Fit invTF (Grid)
  t_inv_grid <- system.time({
    grid <- TwilightFree::makeGrid(lon=c(130, 185), lat=c(-85, -40), cell.size=1.0, mask="sea")
    fit_inv_grid <- TwilightFreeGrid(
      date_time = as.POSIXct(track_data$time, tz="UTC"), 
      light = track_data[[light_col]], 
      grid = grid,
      start_lat = track_data$true_lat[1], start_lon = track_data$true_lon[1],
      end_lat = track_data$true_lat[nrow(track_data)], end_lon = track_data$true_lon[nrow(track_data)],
      step_hours = 4.0,
      diffusion = 60
    )
  })
  
  # Load Caches
  f_flightr <- paste0("cache_flightr", f_prefix, ".rds")
  f_sgat <- paste0("cache_sgat", f_prefix, ".rds")
  
  # Try local and vignettes/
  found_f <- FALSE
  found_s <- FALSE
  if (file.exists(f_flightr)) found_f <- TRUE
  else if (file.exists(paste0("vignettes/", f_flightr))) { f_flightr <- paste0("vignettes/", f_flightr); found_f <- TRUE }
  
  if (file.exists(f_sgat)) found_s <- TRUE
  else if (file.exists(paste0("vignettes/", f_sgat))) { f_sgat <- paste0("vignettes/", f_sgat); found_s <- TRUE }
  
  if (!found_f || !found_s) {
    warning(paste("Caches missing for", scenario))
    return(NULL)
  }
  
  c_flightr <- readRDS(f_flightr)
  c_sgat <- readRDS(f_sgat)
  
  obj_flightr <- find_fit(c_flightr, "Results")
  obj_sgat <- find_fit(c_sgat, "model")
  
  # Extract Paths
  df_inv_guided <- data.frame(
    Time = as.POSIXct(fit_inv_guided$knot_times, origin="1970-01-01", tz="UTC"), 
    Lat = fit_inv_guided$lat, 
    Lon = fit_inv_guided$lon, 
    Method = "invTF (Guided)"
  )
  
  trip_grid <- TwilightFree::trip(fit_inv_grid$fit)
  df_inv_grid <- data.frame(Time = as.POSIXct(trip_grid[, 1], origin="1970-01-01", tz="UTC"), Method = "invTF (Grid)")
  
  # Smart column assignment for Grid results
  c2 <- as.numeric(trip_grid[, 2])
  c3 <- as.numeric(trip_grid[, 3])
  if (all(abs(c2) <= 90, na.rm=TRUE) && !all(abs(c3) <= 90, na.rm=TRUE)) {
    df_inv_grid$Lat <- c2
    df_inv_grid$Lon <- c3
  } else if (all(abs(c3) <= 90, na.rm=TRUE) && !all(abs(c2) <= 90, na.rm=TRUE)) {
    df_inv_grid$Lat <- c3
    df_inv_grid$Lon <- c2
  } else {
    # If both or neither are in [-90, 90], assume standard order but check names
    if ("Lat" %in% colnames(trip_grid)) {
       df_inv_grid$Lat <- as.numeric(trip_grid[, "Lat"])
       df_inv_grid$Lon <- as.numeric(trip_grid[, "Lon"])
    } else {
       df_inv_grid$Lat <- c2
       df_inv_grid$Lon <- c3
    }
  }
  
  df_flightr <- data.frame(Time = obj_flightr$Results$Quantiles$time, Lat = obj_flightr$Results$Quantiles$Medianlat, Lon = obj_flightr$Results$Quantiles$Medianlon, Method = "FLightR")
  
  sgat_time_val <- if (!is.null(obj_sgat$model$time)) obj_sgat$model$time else obj_sgat$model$z
  df_sgat <- data.frame(Time = sgat_time_val, Lat = apply(obj_sgat$x[[1]][,2,], 1, mean), Lon = apply(obj_sgat$x[[1]][,1,], 1, mean), Method = "SGAT")
  
  # Metrics function
  calc_metrics <- function(df, name) {
    if (is.null(df) || nrow(df) == 0) return(NULL)
    t_lat <- approx(as.numeric(track_data$time), track_data$true_lat, xout = as.numeric(df$Time))$y
    t_lon <- approx(as.numeric(track_data$time), track_data$true_lon, xout = as.numeric(df$Time))$y
    
    # Ensure latitudes are within physical bounds for distance calculation
    safe_lat <- pmax(-90, pmin(90, df$Lat))
    errs <- get_dist(df$Lon, safe_lat, t_lon, t_lat)
    
    lat_err <- abs(df$Lat - t_lat)
    lon_err <- abs(df$Lon - t_lon)
    
    data.frame(
      Method = name,
      RMSE = sqrt(mean(errs^2, na.rm=TRUE)),
      Mean_Lat_Err = mean(lat_err, na.rm=TRUE),
      Mean_Lon_Err = mean(lon_err, na.rm=TRUE),
      Time_s = NA
    )
  }
  
  m_guided <- calc_metrics(df_inv_guided, "invTF (Guided)")
  m_guided$Time_s <- t_inv_guided["elapsed"]
  
  m_grid <- calc_metrics(df_inv_grid, "invTF (Grid)")
  m_grid$Time_s <- t_inv_grid["elapsed"]
  
  m_flightr <- calc_metrics(df_flightr, "FLightR")
  m_flightr$Time_s <- c_flightr$time["elapsed"]
  
  m_sgat <- calc_metrics(df_sgat, "SGAT")
  m_sgat$Time_s <- c_sgat$time["elapsed"]
  
  metrics <- bind_rows(m_guided, m_grid, m_flightr, m_sgat)
  
  list(
    paths = bind_rows(df_inv_guided, df_inv_grid, df_flightr, df_sgat) %>% mutate(Scenario = scenario),
    metrics = metrics %>% mutate(Scenario = scenario)
  )
}


## ----load_data, message=FALSE, warning=FALSE----------------------------------
# Path adjustment for vignette rendering
track_path <- if(file.exists("scratch/simulated_seal_light_scenarios.rds")) "scratch/simulated_seal_light_scenarios.rds" else "../scratch/simulated_seal_light_scenarios.rds"
track_data <- readRDS(track_path)
scenarios <- c("ideal", "cloudy", "shaded", "alan")
all_results <- lapply(scenarios, function(s) get_bench_data(s, track_data))
all_paths <- bind_rows(lapply(all_results, function(x) x$paths))
all_metrics <- bind_rows(lapply(all_results, function(x) x$metrics))


## ----summary_table, echo=FALSE------------------------------------------------
summary_tab <- all_metrics %>%
  select(Scenario, Method, RMSE, Time_s) %>%
  mutate(RMSE = round(RMSE, 1), Time_s = round(Time_s, 1)) %>%
  pivot_wider(names_from = Scenario, values_from = c(RMSE, Time_s))

knitr::kable(summary_tab, caption = "Accuracy (RMSE km) and Performance (Time s)")


## ----plot_trajectories, echo=FALSE, fig.height=14-----------------------------
plot_scenario <- function(scen) {
  scen_paths <- all_paths %>% filter(Scenario == scen)
  
  ggplot() +
    geom_path(data = track_data, aes(x = true_lon, y = true_lat), color = "black", linewidth = 2, alpha = 0.15) +
    geom_path(data = scen_paths, aes(x = Lon, y = Lat, color = Method), linewidth = 0.8) +
    facet_wrap(~Method, ncol = 2) +
    theme_minimal() +
    labs(title = paste("Scenario:", toupper(scen)), x = "Longitude", y = "Latitude") +
    theme(legend.position = "none")
}

(plot_scenario("ideal") / plot_scenario("cloudy"))
(plot_scenario("shaded") / plot_scenario("alan"))

