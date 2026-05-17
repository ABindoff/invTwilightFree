## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 12,
  fig.height = 10,
  out.width = "100%"
)


## ----setup, message=FALSE, warning=FALSE--------------------------------------

library(ggplot2)
library(patchwork)
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
  prefix <- paste0("_seal_", scenario)
  light_col <- paste0("light_", scenario)
  if (scenario == "cloudy") light_col <- "light_ideal"
  
  # Fit invTF (Guided)
  t_inv_guided <- system.time({
    fit_inv_guided <- TwilightFreeSMC(
      date_time = track_data$time, light = track_data[[light_col]],
      start_lat = track_data$true_lat[1], start_lon = track_data$true_lon[1],
      end_lat = track_data$true_lat[nrow(track_data)], end_lon = track_data$true_lon[nrow(track_data)],
      method = "guided", n_particles = 1000
    )
  })
  
  # Fit invTF (Grid)
  t_inv_grid <- system.time({
    grid <- TwilightFree::makeGrid(lon=c(130, 185), lat=c(-85, -40), cell.size=2, mask="sea")
    fit_inv_grid <- TwilightFreeGrid(
      date_time = as.POSIXct(track_data$time, tz="UTC"), 
      light = track_data[[light_col]], 
      grid = grid,
      start_lat = track_data$true_lat[1], start_lon = track_data$true_lon[1],
      end_lat = track_data$true_lat[nrow(track_data)], end_lon = track_data$true_lon[nrow(track_data)],
      step_hours = 4.0
    )
  })
  
  # Load Caches
  f_flightr <- paste0("cache_flightr", prefix, ".rds")
  f_sgat <- paste0("cache_sgat", prefix, ".rds")
  f_orig <- paste0("cache_orig_tf", prefix, ".rds")
  
  if (!file.exists(f_flightr) || !file.exists(f_sgat)) return(NULL)
  
  c_flightr <- readRDS(f_flightr)
  c_sgat <- readRDS(f_sgat)
  
  # Unwrap
  obj_flightr <- if ("fit" %in% names(c_flightr$fit)) c_flightr$fit$fit else c_flightr$fit
  obj_sgat <- if ("fit" %in% names(c_sgat$fit)) c_sgat$fit$fit else c_sgat$fit
  
  df_inv_guided <- data.frame(Time = fit_inv_guided$knot_times, Lat = fit_inv_guided$lat, Lon = fit_inv_guided$lon, Method = "invTF (Guided)")
  trip_grid <- TwilightFree::trip(fit_inv_grid$fit)
  df_inv_grid <- data.frame(Time = as.POSIXct(trip_grid[, "Date"], origin="1970-01-01", tz="UTC"), Lat = trip_grid[, "Lat"], Lon = trip_grid[, "Lon"], Method = "invTF (Grid)")
  df_flightr <- data.frame(Time = obj_flightr$Results$Quantiles$time, Lat = obj_flightr$Results$Quantiles$Medianlat, Lon = obj_flightr$Results$Quantiles$Medianlon, Method = "FLightR")
  df_sgat <- data.frame(Time = obj_sgat$model$z, Lat = apply(obj_sgat$x[[1]][,2,], 1, mean), Lon = apply(obj_sgat$x[[1]][,1,], 1, mean), Method = "SGAT")
  df_orig_tf <- NULL
  if (file.exists(f_orig)) {
    obj_orig <- readRDS(f_orig)
    trip_orig <- TwilightFree::trip(obj_orig$fit)
    df_orig_tf <- data.frame(Time = as.POSIXct(trip_orig[, "Date"], origin="1970-01-01", tz="UTC"), Lat = trip_orig[, "Lat"], Lon = trip_orig[, "Lon"], Method = "TwilightFree")
  }
  
  # Metrics
  calc_metrics <- function(df) {
    t_lat <- approx(as.numeric(track_data$time), track_data$true_lat, xout = as.numeric(df$Time))$y
    t_lon <- approx(as.numeric(track_data$time), track_data$true_lon, xout = as.numeric(df$Time))$y
    errs <- get_dist(df$Lon, df$Lat, t_lon, t_lat)
    list(RMSE = sqrt(mean(errs^2, na.rm=TRUE)), HDI = quantile(errs, c(0.025, 0.975), na.rm=TRUE))
  }
  
  m_inv_guided <- calc_metrics(df_inv_guided)
  m_inv_grid <- calc_metrics(df_inv_grid)
  m_flightr <- calc_metrics(df_flightr)
  m_sgat <- calc_metrics(df_sgat)
  m_orig_tf <- if (!is.null(df_orig_tf)) calc_metrics(df_orig_tf) else list(RMSE=NA, HDI=c(NA,NA))
  
  methods_out <- c("invTF (Guided)", "invTF (Grid)", "FLightR", "SGAT", "TwilightFree")
  rmses <- c(m_inv_guided$RMSE, m_inv_grid$RMSE, m_flightr$RMSE, m_sgat$RMSE, m_orig_tf$RMSE)
  
  list(
    inv_guided = df_inv_guided, inv_grid = df_inv_grid, flightr = df_flightr, sgat = df_sgat, orig_tf = df_orig_tf,
    metrics = data.frame(
      Method = methods_out,
      `RMSE (km)` = round(rmses, 1)
    )
  )
}


## ----cloudy, message=FALSE, warning=FALSE-------------------------------------
res_cloudy <- get_paths("cloudy")
if (!is.null(res_cloudy)) {
  ggplot() +
    geom_path(data = track_data, aes(x = true_lon, y = true_lat), color = "black", size = 1.5, alpha = 0.2) +
    geom_path(data = res_cloudy$inv_guided, aes(x = Lon, y = Lat, color = "invTF (Guided)"), size = 1.0) +
    geom_path(data = res_cloudy$inv_grid, aes(x = Lon, y = Lat, color = "invTF (Grid)"), size = 1.0) +
    geom_path(data = res_cloudy$orig_tf, aes(x = Lon, y = Lat, color = "TwilightFree"), size = 1.0) +
    geom_path(data = res_cloudy$flightr, aes(x = Lon, y = Lat, color = "FLightR"), size = 1.0) +
    geom_path(data = res_cloudy$sgat, aes(x = Lon, y = Lat, color = "SGAT"), size = 1.0) +
    theme_minimal() + labs(title = "Cloudy Scenario", color = "Method") +
    scale_color_manual(values = c("invTF (Guided)" = "#E41A1C", "invTF (Grid)" = "#FF7F00", "TwilightFree" = "#984EA3", "FLightR" = "#377EB8", "SGAT" = "#4DAF4A"))
} else {
  print("Caches not yet ready for Cloudy scenario.")
}


## ----shaded, message=FALSE, warning=FALSE-------------------------------------
res_shaded <- get_paths("shaded")
if (!is.null(res_shaded)) {
  ggplot() +
    geom_path(data = track_data, aes(x = true_lon, y = true_lat), color = "black", size = 1.5, alpha = 0.2) +
    geom_path(data = res_shaded$inv_guided, aes(x = Lon, y = Lat, color = "invTF (Guided)"), size = 1.0) +
    geom_path(data = res_shaded$inv_grid, aes(x = Lon, y = Lat, color = "invTF (Grid)"), size = 1.0) +
    geom_path(data = res_shaded$orig_tf, aes(x = Lon, y = Lat, color = "TwilightFree"), size = 1.0) +
    geom_path(data = res_shaded$flightr, aes(x = Lon, y = Lat, color = "FLightR"), size = 1.0) +
    geom_path(data = res_shaded$sgat, aes(x = Lon, y = Lat, color = "SGAT"), size = 1.0) +
    theme_minimal() + labs(title = "Shaded ARS Scenario", color = "Method") +
    scale_color_manual(values = c("invTF (Guided)" = "#E41A1C", "invTF (Grid)" = "#FF7F00", "TwilightFree" = "#984EA3", "FLightR" = "#377EB8", "SGAT" = "#4DAF4A"))
}


## ----alan, message=FALSE, warning=FALSE---------------------------------------
res_alan <- get_paths("alan")
if (!is.null(res_alan)) {
  ggplot() +
    geom_path(data = track_data, aes(x = true_lon, y = true_lat), color = "black", size = 1.5, alpha = 0.2) +
    geom_path(data = res_alan$inv_guided, aes(x = Lon, y = Lat, color = "invTF (Guided)"), size = 1.0) +
    geom_path(data = res_alan$inv_grid, aes(x = Lon, y = Lat, color = "invTF (Grid)"), size = 1.0) +
    geom_path(data = res_alan$orig_tf, aes(x = Lon, y = Lat, color = "TwilightFree"), size = 1.0) +
    geom_path(data = res_alan$flightr, aes(x = Lon, y = Lat, color = "FLightR"), size = 1.0) +
    geom_path(data = res_alan$sgat, aes(x = Lon, y = Lat, color = "SGAT"), size = 1.0) +
    theme_minimal() + labs(title = "ALAN + Shading Scenario", color = "Method") +
    scale_color_manual(values = c("invTF (Guided)" = "#E41A1C", "invTF (Grid)" = "#FF7F00", "TwilightFree" = "#984EA3", "FLightR" = "#377EB8", "SGAT" = "#4DAF4A"))
}


## ----sum_cloudy---------------------------------------------------------------
if(!is.null(res_cloudy)) knitr::kable(res_cloudy$metrics)


## ----sum_shaded---------------------------------------------------------------
if(!is.null(res_shaded)) knitr::kable(res_shaded$metrics)


## ----sum_alan-----------------------------------------------------------------
if(!is.null(res_alan)) knitr::kable(res_alan$metrics)

