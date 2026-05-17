
library(invTwilightFree)
library(ggplot2)
library(dplyr)

# Load the new scenarios
track_data <- readRDS("scratch/simulated_seal_light_scenarios.rds")

# Helper for Distance
get_dist <- function(lon1, lat1, lon2, lat2) {
  rad <- pi/180
  a1 <- lat1 * rad; a2 <- lat2 * rad
  b1 <- lon1 * rad; b2 <- lon2 * rad
  dlon <- b2 - b1
  acos(pmin(1, sin(a1)*sin(a2) + cos(a1)*cos(a2)*cos(dlon))) * 6371.0
}

# 1. Fit invTF (Guided) - Auto-calibration
cat("Fitting invTF (Auto-calibration)...\n")
fit_inv_auto <- TwilightFreeSMC(
  date_time = track_data$time, light = track_data$light_ideal,
  start_lat = track_data$true_lat[1], start_lon = track_data$true_lon[1],
  end_lat = track_data$true_lat[nrow(track_data)], end_lon = track_data$true_lon[nrow(track_data)],
  method = "guided", n_particles = 1000,
  diffusion = 100
)

cat("Auto-calibrated parameters:\n")
cat("  Calibration:", fit_inv_auto$calibration, "\n")
cat("  Likelihood:", fit_inv_auto$likelihood_params, "\n")

# 2. Fit invTF (Guided) - True parameters
cat("\nFitting invTF (True parameters)...\n")
true_cal <- c(558.5, 5.818)
true_lik <- c(1.0 / (64 * 0.5), 64, 0.1)

fit_inv <- TwilightFreeSMC(
  date_time = track_data$time, light = track_data$light_ideal,
  start_lat = track_data$true_lat[1], start_lon = track_data$true_lon[1],
  end_lat = track_data$true_lat[nrow(track_data)], end_lon = track_data$true_lon[nrow(track_data)],
  method = "guided", n_particles = 1000,
  diffusion = 100,
  calibration = true_cal, likelihood_params = true_lik
)

# 2. Calculate Errors
df_inv <- data.frame(
  Time = as.POSIXct(fit_inv$knot_times, origin="1970-01-01", tz="UTC"),
  Lat = fit_inv$lat,
  Lon = fit_inv$lon
)

df_inv <- df_inv %>%
  mutate(
    true_lat = approx(as.numeric(track_data$time), track_data$true_lat, xout = as.numeric(Time))$y,
    true_lon = approx(as.numeric(track_data$time), track_data$true_lon, xout = as.numeric(Time))$y,
    lat_err = Lat - true_lat,
    lon_err = Lon - true_lon,
    dist_err = get_dist(Lon, Lat, true_lon, true_lat)
  )

rmse <- sqrt(mean(df_inv$dist_err^2, na.rm=TRUE))
cat("\n--- Recovery Results (Cloudy) ---\n")
cat("RMSE (km):", round(rmse, 2), "\n")
cat("Mean Lat Error (deg):", round(mean(df_inv$lat_err, na.rm=TRUE), 3), "\n")
cat("Mean Lon Error (deg):", round(mean(df_inv$lon_err, na.rm=TRUE), 3), "\n")

# Plot
p <- ggplot() +
  geom_path(data = track_data, aes(x = true_lon, y = true_lat), color = "black", size = 1.5, alpha = 0.2) +
  geom_path(data = df_inv, aes(x = Lon, y = Lat), color = "red", size = 1) +
  theme_minimal() + labs(title = paste("Recovery Test - RMSE:", round(rmse, 2), "km"))
ggsave("scratch/test_recovery_cloudy.png", p, width = 8, height = 8)
cat("Plot saved to scratch/test_recovery_cloudy.png\n")
