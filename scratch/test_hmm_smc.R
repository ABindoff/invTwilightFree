devtools::load_all(".")
library(geosphere)

get_dist <- function(lon1, lat1, lon2, lat2) {
  distHaversine(cbind(lon1, lat1), cbind(lon2, lat2)) / 1000
}

track_data <- readRDS("scratch/simulated_seal_light_scenarios.rds")
light_col <- "light_ideal"

cat("Running 1-state model...\n")
fit_1 <- TwilightFreeSMC(
  date_time = track_data$time, light = track_data[[light_col]],
  start_lat = track_data$true_lat[1], start_lon = track_data$true_lon[1],
  end_lat = track_data$true_lat[nrow(track_data)], end_lon = track_data$true_lon[nrow(track_data)],
  method = "guided", n_particles = 5000,
  diffusion = 60,
  calibration = c(311.57, 3.19),
  likelihood_params = c(1.0, 64, 0.05)
)

t_lat_1 <- approx(as.numeric(track_data$time), track_data$true_lat, xout = as.numeric(fit_1$knot_times))$y
t_lon_1 <- approx(as.numeric(track_data$time), track_data$true_lon, xout = as.numeric(fit_1$knot_times))$y
rmse_1 <- sqrt(mean(get_dist(fit_1$lon, fit_1$lat, t_lon_1, t_lat_1)^2, na.rm=TRUE))
cat(sprintf("1-State RMSE: %.2f km\n", rmse_1))

cat("Running 2-state model (default sticky transitions)...\n")
fit_2 <- TwilightFreeSMC(
  date_time = track_data$time, light = track_data[[light_col]],
  start_lat = track_data$true_lat[1], start_lon = track_data$true_lon[1],
  end_lat = track_data$true_lat[nrow(track_data)], end_lon = track_data$true_lon[nrow(track_data)],
  method = "guided", n_particles = 5000,
  diffusion = c(10, 80),
  calibration = c(311.57, 3.19),
  likelihood_params = c(1.0, 64, 0.05)
)

rmse_2 <- sqrt(mean(get_dist(fit_2$lon, fit_2$lat, t_lon_1, t_lat_1)^2, na.rm=TRUE))
cat(sprintf("2-State RMSE (default transitions): %.2f km\n", rmse_2))

cat("Running 2-state model (rapid switching)...\n")
trans_rapid <- matrix(c(0.5, 0.5, 0.5, 0.5), nrow=2)
fit_3 <- TwilightFreeSMC(
  date_time = track_data$time, light = track_data[[light_col]],
  start_lat = track_data$true_lat[1], start_lon = track_data$true_lon[1],
  end_lat = track_data$true_lat[nrow(track_data)], end_lon = track_data$true_lon[nrow(track_data)],
  method = "guided", n_particles = 5000,
  diffusion = c(10, 80),
  trans_prob = as.numeric(t(trans_rapid)),
  calibration = c(311.57, 3.19),
  likelihood_params = c(1.0, 64, 0.05)
)

rmse_3 <- sqrt(mean(get_dist(fit_3$lon, fit_3$lat, t_lon_1, t_lat_1)^2, na.rm=TRUE))
cat(sprintf("2-State RMSE (rapid switching): %.2f km\n", rmse_3))

cat("Running 2-state model (highly sticky ARS)...\n")
trans_sticky <- matrix(c(0.99, 0.01, 0.05, 0.95), nrow=2, byrow=TRUE)
fit_4 <- TwilightFreeSMC(
  date_time = track_data$time, light = track_data[[light_col]],
  start_lat = track_data$true_lat[1], start_lon = track_data$true_lon[1],
  end_lat = track_data$true_lat[nrow(track_data)], end_lon = track_data$true_lon[nrow(track_data)],
  method = "guided", n_particles = 5000,
  diffusion = c(10, 80),
  trans_prob = as.numeric(t(trans_sticky)),
  calibration = c(311.57, 3.19),
  likelihood_params = c(1.0, 64, 0.05)
)
rmse_4 <- sqrt(mean(get_dist(fit_4$lon, fit_4$lat, t_lon_1, t_lat_1)^2, na.rm=TRUE))
cat(sprintf("2-State RMSE (highly sticky): %.2f km\n", rmse_4))
