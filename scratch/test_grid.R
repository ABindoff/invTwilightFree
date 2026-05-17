devtools::load_all(".")
library(geosphere)

get_dist <- function(lon1, lat1, lon2, lat2) {
  distHaversine(cbind(lon1, lat1), cbind(lon2, lat2)) / 1000
}

track_data <- readRDS("scratch/simulated_seal_light_scenarios.rds")
light_col <- "light_ideal"

cat("Creating grid...\n")
g <- TwilightFree::makeGrid(lon = c(-180, 180), lat = c(-80, -20), cell.size = 2)

run_test <- function(name, diff_vec, trans_mat = NULL) {
  cat(sprintf("\nRunning Grid HMM: %s...\n", name))
  
  if (!is.null(trans_mat)) trans_mat <- as.numeric(t(trans_mat))
  
  t_fit <- system.time({
    fit_grid <- TwilightFreeGrid(
      date_time = track_data$time, light = track_data[[light_col]],
      grid = g,
      start_lat = track_data$true_lat[1], start_lon = track_data$true_lon[1],
      end_lat = track_data$true_lat[nrow(track_data)], end_lon = track_data$true_lon[nrow(track_data)],
      step_hours = 12.0,
      diffusion = diff_vec,
      trans_prob = trans_mat,
      calibration = c(311.57, 3.19),
      likelihood_params = c(1.0, 64, 0.05)
    )
  })
  
  t_lat <- approx(as.numeric(track_data$time), track_data$true_lat, xout = as.numeric(fit_grid$fit$time))$y
  t_lon <- approx(as.numeric(track_data$time), track_data$true_lon, xout = as.numeric(fit_grid$fit$time))$y
  rmse <- sqrt(mean(get_dist(fit_grid$fit$lon, fit_grid$fit$lat, t_lon, t_lat)^2, na.rm=TRUE))
  
  cat(sprintf("RMSE: %.2f km (Took %.1f seconds)\n", rmse, t_fit["elapsed"]))
}

run_test("1-State Model", 60)

run_test("2-State Model (default sticky transitions)", c(10, 80))

trans_rapid <- matrix(c(0.5, 0.5, 0.5, 0.5), nrow=2)
run_test("2-State Model (rapid switching)", c(10, 80), trans_rapid)

trans_sticky <- matrix(c(0.99, 0.01, 0.05, 0.95), nrow=2, byrow=TRUE)
run_test("2-State Model (highly sticky ARS)", c(10, 80), trans_sticky)
