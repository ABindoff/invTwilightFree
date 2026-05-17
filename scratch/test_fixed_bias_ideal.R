devtools::load_all(".")
library(geosphere)

# Paths
track_path <- "scratch/simulated_seal_light_scenarios.rds"
md_path <- "scratch/fixed_bias_results.md"

track_data <- readRDS(track_path)

get_dist <- function(lon1, lat1, lon2, lat2) {
  distHaversine(cbind(lon1, lat1), cbind(lon2, lat2)) / 1000
}

cat("# Bias-Corrected Evaluation Log\n\n", file = md_path)
cat("| Scenario | RMSE (km) | Time (s) |\n", file = md_path, append = TRUE)
cat("|----------|-----------|----------|\n", file = md_path, append = TRUE)

trans_sticky <- matrix(c(0.9, 0.1, 0.1, 0.9), nrow=2)
trans_flat <- as.numeric(t(trans_sticky))

scenarios <- c("ideal", "cloudy", "shaded", "alan")

for (scen in scenarios) {
  light_col <- if (scen == "cloudy") "light_ideal" else paste0("light_", scen)
  
  cat(sprintf("Evaluating %s with fixed calibration and hyperprior...\n", toupper(scen)))
  
  t_run <- system.time({
    fit <- TwilightFreeSMC(
      date_time = track_data$time, light = track_data[[light_col]],
      start_lat = track_data$true_lat[1], start_lon = track_data$true_lon[1],
      end_lat = track_data$true_lat[nrow(track_data)], end_lon = track_data$true_lon[nrow(track_data)],
      method = "guided", n_particles = 5000,
      diffusion = c(15, 100), 
      trans_prob = trans_flat,
      calibration = c(558.5, 5.818), # THE FIX!
      likelihood_params = c(1.0, 64, 2.0, 20.0) # THE HYPERPRIOR!
    )
  })
  
  t_lat <- approx(as.numeric(track_data$time), track_data$true_lat, xout = as.numeric(fit$knot_times))$y
  t_lon <- approx(as.numeric(track_data$time), track_data$true_lon, xout = as.numeric(fit$knot_times))$y
  safe_lat <- pmax(-90, pmin(90, fit$lat))
  rmse <- sqrt(mean(get_dist(fit$lon, safe_lat, t_lon, t_lat)^2, na.rm=TRUE))
  
  cat(sprintf("| %s | %.2f | %.1f |\n", toupper(scen), rmse, t_run["elapsed"]), file = md_path, append = TRUE)
}

cat("Finished evaluating all scenarios.\n")
