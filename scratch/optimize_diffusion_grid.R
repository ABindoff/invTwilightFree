devtools::load_all(".")
library(geosphere)
library(dplyr)

# Paths
track_path <- "scratch/simulated_seal_light_scenarios.rds"
md_path <- "scratch/diffusion_optimization_grid.md"

track_data <- readRDS(track_path)

get_dist <- function(lon1, lat1, lon2, lat2) {
  distHaversine(cbind(lon1, lat1), cbind(lon2, lat2)) / 1000
}

# Initialise Markdown
cat("# Diffusion Parameter Optimization Log (Grid HMM)\n\n", file = md_path)
cat("## Grid Search on Ideal Scenario\n\n", file = md_path, append = TRUE)
cat("| ARS Diffusion (D1) | Transit Diffusion (D2) | RMSE (km) | Time (s) |\n", file = md_path, append = TRUE)
cat("|--------------------|------------------------|-----------|----------|\n", file = md_path, append = TRUE)

# Parameter Grids
d1_grid <- c(5, 10, 15, 20, 25)
d2_grid <- c(40, 60, 80, 100, 120, 150)
trans_rapid <- as.numeric(t(matrix(c(0.5, 0.5, 0.5, 0.5), nrow=2)))

results <- data.frame(D1 = numeric(), D2 = numeric(), RMSE = numeric(), Time = numeric())

cat("Starting Grid HMM grid search on IDEAL scenario...\n")

g <- TwilightFree::makeGrid(lon=c(130, 185), lat=c(-85, -40), cell.size=1.0, mask="sea")

for (d1 in d1_grid) {
  for (d2 in d2_grid) {
    if (d1 >= d2) next
    
    cat(sprintf("Testing D1 = %d, D2 = %d...\n", d1, d2))
    
    t_run <- system.time({
      fit_inv_grid <- TwilightFreeGrid(
        date_time = as.POSIXct(track_data$time, tz="UTC"), 
        light = track_data$light_ideal, 
        grid = g,
        start_lat = track_data$true_lat[1], start_lon = track_data$true_lon[1],
        end_lat = track_data$true_lat[nrow(track_data)], end_lon = track_data$true_lon[nrow(track_data)],
        step_hours = 4.0,
        diffusion = c(d1, d2),
        trans_prob = trans_rapid,
        calibration = c(311.57, 3.19),
        likelihood_params = c(1.0, 64, 0.05)
      )
    })
    
    t_lat <- approx(as.numeric(track_data$time), track_data$true_lat, xout = as.numeric(fit_inv_grid$fit$time))$y
    t_lon <- approx(as.numeric(track_data$time), track_data$true_lon, xout = as.numeric(fit_inv_grid$fit$time))$y
    
    safe_lat <- pmax(-90, pmin(90, fit_inv_grid$fit$lat))
    errs <- get_dist(fit_inv_grid$fit$lon, safe_lat, t_lon, t_lat)
    rmse <- sqrt(mean(errs^2, na.rm=TRUE))
    
    # Log to Markdown
    cat(sprintf("| %d | %d | %.2f | %.1f |\n", d1, d2, rmse, t_run["elapsed"]), file = md_path, append = TRUE)
    
    results <- rbind(results, data.frame(D1 = d1, D2 = d2, RMSE = rmse, Time = t_run["elapsed"]))
  }
}

# Find best parameters
best_run <- results[which.min(results$RMSE), ]
cat(sprintf("\nBest Parameters Found: D1 = %d, D2 = %d with RMSE %.2f km\n", best_run$D1, best_run$D2, best_run$RMSE))

cat("\n## Evaluation on Realistic Scenarios\n\n", file = md_path, append = TRUE)
cat(sprintf("Using optimized parameters: D1 = %d, D2 = %d\n\n", best_run$D1, best_run$D2), file = md_path, append = TRUE)
cat("| Scenario | RMSE (km) | Time (s) |\n", file = md_path, append = TRUE)
cat("|----------|-----------|----------|\n", file = md_path, append = TRUE)

scenarios <- c("cloudy", "shaded", "alan")

for (scen in scenarios) {
  light_col <- if (scen == "cloudy") "light_ideal" else paste0("light_", scen)
  
  cat(sprintf("Evaluating on %s...\n", scen))
  t_run <- system.time({
    fit_inv_grid <- TwilightFreeGrid(
        date_time = as.POSIXct(track_data$time, tz="UTC"), 
        light = track_data[[light_col]], 
        grid = g,
        start_lat = track_data$true_lat[1], start_lon = track_data$true_lon[1],
        end_lat = track_data$true_lat[nrow(track_data)], end_lon = track_data$true_lon[nrow(track_data)],
        step_hours = 4.0,
        diffusion = c(best_run$D1, best_run$D2),
        trans_prob = trans_rapid,
        calibration = c(311.57, 3.19),
        likelihood_params = c(1.0, 64, 0.05)
      )
  })
  
  t_lat <- approx(as.numeric(track_data$time), track_data$true_lat, xout = as.numeric(fit_inv_grid$fit$time))$y
  t_lon <- approx(as.numeric(track_data$time), track_data$true_lon, xout = as.numeric(fit_inv_grid$fit$time))$y
  
  safe_lat <- pmax(-90, pmin(90, fit_inv_grid$fit$lat))
  errs <- get_dist(fit_inv_grid$fit$lon, safe_lat, t_lon, t_lat)
  rmse <- sqrt(mean(errs^2, na.rm=TRUE))
  
  cat(sprintf("| %s | %.2f | %.1f |\n", toupper(scen), rmse, t_run["elapsed"]), file = md_path, append = TRUE)
}

cat("Optimization and Evaluation Complete.\n")
