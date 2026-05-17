devtools::load_all(".")
library(geosphere)
library(dplyr)

# Paths
track_path <- "scratch/simulated_seal_light_scenarios.rds"
md_path <- "scratch/comprehensive_optimization.md"

track_data <- readRDS(track_path)

get_dist <- function(lon1, lat1, lon2, lat2) {
  distHaversine(cbind(lon1, lat1), cbind(lon2, lat2)) / 1000
}

# Initialise Markdown
cat("# Comprehensive Parameter Optimization Log\n\n", file = md_path)

# =========================================================================
# STAGE 1: Biological Optimization (Diffusion & Stickiness on Ideal Track)
# =========================================================================
cat("## Stage 1: Biological Prior Optimization (Ideal Track)\n\n", file = md_path, append = TRUE)
cat("Optimizing Diffusion parameters (D1, D2) and Behavioural Stickiness (P(stay)).\n\n", file = md_path, append = TRUE)
cat("| D1 (ARS) | D2 (Transit) | Stickiness | RMSE (km) | Time (s) |\n", file = md_path, append = TRUE)
cat("|----------|--------------|------------|-----------|----------|\n", file = md_path, append = TRUE)

d1_grid <- c(5, 10, 15)
d2_grid <- c(60, 80, 100)
stickiness_grid <- c(0.5, 0.7, 0.9)

stage1_results <- data.frame(D1 = numeric(), D2 = numeric(), Stickiness = numeric(), RMSE = numeric())

cat("Starting Stage 1: Biological Optimization...\n")

for (d1 in d1_grid) {
  for (d2 in d2_grid) {
    if (d1 >= d2) next
    for (p in stickiness_grid) {
      trans_mat <- as.numeric(t(matrix(c(p, 1-p, 1-p, p), nrow=2)))
      
      cat(sprintf("Testing D1 = %d, D2 = %d, Stickiness = %.2f...\n", d1, d2, p))
      
      t_run <- system.time({
        fit <- TwilightFreeSMC(
          date_time = track_data$time, light = track_data$light_ideal,
          start_lat = track_data$true_lat[1], start_lon = track_data$true_lon[1],
          end_lat = track_data$true_lat[nrow(track_data)], end_lon = track_data$true_lon[nrow(track_data)],
          method = "guided", n_particles = 5000,
          diffusion = c(d1, d2), trans_prob = trans_mat,
          calibration = c(311.57, 3.19), likelihood_params = c(1.0, 64, 0.05)
        )
      })
      
      t_lat <- approx(as.numeric(track_data$time), track_data$true_lat, xout = as.numeric(fit$knot_times))$y
      t_lon <- approx(as.numeric(track_data$time), track_data$true_lon, xout = as.numeric(fit$knot_times))$y
      
      safe_lat <- pmax(-90, pmin(90, fit$lat))
      rmse <- sqrt(mean(get_dist(fit$lon, safe_lat, t_lon, t_lat)^2, na.rm=TRUE))
      
      cat(sprintf("| %d | %d | %.2f | %.2f | %.1f |\n", d1, d2, p, rmse, t_run["elapsed"]), file = md_path, append = TRUE)
      stage1_results <- rbind(stage1_results, data.frame(D1=d1, D2=d2, Stickiness=p, RMSE=rmse))
    }
  }
}

best_bio <- stage1_results[which.min(stage1_results$RMSE), ]
best_trans <- as.numeric(t(matrix(c(best_bio$Stickiness, 1-best_bio$Stickiness, 1-best_bio$Stickiness, best_bio$Stickiness), nrow=2)))

cat(sprintf("\n**WINNING BIOLOGICAL MODEL:** D1 = %d, D2 = %d, Stickiness = %.2f (RMSE = %.2f km)\n\n", 
            best_bio$D1, best_bio$D2, best_bio$Stickiness, best_bio$RMSE), file = md_path, append = TRUE)


# =========================================================================
# STAGE 2: Likelihood Optimization (Noise Parameters on Cloudy/Shaded/ALAN)
# =========================================================================
cat("## Stage 2: Spike-and-Slab Likelihood Optimization\n\n", file = md_path, append = TRUE)
cat("Optimizing precision (`lambda`) and outlier probability (`prob_slab`) for noisy tracks.\n\n", file = md_path, append = TRUE)

lambda_grid <- c(0.5, 1.0, 2.0)
prob_slab_grid <- c(0.01, 0.05, 0.15)
scenarios <- c("cloudy", "shaded", "alan")

cat("Starting Stage 2: Likelihood Optimization...\n")

for (scen in scenarios) {
  light_col <- if (scen == "cloudy") "light_ideal" else paste0("light_", scen)
  
  cat(sprintf("### Scenario: %s\n\n", toupper(scen)), file = md_path, append = TRUE)
  cat("| Lambda | Prob Slab | RMSE (km) | Time (s) |\n", file = md_path, append = TRUE)
  cat("|--------|-----------|-----------|----------|\n", file = md_path, append = TRUE)
  
  scen_results <- data.frame(Lambda = numeric(), ProbSlab = numeric(), RMSE = numeric())
  
  for (lam in lambda_grid) {
    for (pslab in prob_slab_grid) {
      cat(sprintf("Testing %s -> Lambda = %.1f, ProbSlab = %.2f...\n", scen, lam, pslab))
      
      t_run <- system.time({
        fit <- TwilightFreeSMC(
          date_time = track_data$time, light = track_data[[light_col]],
          start_lat = track_data$true_lat[1], start_lon = track_data$true_lon[1],
          end_lat = track_data$true_lat[nrow(track_data)], end_lon = track_data$true_lon[nrow(track_data)],
          method = "guided", n_particles = 5000,
          diffusion = c(best_bio$D1, best_bio$D2), trans_prob = best_trans,
          calibration = c(311.57, 3.19), likelihood_params = c(lam, 64, pslab)
        )
      })
      
      t_lat <- approx(as.numeric(track_data$time), track_data$true_lat, xout = as.numeric(fit$knot_times))$y
      t_lon <- approx(as.numeric(track_data$time), track_data$true_lon, xout = as.numeric(fit$knot_times))$y
      safe_lat <- pmax(-90, pmin(90, fit$lat))
      rmse <- sqrt(mean(get_dist(fit$lon, safe_lat, t_lon, t_lat)^2, na.rm=TRUE))
      
      cat(sprintf("| %.1f | %.2f | %.2f | %.1f |\n", lam, pslab, rmse, t_run["elapsed"]), file = md_path, append = TRUE)
      scen_results <- rbind(scen_results, data.frame(Lambda=lam, ProbSlab=pslab, RMSE=rmse))
    }
  }
  
  best_like <- scen_results[which.min(scen_results$RMSE), ]
  cat(sprintf("\n**BEST FOR %s:** Lambda = %.1f, Prob Slab = %.2f (RMSE = %.2f km)\n\n", 
              toupper(scen), best_like$Lambda, best_like$ProbSlab, best_like$RMSE), file = md_path, append = TRUE)
}

cat("Optimization Protocol Complete.\n")
