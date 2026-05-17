
library(ggplot2)
library(dplyr)

# Macquarie Island Coordinates
colony_lat <- -54.619
colony_lon <- 158.860
target_lat <- -60.0
target_lon <- 175.0

simulate_crw_seal_4min <- function() {
  set.seed(999)
  
  # Parameters
  duration_days <- 40
  obs_per_day <- 360 # 4-min intervals
  total_obs <- duration_days * obs_per_day
  times <- seq(as.POSIXct("2024-01-15", tz="UTC"), by = "4 mins", length.out = total_obs)
  
  lats <- numeric(total_obs)
  lons <- numeric(total_obs)
  states <- numeric(total_obs)
  
  lats[1] <- colony_lat
  lons[1] <- colony_lon
  states[1] <- 1
  current_angle <- atan2(target_lat - colony_lat, target_lon - colony_lon)
  
  # Speeds (per 4-min step)
  # 80 km/day ~ 0.002 deg/step, 120 km/day ~ 0.003 deg/step
  transit_speed_min <- 0.002
  transit_speed_max <- 0.003
  ars_speed     <- 0.0006 
  
  transit_rho <- 0.99
  ars_rho     <- 0.5
  
  for (i in 2:total_obs) {
    day <- i / obs_per_day
    
    # Residency period (Days 0-7)
    if (day <= 7) {
      lats[i] <- colony_lat
      lons[i] <- colony_lon
      states[i] <- 0 # Stationary/Calibrating
      next
    }
    
    # State switching logic (Transit vs ARS)
    if (day <= 15 || day >= 30) {
      states[i] <- 1 # Transit
    } else {
      if (states[i-1] == 2) {
        states[i] <- if (runif(1) < 0.995) 2 else 1
      } else {
        states[i] <- if (runif(1) < 0.95) 1 else 2
      }
    }
    
    bias_angle <- if (day < 30) atan2(target_lat - lats[i-1], target_lon - lons[i-1]) 
                  else atan2(colony_lat - lats[i-1], colony_lon - lons[i-1])
    
    rho <- if (states[i] == 1) transit_rho else ars_rho
    current_angle <- current_angle + (1 - rho) * (bias_angle - current_angle) + rnorm(1, 0, (1 - rho) * 2)
    
    speed <- if (states[i] == 1) runif(1, transit_speed_min, transit_speed_max) else ars_speed
    speed <- speed * runif(1, 0.9, 1.1)
    
    lats[i] <- lats[i-1] + speed * sin(current_angle)
    lons[i] <- lons[i-1] + speed * cos(current_angle) / cos(lats[i-1] * pi / 180)
    
    if (lats[i] < -80) lats[i] <- -80
    if (lons[i] > 180) lons[i] <- 180 - (lons[i] - 180)
  }
  
  data.frame(time = times, true_lat = lats, true_lon = lons, state = states)
}

seal_track_4min <- simulate_crw_seal_4min()
saveRDS(seal_track_4min, "scratch/simulated_seal_track_crw.rds")
cat("4-minute CRW simulation complete.\n")
