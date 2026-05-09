# Simulate ground-truth tracks for invTwilightFree

library(usethis)
library(SGAT) # for zenith calculations

simulate_track <- function(name, start_date, duration_days, waypoints, noise_prob = 0.05) {
  # waypoints: list of c(lat, lon, days_at_location)
  
  obs_interval_mins <- 10
  n_obs <- duration_days * 24 * (60 / obs_interval_mins)
  times <- seq(as.POSIXct(start_date, tz="UTC"), 
               by = paste(obs_interval_mins, "mins"), 
               length.out = n_obs)
               
  # Interpolate path
  total_duration <- sum(sapply(waypoints, function(x) x[3]))
  if (total_duration != duration_days) {
    # Scale to fit
    scale <- duration_days / total_duration
    for(i in seq_along(waypoints)) waypoints[[i]][3] <- waypoints[[i]][3] * scale
  }
  
  lats <- numeric(n_obs)
  lons <- numeric(n_obs)
  
  idx <- 1
  for(i in 1:(length(waypoints)-1)) {
    p1 <- waypoints[[i]]
    p2 <- waypoints[[i+1]]
    
    steps <- floor((p1[3] * 24 * 60) / obs_interval_mins)
    if(i == length(waypoints)-1) steps <- n_obs - idx + 1
    
    # Linear interpolation in lat/lon for simplicity of simulation
    lats[idx:(idx+steps-1)] <- seq(p1[1], p2[1], length.out = steps)
    lons[idx:(idx+steps-1)] <- seq(p1[2], p2[2], length.out = steps)
    idx <- idx + steps
  }
  
  get_solar_zenith <- function(t, lon, lat) {
    phi <- lat * pi / 180
    d <- (as.numeric(t) - 946728000) / 86400
    g <- (357.529 + 0.98560028 * d) * pi / 180
    q <- 280.459 + 0.98564736 * d
    l_sun <- (q + 1.915 * sin(g) + 0.020 * sin(2 * g)) * pi / 180
    e <- (23.439 - 0.00000036 * d) * pi / 180
    delta <- asin(sin(e) * sin(l_sun))
    ra <- atan2(cos(e) * sin(l_sun), cos(l_sun))
    gmst <- (18.697374558 + 24.06570982441908 * d) * 15 * pi / 180
    lmst <- gmst + lon * pi / 180
    h <- lmst - ra
    cos_z <- sin(phi) * sin(delta) + cos(phi) * cos(delta) * cos(h)
    return(acos(cos_z) * 180 / pi)
  }
  
  light <- numeric(n_obs)
  for(i in 1:n_obs) {
    z <- get_solar_zenith(times[i], lons[i], lats[i])
    l <- 64 - 0.7 * z
    l <- min(max(l, 0), 64)
    light[i] <- l
  }
  
  # Add noise (shading / censoring)
  is_noise <- runif(n_obs) < noise_prob
  light[is_noise] <- pmax(0, light[is_noise] - runif(sum(is_noise), 10, 40))
  
  # Add complete block censoring
  censor_blocks <- sample(1:n_obs, 3)
  for(cb in censor_blocks) {
    end_cb <- min(n_obs, cb + sample(20:100, 1))
    light[cb:end_cb] <- 0
  }
  
  df <- data.frame(
    time = times,
    true_lat = lats,
    true_lon = lons,
    light = light
  )
  
  return(df)
}

set.seed(42)

# 1. Short Track (2 weeks)
sim_short <- simulate_track(
  "short", "2024-01-01", 14,
  list(c(-45, 140, 7), c(-48, 145, 7))
)

# 2. Polar Foray (Alaskan Summer - 4 weeks)
sim_polar <- simulate_track(
  "polar", "2024-06-01", 28,
  list(c(55, -135, 7), c(70, -145, 14), c(55, -135, 7))
)

# 3. Equator Crossing at Equinox
sim_equinox <- simulate_track(
  "equinox", "2024-09-10", 30,
  list(c(20, -30, 10), c(-20, -30, 10), c(-40, -30, 10))
)

# 4. Long Track (12 months)
sim_long <- simulate_track(
  "long", "2023-01-01", 365,
  list(c(-45, 140, 60), c(-60, 110, 120), c(-65, 80, 60), c(-45, 140, 125))
)

use_data(sim_short, overwrite = TRUE)
use_data(sim_polar, overwrite = TRUE)
use_data(sim_equinox, overwrite = TRUE)
use_data(sim_long, overwrite = TRUE)
