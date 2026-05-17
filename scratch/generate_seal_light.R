
library(dplyr)

# Load the 4-min CRW track
track <- readRDS("scratch/simulated_seal_track_crw.rds")
n_obs <- nrow(track)

# Solar Zenith Function
get_z <- function(t, lon, lat) {
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
  phi <- lat * pi / 180
  cos_z <- sin(phi) * sin(delta) + cos(phi) * cos(delta) * cos(h)
  acos(pmax(-1, pmin(1, cos_z))) * 180 / pi
}

# 1. Base Light (4-min resolution)
track$z <- as.numeric(mapply(get_z, track$time, track$true_lon, track$true_lat))
track$light_base <- as.numeric(pmin(pmax(558.5 - 5.818 * track$z, 0), 64))

# 2. Add "Tag Sampling Error" (Max-of-4 approximation)
set.seed(555)
track$light_sampled <- as.numeric(sapply(track$light_base, function(l) {
  samples <- l + rnorm(4, 0, 0.5) 
  max(pmax(0, pmin(64, samples)))
}))

# 3. Weather Simulation
weather_state <- 1
weather_seq <- numeric(n_obs)
for(i in 1:n_obs) {
  if (i %% 360 == 0) {
    weather_state <- if (weather_state == 1) (if(runif(1) < 0.25) 2 else 1) else (if(runif(1) < 0.2) 1 else 2)
  }
  weather_seq[i] <- weather_state
}
track$cloud_mult <- as.numeric(ifelse(weather_seq == 2, runif(n_obs, 0.4, 0.7), runif(n_obs, 0.98, 1.02)))

# Scenario A: Ideal with Cloud Cover
track$light_ideal <- as.numeric(track$light_sampled * track$cloud_mult)
track$light_ideal <- pmax(0, pmin(64, track$light_ideal))

# Scenario B: Shaded during ARS (Diving)
track$light_shaded <- as.numeric(track$light_ideal)
ars_idx <- which(track$state == 2)
shade_events <- runif(length(ars_idx)) < 0.2
track$light_shaded[ars_idx[shade_events]] <- track$light_shaded[ars_idx[shade_events]] * runif(sum(shade_events), 0, 0.1)

# Scenario C: Shaded + ALAN near Macca
track$light_alan <- as.numeric(track$light_shaded)
colony_lat <- -54.619
colony_lon <- 158.860
dist_to_macca <- sqrt((track$true_lat - colony_lat)^2 + (track$true_lon - colony_lon)^2) * 111
near_colony <- which(dist_to_macca < 100)
is_night <- which(track$z > 96)
alan_pool <- intersect(near_colony, is_night)

if (length(alan_pool) > 0) {
  n_alan <- min(length(alan_pool), 20)
  alan_events <- sample(alan_pool, size = n_alan)
  track$light_alan[alan_events] <- runif(n_alan, 5, 25)
}

# Sanity Check
cat("Light Range (Ideal):", range(track$light_ideal), "\n")

saveRDS(track, "scratch/simulated_seal_light_scenarios.rds")
cat("4-minute Light synthesis complete.\n")

