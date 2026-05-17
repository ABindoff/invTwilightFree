
library(invTwilightFree)
library(dplyr)

track_data <- readRDS("scratch/simulated_seal_light_scenarios.rds")

# Make sure time is POSIXct
track_data$time <- as.POSIXct(track_data$time, tz="UTC")

true_cal <- c(558.5, 5.818)
true_lik <- c(1.0 / (64 * 0.5), 64, 0.1)

cat("Fitting invTF...\n")
fit_inv <- TwilightFreeSMC(
  date_time = track_data$time, 
  light = track_data$light_ideal, 
  start_lat = track_data$true_lat[1], 
  start_lon = track_data$true_lon[1], 
  end_lat = track_data$true_lat[nrow(track_data)], 
  end_lon = track_data$true_lon[nrow(track_data)], 
  method = 'ffbs', 
  n_particles = 1000, 
  diffusion = 100, 
  calibration = true_cal, 
  likelihood_params = true_lik
)

df <- data.frame(
  Time = as.POSIXct(fit_inv$knot_times, origin='1970-01-01', tz='UTC'), 
  Lat = fit_inv$lat, 
  Lon = fit_inv$lon
)

df <- df %>% 
  mutate(
    true_lon = approx(as.numeric(track_data$time), track_data$true_lon, xout = as.numeric(Time))$y, 
    true_lat = approx(as.numeric(track_data$time), track_data$true_lat, xout = as.numeric(Time))$y,
    lon_err = Lon - true_lon,
    lat_err = Lat - true_lat
  ) %>%
  mutate(
    lon_err = ifelse(lon_err > 180, lon_err - 360, ifelse(lon_err < -180, lon_err + 360, lon_err))
  )

cat("\n--- First 6 records ---\n")
print(head(df))
cat("\n--- Last 6 records ---\n")
print(tail(df))
cat("\n--- Summary of Longitude Error ---\n")
print(summary(df$lon_err))
