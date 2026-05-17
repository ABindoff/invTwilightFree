## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 10,
  fig.height = 8,
  out.width = "100%"
)

library(ggplot2)
library(dplyr)
library(geosphere)

# Distance helper
get_dist <- function(lon1, lat1, lon2, lat2) {
  distHaversine(cbind(lon1, lat1), cbind(lon2, lat2)) / 1000
}


## ----generate_easy_track------------------------------------------------------
# Generate a simple 20-day track
set.seed(42)
times <- seq(as.POSIXct("2024-11-01", tz="UTC"), by="5 min", length.out = 20 * 24 * 12)
n <- length(times)

# Hobart to Adelaide-ish
true_lat <- seq(-42.8, -35.0, length.out = n)
true_lon <- seq(147.3, 138.6, length.out = n)

# Generate ideal light data using the package's own solar model
# (This ensures we aren't testing cross-package solar math differences yet)
light_ideal <- invTwilightFree::solar_zenith(as.numeric(times), true_lon, true_lat)
# Convert zenith to a "light" value (e.g. max 64, transition at 96)
light_ideal <- pmax(0, 64 * (96 - light_ideal) / (96 - 85))
light_ideal[light_ideal > 64] <- 64

easy_track <- data.frame(
  time = times,
  true_lat = true_lat,
  true_lon = true_lon,
  light = light_ideal
)


## ----tuning_invtf-------------------------------------------------------------
fit_guided <- TwilightFreeSMC(
  date_time = easy_track$time, 
  light = easy_track$light,
  start_lat = -42.8, start_lon = 147.3,
  end_lat = -35.0, end_lon = 138.6,
  method = "guided", n_particles = 1000,
  diffusion = 50 # km/day
)

# Extract metrics
t_lat <- approx(as.numeric(easy_track$time), easy_track$true_lat, xout = as.numeric(fit_guided$knot_times))$y
t_lon <- approx(as.numeric(easy_track$time), easy_track$true_lon, xout = as.numeric(fit_guided$knot_times))$y
errs <- get_dist(fit_guided$lon, fit_guided$lat, t_lon, t_lat)
rmse_guided <- sqrt(mean(errs^2, na.rm=TRUE))

print(paste("Guided RMSE:", round(rmse_guided, 2), "km"))


## ----tuning_grid--------------------------------------------------------------
grid <- TwilightFree::makeGrid(lon=c(130, 155), lat=c(-50, -30), cell.size=1, mask="sea")
fit_grid <- TwilightFreeGrid(
  date_time = easy_track$time, 
  light = easy_track$light,
  grid = grid,
  start_lat = -42.8, start_lon = 147.3,
  end_lat = -35.0, end_lon = 138.6,
  step_hours = 6.0,
  diffusion = 50
)

trip_mat <- TwilightFree::trip(fit_grid$fit)
# Correct indices: Date, Lon, Lat (based on previous findings)
df_grid <- data.frame(
  Time = as.POSIXct(trip_mat[, 1], origin="1970-01-01", tz="UTC"),
  Lon = trip_mat[, 2],
  Lat = trip_mat[, 3]
)

t_lat_g <- approx(as.numeric(easy_track$time), easy_track$true_lat, xout = as.numeric(df_grid$Time))$y
t_lon_g <- approx(as.numeric(easy_track$time), easy_track$true_lon, xout = as.numeric(df_grid$Time))$y
errs_g <- get_dist(df_grid$Lon, df_grid$Lat, t_lon_g, t_lat_g)
rmse_grid <- sqrt(mean(errs_g^2, na.rm=TRUE))

print(paste("Grid RMSE:", round(rmse_grid, 2), "km"))


## ----plots--------------------------------------------------------------------
fit_guided_df <- data.frame(
  Time = as.POSIXct(fit_guided$knot_times, origin="1970-01-01", tz="UTC"),
  Lat = fit_guided$lat,
  Lon = fit_guided$lon
)

ggplot() +
  geom_path(data = easy_track, aes(x = true_lon, y = true_lat), color = "black", linewidth = 2, alpha = 0.2) +
  geom_path(data = fit_guided_df, aes(x = Lon, y = Lat), color = "cyan", linewidth = 1) +
  geom_path(data = df_grid, aes(x = Lon, y = Lat), color = "green", linewidth = 1) +
  theme_minimal() +
  labs(title = "Easy Track Recovery", subtitle = "Cyan: Guided, Green: Grid, Grey: Truth")

