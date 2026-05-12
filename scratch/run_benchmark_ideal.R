
library(invTwilightFree)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(MASS)
library(SGAT)
library(TwGeos)
library(FLightR)

# Simulation - IDEAL CONDITIONS
set.seed(42) # Different seed for a fresh comparison
simulate_benchmark_ideal <- function(duration_days = 180) {
  waypoints <- list(c(-45, 140, 30), c(-60, 110, 60), c(-65, 80, 60), c(-45, 140, 30))
  obs_interval_mins <- 10
  n_obs <- duration_days * 24 * (60 / obs_interval_mins)
  times <- seq(as.POSIXct("2024-01-01 00:00:00", tz="UTC"), by = "10 mins", length.out = n_obs)
  lats <- approx(cumsum(c(0, 30, 60, 60, 30)), c(-45, -60, -65, -45, -45), xout = seq(0, 180, length.out = n_obs))$y
  lons <- approx(cumsum(c(0, 30, 60, 60, 30)), c(140, 110, 80, 140, 140), xout = seq(0, 180, length.out = n_obs))$y
  
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
  
  light <- numeric(n_obs)
  for(i in 1:n_obs) {
    z <- get_z(times[i], lons[i], lats[i])
    l <- pmin(pmax(558.5 - 5.818 * z, 0), 64)
    # NO SHADING Spikes
    # Subtle cloud variation (±5% on intensity)
    l <- l * runif(1, 0.95, 1.05)
    # Very low background noise
    l <- l + rnorm(1, 0, 0.1)
    light[i] <- pmax(0, pmin(64, l))
  }
  data.frame(time = times, true_lat = lats, true_lon = lons, light = light)
}

track <- simulate_benchmark_ideal()

# --- invTF ---
cat("Running invTF (Ideal)...\n")
fit_inv_guided <- TwilightFreeSMC(
    date_time = track$time, light = track$light,
    start_lat = track$true_lat[1], start_lon = track$true_lon[1],
    end_lat = track$true_lat[nrow(track)], end_lon = track$true_lon[nrow(track)],
    method = "guided", n_particles = 1000
)

# --- Twilights (Perfect conditions) ---
twl_data <- data.frame(Date = track$time, Light = track$light)
twl <- TwGeos::findTwilights(twl_data, threshold = 2, include = twl_data$Date)

# --- FLightR ---
cat("Running FLightR (Ideal)...\n")
flightr_data <- data.frame(
  datetime = format(track$time, "%Y-%m-%dT%H:%M:%SZ"),
  light = track$light,
  twilight = 0, interp = FALSE, excluded = FALSE
)
for (i in 1:nrow(twl)) {
  idx <- which.min(abs(as.numeric(track$time) - as.numeric(twl$Twilight[i])))
  flightr_data$twilight[idx] <- ifelse(twl$Rise[i], 1, 2)
}
flightr_data$interp[1] <- TRUE 
tmp_file <- "scratch/benchmark_flightr_input_ideal.csv"
write.csv(flightr_data, tmp_file, row.names = FALSE, quote = FALSE)
proc_data <- FLightR::get.tags.data(tmp_file, log.light.borders = c(0, 64), measurement.period = 600)
Grid <- FLightR::make.grid(left = 70, right = 180, bottom = -80, top = 0, distance.from.land.allowed.to.use = c(-Inf, Inf))
calib_df <- data.frame(calibration.start = track$time[1], calibration.stop = track$time[1] + 3*86400, lon = track$true_lon[1], lat = track$true_lat[1])
calib_flightr <- FLightR::make.calibration(proc_data, Calibration.periods = calib_df)
prerun <- FLightR::make.prerun.object(proc_data, Grid, start = c(track$true_lon[1], track$true_lat[1]), Calibration = calib_flightr)
fit_flightr <- FLightR::run.particle.filter(prerun, threads = 1, nParticles = 50000, plot = FALSE)
saveRDS(list(fit = fit_flightr, time = 0), "vignettes/cache_flightr_ideal.rds")
cat("FLightR complete.\n")

flightr_lon <- fit_flightr$Results$Quantiles$Medianlon
flightr_lat <- fit_flightr$Results$Quantiles$Medianlat
flightr_time <- fit_flightr$Results$Quantiles$time

# --- SGAT ---
cat("Running SGAT (Ideal)...\n")
# Start with truth as this is "Ideal"
true_lat_sgat_init <- approx(as.numeric(track$time), track$true_lat, xout = as.numeric(twl$Twilight))$y
true_lon_sgat_init <- approx(as.numeric(track$time), track$true_lon, xout = as.numeric(twl$Twilight))$y
x0 <- cbind(true_lon_sgat_init, true_lat_sgat_init)

z0 <- twl$Twilight
beta <- c(2.2, 0.08)
logp.x <- function(x) dnorm(x[,2], 0, 90, log = TRUE)
alpha <- c(0, 5) # Much tighter error SD for ideal conditions
model <- SGAT::thresholdModel(twl$Twilight, twl$Rise, twilight.model = "Normal", alpha = alpha, beta = beta, logp.x = logp.x, x0 = x0, z0 = z0)
proposal <- SGAT::mvnorm(S = diag(c(0.05, 0.05)^2), n = nrow(x0))
fit_sgat <- SGAT::stellaMetropolis(model, proposal, x0 = x0, iters = 100000, thin = 100)
saveRDS(list(fit = fit_sgat, time = 0), "vignettes/cache_sgat_ideal.rds")
cat("SGAT complete.\n")

sgat_chain <- fit_sgat$x[[1]]
sgat_lon <- apply(sgat_chain[,1,], 1, mean)
sgat_lat <- apply(sgat_chain[,2,], 1, mean)
sgat_time <- twl$Twilight

# --- Accuracy & Plotting ---
cat("Finalizing plots...\n")
true_lat_k <- approx(as.numeric(track$time), track$true_lat, xout = fit_inv_guided$knot_times)$y
true_lon_k <- approx(as.numeric(track$time), track$true_lon, xout = fit_inv_guided$knot_times)$y
true_lat_flightr <- approx(as.numeric(track$time), track$true_lat, xout = as.numeric(flightr_time))$y
true_lon_flightr <- approx(as.numeric(track$time), track$true_lon, xout = as.numeric(flightr_time))$y
true_lat_sgat <- approx(as.numeric(track$time), track$true_lat, xout = as.numeric(sgat_time))$y
true_lon_sgat <- approx(as.numeric(track$time), track$true_lon, xout = as.numeric(sgat_time))$y

df_err <- rbind(
  data.frame(Time = fit_inv_guided$knot_times, Lat_Error = fit_inv_guided$lat - true_lat_k, Lon_Error = fit_inv_guided$lon - true_lon_k, Method = "invTF (Guided)"),
  data.frame(Time = flightr_time, Lat_Error = flightr_lat - true_lat_flightr, Lon_Error = flightr_lon - true_lon_flightr, Method = "FLightR"),
  data.frame(Time = sgat_time, Lat_Error = sgat_lat - true_lat_sgat, Lon_Error = sgat_lon - true_lon_sgat, Method = "SGAT")
)

p1 <- ggplot(df_err, aes(x = Time, y = Lat_Error, color = Method)) +
  geom_line(alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_minimal() +
  coord_cartesian(ylim = c(-5, 5)) + # Tighter limits for ideal case
  labs(title = "Latitude Error (Ideal Conditions)", y = "Degrees")

p2 <- ggplot(df_err, aes(x = Time, y = Lon_Error, color = Method)) +
  geom_line(alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_minimal() +
  coord_cartesian(ylim = c(-5, 5)) +
  labs(title = "Longitude Error (Ideal Conditions)", y = "Degrees")

final_plot_ideal <- p1 / p2 + plot_layout(guides = "collect")
ggsave("scratch/benchmark_ideal_highres.png", final_plot_ideal, width = 12, height = 8)

cat("Ideal benchmark complete. Plots saved.\n")
