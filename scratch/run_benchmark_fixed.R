
library(invTwilightFree)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(MASS)
library(SGAT)
library(TwGeos)

# Helper for Distance
get_dist <- function(lon1, lat1, lon2, lat2) {
  rad <- pi/180
  a1 <- lat1 * rad; a2 <- lat2 * rad
  b1 <- lon1 * rad; b2 <- lon2 * rad
  dlon <- b2 - b1
  acos(pmin(1, sin(a1)*sin(a2) + cos(a1)*cos(a2)*cos(dlon))) * 6371.0
}

# Simulation
set.seed(123)
simulate_benchmark <- function(duration_days = 180) {
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
    l <- l + rnorm(1, 0, 0.5)
    if(runif(1) < 0.05) l <- pmax(0, l - runif(1, 10, 40))
    light[i] <- pmax(0, pmin(64, l))
  }
  data.frame(time = times, true_lat = lats, true_lon = lons, light = light)
}

track <- simulate_benchmark()

# invTF
fit_inv_guided <- TwilightFreeSMC(
    date_time = track$time, light = track$light,
    start_lat = track$true_lat[1], start_lon = track$true_lon[1],
    end_lat = track$true_lat[nrow(track)], end_lon = track$true_lon[nrow(track)],
    method = "guided", n_particles = 1000
)

# Preprocess for SGAT
twl_data <- data.frame(Date = track$time, Light = track$light)
twl <- TwGeos::findTwilights(twl_data, threshold = 5, include = twl_data$Date)

# SGAT
path <- SGAT::thresholdPath(twl$Twilight, twl$Rise, unfold = FALSE)
x0 <- path$x
if(any(is.na(x0[,2]))) {
    x0[,2] <- approx(x = which(!is.na(x0[,2])), y = x0[!is.na(x0[,2]), 2], xout = 1:nrow(x0), rule = 2)$y
}
z0 <- twl$Twilight
beta <- c(2, 2)
logp.x <- function(x) dnorm(x[,2], 0, 90, log = TRUE)
alpha_mat <- matrix(c(0, 20), nrow = nrow(twl), ncol = 2, byrow = TRUE)
model <- SGAT::thresholdModel(twl$Twilight, twl$Rise, 
                              twilight.model = "Normal",
                              alpha = alpha_mat, beta = beta,
                              logp.x = logp.x, x0 = x0, z0 = z0)
proposal <- SGAT::mvnorm(S = diag(c(0.1, 0.1)^2), n = nrow(x0))
fit_sgat <- SGAT::stellaMetropolis(model, proposal, x0 = x0, iters = 50000, thin = 50)

sgat_chain <- fit_sgat$x[[1]]
sgat_lon <- apply(sgat_chain[,1,], 1, mean)
sgat_lat <- apply(sgat_chain[,2,], 1, mean)
sgat_time <- twl$Twilight

# Prep for Plot
true_lat_k <- approx(as.numeric(track$time), track$true_lat, xout = fit_inv_guided$knot_times)$y
true_lon_k <- approx(as.numeric(track$time), track$true_lon, xout = fit_inv_guided$knot_times)$y
true_lat_sgat <- approx(as.numeric(track$time), track$true_lat, xout = as.numeric(sgat_time))$y
true_lon_sgat <- approx(as.numeric(track$time), track$true_lon, xout = as.numeric(sgat_time))$y

df_err <- rbind(
  data.frame(Time = fit_inv_guided$knot_times, Lat_Error = fit_inv_guided$lat - true_lat_k, Lon_Error = fit_inv_guided$lon - true_lon_k, Method = "invTF (Guided)"),
  data.frame(Time = sgat_time, Lat_Error = sgat_lat - true_lat_sgat, Lon_Error = sgat_lon - true_lon_sgat, Method = "SGAT")
)

p1 <- ggplot(df_err, aes(x = Time, y = Lat_Error, color = Method)) +
  geom_line(alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(title = "Latitude Error (Estimated - True)", y = "Degrees")

p2 <- ggplot(df_err, aes(x = Time, y = Lon_Error, color = Method)) +
  geom_line(alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(title = "Longitude Error (Estimated - True)", y = "Degrees")

final_plot <- p1 / p2 + plot_layout(guides = "collect")
ggsave("scratch/benchmark_fixed.png", final_plot, width = 12, height = 8)
cat("Plot saved to scratch/benchmark_fixed.png\n")

# Save to cache for vignette
saveRDS(list(fit = fit_sgat, time = 0), "vignettes/cache_sgat.rds")
cat("Cache saved to vignettes/cache_sgat.rds\n")
