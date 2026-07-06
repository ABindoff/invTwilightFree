# H1: simulate a panel of N individuals for the hierarchical sampler.
# Each individual is a correlated random walk with its OWN per-step movement
# variance sig2_i (km^2), drawn around a population level. Short mid-latitude
# tracks, away from equinox, so latitude is identified (the hard equinox case is
# a separate axis; here we validate the POOLING). Forward zenith->light model
# matches data-raw/simulate_tracks.R.
suppressPackageStartupMessages(library(invTwilightFree))

# forward solar zenith (same formula as data-raw/simulate_tracks.R / the engine)
zenith_R <- function(t, lon, lat) {
  as.numeric(solar_zenith(as.numeric(t), lon, lat))
}

simulate_panel <- function(N = 6, days = 14, step_min = 10,
                           start_date = "2024-05-01",
                           pop_km_per_day = 45,      # population movement level
                           lat0_range = c(-52, -40), lon0_range = c(135, 150),
                           noise_prob = 0.05, seed = 20240501) {
  set.seed(seed)
  days_i <- if (length(days) == 1L) rep(days, N) else days   # per-individual track length
  stopifnot(length(days_i) == N)
  dt_day <- step_min / (60 * 24)                 # obs interval in days
  km_lat <- 111.0
  # per-individual daily movement sd (km/day), heterogeneous around the pop level
  sd_day_i <- pop_km_per_day * exp(rnorm(N, 0, 0.35))
  ind <- vector("list", N)
  for (i in seq_len(N)) {
    n_obs <- round(days_i[i] * 24 * (60 / step_min))
    times <- seq(as.POSIXct(start_date, tz = "UTC"), by = paste(step_min, "mins"), length.out = n_obs)
    lat0 <- runif(1, lat0_range[1], lat0_range[2]); lon0 <- runif(1, lon0_range[1], lon0_range[2])
    km_lon <- 111.0 * cos(lat0 * pi / 180)
    step_sd_km <- sd_day_i[i] * sqrt(dt_day)     # per-obs-step km sd
    lat <- numeric(n_obs); lon <- numeric(n_obs); lat[1] <- lat0; lon[1] <- lon0
    for (t in 2:n_obs) {
      lat[t] <- lat[t-1] + rnorm(1, 0, step_sd_km / km_lat)
      lon[t] <- lon[t-1] + rnorm(1, 0, step_sd_km / km_lon)
    }
    z <- zenith_R(times, lon, lat)
    light <- pmin(pmax(558.5 - 5.818 * z, 0), 64)
    isn <- runif(n_obs) < noise_prob
    light[isn] <- pmax(0, light[isn] - runif(sum(isn), 10, 40))
    ind[[i]] <- data.frame(time = times, light = light, true_lat = lat, true_lon = lon)
  }
  list(ind = ind, sd_day_true = sd_day_i, pop_km_per_day = pop_km_per_day, days = days_i,
       start = data.frame(lat = sapply(ind, function(d) d$true_lat[1]),
                          lon = sapply(ind, function(d) d$true_lon[1])))
}

if (sys.nframe() == 0L) {
  P <- simulate_panel()
  cat(sprintf("panel: N=%d, obs/ind=%d, pop=%.0f km/day\n",
              length(P$ind), nrow(P$ind[[1]]), P$pop_km_per_day))
  cat("true per-individual movement (km/day):", paste(round(P$sd_day_true, 1), collapse = ", "), "\n")
  cat("start lat range:", paste(round(range(P$start$lat), 1), collapse = " to "),
      "| lon range:", paste(round(range(P$start$lon), 1), collapse = " to "), "\n")
  # quick sanity plot of the tracks + one light series
  grDevices::png("panel_check.png", width = 1100, height = 460)
  graphics::par(mfrow = c(1, 2), mar = c(4,4,2,1))
  plot(NA, xlim = range(sapply(P$ind, function(d) range(d$true_lon))),
       ylim = range(sapply(P$ind, function(d) range(d$true_lat))),
       xlab = "lon", ylab = "lat", main = "simulated panel tracks")
  for (i in seq_along(P$ind)) graphics::lines(P$ind[[i]]$true_lon, P$ind[[i]]$true_lat, col = i)
  plot(P$ind[[1]]$time, P$ind[[1]]$light, type = "l", xlab = "time", ylab = "light",
       main = "individual 1 light series")
  grDevices::dev.off()
  cat("saved panel_check.png\n")
}
