test_that("SMC recovers sim_short track RMSE < 50 km", {
  data("sim_short", package = "invTwilightFree", envir = environment())

  fit <- TwilightFreeSMC(
    date_time   = sim_short$time,
    light       = sim_short$light,
    start_lat   = sim_short$true_lat[1],
    start_lon   = sim_short$true_lon[1],
    end_lat     = sim_short$true_lat[nrow(sim_short)],
    end_lon     = sim_short$true_lon[nrow(sim_short)],
    n_particles = 500,
    step_hours  = 12,
    diffusion   = 50,
    method      = "ffbs",
    seed        = 42
  )

  # Interpolate truth to knot times
  knot_t    <- fit$knot_times
  truth_lat <- approx(as.numeric(sim_short$time), sim_short$true_lat, xout = knot_t)$y
  truth_lon <- approx(as.numeric(sim_short$time), sim_short$true_lon, xout = knot_t)$y

  # Great-circle distance in km using Haversine
  earth_r <- 6371
  lat1 <- fit$lat * pi / 180
  lat2 <- truth_lat * pi / 180
  lon1 <- fit$lon * pi / 180
  lon2 <- truth_lon * pi / 180
  a <- sin((lat2 - lat1) / 2)^2 + cos(lat1) * cos(lat2) * sin((lon2 - lon1) / 2)^2
  dist_km <- 2 * earth_r * asin(sqrt(a))

  rmse_km <- sqrt(mean(dist_km^2, na.rm = TRUE))

  # Calibration vignette uses < 50 km as the passing bar
  expect_lt(rmse_km, 50,
            label = sprintf("RMSE = %.1f km (target < 50 km)", rmse_km))
})
