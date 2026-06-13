# Verify the adaptive diffusion heuristic is properly parameterised.
# Strategy: use a constant-light series (flat; max_obs - min_obs = 0 < threshold).
# With flat_light_scale=0.1 (default), particles should stay near the start.
# With flat_light_scale=1.0 (disabled), particles are free to wander widely.

make_flat_smc <- function(flat_light_scale, n_particles = 300, seed = 7) {
  times <- seq(
    as.POSIXct("2020-06-01", tz = "UTC"),
    as.POSIXct("2020-06-04", tz = "UTC"),
    by = "10 min"
  )
  # Constant light = perfectly flat; knot range is always 0 < default threshold 10.
  light <- rep(32.0, length(times))

  TwilightFreeSMC(
    date_time          = times,
    light              = light,
    start_lat          = -43.0,
    start_lon          = 147.0,
    n_particles        = n_particles,
    step_hours         = 12,
    diffusion          = 500,       # large diffusion to amplify the contrast
    calibration        = c(64, 64 / 90),
    likelihood_params  = c(0.5, 64, 0.05),
    method             = "forward",
    seed               = seed,
    flat_light_scale   = flat_light_scale,
    flat_light_threshold = 10.0
  )
}

test_that("flat_light_scale=0.1 constrains diffusion on flat-light knots", {
  r_constrained <- make_flat_smc(flat_light_scale = 0.1)
  r_free        <- make_flat_smc(flat_light_scale = 1.0)

  sd_constrained <- sd(r_constrained$lat)
  sd_free        <- sd(r_free$lat)

  # Constrained run should have much less latitude spread
  expect_lt(sd_constrained, sd_free,
            label = "constrained lat SD < free lat SD on flat-light series")
})

test_that("flat_light_threshold=0 disables the heuristic (same as scale=1.0)", {
  # threshold=0 means max_obs-min_obs is never < 0, so diff_scale is always 1.0
  r_threshold0 <- make_flat_smc(flat_light_scale = 0.1)  # scale won't matter
  # With threshold=0 the heuristic is always skipped regardless of scale
  r_disabled   <- TwilightFreeSMC(
    date_time          = seq(as.POSIXct("2020-06-01", tz = "UTC"),
                             as.POSIXct("2020-06-04", tz = "UTC"), by = "10 min"),
    light              = rep(32.0, length(seq(as.POSIXct("2020-06-01", tz = "UTC"),
                                              as.POSIXct("2020-06-04", tz = "UTC"), by = "10 min"))),
    start_lat = -43, start_lon = 147,
    n_particles = 300, step_hours = 12, diffusion = 500,
    calibration = c(64, 64/90), likelihood_params = c(0.5, 64, 0.05),
    method = "forward", seed = 7,
    flat_light_threshold = 0.0,   # disable heuristic
    flat_light_scale     = 0.1    # scale value is irrelevant when threshold=0
  )
  r_free <- make_flat_smc(flat_light_scale = 1.0)   # explicitly disabled via scale

  # Both should wander widely; disabled-threshold run should have similar spread to scale=1
  expect_gt(sd(r_disabled$lat), sd(make_flat_smc(0.1)$lat),
            label = "threshold=0 gives larger spread than constrained run")
})
