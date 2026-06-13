# End-to-end integration test: bathy floor_rule removes track excursions over
# implausibly shallow water.
#
# Physical setup:
#   Light observations are generated at lon=-20, lat=30 (western half of the
#   grid) on June 1-3 (northern solstice, unambiguous lat). The Grid HMM light-
#   only MAP therefore sits in the western half. The same animal recorded 500 m
#   max dive depths at every knot — impossible in water only 50 m deep. Adding a
#   floor_rule() bathy term assigns -Inf to every western cell (depth 50 m) and
#   0 to every eastern cell (depth 2000 m), flipping the MAP entirely to the
#   eastern (deep) half.
#
# Regression value: ensures build_aux_matrix -> as.numeric(t()) -> Rust indexing
# round-trip is correct for a depth raster, and that floor_rule integrates with
# TwilightFreeGrid end-to-end.

test_that("bathy floor_rule removes MAP excursions over shallow water", {
  skip_if_not_installed("raster")

  # --- Noiseless obs from lon=-20, lat=30 (western/shallow side) ---
  # June solstice: latitude is unambiguous. Longitude is encoded in twilight
  # timing (40 deg = 2h40min offset between west and east halves of the grid).
  # With noiseless obs the light-only MAP is tightly pinned to the west.
  times  <- seq(
    as.POSIXct("2024-06-01", tz = "UTC"),
    as.POSIXct("2024-06-03", tz = "UTC"),
    by = "10 min"
  )
  t_unix <- as.numeric(times)
  cal    <- c(64, 64 / 90)
  lp     <- c(0.5, 64, 0.05)

  zenith <- solar_zenith(t_unix, rep(-20, length(t_unix)), rep(30, length(t_unix)))
  obs    <- pmax(0, pmin(64, cal[1] - cal[2] * zenith))

  # Grid: equal extent west (lon < 0) and east (lon >= 0) of the prime meridian
  grid <- makeGrid(lon = c(-40, 40), lat = c(25, 35), cell.size = 5)

  # 1x2 bathy raster: west cell = 50 m (shallow), east cell = 2000 m (deep).
  # raster::extract() assigns any lon < 0 to the west cell (value 50) and any
  # lon >= 0 to the east cell (value 2000), matching grid cell centres exactly.
  bathy_r <- raster::raster(
    matrix(c(50, 2000), nrow = 1),
    xmn = -40, xmx = 40, ymn = 25, ymx = 35
  )

  # ---- Light only ----
  fit_lo <- TwilightFreeGrid(
    date_time         = times,
    light             = obs,
    grid              = grid,
    calibration       = cal,
    likelihood_params = lp,
    step_hours        = 12,
    diffusion         = 200
  )

  # ---- Light + bathy floor (500 m max dive, hard constraint) ----
  # Dense dive tag: 500 m reading at every light observation time so that
  # reduce="max" returns 500 for every knot interval.
  dive_df    <- data.frame(time = times, value = 500)
  term_bathy <- location_term(
    "bathy",
    tag    = dive_df,
    reduce = "max",
    source = bathy_source(rast = bathy_r),
    rule   = floor_rule()
  )
  fit_bathy <- TwilightFreeGrid(
    date_time         = times,
    light             = obs,
    grid              = grid,
    calibration       = cal,
    likelihood_params = lp,
    step_hours        = 12,
    diffusion         = 200,
    terms             = list(term_bathy)
  )

  lon_lo    <- fit_lo$fit$lon
  lon_bathy <- fit_bathy$fit$lon

  # Light-only MAP is pinned to the western (shallow) half
  expect_true(all(lon_lo < 0),
    label = paste("light-only lon should all be <0; got",
                  paste(round(lon_lo, 1), collapse = " ")))

  # Bathy floor forces every knot into the eastern (deep) half
  expect_true(all(lon_bathy >= 0),
    label = paste("bathy-constrained lon should all be >=0; got",
                  paste(round(lon_bathy, 1), collapse = " ")))

  # The westward-to-eastward shift must be substantial (not numerical noise).
  # The exact value depends on grid resolution; 15 is a conservative lower bound.
  expect_gt(mean(lon_bathy) - mean(lon_lo), 15)
})
