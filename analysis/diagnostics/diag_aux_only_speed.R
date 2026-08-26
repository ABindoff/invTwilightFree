# IS AN AUX-ONLY GRID HMM A FAST, CHEAP TRACK?
#
# The whole case for wiring the edge estimator into this engine rests on the HMM
# being affordable once the light term is gone. The Rust says it should be:
#
#   "This runs BEFORE the light loop so a hard-masked cell can skip the
#    per-observation light likelihood, which is the dominant cost."
#
# The light term costs n_cells x n_obs zenith-and-density evaluations: for a
# 240-day deployment that is about 5000 x 11000 = 5.5e7. An aux-only run replaces
# all of it with a precomputed k x n matrix lookup, leaving only the forward-backward
# and the movement kernel. If the comment is right the difference should be large.
#
# Measured here on the SAME window and grid, calling run_grid_hmm directly so the
# light vectors can be emptied, which the R wrapper has no path for. knot_times and
# obs_times are separate arguments in the Rust, so this needs no engine change --
# only the wrapper does, and only for ergonomics.
suppressMessages({ library(data.table) })
source("analysis/diagnostics/nes_common.R")
suppressMessages(devtools::load_all(".", quiet = TRUE))

P <- readRDS("analysis/cache/nes/panel_v1.rds")
TG <- "2021023"; DAYS <- 40; CELL <- 1; STEP_H <- 12
g <- P$tags[[TG]]; r <- P$responses[[TG]]
L <- as.data.table(g$light)
t0 <- min(as.numeric(L$time)) + 60 * 86400
W <- L[as.numeric(time) >= t0 & as.numeric(time) < t0 + DAYS * 86400]

grid <- rast(xmin = 150, xmax = 250, ymin = 20, ymax = 70, resolution = CELL,
             crs = "EPSG:4326")
values(grid) <- 1
xy <- terra::xyFromCell(grid, seq_len(terra::ncell(grid)))
lon_vec <- xy[, 1]; lat_vec <- xy[, 2]
n <- length(lon_vec)
knot_times <- seq(min(as.numeric(W$time)), max(as.numeric(W$time)), by = STEP_H * 3600)
K <- length(knot_times)
cat(sprintf("tag %s | %d days | %d light obs | %d knots | %d cells\n",
            TG, DAYS, nrow(W), K, n))
cat(sprintf("light term cost ~ cells x obs = %.2e evaluations\n\n", n * nrow(W)))

diffusion <- 110; trans_prob <- 1
obs <- pmax(0, W$light - r$baseline)
lp <- c(1.0 / (r$max_light * 0.5), r$max_light, 0.10)

run <- function(with_light, aux) {
  run_grid_hmm(
    lon = lon_vec, lat = lat_vec, knot_times = knot_times,
    obs_times = if (with_light) as.numeric(W$time) else numeric(0),
    obs_light = if (with_light) as.numeric(obs) else numeric(0),
    fixed_idx = integer(0), fixed_lon = numeric(0), fixed_lat = numeric(0),
    diffusion = diffusion, diffusion_lon = numeric(0), trans_prob = trans_prob,
    calibration = as.numeric(r$calibration), likelihood_params = lp,
    shade_ratio = 2, aux_logl = aux, area_correction = TRUE,
    lambda_scale = numeric(0), drift_correction = FALSE)
}

# A stand-in for the edge estimator's output: a per-day Gaussian day-length
# likelihood on the latitude grid. Shape and cost are what matter here, not
# accuracy, so this uses the TRUE day length plus noise rather than a real detector.
d2r <- pi/180
A <- as.data.table(g$argos)[is.finite(lat)]
lat_at <- approxfun(as.numeric(A$time), A$lat, rule = 2)
dayl <- function(phi, dec, z0 = 91.92) {
  x <- (cos(z0*d2r) - sin(phi*d2r)*sin(dec*d2r)) / (cos(phi*d2r)*cos(dec*d2r))
  2*acos(pmin(1, pmax(-1, x)))/d2r/15
}
aux <- matrix(0.0, nrow = K, ncol = n)
set.seed(1)
for (k in seq_len(K)) {
  tk <- as.POSIXct(knot_times[k], origin = "1970-01-01", tz = "UTC")
  dec <- solar_declination(tk)
  obs_dl <- dayl(lat_at(knot_times[k]), dec) + rnorm(1, 0, 12.4/60)   # 12.4 min SD
  pred <- dayl(lat_vec, dec)
  ll <- -0.5 * ((obs_dl - pred) * 60 / 12.4)^2
  ll[!is.finite(ll)] <- -Inf
  aux[k, ] <- ll
}
aux_flat <- as.numeric(t(aux))

cat("=== timing, same grid and knots ===\n")
t1 <- system.time(f1 <- run(TRUE, numeric(0)))[["elapsed"]]
cat(sprintf("  light only        : %7.2f s   log_z %10.1f\n", t1, f1$log_z))
t2 <- system.time(f2 <- run(FALSE, aux_flat))[["elapsed"]]
cat(sprintf("  aux only          : %7.2f s   log_z %10.1f\n", t2, f2$log_z))
t3 <- system.time(f3 <- run(TRUE, aux_flat))[["elapsed"]]
cat(sprintf("  light + aux       : %7.2f s   log_z %10.1f\n", t3, f3$log_z))
cat(sprintf("\n  speedup, aux-only vs light-only: %.1fx\n", t1 / t2))
cat(sprintf("  extrapolated to a 240-day deployment: light %.1f min -> aux %.1f s\n",
            t1 * (240/DAYS) / 60, t2 * (240/DAYS)))
cat(sprintf("  whole 29-tag panel, aux only: about %.1f min\n", t2 * (240/DAYS) * 29 / 60))
cat("\n  (log_z is NOT comparable across rows: different observation sets)\n")
