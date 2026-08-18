# THE (lambda, temper) PLANE. Targeting COVERAGE, the last unexplained failure.
#
# WHY A PLANE AND NOT ANOTHER AXIS. Two sweeps already exist and were never put side
# by side. They move bias in OPPOSITE directions:
#
#   lambda at temper 1 (stage 2, 29 tags)      temper at lambda 1 (8 tags, partial)
#     4.00  cover 0.421  bias -0.65              1.00  cover 0.679  bias +0.05
#     1.00  cover 0.656  bias -1.54              0.50  cover 0.729  bias +1.63
#     0.50  cover 0.779  bias -1.69              0.30  cover 0.654  bias +3.69
#     0.25  cover 0.857  bias -1.68              0.15  cover 0.454  bias +7.87
#
# Loosening lambda drives bias SOUTH and improves coverage; tempering drives bias
# NORTH, hard. So they are not degenerate, and there should be a point in the plane
# where the bias contributions cancel while the coverage gain survives. No
# single-axis sweep can find it, which is why weeks of them did not.
#
# WHY IT IS PRINCIPLED, not curve-fitting. `lambda` is currently doing two jobs: the
# per-observation noise scale -- measured at 4 to 9 times TIGHTER than shipped
# (`diag_emission_spread.R`, on the pipeline's own calibration, clamped bands
# excluded) -- and the independence correction, which the thinning experiment puts at
# an effective sample size near a quarter of the count. One parameter can only
# compromise between them. `temper` exists to carry the second and was added with
# exactly that rationale.
#
# PREDICTION, WRITTEN BEFORE RUNNING. Reading increments off the two tables above,
# lambda 0.25 with temper 0.5 should land bias near -0.1 with coverage above 0.8.
# That would be the first setting to get both. The additivity is approximate -- the
# two sweeps used different tag subsets -- so treat it as an estimate, and treat a
# coverage above 0.8 at |bias| under 0.5 as the success criterion rather than the
# point value.
#
# CAVEAT I CHECKED: to leading order both the likelihood curvature and the spike
# normaliser's tilt scale as lambda*temper, so the pure diagonal is close to inert.
# The plane is worth running because the EMPIRICAL bias responses have opposite sign,
# which means higher-order structure is doing the work.
#
# HOW TEMPER IS PASSED. It is the 7th `likelihood_params` entry and requires the
# 6-entry darkness form first. Setting `shade_ratio_dark` = `shade_ratio` = 2 and
# `prob_slab_dark` = `prob_slab` = 0.10 makes the darkness regime INERT, so the only
# active change is the temper. That is the trick fit_temper_sweep.R used.
#
# CONTROL: the (lambda 1, temper 1) cell must reproduce rescore29's `tangent_areaON`.
# It is run FIRST so a failure stops the run being interpreted.
#
# Cells are ordered so partial results are already useful: control, the prediction,
# then the axes, then the rest.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))
suppressMessages(devtools::load_all(".", quiet = TRUE))

P <- readRDS("analysis/cache/nes/panel_v1.rds")
tags <- P$tags; responses <- P$responses
OUT <- "scratch/nes_calibration/lambda_temper_plane.csv"
STEP_HOURS <- 12; DIFFUSION <- 110; CELL <- 1
SUBSET <- c("2021033", "2021027", "2023037", "2023032", "2021032", "2022041")

CELLS <- rbind(
  data.table(lam = 1.00, tmp = 1.0),   # control
  data.table(lam = 0.25, tmp = 0.5),   # the prediction
  data.table(lam = 0.25, tmp = 1.0),
  data.table(lam = 1.00, tmp = 0.5),
  data.table(lam = 0.50, tmp = 0.5),
  data.table(lam = 0.50, tmp = 1.0),
  data.table(lam = 0.25, tmp = 0.3),
  data.table(lam = 0.50, tmp = 0.3),
  data.table(lam = 1.00, tmp = 0.3))

grid <- rast(xmin = 150, xmax = 250, ymin = 20, ymax = 70,
             resolution = CELL, crs = "EPSG:4326")
values(grid) <- 1

done <- if (file.exists(OUT)) {
  d0 <- fread(OUT, colClasses = list(character = "id")); paste(d0$id, d0$arm)
} else character(0)
cat(sprintf("%d cells x %d tags = %d fits; %d done\n\n",
            nrow(CELLS), length(SUBSET), nrow(CELLS) * length(SUBSET), length(done)))

for (ci in seq_len(nrow(CELLS))) for (tg in intersect(SUBSET, names(tags))) {
  arm <- sprintf("lam%0.2f_t%0.2f", CELLS$lam[ci], CELLS$tmp[ci])
  if (paste(tg, arm) %in% done) next
  g <- tags[[tg]]; r <- responses[[tg]]
  if (is.null(r)) next
  lam <- CELLS$lam[ci] / (r$max_light * 0.5)
  # 7-entry form; darkness regime made inert by matching the light-regime values
  lp <- c(lam, r$max_light, 0.10, 0.35, 2, 0.10, CELLS$tmp[ci])
  message("fit: ", tg, " / ", arm)
  invisible(capture.output(
    f <- TwilightFreeGrid(g$light$time, pmax(0, g$light$light - r$baseline), grid = grid,
      start_lon = unname(g$p0[1]), start_lat = unname(g$p0[2]),
      end_lon = unname(g$p1[1]), end_lat = unname(g$p1[2]),
      step_hours = STEP_HOURS, diffusion = DIFFUSION,
      calibration = r$calibration, likelihood_params = lp)))
  gp <- grid_posterior(f)
  sdp <- numeric(nrow(gp$P))
  for (i in seq_len(nrow(gp$P))) {
    w <- gp$P[i, ]; s <- sum(w)
    if (!is.finite(s) || s <= 0) { sdp[i] <- NA; next }
    w <- w / s; sdp[i] <- sqrt(sum(w * (gp$lat - sum(w * gp$lat))^2))
  }
  tm <- as.POSIXct(f$fit$time, origin = "1970-01-01", tz = "UTC")
  tr <- argos_at(as.data.table(g$argos), tm)
  ep <- f$fit$lat - tr$lat
  e  <- gc_km(lon360(f$fit$lon), f$fit$lat, tr$lon, tr$lat)
  ok <- is.finite(ep) & is.finite(sdp) & sdp > 0
  fwrite(data.table(id = tg, arm = arm, lam = CELLS$lam[ci], tmp = CELLS$tmp[ci],
                    n = sum(ok), median_km = median(e[ok]),
                    bias_lat = mean(ep[ok]), rmse_lat = sqrt(mean(ep[ok]^2)),
                    lat_sd = mean(sdp[ok]),
                    cover_lat = mean(abs(ep[ok]) <= 1.96 * sdp[ok]),
                    log_z = f$log_z), OUT, append = file.exists(OUT))
}

D <- fread(OUT, colClasses = list(character = "id"))
S <- D[, .(n = .N, km = round(median(median_km)), bias = round(mean(bias_lat), 3),
           abs_bias = round(mean(abs(bias_lat)), 3), rmse = round(mean(rmse_lat), 2),
           lat_sd = round(mean(lat_sd), 2), cover = round(mean(cover_lat), 3)),
       by = .(lam, tmp)][order(-lam, -tmp)]
cat("\n=== the (lambda, temper) plane ===\n")
print(as.data.frame(S), row.names = FALSE)

ctl <- S[lam == 1 & tmp == 1]
if (nrow(ctl)) cat(sprintf("\nCONTROL (lam 1, temper 1): km %d, cover %.3f -- rescore29 gave 250 km, 0.593\n",
                           ctl$km, ctl$cover))
best <- S[abs_bias < 0.5][order(-cover)]
cat("\n=== cells meeting |bias| < 0.5, ranked by coverage ===\n")
if (nrow(best)) print(as.data.frame(best), row.names = FALSE) else
  cat("  none yet\n")
cat("\n  success criterion, pre-registered: coverage > 0.8 at |bias| < 0.5.\n")
cat("  No single-axis sweep has reached both (best on lambda: 0.857 at bias -1.68;\n")
cat("  best on temper: 0.729 at bias +1.63).\n")
