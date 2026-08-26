# THE REPORTING RUN. All 29 deployments, four arms, each adding ONE change.
#
# Every change below is measured, not guessed, and each arm isolates one of them so
# the manuscript can report a decomposition rather than a single number.
#
#   A  shipped            REPRODUCTION CONTROL. Must return `rescore29_results.csv`
#                         arm `tangent_areaON`. If it does not, nothing else here is
#                         readable.
#   B  A + depth gate     Drop observations with `depth_min` > 10 m. `depth_min` is a
#                         GATE, not a gradient: median residual is flat at +0.028 /
#                         +0.022 for <2.5 m and 2.5-10 m, then falls to -0.103 /
#                         -0.246 / -0.492 for 10-25 / 25-100 / >100 m. 95.1% of
#                         slab-scale events sit below it. Gating collapses the
#                         SEASONAL emission drift by 92% (swing 0.910 -> 0.076
#                         zenith-equivalent deg) and twilight tail mass to ~0.
#                         Costs almost nothing: 96.3% of observations kept, zero
#                         knots emptied, 5th-pct knot still has 20 observations.
#   C  B + prob_slab 0.01 Measured tail mass on gated light is 0.000-0.001; the
#                         shipped 0.10 is a ~100x overstatement. On gated light the
#                         slab is not compensating for anything, because what it was
#                         absorbing was dive shading.
#   D  C + lambda 0.5     Stage 2 (29 tags) found lam 0.5 beats the shipped 1.0 on
#                         BOTH accuracy and coverage (266 km / 0.779 vs 276 / 0.656).
#                         Validated but never shipped.
#
# PRE-REGISTERED, and two of these could go the wrong way:
#   B: seasonal structure in per-knot bias should shrink. Overall bias may barely
#      move -- the LATITUDE drift is a property of the light itself and gating does
#      not touch it (slope -0.0850 ungated vs -0.0707 gated, only 17%).
#   C: should IMPROVE accuracy (the slab dilutes real information) but may WORSEN
#      coverage, because a lighter slab sharpens the likelihood and narrows the
#      posterior. If both worsen, B is the recommended configuration, not C.
#   D: should improve coverage and accuracy together as it did in stage 2 -- but
#      stage 2 was on UNGATED, heavy-slab light, so the optimum may have moved.
#
# NOT VARIED HERE, deliberately: the gate applies to the ENGINE INPUT only, not to
# the light used to FIT the response. Gating the calibration too would change
# `baseline` (q05) and `max_light` at the same time, breaking the one-change-per-arm
# design. That interaction is a separate question and is untested.
#
# Also not included: pooled-evidence z50 profiling. Validated (recovers a known truth
# exactly; worth ~125 km per degree against a recorded 5.6 deg spread across tags)
# but it needs ~7 fits per tag to profile, which is its own run.
#
# Per-knot errors are written as well as per-tag summaries, so the seasonal and
# latitude structure can be reported rather than just the headline.
#
# RESUMABLE: checkpointed on (id, arm).
SCRATCH <- "analysis/diagnostics"
source(file.path(SCRATCH, "nes_common.R"))
suppressMessages(devtools::load_all(".", quiet = TRUE))

P <- readRDS("analysis/cache/nes/panel_v1.rds")
tags <- P$tags; responses <- P$responses
OUT_TAG  <- "scratch/nes_calibration/final29_tags.csv"
OUT_KNOT <- "scratch/nes_calibration/final29_knots.csv"
STEP_HOURS <- 12; DIFFUSION <- 110; CELL <- 1; GATE_M <- 10

ARMS <- data.table(
  arm   = c("A_shipped", "B_gate", "C_gate_slab", "D_gate_slab_lam"),
  gate  = c(FALSE, TRUE, TRUE, TRUE),
  pslab = c(0.10, 0.10, 0.01, 0.01),
  lam   = c(1.00, 1.00, 1.00, 0.50))

# Optional arm filter, e.g. `Rscript fit_final_29.R B_gate`. No arguments keeps the
# original behaviour of running all four. Added because arm B is now measured as
# HARMFUL (worse accuracy, bias and coverage, dose-dependent, 21 of 23 tags), and
# arms C and D are both built on top of it, so they can only quantify how harmful
# the combination is. Finishing B completes a reportable decomposition; C and D are
# not worth 58 fits. The ungated versions of their changes are the better question.
SEL_ARMS <- commandArgs(trailingOnly = TRUE)
if (length(SEL_ARMS)) {
  ARMS <- ARMS[arm %in% SEL_ARMS]
  if (!nrow(ARMS)) stop("no arm matched: ", paste(SEL_ARMS, collapse = ", "))
  cat("arm filter active:", paste(ARMS$arm, collapse = ", "), "\n")
}

grid <- rast(xmin = 150, xmax = 250, ymin = 20, ymax = 70,
             resolution = CELL, crs = "EPSG:4326")
values(grid) <- 1

done <- if (file.exists(OUT_TAG)) {
  d0 <- fread(OUT_TAG, colClasses = list(character = "id")); paste(d0$id, d0$arm)
} else character(0)
cat(sprintf("%d tags x %d arms = %d fits; %d done\n\n",
            length(tags), nrow(ARMS), length(tags) * nrow(ARMS), length(done)))

for (ai in seq_len(nrow(ARMS))) for (tg in names(tags)) {
  a <- ARMS[ai]
  if (paste(tg, a$arm) %in% done) next
  g <- tags[[tg]]; r <- responses[[tg]]
  if (is.null(r)) { message(tg, ": no response"); next }
  L <- as.data.table(g$light)
  if (a$gate) {
    if (!"depth_min" %in% names(L)) { message(tg, ": no depth channel, skipped"); next }
    L <- L[is.finite(depth_min) & depth_min <= GATE_M]
    if (nrow(L) < 1000) { message(tg, ": too few observations after gating"); next }
  }
  message("fit: ", tg, " / ", a$arm)
  invisible(capture.output(
    f <- TwilightFreeGrid(L$time, pmax(0, L$light - r$baseline), grid = grid,
      start_lon = unname(g$p0[1]), start_lat = unname(g$p0[2]),
      end_lon = unname(g$p1[1]), end_lat = unname(g$p1[2]),
      step_hours = STEP_HOURS, diffusion = DIFFUSION,
      calibration = r$calibration,
      likelihood_params = c(a$lam / (r$max_light * 0.5), r$max_light, a$pslab))))
  gp <- grid_posterior(f)
  sdp <- numeric(nrow(gp$P))
  for (i in seq_len(nrow(gp$P))) {
    w <- gp$P[i, ]; s <- sum(w)
    if (!is.finite(s) || s <= 0) { sdp[i] <- NA; next }
    w <- w / s; sdp[i] <- sqrt(sum(w * (gp$lat - sum(w * gp$lat))^2))
  }
  tm <- as.POSIXct(f$fit$time, origin = "1970-01-01", tz = "UTC")
  tr <- argos_at(as.data.table(g$argos), tm)      # strict scorer, NA across >24 h gaps
  ep <- f$fit$lat - tr$lat
  e  <- gc_km(lon360(f$fit$lon), f$fit$lat, tr$lon, tr$lat)
  ok <- is.finite(ep) & is.finite(sdp) & sdp > 0
  fwrite(data.table(id = tg, arm = a$arm, season = g$season, family = g$family,
                    n = sum(ok), median_km = median(e[ok]),
                    bias_lat = mean(ep[ok]), rmse_lat = sqrt(mean(ep[ok]^2)),
                    lat_sd = mean(sdp[ok]),
                    cover_lat = mean(abs(ep[ok]) <= 1.96 * sdp[ok]),
                    log_z = f$log_z), OUT_TAG, append = file.exists(OUT_TAG))
  fwrite(data.table(id = tg, arm = a$arm, time = tm[ok], true_lat = tr$lat[ok],
                    err_lat = ep[ok], err_km = e[ok], lat_sd = sdp[ok],
                    decl = solar_declination(tm[ok])),
         OUT_KNOT, append = file.exists(OUT_KNOT))
}

D <- fread(OUT_TAG, colClasses = list(character = "id"))
D[, arm := factor(arm, levels = ARMS$arm)]
cat(sprintf("\n%d fits, %d tags\n\n=== all deployments ===\n", nrow(D), uniqueN(D$id)))
print(as.data.frame(D[, .(n = .N, km = round(median(median_km)),
                          bias = round(mean(bias_lat), 3),
                          abs_bias = round(mean(abs(bias_lat)), 3),
                          rmse_lat = round(mean(rmse_lat), 2),
                          lat_sd = round(mean(lat_sd), 2),
                          cover = round(mean(cover_lat), 3)),
                      by = arm][order(arm)]), row.names = FALSE)

ref <- fread("scratch/nes_calibration/rescore29_results.csv",
             colClasses = list(character = "id"))[arm == "tangent_areaON"]
m <- merge(D[arm == "A_shipped", .(id, mine = median_km)], ref[, .(id, ref = median_km)],
           by = "id")
if (nrow(m))
  cat(sprintf("\nCONTROL: A_shipped median %.0f km vs rescore29 %.0f km; max |ratio-1| %.3f\n",
              median(m$mine), median(m$ref), max(abs(m$mine / m$ref - 1))))

for (x in setdiff(ARMS$arm, "A_shipped")) {
  w <- merge(D[arm == "A_shipped", .(id, k0 = median_km, b0 = bias_lat, c0 = cover_lat)],
             D[arm == x, .(id, k1 = median_km, b1 = bias_lat, c1 = cover_lat)], by = "id")
  if (!nrow(w)) next
  cat(sprintf("\n  %-16s vs shipped (n=%2d): d km %+5.0f (better %2d/%2d) | d|bias| %+0.3f | d cover %+0.3f\n",
              x, nrow(w), mean(w$k1 - w$k0), sum(w$k1 < w$k0), nrow(w),
              mean(abs(w$b1)) - mean(abs(w$b0)), mean(w$c1 - w$c0)))
}
cat("\nPer-knot errors in final29_knots.csv for the seasonal and latitude structure.\n")
