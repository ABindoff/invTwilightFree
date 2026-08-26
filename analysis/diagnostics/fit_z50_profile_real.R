# PROFILE z50 ON REAL TAGS, by pooled marginal likelihood.
#
# WHY, and what this can and cannot do -----------------------------------------
# diag_profile_recovery.R already ran this on the synthetic battery where truth is
# known. Results, which set the expectations here:
#
#   * pooled argmax recovered the true z50 EXACTLY (+0.00), while per-tag argmaxes
#     scattered with SD 1.03 deg. So pool; do not trust per-tag maxima.
#   * one tag of six had its profile run to the search bound, identifying nothing.
#     Expect some real tags to do the same and exclude them explicitly.
#   * d(latitude bias)/d(z50) came out at 0.07 deg per deg, NOT the 1.4-4.4 that
#     the Fisher information suggested, and pooled bias never crossed zero over a
#     4-degree swing (-0.84 to -0.49). The within-day exchange rate changes sign
#     across the year, so it cancels in the deployment mean.
#
# THEREFORE: this run will NOT remove the latitude bias. Anyone reading the output
# expecting that will misread it. What it CAN do is catch a common z50 that is
# simply in the wrong place, which costs real accuracy: on the synthetic grid a
# 2-degree error cost 160 km (350 against 189).
#
# WHY THAT IS PLAUSIBLE HERE. These are double-tagging TDRs deployed to measure
# depth; the light channel was never intended for geolocation and there was no
# field calibration. `fit_light_response()` therefore estimates the response from a
# hauled-out animal at the colony, frequently under other seals and their filth,
# and that value is then applied to clean open-water surfacings for 240 days. The
# calibration and the movement are in different optical regimes, and nothing in the
# pipeline measures the difference. Profiling z50 from the AT-SEA evidence bypasses
# the haul-out fit entirely, so the gap between the two IS that mismatch.
#
# DESIGN. The pipeline already POOLS z50, so all 25 Mk9_219 tags share 91.92 and
# the 4 F18A share 91.74. The live question is whether the common value is right,
# not whether tags differ, so dz is applied relative to each tag's current value
# and the pooled maximiser is the estimate. Mk9_219 only: F18A pools from a single
# deployment and cannot support its own profile.
#
# COST. 12 tags x 7 dz = 84 fits at roughly 9 min each, so about 13 hours. The
# synthetic test reached the exact answer pooling 6 tags, so 12 is deliberate
# headroom rather than a minimum. Checkpointed on (id, dz).
SCRATCH <- "analysis/diagnostics"
source(file.path(SCRATCH, "nes_common.R"))
suppressMessages(devtools::load_all(".", quiet = TRUE))

P <- readRDS("analysis/cache/nes/panel_v1.rds")
OUT <- "scratch/nes_calibration/z50_profile_real.csv"
STEP_HOURS <- 12; DIFFUSION <- 110; CELL <- 1; PSLAB <- 0.10; LAM <- 1.00
# Grid extended UPWARD and the loop reordered, both for reasons found after v1
# started:
#
#   * dz = -2 costs 541 km against 243 at dz = 0, so the optimum is at or above
#     the current value and a symmetric window wastes half its fits below it.
#   * an independent edge-based estimator fits the steepest point of the response
#     at 97.8 deg for this sensor family. That is NOT the same parameter as z50
#     here -- the engine's panel calibration is a two-parameter CLAMPED LINEAR
#     ramp (78.0 to 105.8 deg, midpoint 91.92), which has no distinguished
#     steepest point -- so it is not evidence of a 6 deg error. But it is a reason
#     to look above the current value rather than symmetrically around it.
DZ <- c(-2, -1, 0, 1, 2, 3, 4)
N_TAGS <- 12

fam <- vapply(P$tags, function(g) g$family, "")
cand <- names(P$tags)[fam == "Mk9_219"]
# spread the selection over latitude and season rather than taking the first 12,
# because the seasonal contrast is what identifies z50 in the first place
info <- rbindlist(lapply(cand, function(tg) {
  A <- as.data.table(P$tags[[tg]]$argos)
  data.table(id = tg, season = P$tags[[tg]]$season,
             lat = median(A$lat, na.rm = TRUE), n = nrow(P$tags[[tg]]$light))
}))
setorder(info, season, lat)
sel <- info$id[round(seq(1, nrow(info), length.out = N_TAGS))]
cat(sprintf("Mk9_219 candidates %d, selected %d:\n  %s\n\n",
            nrow(info), length(sel), paste(sel, collapse = ", ")))

grid <- rast(xmin = 150, xmax = 250, ymin = 20, ymax = 70,
             resolution = CELL, crs = "EPSG:4326")
values(grid) <- 1

done <- if (file.exists(OUT)) {
  d0 <- fread(OUT, colClasses = list(character = "id")); paste(d0$id, d0$dz)
} else character(0)
cat(sprintf("%d of %d fits done\n\n", length(done), length(DZ) * length(sel)))

# TAGS OUTER, dz inner. The reverse order finishes every tag at one dz before
# starting the next, so nothing is learnable until the whole 13-hour run is done.
# This way a complete profile exists after 7 fits, about an hour, and the grid can
# be corrected early if it turns out to be in the wrong place.
for (tg in sel) for (dz in DZ) {
  if (paste(tg, dz) %in% done) next
  g <- P$tags[[tg]]; r <- P$responses[[tg]]
  if (is.null(r)) next
  L <- as.data.table(g$light)
  # calibration is c(slope * (z50 + 2*scale), slope); shifting z50 by dz shifts
  # only the intercept, by slope*dz. Slope and scale are held fixed so the
  # comparison is over the threshold alone.
  cal <- c(r$calibration[1] + r$slope * dz, r$calibration[2])
  message("fit: ", tg, " dz ", dz)
  ok_fit <- try({
    invisible(capture.output(
      f <- TwilightFreeGrid(L$time, pmax(0, L$light - r$baseline), grid = grid,
        start_lon = unname(g$p0[1]), start_lat = unname(g$p0[2]),
        end_lon = unname(g$p1[1]), end_lat = unname(g$p1[2]),
        step_hours = STEP_HOURS, diffusion = DIFFUSION, calibration = cal,
        likelihood_params = c(LAM / (r$max_light * 0.5), r$max_light, PSLAB))))
    TRUE }, silent = TRUE)
  if (inherits(ok_fit, "try-error")) {
    message("  FAILED: ", conditionMessage(attr(ok_fit, "condition"))); next
  }
  gp <- grid_posterior(f)
  sdp <- numeric(nrow(gp$P)); mu <- numeric(nrow(gp$P))
  for (i in seq_len(nrow(gp$P))) {
    w <- gp$P[i, ]; s <- sum(w)
    if (!is.finite(s) || s <= 0) { sdp[i] <- NA; mu[i] <- NA; next }
    w <- w / s; mu[i] <- sum(w * gp$lat)
    sdp[i] <- sqrt(sum(w * (gp$lat - mu[i])^2))
  }
  tm <- as.POSIXct(f$fit$time, origin = "1970-01-01", tz = "UTC")
  tr <- argos_at(as.data.table(g$argos), tm)
  ep <- f$fit$lat - tr$lat
  e  <- gc_km(lon360(f$fit$lon), f$fit$lat, tr$lon, tr$lat)
  keep <- is.finite(ep) & is.finite(sdp) & sdp > 0
  fwrite(data.table(id = tg, dz = dz, z50_base = r$z50, z50 = r$z50 + dz,
                    season = g$season, n = sum(keep), log_z = f$log_z,
                    median_km = median(e[keep]),
                    bias_lat = mean(ep[keep]),
                    bias_mean = mean((mu - tr$lat)[keep], na.rm = TRUE),
                    rmse_lat = sqrt(mean(ep[keep]^2)),
                    lat_sd = mean(sdp[keep]),
                    cover_lat = mean(abs(ep[keep]) <= 1.96 * sdp[keep])),
         OUT, append = file.exists(OUT))
}
cat("\ndone\n")
