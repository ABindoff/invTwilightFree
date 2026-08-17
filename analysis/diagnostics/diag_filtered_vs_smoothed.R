# Does the latitude bias appear in the FORWARD FILTER?
#
# The Ito correction was exact and inert (notes/latitude_bias_investigation.md 4a).
# The explanation offered was that a forward-backward smoother cancels a drift
# entering from both time directions. That explanation is itself a hypothesis, and
# this is the experiment that tests it -- because if it is right, the drift must be
# VISIBLE in the one-sided forward filter, where nothing cancels it.
#
# THE DISCRIMINATOR IS THE SHAPE IN KNOT INDEX, not the average:
#
#   DRIFT  filtered bias GROWS with k (each step adds -tan(phi)sigma^2/2R^2 until
#          the likelihood balances it), and the smoothed bias is much smaller.
#          Turning drift_correction ON should then flatten the FILTERED curve even
#          though it left the smoothed one unmoved -- that is the decisive
#          signature, because it would show the correction works and the smoother
#          hides it.
#
#   TILT   filtered bias is FLAT in k and roughly equal to the smoothed bias. Then
#          nothing accumulates, the cancellation story is wrong or irrelevant, and
#          the mechanism is a per-knot reweighting of the marginal.
#
#   NEITHER  if the filtered bias is flat AND drift_correction still does nothing
#          to it, the movement prior is exonerated altogether and the residual is
#          in the emission or the knot structure.
#
# Note the two marginals agree at the FINAL knot by construction (beta is uniform
# there), which doubles as an internal check on the new `filtered` output.
#
# Six tracks, offset 0, perfect light, 1 degree, area ON -- the same configuration
# as every other experiment in this series so the numbers are comparable.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))
suppressMessages(devtools::load_all(".", quiet = TRUE))
BAT <- readRDS("scratch/nes_calibration/light_battery.rds")
OUT <- "scratch/nes_calibration/filtered_vs_smoothed.csv"
STEP_H <- 12; CELL <- 1; PSLAB <- 0.10; DIFF <- 110

REC <- Filter(function(b) b$offset == 0, BAT)

mean_lat <- function(P, lat) {
  apply(P, 1, function(w) { s <- sum(w); if (!is.finite(s) || s <= 0) NA else sum(w / s * lat) })
}

fit_one <- function(b, drift) {
  gg <- rast(xmin = floor(min(b$lon)) - 12, xmax = ceiling(max(b$lon)) + 12,
             ymin = max(-88, floor(min(b$lat)) - 12),
             ymax = min(88, ceiling(max(b$lat)) + 12),
             resolution = CELL, crs = "EPSG:4326")
  values(gg) <- 1
  lam <- 1 / (b$max_light * 0.5)
  invisible(capture.output(
    f <- TwilightFreeGrid(b$time, b$perfect, grid = gg,
      start_lon = b$lon[1], start_lat = b$lat[1],
      end_lon = b$lon[length(b$lon)], end_lat = b$lat[length(b$lat)],
      step_hours = STEP_H, diffusion = DIFF, calibration = b$response,
      likelihood_params = c(lam, b$max_light, PSLAB),
      area_correction = TRUE, drift_correction = drift)))

  sm <- grid_posterior(f, filtered = FALSE)
  fl <- grid_posterior(f, filtered = TRUE)
  mu_s <- mean_lat(sm$P, sm$lat)
  mu_f <- mean_lat(fl$P, fl$lat)

  tk <- as.numeric(f$fit$time); tt <- as.numeric(b$time)
  tlat <- approx(tt, b$lat, tk, rule = 2)$y
  sup  <- approx(tt, as.numeric(b$supported), tk, rule = 2)$y > 0.999
  K <- length(tk)

  data.table(id = b$id, drift = drift, k = seq_len(K), K = K,
             frac = seq_len(K) / K, sup = sup,
             truth = tlat, smoothed = mu_s, filt = mu_f)
}

R <- list()
for (dr in c(FALSE, TRUE)) for (b in REC) {
  r <- try(fit_one(b, dr), silent = TRUE)
  if (inherits(r, "try-error")) { message("FAILED ", b$id, " drift=", dr,
                                          ": ", conditionMessage(attr(r, "condition"))); next }
  R[[length(R) + 1]] <- r
  message(b$id, " drift=", dr, " done")
}
D <- rbindlist(R); fwrite(D, OUT)

D[, e_s := smoothed - truth]
D[, e_f := filt - truth]
S <- D[sup == TRUE & is.finite(e_s) & is.finite(e_f)]

cat("\n=== internal check: filtered == smoothed at the FINAL knot ===\n")
last <- D[k == K & is.finite(e_s) & is.finite(e_f)]
cat(sprintf("  max |filtered - smoothed| at k = K: %.6f deg (should be ~0)\n",
            max(abs(last$filt - last$smoothed))))

cat("\n=== overall latitude bias, by track then averaged ===\n")
for (dr in c(FALSE, TRUE)) {
  x <- S[drift == dr, .(s = mean(e_s), f = mean(e_f)), by = id]
  cat(sprintf("  drift %-5s : smoothed %+0.3f | FILTERED %+0.3f\n",
              dr, mean(x$s), mean(x$f)))
}

cat("\n=== THE DISCRIMINATOR: filtered bias by position along the track ===\n")
cat(sprintf("%-7s %-9s %8s %8s %8s %8s %8s\n",
            "drift", "series", "0-20%", "20-40%", "40-60%", "60-80%", "80-100%"))
for (dr in c(FALSE, TRUE)) {
  x <- S[drift == dr]
  x[, bin := cut(frac, seq(0, 1, 0.2), include.lowest = TRUE, labels = FALSE)]
  for (v in c("e_f", "e_s")) {
    m <- x[, .(b = mean(get(v))), by = bin][order(bin)]
    cat(sprintf("%-7s %-9s %8.3f %8.3f %8.3f %8.3f %8.3f\n",
                dr, if (v == "e_f") "filtered" else "smoothed",
                m$b[1], m$b[2], m$b[3], m$b[4], m$b[5]))
  }
}

cat("\n=== does the filtered bias TREND with knot index? ===\n")
for (dr in c(FALSE, TRUE)) {
  x <- S[drift == dr]
  fitl <- lm(e_f ~ frac, data = x)
  cat(sprintf("  drift %-5s : slope %+0.3f deg across the track (p = %.2g)\n",
              dr, coef(fitl)[2], summary(fitl)$coefficients[2, 4]))
}

cat("\nREADING\n")
cat("  Filtered bias growing with k and flattened by drift_correction => the\n")
cat("  drift is real and the SMOOTHER was hiding it.\n")
cat("  Filtered bias flat and ~= smoothed => a per-knot TILT, and the movement\n")
cat("  prior is not the mechanism at all.\n")
