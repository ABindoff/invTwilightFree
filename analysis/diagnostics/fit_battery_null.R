# THE NULL GATE. Fit the battery's PERFECT light with the response that generated it.
#
# If the estimator cannot recover a known track from noiseless, correctly-modelled
# light, then the latitude bias is in the estimator and everything downstream of it --
# the lambda sweeps, the shading story, the declination schedule -- is measuring
# something else. This is 24 fits and it can invalidate the design of the 648 that
# follow, so it runs first and alone.
#
# Also fits the NOISY level (clear sky plus the engine's own spike-and-slab, no
# attenuation) at lambda multiplier 1, which is a second null: there the model is
# exactly correct, so any bias is again the estimator's rather than a misspecification.
#
# Response is the LOGISTIC one used to generate the light. A failure here is therefore
# unambiguous -- it cannot be blamed on response mismatch.
#
# Scored ONLY on samples whose truth is supported (within 12 h of an Argos fix).
# Statistics cluster by TRACK, not by record: the four date offsets share an
# underlying track and are not independent.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))
suppressMessages(devtools::load_all(".", quiet = TRUE))
BAT <- readRDS("scratch/nes_calibration/light_battery.rds")
OUT <- "scratch/nes_calibration/battery_null.csv"
STEP_H <- 12; DIFF <- 110; CELL <- 1; PSLAB <- 0.10

fit_one <- function(b, which_light, lam_mult) {
  y <- b[[which_light]]
  # knots are 12-hourly; keep only observations we can score against
  gg <- rast(xmin = floor(min(b$lon)) - 12, xmax = ceiling(max(b$lon)) + 12,
             ymin = max(-88, floor(min(b$lat)) - 12),
             ymax = min(88, ceiling(max(b$lat)) + 12),
             resolution = CELL, crs = "EPSG:4326")
  values(gg) <- 1
  lam <- lam_mult / (b$max_light * 0.5)
  invisible(capture.output(
    f <- TwilightFreeGrid(b$time, y, grid = gg,
      start_lon = b$lon[1], start_lat = b$lat[1],
      end_lon = b$lon[length(b$lon)], end_lat = b$lat[length(b$lat)],
      step_hours = STEP_H, diffusion = DIFF,
      calibration = b$response,                 # the LOGISTIC response that generated it
      likelihood_params = c(lam, b$max_light, PSLAB))))
  # truth and support at the knot times
  tk <- as.numeric(f$fit$time); tt <- as.numeric(b$time)
  tlat <- approx(tt, b$lat, tk, rule = 2)$y
  tlon <- approx(tt, lon360(b$lon), tk, rule = 2)$y
  sup  <- approx(tt, as.numeric(b$supported), tk, rule = 2)$y > 0.999
  gp <- grid_posterior(f); sdp <- numeric(nrow(gp$P))
  for (k in seq_len(nrow(gp$P))) {
    w <- gp$P[k, ]; s <- sum(w)
    if (!is.finite(s) || s <= 0) { sdp[k] <- NA; next }
    w <- w / s; mu <- sum(w * gp$lat); sdp[k] <- sqrt(sum(w * (gp$lat - mu)^2))
  }
  e <- f$fit$lat - tlat
  ok <- sup & is.finite(e) & is.finite(sdp) & sdp > 0
  data.table(id = b$id, offset = b$offset, light = which_light, lam_mult = lam_mult,
             n = sum(ok), frac_sup = mean(sup),
             bias = mean(e[ok]), rmse = sqrt(mean(e[ok]^2)),
             cover = mean(abs(e[ok]) <= 1.96 * sdp[ok]),
             km = median(gc_km(lon360(f$fit$lon)[ok], f$fit$lat[ok], tlon[ok], tlat[ok])))
}

# CHECKPOINTED per fit: the first attempt at this lost five completed fits to an
# overnight restart because it only wrote at the end. Key is (id, offset, light);
# rows already on disk are skipped, so re-running resumes.
done <- if (file.exists(OUT)) {
  d0 <- fread(OUT, colClasses = list(character = "id"))
  paste(d0$id, d0$offset, d0$light)
} else character(0)
n_want <- length(BAT) * 2
cat(sprintf("%d of %d fits already done\n\n", length(done), n_want))

for (b in BAT) {
  for (cfg in list(c("perfect", 1), c("noisy", 1))) {
    key <- paste(b$id, b$offset, cfg[1])
    if (key %in% done) next
    r <- try(fit_one(b, cfg[1], as.numeric(cfg[2])), silent = TRUE)
    if (inherits(r, "try-error")) {
      # do not swallow: a systematic failure must not masquerade as a short table
      message("FAILED ", key, ": ", conditionMessage(attr(r, "condition")))
      next
    }
    fwrite(r, OUT, append = file.exists(OUT))
    message(key, " done")
  }
}
D <- fread(OUT, colClasses = list(character = "id"))

cat(sprintf("\n%d of %d fits\n", nrow(D), n_want))
if (nrow(D) < n_want)
  cat("*** INCOMPLETE -- the verdict below is provisional; re-run to finish.\n")
cat("\n")
for (lv in unique(D$light)) {
  x <- D[light == lv]
  # cluster by TRACK: the four offsets share a track and are not independent
  cl <- x[, .(bias = mean(bias), rmse = mean(rmse), cover = mean(cover),
              km = mean(km)), by = id]
  cat(sprintf("=== %s light (lambda x1, logistic response) ===\n", toupper(lv)))
  cat(sprintf("  by record  (n=%2d): bias %+0.3f | rmse %.3f | cover %.3f | km %.0f\n",
              nrow(x), mean(x$bias), mean(x$rmse), mean(x$cover), median(x$km)))
  cat(sprintf("  by track   (n=%2d): bias %+0.3f | rmse %.3f | cover %.3f | km %.0f\n",
              nrow(cl), mean(cl$bias), mean(cl$rmse), mean(cl$cover), median(cl$km)))
  p <- suppressWarnings(wilcox.test(cl$bias)$p.value)
  cat(sprintf("  bias differs from 0 (by track): p = %.4f\n", p))
  cat(sprintf("  per-track bias: %s\n\n",
              paste(sprintf("%s %+0.2f", cl$id, cl$bias), collapse = "  ")))
}
cat("VERDICT\n")
pf <- D[light == "perfect", .(b = mean(bias)), by = id]
cat(sprintf("  perfect light, mean |bias| by track = %.3f deg\n", mean(abs(pf$b))))
cat(if (mean(abs(pf$b)) < 0.25)
      "  GATE PASSED: the estimator recovers a known track from correct light.\n  The lambda sweep will therefore be measuring shading, not an estimator defect.\n"
    else
      "  *** GATE FAILED: the estimator is biased on perfect, correctly-modelled light.\n  The latitude bias is not (only) about shading, and the sweep as designed would\n  be measuring the wrong thing.\n")
cat("\nby date offset (does season matter even with perfect light?):\n")
print(as.data.frame(D[, .(bias = round(mean(bias), 3), rmse = round(mean(rmse), 3),
                          cover = round(mean(cover), 3)), by = .(light, offset)]),
      row.names = FALSE)
