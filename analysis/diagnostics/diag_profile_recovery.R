# PROFILE-RECOVERY TEST: can the marginal likelihood identify the calibration?
#
# THE DIAGNOSIS THIS TESTS. The emission is sloppy (condition number 3.8e4) and
# latitude is essentially the soft direction within a day (weight 0.881), trading
# against the threshold z50 at 1.4-4.4 deg of latitude per degree of z50, with the
# exchange rate CHANGING SIGN across the year. Measured gauge inflation falls from
# 1.81 on a 5-day window to 1.08 over a full 240-day deployment, and to 1.01 pooling
# six tags, with se(z50) shrinking as 1/sqrt(N) to 0.045 deg.
#
# So the seasonal contrast over a full deployment IDENTIFIES the calibration -- but
# `fit_light_response()` estimates z50 on a 15-day haul-out window spanning under a
# degree of declination, where inflation is 1.7 and z50 is worst identified, then
# freezes it for a deployment spanning 47 degrees. The information that breaks the
# gauge is in the data and is discarded by construction.
#
# THE PROPOSED FIX NEEDS NO ENGINE CHANGE. `log_z` is already the marginal likelihood
# of the whole track. Profile it over z50, summed across tags, and take the
# maximiser: empirical Bayes on the shared systematic, using the full-deployment
# contrast. Comparable across the sweep because grid, area_correction and the support
# [0, max_light] are all held fixed; only the model's calibration varies.
#
# THE TEST. The battery generates from a KNOWN response, so the true z50 is known and
# the answer is checkable.
#
# PRE-REGISTERED PREDICTIONS:
#   1. Each tag's log_z profile has an INTERIOR maximum. (Real risk of failure: the
#      clock estimator's profile ran to its search bound on all 10 tags and every
#      estimate had to be discarded. A profile that hits the edge identifies nothing.)
#   2. The per-tag maximisers scatter around dz50 = 0; the POOLED maximiser is closer
#      to 0 than the typical per-tag one.
#   3. Latitude bias varies strongly with dz50 -- of order 1.4-4.4 deg per degree --
#      confirming the exchange rate measured from the Fisher information.
#   4. THE SHARP ONE: the dz50 that MAXIMISES EVIDENCE coincides with the dz50 that
#      ZEROES LATITUDE BIAS. If they diverge, maximising the marginal likelihood is
#      not a route to removing the bias, and this whole approach fails however
#      well-identified the parameter is.
#
# Noisy arm, 6 tracks, offset 0, 7 values of dz50 = 42 fits. Checkpointed on (id, dz).
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))
suppressMessages(devtools::load_all(".", quiet = TRUE))
BAT <- readRDS("scratch/nes_calibration/light_battery.rds")
REC <- Filter(function(z) z$offset == 0, BAT)
OUT <- "scratch/nes_calibration/profile_recovery.csv"
STEP_H <- 12; CELL <- 1; PSLAB <- 0.10; DIFF <- 110
DZ <- c(-2, -1, -0.5, 0, 0.5, 1, 2)

fit_one <- function(b, dz) {
  r <- b$response; r[3] <- r[3] + dz          # perturb the threshold only
  gg <- rast(xmin = floor(min(b$lon)) - 12, xmax = ceiling(max(b$lon)) + 12,
             ymin = max(-88, floor(min(b$lat)) - 12),
             ymax = min(88, ceiling(max(b$lat)) + 12),
             resolution = CELL, crs = "EPSG:4326")
  values(gg) <- 1
  lam <- 1 / (b$max_light * 0.5)
  invisible(capture.output(
    f <- TwilightFreeGrid(b$time, b$noisy, grid = gg,
      start_lon = b$lon[1], start_lat = b$lat[1],
      end_lon = b$lon[length(b$lon)], end_lat = b$lat[length(b$lat)],
      step_hours = STEP_H, diffusion = DIFF, calibration = r,
      likelihood_params = c(lam, b$max_light, PSLAB),
      area_correction = TRUE)))

  gp <- grid_posterior(f); K <- nrow(gp$P)
  mu <- rep(NA_real_, K)
  for (k in seq_len(K)) {
    w <- gp$P[k, ]; s <- sum(w)
    if (is.finite(s) && s > 0) mu[k] <- sum(w / s * gp$lat)
  }
  tk <- as.numeric(f$fit$time); tt <- as.numeric(b$time)
  tlat <- approx(tt, b$lat, tk, rule = 2)$y
  tlon <- approx(tt, lon360(b$lon), tk, rule = 2)$y
  sup  <- approx(tt, as.numeric(b$supported), tk, rule = 2)$y > 0.999
  ok <- sup & is.finite(mu)
  data.table(id = b$id, dz = dz, n = sum(ok), log_z = f$log_z,
             bias = mean((mu - tlat)[ok]),
             bias_mode = mean((f$fit$lat - tlat)[ok]),
             km = median(gc_km(lon360(f$fit$lon)[ok], f$fit$lat[ok],
                               tlon[ok], tlat[ok])))
}

done <- if (file.exists(OUT)) {
  d0 <- fread(OUT, colClasses = list(character = "id"))
  paste(d0$id, d0$dz)
} else character(0)
cat(sprintf("%d of %d fits already done\n\n", length(done), length(DZ) * length(REC)))

for (dz in DZ) for (b in REC) {
  key <- paste(b$id, dz)
  if (key %in% done) next
  r <- try(fit_one(b, dz), silent = TRUE)
  if (inherits(r, "try-error")) {
    message("FAILED ", key, ": ", conditionMessage(attr(r, "condition"))); next
  }
  fwrite(r, OUT, append = file.exists(OUT))
  message(key, " done")
}

D <- fread(OUT, colClasses = list(character = "id"))
cat(sprintf("\n%d fits\n", nrow(D)))

cat("\n=== per tag: log_z profile (relative to its own maximum) ===\n")
D[, rel := log_z - max(log_z), by = id]
print(as.data.frame(dcast(D, id ~ dz, value.var = "rel")[, lapply(.SD, function(x)
  if (is.numeric(x)) round(x, 1) else x)]), row.names = FALSE)

pk <- D[, .(dz_hat = dz[which.max(log_z)],
            interior = dz[which.max(log_z)] > min(DZ) & dz[which.max(log_z)] < max(DZ)),
        by = id]
cat("\n  per-tag argmax dz50: ", paste(sprintf("%s %+0.1f", pk$id, pk$dz_hat),
                                       collapse = "  "), "\n")
cat(sprintf("  PREDICTION 1 (interior maxima): %d of %d interior -> %s\n",
            sum(pk$interior), nrow(pk),
            if (all(pk$interior)) "HOLDS" else "*** FAILS for some tags"))

P <- D[, .(log_z = sum(log_z), bias = mean(bias), bias_mode = mean(bias_mode),
           km = median(km)), by = dz][order(dz)]
P[, rel := log_z - max(log_z)]
cat("\n=== pooled across tags ===\n")
print(as.data.frame(P[, .(dz, rel_logz = round(rel, 1), bias = round(bias, 3),
                          bias_mode = round(bias_mode, 3), km = round(km))]),
      row.names = FALSE)

dz_ev <- P$dz[which.max(P$log_z)]
# where does bias cross zero? linear interpolation on the pooled profile
f_bias <- approxfun(P$dz, P$bias)
rt <- try(uniroot(f_bias, range(P$dz))$root, silent = TRUE)
dz_bias <- if (inherits(rt, "try-error")) NA_real_ else rt

cat(sprintf("\n  PREDICTION 2: pooled argmax dz50 = %+0.2f (truth 0; per-tag spread %.2f)\n",
            dz_ev, sd(pk$dz_hat)))
slope <- coef(lm(bias ~ dz, P))[2]
cat(sprintf("  PREDICTION 3: d(bias)/d(z50) = %.2f deg latitude per deg z50\n", slope))
cat(sprintf("               (Fisher exchange rate was 1.4 - 4.4) -> %s\n",
            if (abs(slope) >= 1.0 && abs(slope) <= 5.5) "CONSISTENT" else "INCONSISTENT"))
cat(sprintf("  PREDICTION 4: evidence peak %+0.2f vs bias-zero %+0.2f  -> %s\n",
            dz_ev, dz_bias,
            if (!is.na(dz_bias) && abs(dz_ev - dz_bias) <= 0.5)
              "COINCIDE: maximising evidence removes the bias" else
              "*** DIVERGE: evidence maximisation does NOT remove the bias"))
