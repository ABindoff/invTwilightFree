# TWO MEASUREMENTS ON REAL DATA. No fitting; facts before experiments.
#
# (A) MOVEMENT ANISOTROPY. The movement prior is isotropic, but a central-place
#     forager runs out and back, mostly east-west. An isotropic prior cannot bias
#     latitude by itself -- it is symmetric -- but it interacts with the cell-area
#     tilt, which is NOT symmetric and which we have just measured at +5.45 deg when
#     removed. Bias from a tilt scales with posterior variance, so a latitude prior
#     wider than the animal's true north-south movement amplifies it.
#
# (B) THE SPECTRAL FILTER. The tags carry an underdescribed blue-weighted filter, so
#     the light-vs-zenith curve in real data need not be the shape a geometric
#     calculation implies. Rayleigh scattering is wavelength dependent, so a
#     blue-weighted sensor's response is not a rescaled broadband one, and the
#     departure grows with AIR MASS -- largest at high zenith, exactly where latitude
#     is read. The synthetic battery cannot see this: it GENERATES from the fitted
#     logistic and is self-consistent by construction.
#
#     A constant offset is harmless -- the fitted z50 absorbs it. What is NOT
#     harmless is STRUCTURE: a dawn/dusk difference at equal zenith (impossible from
#     solar geometry alone), or a curve that MOVES WITH SEASON at fixed zenith (also
#     impossible for a pure function of zenith). Either distorts apparent day length,
#     and apparent day length is latitude.
#
# Truth comes from the battery's offset-0 records, whose positions are Argos
# interpolated and already clipped per deployment -- avoiding the PTT-reuse trap in
# the raw Argos files. Real light comes from the decimated archives at the same
# timestamps.
#
# Aggregates only. No per-observation output.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))
suppressMessages(devtools::load_all(".", quiet = TRUE))
BAT  <- readRDS("scratch/nes_calibration/light_battery.rds")
REC  <- Filter(function(z) z$offset == 0, BAT)
arch <- lapply(readRDS("analysis/cache/nes/archives_v3_30min.rds"), `[[`, "main")

# ---------------------------------------------------------------- (A) anisotropy
cat("=== (A) TRUE MOVEMENT ANISOTROPY (Argos truth, 12 h steps) ===\n")
A <- rbindlist(lapply(REC, function(b) {
  tt <- as.numeric(b$time)
  g  <- seq(min(tt), max(tt), by = 12 * 3600)
  gl <- approx(tt, b$lat, g)$y
  gn <- approx(tt, lon360(b$lon), g)$y
  dlat <- diff(gl) * 111.32
  dlon <- diff(gn); dlon[dlon >  180] <- dlon[dlon >  180] - 360
  dlon[dlon < -180] <- dlon[dlon < -180] + 360
  dlon <- dlon * 111.32 * cos(head(gl, -1) * pi / 180)
  data.table(id = b$id, n = length(dlat), sd_ns = sd(dlat), sd_ew = sd(dlon),
             ratio = sd(dlon) / sd(dlat))
}))
print(as.data.frame(A[, .(id, n, sd_ns = round(sd_ns, 1), sd_ew = round(sd_ew, 1),
                          ratio = round(ratio, 2))]), row.names = FALSE)
prior_step <- 110 * sqrt(0.5)
cat(sprintf("\n  median NS %.1f km, EW %.1f km, EW/NS ratio %.2f\n",
            median(A$sd_ns), median(A$sd_ew), median(A$ratio)))
cat(sprintf("  the prior assumes ratio 1.00 and %.1f km per 12 h ON BOTH AXES\n",
            prior_step))
cat(sprintf("  => latitude prior is %.2fx the true NS step; longitude %.2fx the EW step\n",
            prior_step / median(A$sd_ns), prior_step / median(A$sd_ew)))

# ------------------------------------------------------------------- (B) filter
cat("\n=== (B) REAL LIGHT vs GEOMETRIC ZENITH, AT KNOWN POSITIONS ===\n")
res <- list()
for (b in REC) {
  g <- arch[[b$id]]
  if (is.null(g) || !nrow(g)) next
  gt <- as.numeric(g$time); bt <- as.numeric(b$time)
  # real light at the battery's timestamps (which are the real dates at offset 0)
  y <- approx(gt, as.numeric(g$light), bt, rule = 2)$y
  ok <- b$supported & is.finite(y)
  if (sum(ok) < 2000) next
  tt <- bt[ok]; la <- b$lat[ok]; lo <- b$lon[ok]
  z  <- solar_zenith(tt, lo, la)
  z2 <- solar_zenith(tt + 900, lo, la)
  res[[b$id]] <- data.table(id = b$id, z = z, y = y[ok],
                            limb = ifelse(z2 < z, "dawn", "dusk"),
                            decl = solar_declination(as.POSIXct(tt, origin = "1970-01-01",
                                                                tz = "UTC")))
}
D <- rbindlist(res)
cat(sprintf("  %d tags, %d supported observations\n", uniqueN(D$id), nrow(D)))
D[, rel := (y - quantile(y, 0.05)) / (quantile(y, 0.95) - quantile(y, 0.05)), by = id]
D[, rel := pmin(pmax(rel, 0), 1)]
D[, zb := cut(z, seq(70, 120, 2.5))]
P <- D[!is.na(zb), .(n = .N, z = round(mean(z), 1), rel = round(median(rel), 3)),
       by = zb][order(zb)][n > 500]

cat("\n  median relative light by zenith (pooled):\n")
print(as.data.frame(P), row.names = FALSE)

cat("\n  DAWN vs DUSK at equal zenith -- solar geometry says these must be EQUAL:\n")
L <- merge(dcast(D[!is.na(zb), .(m = median(rel)), by = .(zb, limb)], zb ~ limb,
                 value.var = "m"), P[, .(zb, z)], by = "zb")
L[, d := dusk - dawn]
print(as.data.frame(L[, .(z, dawn = round(dawn, 3), dusk = round(dusk, 3),
                          diff = round(d, 3))]), row.names = FALSE)
tw <- L[z >= 85 & z <= 100]
cat(sprintf("  mean |dawn - dusk| across the twilight band: %.3f of range\n",
            mean(abs(tw$d), na.rm = TRUE)))

cat("\n  SEASONAL at fixed zenith -- a pure function of zenith CANNOT move here:\n")
D[, dband := cut(decl, c(-24, -10, 10, 24), labels = c("winter", "equinox", "summer"))]
S <- merge(dcast(D[!is.na(zb) & !is.na(dband), .(m = median(rel)), by = .(zb, dband)],
                 zb ~ dband, value.var = "m"), P[, .(zb, z)], by = "zb")
S[, swing := pmax(winter, equinox, summer, na.rm = TRUE) -
             pmin(winter, equinox, summer, na.rm = TRUE)]
print(as.data.frame(S[, .(z, winter = round(winter, 3), equinox = round(equinox, 3),
                          summer = round(summer, 3), swing = round(swing, 3))]),
      row.names = FALSE)
cat(sprintf("\n  mean seasonal swing in the twilight band (85-100): %.3f of range\n",
            mean(S[z >= 85 & z <= 100]$swing, na.rm = TRUE)))
cat("  The twilight band spans about 1.0 of relative range in total, so a swing of\n")
cat("  0.10 is a tenth of the entire signal that latitude is read from.\n")
