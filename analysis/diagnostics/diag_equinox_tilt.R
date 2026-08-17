# Is the latitude bias SEASONAL, and is it the area tilt?
#
# Established so far: the bias is not accumulating (it is near zero for 40% of the
# track, dives to -1.9 deg mid-track, partly recovers), it is not latitude
# (quintiles 2-4 sit at 47.2/47.4/47.6 N with biases -0.03/-1.65/-1.17), and it is
# not the movement prior (the exact Ito correction moves the FILTERED estimate by
# 0.023 deg, where nothing cancels it). Position along the track predicts it
# 7x better than tan(latitude): R2 0.312 vs 0.045.
#
# At offset 0 these records run roughly June to January, so mid-track is the
# September equinox -- where E2, the even-harmonic content of the light, nulls and
# latitude information vanishes. Hypothesis: wherever the light stops constraining
# latitude, something pulls the estimate equatorward.
#
# TWO PREDICTIONS, WRITTEN DOWN BEFORE RUNNING.
#
# (1) SEASON, not track position. The battery holds the SAME movement at four date
#     offsets, which is exactly the instrument for this. A ~240-day track shifted by
#     0/90/180/270 days puts the equinoxes at different fractions along it:
#
#       offset   0  (Jun-Jan): ONE equinox near frac 0.47   -> single mid-track peak
#       offset  90 (Aug-Apr):  TWO, near frac 0.10 and 0.84 -> peaks near BOTH ends
#       offset 180 (Nov-Jul):  ONE near frac 0.47           -> single mid-track peak
#       offset 270 (Feb-Oct):  TWO, near frac 0.09 and 0.87 -> peaks near BOTH ends
#
#     If the bias is seasonal, offsets 90 and 270 must show a DOUBLE-humped profile
#     with the middle relatively clean -- the opposite of offset 0. If instead the
#     bias stays mid-track at every offset, season is wrong and something about
#     track position (endpoint pinning, say) is responsible.
#
# (2) THE AREA TILT, quantitatively. A cos(latitude) tilt on a knot's marginal
#     shifts its mean by about -sigma_post^2 * tan(lat) (sigma in radians). Tested
#     against the per-knot posterior sd, which is saved here for the first time.
#
#     HONEST ARITHMETIC UP FRONT: at the observed sigma_post of ~1.1 deg this
#     predicts -0.021 deg, roughly EIGHTY times too small to explain -1.6. So the
#     simple tilt is very likely already dead, and reproducing the -1.6 would need
#     sigma_post near 9-10 deg. What is still worth measuring is the SHAPE: if the
#     bias tracks sigma_post^2 * tan(lat) up to a large constant factor, the tilt is
#     acting through the chain rather than knot-by-knot. If it does not track it at
#     all, the tilt is out and the equatorward pull is something else.
#
# 24 fits, checkpointed per (id, offset). Perfect light, 1 degree, area ON,
# drift OFF -- the same configuration as everything else in this series.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))
suppressMessages(devtools::load_all(".", quiet = TRUE))
BAT <- readRDS("scratch/nes_calibration/light_battery.rds")
OUT <- "scratch/nes_calibration/equinox_tilt.csv"
STEP_H <- 12; CELL <- 1; PSLAB <- 0.10; DIFF <- 110

fit_one <- function(b) {
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
      area_correction = TRUE)))

  gp <- grid_posterior(f)
  K <- nrow(gp$P)
  mu <- sd <- rep(NA_real_, K)
  for (k in seq_len(K)) {
    w <- gp$P[k, ]; s <- sum(w)
    if (!is.finite(s) || s <= 0) next
    w <- w / s
    mu[k] <- sum(w * gp$lat)
    sd[k] <- sqrt(sum(w * (gp$lat - mu[k])^2))
  }

  tk <- as.numeric(f$fit$time); tt <- as.numeric(b$time)
  tlat <- approx(tt, b$lat, tk, rule = 2)$y
  sup  <- approx(tt, as.numeric(b$supported), tk, rule = 2)$y > 0.999
  ktime <- as.POSIXct(tk, origin = "1970-01-01", tz = "UTC")

  data.table(id = b$id, offset = b$offset, k = seq_len(K), K = K,
             frac = seq_len(K) / K, sup = sup,
             time = ktime, decl = solar_declination(ktime),
             truth = tlat, mu = mu, sd = sd, mode = f$fit$lat)
}

done <- if (file.exists(OUT)) {
  d0 <- fread(OUT, colClasses = list(character = "id"))
  unique(paste(d0$id, d0$offset))
} else character(0)
cat(sprintf("%d of %d fits already done\n\n", length(done), length(BAT)))

for (b in BAT) {
  key <- paste(b$id, b$offset)
  if (key %in% done) next
  r <- try(fit_one(b), silent = TRUE)
  if (inherits(r, "try-error")) {
    message("FAILED ", key, ": ", conditionMessage(attr(r, "condition"))); next
  }
  fwrite(r, OUT, append = file.exists(OUT))
  message(key, " done")
}

D <- fread(OUT, colClasses = list(character = "id"))
cat(sprintf("\n%d fits, %d knots\n", uniqueN(paste(D$id, D$offset)), nrow(D)))
D <- D[sup == TRUE & is.finite(mu) & is.finite(sd) & sd > 0]
D[, e := mu - truth]
D[, adecl := abs(decl)]

cat("\n=== PREDICTION 1: bias by position along track, SPLIT BY OFFSET ===\n")
D[, q := cut(frac, seq(0, 1, 0.2), include.lowest = TRUE, labels = FALSE)]
tab <- dcast(D[, .(b = round(mean(e), 2)), by = .(offset, q)], offset ~ q,
             value.var = "b")
setnames(tab, c("offset", "0-20%", "20-40%", "40-60%", "60-80%", "80-100%"))
print(as.data.frame(tab), row.names = FALSE)
cat("  seasonal => offsets 90 and 270 humped at the ENDS, 0 and 180 in the MIDDLE\n")

cat("\n=== the same thing keyed on |solar declination| instead of position ===\n")
D[, dband := cut(adecl, c(-Inf, 5, 10, 15, 20, Inf),
                 labels = c("<5 (equinox)", "5-10", "10-15", "15-20", ">20 (solstice)"))]
print(as.data.frame(D[, .(n = .N, bias = round(mean(e), 3),
                          sd_post = round(mean(sd), 2)), by = dband][order(dband)]),
      row.names = FALSE)

cat("\n  per offset, so this is not one season doing all the work:\n")
print(as.data.frame(dcast(D[, .(b = round(mean(e), 2)), by = .(offset, dband)],
                          offset ~ dband, value.var = "b")), row.names = FALSE)

cat("\n=== PREDICTION 2: does bias track -sigma_post^2 * tan(lat)? ===\n")
D[, tilt := -(sd * pi / 180)^2 * tan(truth * pi / 180) * 180 / pi]
cat(sprintf("  mean predicted tilt %+0.4f deg vs mean observed bias %+0.4f deg\n",
            mean(D$tilt), mean(D$e)))
cat(sprintf("  ratio observed/predicted: %.1fx\n", mean(D$e) / mean(D$tilt)))
m <- lm(e ~ tilt, data = D)
cat(sprintf("  e ~ tilt : R2 %.4f, slope %+0.2f (1.0 would be exact)\n",
            summary(m)$r.squared, coef(m)[2]))
cat(sprintf("  Spearman(e, tilt) = %+0.3f\n",
            suppressWarnings(cor(D$e, D$tilt, method = "spearman"))))
cat(sprintf("  Spearman(e, sd)   = %+0.3f   [is a wide posterior a biased one?]\n",
            suppressWarnings(cor(D$e, D$sd, method = "spearman"))))
cat(sprintf("  Spearman(sd, |declination|) = %+0.3f  [does the light go blind at the equinox?]\n",
            suppressWarnings(cor(D$sd, D$adecl, method = "spearman"))))

cat("\n=== model comparison on the same rows ===\n")
for (f in c("e ~ tilt", "e ~ sd", "e ~ adecl", "e ~ frac", "e ~ sd + adecl")) {
  r2 <- summary(lm(as.formula(f), data = D))$r.squared
  cat(sprintf("  %-16s R2 %.4f\n", f, r2))
}
