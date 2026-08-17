# RE-RUN THE KEY DIAGNOSTICS ON THE HONEST ARM.
#
# WHY THIS EXISTS. Everything diagnosed today was measured on the battery's
# `perfect` (noiseless) arm. That arm is NOT a valid null: noiseless light is the
# EXPECTATION of the emission, not a draw from it, and evaluating a likelihood at
# E[y] is a different object from evaluating it at a sample. Demonstrated, not
# argued: the per-knot latitude tilt is -0.113 nats/deg on the perfect arm and
# -0.0024 on the noisy arm, and every "fix" that helped the perfect arm made the
# noisy arm WORSE (censoring 20x the seasonal swing, censoring + symmetric arms
# 17x). So the emission tilt is an artefact of the test, and the conclusions that
# rested on the perfect arm have to be re-established or dropped.
#
# The `noisy` arm is clear-sky light plus a draw from the engine's OWN spike-and-slab
# mixture, verified by KS test against `LightMix` at the fitted parameters. So the
# EMISSION is correctly specified there. The TRACK is still a real Argos track rather
# than a draw from the movement prior, which is exactly the situation on real data,
# and is what makes this the honest proxy rather than a circular self-consistency
# check.
#
# TWO QUESTIONS, both re-asked on the honest arm:
#
#   1. IS THE BIAS SEASONAL? The signature that drove the whole diagnosis
#      (-2.156 deg at declination -7, +0.101 POLEWARD at +22) was measured on the
#      perfect arm. If it is absent here, it was an artefact and the manuscript
#      should say nothing about seasonal latitude bias.
#
#   2. DOES THE AREA CORRECTION STILL MATTER? It is the one surviving candidate,
#      and unlike the emission tilt it is DATA-INDEPENDENT -- `log_cell_area(lat_i)`
#      is added to alpha regardless of what was observed -- so it cannot be an
#      artefact of off-model data. Its fit-level effect (+0.254 deg, 6/6 tracks) was
#      measured on the perfect arm and needs re-measuring here.
#
# REPRODUCTION CONTROL: area ON, offset 0, MODE bias must return -0.446, the null
# gate's noisy/offset-0 figure over these same six tracks.
#
# 6 tracks x 4 offsets x {area ON, OFF} = 48 fits, checkpointed on (id, offset, area).
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))
suppressMessages(devtools::load_all(".", quiet = TRUE))
BAT <- readRDS("scratch/nes_calibration/light_battery.rds")
OUT <- "scratch/nes_calibration/noisy_arm_fits.csv"
STEP_H <- 12; CELL <- 1; PSLAB <- 0.10; DIFF <- 110

fit_one <- function(b, area) {
  gg <- rast(xmin = floor(min(b$lon)) - 12, xmax = ceiling(max(b$lon)) + 12,
             ymin = max(-88, floor(min(b$lat)) - 12),
             ymax = min(88, ceiling(max(b$lat)) + 12),
             resolution = CELL, crs = "EPSG:4326")
  values(gg) <- 1
  lam <- 1 / (b$max_light * 0.5)
  invisible(capture.output(
    f <- TwilightFreeGrid(b$time, b$noisy, grid = gg,          # <- the HONEST arm
      start_lon = b$lon[1], start_lat = b$lat[1],
      end_lon = b$lon[length(b$lon)], end_lat = b$lat[length(b$lat)],
      step_hours = STEP_H, diffusion = DIFF, calibration = b$response,
      likelihood_params = c(lam, b$max_light, PSLAB),
      area_correction = area)))

  gp <- grid_posterior(f); K <- nrow(gp$P)
  mu <- sd <- rep(NA_real_, K)
  for (k in seq_len(K)) {
    w <- gp$P[k, ]; s <- sum(w)
    if (!is.finite(s) || s <= 0) next
    w <- w / s
    mu[k] <- sum(w * gp$lat); sd[k] <- sqrt(sum(w * (gp$lat - mu[k])^2))
  }
  tk <- as.numeric(f$fit$time); tt <- as.numeric(b$time)
  tlat <- approx(tt, b$lat, tk, rule = 2)$y
  sup  <- approx(tt, as.numeric(b$supported), tk, rule = 2)$y > 0.999
  ktime <- as.POSIXct(tk, origin = "1970-01-01", tz = "UTC")

  data.table(id = b$id, offset = b$offset, area = area,
             k = seq_len(K), frac = seq_len(K) / K, sup = sup,
             decl = solar_declination(ktime), truth = tlat,
             mu = mu, mode = f$fit$lat, sd = sd)
}

done <- if (file.exists(OUT)) {
  d0 <- fread(OUT, colClasses = list(character = "id"))
  unique(paste(d0$id, d0$offset, d0$area))
} else character(0)
cat(sprintf("%d of %d fits already done\n\n", length(done), 2 * length(BAT)))

for (ar in c(TRUE, FALSE)) for (b in BAT) {
  key <- paste(b$id, b$offset, ar)
  if (key %in% done) next
  r <- try(fit_one(b, ar), silent = TRUE)
  if (inherits(r, "try-error")) {
    message("FAILED ", key, ": ", conditionMessage(attr(r, "condition"))); next
  }
  fwrite(r, OUT, append = file.exists(OUT))
  message(key, " done")
}

D <- fread(OUT, colClasses = list(character = "id"))
cat(sprintf("\n%d fits, %d knots\n", uniqueN(paste(D$id, D$offset, D$area)), nrow(D)))
S <- D[sup == TRUE & is.finite(mu) & is.finite(sd) & sd > 0]
S[, `:=`(e = mu - truth, e_mode = mode - truth)]

ctl <- S[area == TRUE & offset == 0, .(b = mean(e_mode)), by = id]
cat(sprintf("\nREPRODUCTION CONTROL (area ON, offset 0, MODE): %+0.3f (gate: -0.446) %s\n",
            mean(ctl$b), if (abs(mean(ctl$b) + 0.446) < 0.08) "OK" else "*** DIFFERS"))

cat("\n=== Q1: IS THE BIAS SEASONAL ON THE HONEST ARM? ===\n")
S[, sband := cut(decl, c(-24, -15, -7, 7, 15, 24), include.lowest = TRUE,
                 labels = c("NH winter", "-15..-7", "equinox", "+7..+15", "NH summer"))]
print(as.data.frame(S[area == TRUE, .(n = .N, bias = round(mean(e), 3),
                                      post_sd = round(mean(sd), 2)),
                      by = sband][order(sband)]), row.names = FALSE)
sw <- S[area == TRUE, .(m = mean(e)), by = sband]
cat(sprintf("  seasonal swing: %.3f deg   (perfect arm was 2.26: -2.156 to +0.101)\n",
            diff(range(sw$m))))
cat(sprintf("  sign flip across season? %s\n",
            if (any(sw$m > 0) && any(sw$m < 0)) "YES" else "NO"))

cat("\n  same table, per offset (so one season is not doing all the work):\n")
print(as.data.frame(dcast(S[area == TRUE, .(b = round(mean(e), 2)), by = .(offset, sband)],
                          offset ~ sband, value.var = "b")), row.names = FALSE)

cat("\n=== Q2: DOES THE AREA CORRECTION STILL MATTER? ===\n")
a <- S[, .(bias_mean = round(mean(e), 3), bias_mode = round(mean(e_mode), 3),
           post_sd = round(mean(sd), 3)), by = area]
print(as.data.frame(a), row.names = FALSE)
w <- dcast(S[, .(b = mean(e)), by = .(id, offset, area)], id + offset ~ area,
           value.var = "b")
w <- w[complete.cases(w)]
d <- w[["FALSE"]] - w[["TRUE"]]
cat(sprintf("\n  area OFF minus ON: %+0.3f deg, north on %d/%d, p = %.4f\n",
            mean(d), sum(d > 0), length(d),
            suppressWarnings(wilcox.test(d)$p.value)))
cat("  (perfect arm gave +0.254, 6/6 tracks, p = 0.031)\n")
