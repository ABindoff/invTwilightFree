# Does a correctly-scaled, per-axis movement prior reduce the latitude bias?
#
# MEASURED (diag_axis_scales.R, Argos truth, per-axis D(lag) = sd(disp)/sqrt(lag)):
#   lag 0.5 d : NS 23.0  EW  31.1
#   lag 5   d : NS 47.6  EW  78.6
#   lag 20  d : NS 77.3  EW 148.5
# Both axes are DIRECTED (D grows with lag) but east-west far more so: growth
# 0.5->20 d is 3.36x north-south and 4.77x east-west, and NS saturates near 78 while
# EW keeps extending. The shipped isotropic 110 is therefore wrong in OPPOSITE
# directions on the two axes -- at the 20 d lag that governs marginal intervals it is
# 1.42x too WIDE in latitude and 0.74x too NARROW in longitude. That is exactly the
# recorded coverage failure (latitude 0.64, longitude 0.99) from one isotropic number
# split between two different axes.
#
# WHY IT SHOULD MOVE THE BIAS: bias ~ tilt x sum-of-covariances, and the cell-area
# tilt is what is being amplified. Removing the area factor is worth +5.45 deg on
# this arm, so the tilt is large and its amplification by an over-wide latitude prior
# is the mechanism under test. Direct evidence it scales: perfect arm gave -0.567 at
# D=110 and -1.670 at D=440.
#
# PREDICTIONS, WRITTEN BEFORE RUNNING:
#   1. iso48 and aniso48/79 both reduce |latitude bias| relative to the shipped 110.
#   2. aniso48/79 beats iso48 on LONGITUDE (rmse and coverage), because it does not
#      also squeeze the axis that needed to be wider.
#   3. aniso77/148 has the best COVERAGE on both axes and worse bias than aniso48/79.
#      That is the bias/coverage trade, now split per axis instead of shared.
#   4. If aniso48/79 gives lower latitude bias AND latitude coverage no worse than
#      the shipped setting, the isotropic prior was costing both at once.
#
# THE HONEST ARM: `noisy` light, where the emission is correctly specified. The
# perfect arm is not a valid null (its tilt is -0.113 nats/deg against -0.002 here),
# so nothing measured on it is trusted.
#
# Control: iso110 must return latitude bias -0.487 (mean) / -0.453 (mode).
# 4 arms x 6 tracks = 24 fits, checkpointed on (id, arm).
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))
suppressMessages(devtools::load_all(".", quiet = TRUE))
BAT <- readRDS("scratch/nes_calibration/light_battery.rds")
REC <- Filter(function(z) z$offset == 0, BAT)
OUT <- "scratch/nes_calibration/axis_prior_fits.csv"
STEP_H <- 12; CELL <- 1; PSLAB <- 0.10

ARMS <- data.table(
  arm  = c("iso110", "iso48", "aniso48_79", "aniso77_148"),
  d_ns = c(110, 48, 48, 77),
  d_ew = c(NA, NA, 79, 148))

fit_one <- function(b, a) {
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
      step_hours = STEP_H, diffusion = a$d_ns,
      diffusion_lon = if (is.na(a$d_ew)) NULL else a$d_ew,
      calibration = b$response, likelihood_params = c(lam, b$max_light, PSLAB),
      area_correction = TRUE)))

  gp <- grid_posterior(f); K <- nrow(gp$P)
  mu <- sdl <- mun <- sdn <- rep(NA_real_, K)
  for (k in seq_len(K)) {
    w <- gp$P[k, ]; s <- sum(w)
    if (!is.finite(s) || s <= 0) next
    w <- w / s
    mu[k]  <- sum(w * gp$lat); sdl[k] <- sqrt(sum(w * (gp$lat - mu[k])^2))
    mun[k] <- sum(w * gp$lon); sdn[k] <- sqrt(sum(w * (gp$lon - mun[k])^2))
  }
  tk <- as.numeric(f$fit$time); tt <- as.numeric(b$time)
  tlat <- approx(tt, b$lat, tk, rule = 2)$y
  tlon <- approx(tt, lon360(b$lon), tk, rule = 2)$y
  sup  <- approx(tt, as.numeric(b$supported), tk, rule = 2)$y > 0.999
  el <- mu - tlat
  eo <- dlon(lon360(mun), tlon)
  ok <- sup & is.finite(el) & is.finite(sdl) & sdl > 0 & is.finite(sdn) & sdn > 0

  data.table(id = b$id, arm = a$arm, n = sum(ok),
             bias_lat  = mean(el[ok]),
             bias_mode = mean((f$fit$lat - tlat)[ok]),
             rmse_lat  = sqrt(mean(el[ok]^2)),
             rmse_lon  = sqrt(mean(eo[ok]^2)),
             sd_lat    = mean(sdl[ok]), sd_lon = mean(sdn[ok]),
             cov_lat   = mean(abs(el[ok]) <= 1.96 * sdl[ok]),
             cov_lon   = mean(abs(eo[ok]) <= 1.96 * sdn[ok]),
             km        = median(gc_km(lon360(f$fit$lon)[ok], f$fit$lat[ok],
                                      tlon[ok], tlat[ok])))
}

done <- if (file.exists(OUT)) {
  d0 <- fread(OUT, colClasses = list(character = "id"))
  paste(d0$id, d0$arm)
} else character(0)
cat(sprintf("%d of %d fits already done\n\n", length(done), nrow(ARMS) * length(REC)))

for (i in seq_len(nrow(ARMS))) for (b in REC) {
  key <- paste(b$id, ARMS$arm[i])
  if (key %in% done) next
  r <- try(fit_one(b, ARMS[i]), silent = TRUE)
  if (inherits(r, "try-error")) {
    message("FAILED ", key, ": ", conditionMessage(attr(r, "condition"))); next
  }
  fwrite(r, OUT, append = file.exists(OUT))
  message(key, " done")
}

D <- fread(OUT, colClasses = list(character = "id"))
D[, arm := factor(arm, levels = ARMS$arm)]
cat(sprintf("\n%d fits\n\n=== per-axis prior, noisy arm ===\n", nrow(D)))
print(as.data.frame(D[, .(n = .N,
                          bias_lat = round(mean(bias_lat), 3),
                          rmse_lat = round(mean(rmse_lat), 3),
                          rmse_lon = round(mean(rmse_lon), 3),
                          sd_lat = round(mean(sd_lat), 2),
                          sd_lon = round(mean(sd_lon), 2),
                          cov_lat = round(mean(cov_lat), 3),
                          cov_lon = round(mean(cov_lon), 3),
                          km = round(median(km))), by = arm][order(arm)]),
      row.names = FALSE)

ctl <- D[arm == "iso110"]
cat(sprintf("\nCONTROL iso110: bias_lat %+0.3f (expected -0.487) %s\n",
            mean(ctl$bias_lat),
            if (abs(mean(ctl$bias_lat) + 0.487) < 0.08) "OK" else "*** DIFFERS"))

cat("\n=== paired against the shipped setting, within track ===\n")
for (a in setdiff(ARMS$arm, "iso110")) {
  w <- merge(D[arm == "iso110", .(id, b0 = bias_lat, k0 = km, c0 = cov_lat)],
             D[arm == a, .(id, b1 = bias_lat, k1 = km, c1 = cov_lat)], by = "id")
  cat(sprintf("  %-12s : d|bias| %+0.3f (better %d/%d)  d km %+0.0f  d cov_lat %+0.3f\n",
              a, mean(abs(w$b1)) - mean(abs(w$b0)),
              sum(abs(w$b1) < abs(w$b0)), nrow(w),
              mean(w$k1) - mean(w$k0), mean(w$c1) - mean(w$c0)))
}
