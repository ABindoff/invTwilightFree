# STAGE 2. The null gate's latitude bias is in the POSTERIOR, not the readout.
# Which term in the geometry puts it there?
#
# Stage 1 settled that it is not how the estimate is read off: mode -0.577,
# marginal mode -0.618, median -0.685, MEAN -0.729. The mean is MORE biased than
# the mode, and the ordering is exactly what the measured skewness (-0.135)
# predicts, so the readouts differ only by skew and the posterior itself sits
# ~0.6-0.7 deg south. Zenith and Argos were already excluded by construction.
#
# With a correctly-specified emission, an ISOTROPIC movement prior cannot displace
# latitude -- there is nothing in it that distinguishes north from south. The one
# term that IS asymmetric in latitude is the CELL-AREA factor: weighting cells by
# cos(latitude) down-weights poleward cells and pushes estimates equatorward,
# which is the observed sign.
#
# That factor is not hypothetical here. Omitting it was measured at 0.204 nats per
# step between 38 N and 50 N, and turning the correction ON moved bias south by
# 1.0-1.4 deg on 16 of 16 tags. A term that large is worth checking for OVERSHOOT.
#
# THE 2x2:
#   area TRUE  / D 110  -- REPRODUCTION CONTROL. The gate scored -0.567 on these
#                          six records (perfect light, offset 0). If this arm does
#                          not return that, the harness differs and nothing here
#                          is readable.
#   area FALSE / D 110  -- the correction removed
#   area TRUE  / D 440  -- does the displacement track the movement prior at all?
#                          D 440 is sigma_disp, the value that gave latitude
#                          coverage 0.93 on real data.
#   area FALSE / D 440
#
# How to read it:
#   area off -> bias ~0        the correction overshoots
#   area off -> bias positive  the two bracket zero; the true factor is between
#   bias flat in both          it is in the EMISSION geometry, and only then does
#                              a resolution sweep (16x the cost) earn its keep
#
# Mode AND posterior mean are both recorded, so an area effect that moves one but
# not the other would be visible rather than assumed.
#
# Offset 0 only, 6 tracks: the gate's effect was consistent across all four
# offsets (-0.45 to -1.02) and all six tracks, so season is not what is being
# tested here and paying 4x for it buys nothing.
#
# Checkpointed per fit, keyed on (id, arm).
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))
suppressMessages(devtools::load_all(".", quiet = TRUE))
BAT <- readRDS("scratch/nes_calibration/light_battery.rds")
OUT <- "scratch/nes_calibration/area_diffusion.csv"
STEP_H <- 12; CELL <- 1; PSLAB <- 0.10

ARMS <- expand.grid(area = c(TRUE, FALSE), diff = c(110, 440),
                    stringsAsFactors = FALSE)
ARMS$name <- sprintf("area%d_D%d", as.integer(ARMS$area), ARMS$diff)

REC <- Filter(function(b) b$offset == 0, BAT)

fit_one <- function(b, area, dif) {
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
      step_hours = STEP_H, diffusion = dif,
      calibration = b$response,
      likelihood_params = c(lam, b$max_light, PSLAB),
      area_correction = area)))

  tk <- as.numeric(f$fit$time); tt <- as.numeric(b$time)
  tlat <- approx(tt, b$lat, tk, rule = 2)$y
  tlon <- approx(tt, lon360(b$lon), tk, rule = 2)$y
  sup  <- approx(tt, as.numeric(b$supported), tk, rule = 2)$y > 0.999

  gp <- grid_posterior(f); K <- nrow(gp$P)
  mean_lat <- sdp <- rep(NA_real_, K)
  for (k in seq_len(K)) {
    w <- gp$P[k, ]; s <- sum(w)
    if (!is.finite(s) || s <= 0) next
    w <- w / s
    mu <- sum(w * gp$lat)
    mean_lat[k] <- mu
    sdp[k] <- sqrt(sum(w * (gp$lat - mu)^2))
  }
  e_mode <- f$fit$lat - tlat
  e_mean <- mean_lat  - tlat
  ok <- sup & is.finite(e_mode) & is.finite(sdp) & sdp > 0

  data.table(id = b$id, arm = sprintf("area%d_D%d", as.integer(area), dif),
             area = area, diff = dif, n = sum(ok),
             bias_mode = mean(e_mode[ok]),
             bias_mean = mean(e_mean[ok]),
             rmse_mode = sqrt(mean(e_mode[ok]^2)),
             lat_sd    = mean(sdp[ok]),
             cover     = mean(abs(e_mode[ok]) <= 1.96 * sdp[ok]),
             km        = median(gc_km(lon360(f$fit$lon)[ok], f$fit$lat[ok],
                                      tlon[ok], tlat[ok])))
}

done <- if (file.exists(OUT)) {
  d0 <- fread(OUT, colClasses = list(character = "id"))
  paste(d0$id, d0$arm)
} else character(0)
n_want <- length(REC) * nrow(ARMS)
cat(sprintf("%d records x %d arms = %d fits; %d already done\n\n",
            length(REC), nrow(ARMS), n_want, length(done)))

for (ai in seq_len(nrow(ARMS))) for (b in REC) {
  key <- paste(b$id, ARMS$name[ai])
  if (key %in% done) next
  r <- try(fit_one(b, ARMS$area[ai], ARMS$diff[ai]), silent = TRUE)
  if (inherits(r, "try-error")) {
    message("FAILED ", key, ": ", conditionMessage(attr(r, "condition")))
    next
  }
  fwrite(r, OUT, append = file.exists(OUT))
  message(key, " done")
}

D <- fread(OUT, colClasses = list(character = "id"))
cat(sprintf("\n%d of %d fits\n", nrow(D), n_want))
if (nrow(D) < n_want) cat("*** INCOMPLETE -- provisional.\n")

a <- D[, .(n = .N, bias_mode = mean(bias_mode), bias_mean = mean(bias_mean),
           rmse = mean(rmse_mode), lat_sd = mean(lat_sd),
           cover = mean(cover), km = median(km)), by = arm][order(arm)]
cat("\n=== 2x2 ===\n")
print(as.data.frame(a[, lapply(.SD, function(x)
        if (is.numeric(x)) round(x, 3) else x)]), row.names = FALSE)

ctl <- D[arm == "area1_D110"]
if (nrow(ctl) > 0)
  cat(sprintf("\nREPRODUCTION CONTROL area1_D110: bias_mode %+0.3f (gate: -0.567)%s\n",
              mean(ctl$bias_mode),
              if (abs(mean(ctl$bias_mode) + 0.567) < 0.1) "  OK" else "  *** DIFFERS"))

# paired contrasts within track, at each diffusion
for (dv in unique(D$diff)) {
  x <- D[diff == dv]
  w <- dcast(x, id ~ area, value.var = "bias_mode")
  if (!all(c("TRUE", "FALSE") %in% names(w))) next
  d <- w[["FALSE"]] - w[["TRUE"]]
  cat(sprintf("\nD=%d, area OFF minus ON: %+0.3f deg, same sign %d/%d, p = %.4f\n",
              dv, mean(d), sum(d > 0), length(d),
              suppressWarnings(wilcox.test(d)$p.value)))
}
# paired contrast on diffusion, area ON
x <- D[area == TRUE]
w <- dcast(x, id ~ diff, value.var = "bias_mode")
if (all(c("110", "440") %in% names(w))) {
  d <- w[["440"]] - w[["110"]]
  cat(sprintf("\narea ON, D 440 minus 110: %+0.3f deg, same sign %d/%d, p = %.4f\n",
              mean(d), sum(d > 0), length(d),
              suppressWarnings(wilcox.test(d)$p.value)))
}
