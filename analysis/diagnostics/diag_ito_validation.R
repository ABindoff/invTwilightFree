# Does the Ito correction remove the equatorward drift at FIT level?
#
# The kernel diagnostic showed the drift is real, matches -tan(phi) sigma^2/(2R^2)
# to correlation 0.99994, and scales as sigma^2. `drift_correction = TRUE` now
# re-centres the movement kernel so latitude is a martingale under the prior. This
# is the test that it actually works on fits, not just on the kernel in isolation.
#
# THE PREDICTION, stated before running (both halves must hold):
#
#   1. bias at D = 110 moves from -0.567 toward zero.
#   2. THE DIAGNOSTIC ONE: the D = 440 arm stops being ~3x worse than D = 110.
#      The drift scales as sigma^2, so if that is what was being corrected, the
#      gap between the two diffusions should largely close. If bias improves at
#      D = 110 but D = 440 is still ~3x worse, the correction is being absorbed
#      somewhere rather than removing the mechanism, and the improvement is a
#      coincidence of magnitude.
#
# Design mirrors diag_area_diffusion.R exactly so the two are directly comparable:
# same 6 tracks, same offset 0, same grid, same perfect light, area_correction ON
# throughout. `drift0_D110` is the REPRODUCTION CONTROL and must return -0.567.
#
# Checkpointed per fit, keyed on (id, arm).
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))
suppressMessages(devtools::load_all(".", quiet = TRUE))
BAT <- readRDS("scratch/nes_calibration/light_battery.rds")
OUT <- "scratch/nes_calibration/ito_validation.csv"
STEP_H <- 12; CELL <- 1; PSLAB <- 0.10

ARMS <- expand.grid(drift = c(FALSE, TRUE), diff = c(110, 440))
ARMS$name <- sprintf("drift%d_D%d", as.integer(ARMS$drift), ARMS$diff)

REC <- Filter(function(b) b$offset == 0, BAT)

fit_one <- function(b, drift, dif) {
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
      area_correction = TRUE, drift_correction = drift)))

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

  data.table(id = b$id, arm = sprintf("drift%d_D%d", as.integer(drift), dif),
             drift = drift, diff = dif, n = sum(ok),
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
  r <- try(fit_one(b, ARMS$drift[ai], ARMS$diff[ai]), silent = TRUE)
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

a <- D[, .(n = .N, bias_mode = round(mean(bias_mode), 3),
           bias_mean = round(mean(bias_mean), 3),
           rmse = round(mean(rmse_mode), 3), lat_sd = round(mean(lat_sd), 3),
           cover = round(mean(cover), 3), km = round(median(km))), by = arm][order(arm)]
cat("\n=== 2x2 ===\n")
print(as.data.frame(a), row.names = FALSE)

ctl <- D[arm == "drift0_D110"]
if (nrow(ctl) > 0)
  cat(sprintf("\nREPRODUCTION CONTROL drift0_D110: %+0.3f (expected -0.567) %s\n",
              mean(ctl$bias_mode),
              if (abs(mean(ctl$bias_mode) + 0.567) < 0.1) "OK" else "*** DIFFERS"))

for (dv in sort(unique(D$diff))) {
  w <- dcast(D[diff == dv], id ~ drift, value.var = "bias_mode")
  if (!all(c("TRUE", "FALSE") %in% names(w))) next
  w <- w[complete.cases(w)]
  d <- w[["TRUE"]] - w[["FALSE"]]
  cat(sprintf("\nD=%d  drift ON - OFF: %+0.3f deg, north on %d/%d, p=%.4f\n",
              dv, mean(d), sum(d > 0), length(d),
              suppressWarnings(wilcox.test(d)$p.value)))
}

cat("\n=== PREDICTION 2: does the D=440 penalty close? ===\n")
for (dr in c(FALSE, TRUE)) {
  w <- dcast(D[drift == dr], id ~ diff, value.var = "bias_mode")
  if (!all(c("110", "440") %in% names(w))) next
  w <- w[complete.cases(w)]
  cat(sprintf("  drift %-5s : D110 %+0.3f, D440 %+0.3f, gap %+0.3f (ratio %.2f)\n",
              dr, mean(w[["110"]]), mean(w[["440"]]),
              mean(w[["440"]]) - mean(w[["110"]]),
              mean(w[["440"]]) / mean(w[["110"]])))
}
cat("\n  Before the correction the gap was -1.103 deg (ratio 2.95).\n")
cat("  If the ratio collapses toward 1, the sigma^2 mechanism is what was fixed.\n")
