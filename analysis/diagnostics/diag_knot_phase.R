# INVARIANCE TEST: does the answer depend on where the knot grid happens to start?
#
# `knot_times <- t_start + (0:(k_steps-1)) * t_step` (R/TwilightFreeGrid.R:220) is a
# rigid grid anchored on the FIRST OBSERVATION TIMESTAMP, in UTC. Each knot owns the
# preceding half-open window (`obs_times[j] > t_prev && <= t_curr`) and the position
# is held constant across it. Nothing aligns that grid to solar time.
#
# So a 12-hour bin may contain one twilight, both, or neither, purely as a function
# of an arbitrary origin; and if a boundary BISECTS a twilight, that twilight is
# split across two knots fitted at different positions, each seeing a truncated
# ramp. Latitude is read from the twilight slope, and a truncated slope is a BIASED
# estimator of it, not merely a noisier one.
#
# WHY THIS TEST IS DIFFERENT FROM EVERYTHING TRIED TODAY. It is an invariance check,
# not an explanation. The data are the same; only the arbitrary binning origin
# moves. The answer MUST NOT depend on it. A positive result cannot be argued away
# as "correct probability, wrong biology" the way the spherical drift was -- it
# would simply be wrong. A null result exonerates the binning outright.
#
# METHOD. The knot phase is shifted by trimming the first `phase` hours of
# observations, which moves `t_start` and hence the whole grid. That is a real, if
# tiny, change to the data: at most 9 hours out of ~5760 (0.16%), all of it at the
# START, where the position is pinned to truth anyway. To keep the arms strictly
# comparable, every arm is SCORED only on knots inside a COMMON window (2 days in
# from each end), so no arm is credited or debited for knots the others do not have.
#
# MECHANISTIC COVARIATE. For each arm we also record the fraction of knot BOUNDARIES
# that land in twilight (solar zenith 85-95 at the true position). That is the
# direct measure of "how often does a boundary bisect a twilight", and if the
# binning is the mechanism, bias should track it.
#
# 4 phases x 6 tracks = 24 fits, checkpointed on (id, phase).
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))
suppressMessages(devtools::load_all(".", quiet = TRUE))
BAT <- readRDS("scratch/nes_calibration/light_battery.rds")
OUT <- "scratch/nes_calibration/knot_phase.csv"
STEP_H <- 12; CELL <- 1; PSLAB <- 0.10; DIFF <- 110
PHASES <- c(0, 3, 6, 9)          # hours
TRIM_DAYS <- 2                    # common scoring window, in from each end

REC <- Filter(function(b) b$offset == 0, BAT)

fit_one <- function(b, phase) {
  t0 <- min(as.numeric(b$time))
  keep <- as.numeric(b$time) >= t0 + phase * 3600
  tt <- b$time[keep]; yy <- b$perfect[keep]
  la <- b$lat[keep];  lo <- b$lon[keep]

  gg <- rast(xmin = floor(min(b$lon)) - 12, xmax = ceiling(max(b$lon)) + 12,
             ymin = max(-88, floor(min(b$lat)) - 12),
             ymax = min(88, ceiling(max(b$lat)) + 12),
             resolution = CELL, crs = "EPSG:4326")
  values(gg) <- 1
  lam <- 1 / (b$max_light * 0.5)
  invisible(capture.output(
    f <- TwilightFreeGrid(tt, yy, grid = gg,
      start_lon = lo[1], start_lat = la[1],
      end_lon = lo[length(lo)], end_lat = la[length(la)],
      step_hours = STEP_H, diffusion = DIFF, calibration = b$response,
      likelihood_params = c(lam, b$max_light, PSLAB),
      area_correction = TRUE)))

  gp <- grid_posterior(f); K <- nrow(gp$P)
  mu <- sd <- rep(NA_real_, K)
  for (k in seq_len(K)) {
    w <- gp$P[k, ]; s <- sum(w)
    if (!is.finite(s) || s <= 0) next
    w <- w / s
    mu[k] <- sum(w * gp$lat)
    sd[k] <- sqrt(sum(w * (gp$lat - mu[k])^2))
  }

  tk   <- as.numeric(f$fit$time)
  ttn  <- as.numeric(b$time)
  tlat <- approx(ttn, b$lat, tk, rule = 2)$y
  tlon <- approx(ttn, b$lon, tk, rule = 2)$y
  sup  <- approx(ttn, as.numeric(b$supported), tk, rule = 2)$y > 0.999

  # how often does a knot BOUNDARY land in twilight, at the true position?
  zb <- solar_zenith(tk, tlon, tlat)
  in_twilight <- mean(zb >= 85 & zb <= 95, na.rm = TRUE)

  data.table(id = b$id, phase = phase, k = seq_len(K), time = tk,
             sup = sup, truth = tlat, mu = mu, sd = sd, mode = f$fit$lat,
             zbound = zb, twi_frac = in_twilight)
}

done <- if (file.exists(OUT)) {
  d0 <- fread(OUT, colClasses = list(character = "id"))
  unique(paste(d0$id, d0$phase))
} else character(0)
cat(sprintf("%d of %d fits already done\n\n", length(done), length(REC) * length(PHASES)))

for (p in PHASES) for (b in REC) {
  key <- paste(b$id, p)
  if (key %in% done) next
  r <- try(fit_one(b, p), silent = TRUE)
  if (inherits(r, "try-error")) {
    message("FAILED ", key, ": ", conditionMessage(attr(r, "condition"))); next
  }
  fwrite(r, OUT, append = file.exists(OUT))
  message(key, " done")
}

D <- fread(OUT, colClasses = list(character = "id"))
cat(sprintf("\n%d fits, %d knots\n", uniqueN(paste(D$id, D$phase)), nrow(D)))

# COMMON scoring window per track, so no arm is scored on knots the others lack
D[, `:=`(tmin = min(time), tmax = max(time)), by = id]
W <- D[sup == TRUE & is.finite(mu) &
       time >= tmin + TRIM_DAYS * 86400 & time <= tmax - TRIM_DAYS * 86400]
W[, e := mu - truth]
W[, frac := (time - min(time)) / (max(time) - min(time)), by = .(id, phase)]

cat(sprintf("scored on the common window: %d knots per arm (mean)\n",
            round(nrow(W) / uniqueN(paste(W$id, W$phase)))))

cat("\n=== THE INVARIANCE TEST: bias by knot-grid phase ===\n")
a <- W[, .(n = .N, bias = round(mean(e), 3), sd_e = round(sd(e), 3),
           post_sd = round(mean(sd), 3),
           twi_bound = round(mean(twi_frac), 3)), by = phase][order(phase)]
print(as.data.frame(a), row.names = FALSE)
cat("  the data are the same; only the binning origin moved.\n")

cat("\n=== per track, so this is not one animal driving it ===\n")
w <- dcast(W[, .(b = round(mean(e), 2)), by = .(id, phase)], id ~ phase,
           value.var = "b")
print(as.data.frame(w), row.names = FALSE)
rng <- W[, .(b = mean(e)), by = .(id, phase)][, .(spread = max(b) - min(b)), by = id]
cat(sprintf("\n  within-track spread across phase: min %.2f, median %.2f, max %.2f deg\n",
            min(rng$spread), median(rng$spread), max(rng$spread)))

cat("\n=== does bias track boundaries landing in twilight? ===\n")
bt <- W[, .(b = mean(e), tw = mean(twi_frac)), by = .(id, phase)]
cat(sprintf("  Spearman(bias, twilight-boundary fraction) = %+0.3f over %d arms\n",
            suppressWarnings(cor(bt$b, bt$tw, method = "spearman")), nrow(bt)))

cat("\n=== does the mid-track hump move with phase? ===\n")
W[, q := cut(frac, seq(0, 1, 0.2), include.lowest = TRUE, labels = FALSE)]
tab <- dcast(W[, .(b = round(mean(e), 2)), by = .(phase, q)], phase ~ q, value.var = "b")
setnames(tab, c("phase", "0-20%", "20-40%", "40-60%", "60-80%", "80-100%"))
print(as.data.frame(tab), row.names = FALSE)

cat("\nREADING\n")
cat("  Bias moving materially with phase => the binning is an artefact and the\n")
cat("  knot grid needs aligning to solar time rather than to an arbitrary origin.\n")
cat("  Bias flat across phase => binning alignment is exonerated; keep the\n")
cat("  seasonal/equinox line of enquiry.\n")
