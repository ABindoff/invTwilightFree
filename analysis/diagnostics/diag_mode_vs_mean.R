# Is the null gate's -0.64 deg a POSTERIOR bias or a READOUT artefact?
#
# The gate showed -0.643 deg of latitude bias on perfect, correctly-modelled light.
# Zenith is excluded (the battery generates with the same `solar_zenith()` the
# engine fits with, so a shared error cancels and truth stays the exact maximiser)
# and Argos is excluded (light is generated AT the Argos positions and scored
# against them). What is left is how the point estimate is READ OFF the posterior.
#
# `f$fit$lat` is the JOINT posterior MODE: the single highest-probability cell,
# `best_lat[k] = lat[best_i]` in run_grid_hmm. On a 1-degree grid that is quantised
# to cell centres, and the mode of a skewed density is biased even when the density
# itself is exactly right.
#
# STAGE 1 (this script). Re-fit the SAME records at the SAME resolution as the gate,
# but record the posterior MEAN and MEDIAN alongside the mode. This costs nothing
# beyond the fits themselves and is the primary discriminator:
#
#   mean ~unbiased, mode -0.64        -> readout artefact. Stage 2 unnecessary.
#   mean also ~-0.6                   -> the POSTERIOR is genuinely displaced, and
#                                        quantisation is not the story.
#
# Staging it this way is deliberate. A resolution sweep is the obvious test but it
# is expensive -- the grid HMM is O(K n^2), so halving the cell size is ~16x the
# work, not 4x -- and this campaign has repeatedly paid for experiments that a
# cheaper gating test would have redirected.
#
# Perfect light only: it is the cleanest null, and the noisy level adds sampling
# noise without adding a distinct mechanism.
#
# Checkpointed per fit, keyed on (id, offset).
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))
suppressMessages(devtools::load_all(".", quiet = TRUE))
BAT <- readRDS("scratch/nes_calibration/light_battery.rds")
OUT <- "scratch/nes_calibration/mode_vs_mean.csv"
STEP_H <- 12; DIFF <- 110; CELL <- 1; PSLAB <- 0.10

# weighted quantile on the marginal latitude distribution
wquant <- function(lat, w, p = 0.5) {
  o <- order(lat); lat <- lat[o]; w <- w[o]
  cw <- cumsum(w) / sum(w)
  lat[which(cw >= p)[1]]
}

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
      step_hours = STEP_H, diffusion = DIFF,
      calibration = b$response,
      likelihood_params = c(lam, b$max_light, PSLAB))))

  tk <- as.numeric(f$fit$time); tt <- as.numeric(b$time)
  tlat <- approx(tt, b$lat, tk, rule = 2)$y
  sup  <- approx(tt, as.numeric(b$supported), tk, rule = 2)$y > 0.999

  gp <- grid_posterior(f)
  K <- nrow(gp$P)
  mean_lat <- med_lat <- marg_mode <- skew <- rep(NA_real_, K)
  for (k in seq_len(K)) {
    w <- gp$P[k, ]; s <- sum(w)
    if (!is.finite(s) || s <= 0) next
    w <- w / s
    mu <- sum(w * gp$lat)
    sd <- sqrt(sum(w * (gp$lat - mu)^2))
    mean_lat[k] <- mu
    med_lat[k]  <- wquant(gp$lat, w)
    # marginal-latitude mode: collapse longitude first, then take the peak
    ag <- tapply(w, gp$lat, sum)
    marg_mode[k] <- as.numeric(names(ag))[which.max(ag)]
    skew[k] <- if (sd > 0) sum(w * ((gp$lat - mu) / sd)^3) else NA_real_
  }

  e_mode <- f$fit$lat - tlat      # what the gate scored: JOINT mode
  e_mean <- mean_lat  - tlat
  e_med  <- med_lat   - tlat
  e_mmod <- marg_mode - tlat
  ok <- sup & is.finite(e_mode) & is.finite(e_mean)

  data.table(id = b$id, offset = b$offset, n = sum(ok),
             bias_mode      = mean(e_mode[ok]),
             bias_mean      = mean(e_mean[ok]),
             bias_median    = mean(e_med[ok]),
             bias_margmode  = mean(e_mmod[ok]),
             rmse_mode      = sqrt(mean(e_mode[ok]^2)),
             rmse_mean      = sqrt(mean(e_mean[ok]^2)),
             skew           = mean(skew[ok], na.rm = TRUE))
}

done <- if (file.exists(OUT)) {
  d0 <- fread(OUT, colClasses = list(character = "id"))
  paste(d0$id, d0$offset)
} else character(0)
cat(sprintf("%d of %d fits already done\n\n", length(done), length(BAT)))

for (b in BAT) {
  key <- paste(b$id, b$offset)
  if (key %in% done) next
  r <- try(fit_one(b), silent = TRUE)
  if (inherits(r, "try-error")) {
    message("FAILED ", key, ": ", conditionMessage(attr(r, "condition")))
    next
  }
  fwrite(r, OUT, append = file.exists(OUT))
  message(key, " done")
}

D <- fread(OUT, colClasses = list(character = "id"))
cat(sprintf("\n%d of %d fits\n", nrow(D), length(BAT)))
if (nrow(D) < length(BAT))
  cat("*** INCOMPLETE -- provisional.\n")

cl <- D[, .(mode = mean(bias_mode), mean = mean(bias_mean),
            median = mean(bias_median), margmode = mean(bias_margmode),
            rmse_mode = mean(rmse_mode), rmse_mean = mean(rmse_mean),
            skew = mean(skew)), by = id]

cat("\n=== latitude bias by READOUT (degrees, by track) ===\n")
cat(sprintf("  joint mode (what the gate scored) : %+0.3f\n", mean(cl$mode)))
cat(sprintf("  marginal-latitude mode            : %+0.3f\n", mean(cl$margmode)))
cat(sprintf("  posterior median                  : %+0.3f\n", mean(cl$median)))
cat(sprintf("  posterior MEAN                    : %+0.3f\n", mean(cl$mean)))
cat(sprintf("\n  rmse: mode %.3f | mean %.3f\n", mean(cl$rmse_mode), mean(cl$rmse_mean)))
cat(sprintf("  mean posterior skewness           : %+0.3f\n", mean(cl$skew)))

cat("\nper track (mode -> mean):\n")
for (i in seq_len(nrow(cl)))
  cat(sprintf("  %-9s %+0.2f -> %+0.2f\n", cl$id[i], cl$mode[i], cl$mean[i]))

d <- cl$mode - cl$mean
cat(sprintf("\npaired mode - mean: %+0.3f deg, same sign on %d/%d, p = %.4f\n",
            mean(d), sum(d < 0), length(d),
            suppressWarnings(wilcox.test(d)$p.value)))

cat("\nVERDICT\n")
# braces required: a bare `else` starting a new line at top level is a parse error
if (abs(mean(cl$mean)) < 0.25 && abs(mean(cl$mode)) > 0.4) {
  cat("  READOUT ARTEFACT: the posterior is close to unbiased; the joint mode is not.\n  Stage 2 (resolution sweep) is unnecessary -- report the mean, or refine the grid.\n")
} else if (abs(mean(cl$mean)) > 0.4) {
  cat("  GENUINE POSTERIOR BIAS: the mean is displaced too, so this is not quantisation.\n  Stage 2 should sweep resolution to check discretisation of the LIKELIHOOD, not the readout.\n")
} else {
  cat("  INTERMEDIATE: neither clean. Read the per-track numbers before designing stage 2.\n")
}
