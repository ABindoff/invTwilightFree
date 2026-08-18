# TANGENT vs LOGISTIC ON REAL LIGHT. Does the response form matter now?
#
# WHAT PROMPTED THIS (diag_real_residual.R, 6 tags, 68k supported observations,
# residuals against each form on the baselined scale the engine sees):
#
#   zenith    P(below), tangent   P(below), logistic   [model assumes 0.667]
#   day <80        0.774                0.754
#   88-92          0.110                0.296
#   92-96          0.071                0.318
#
# In the twilight band -- where latitude is read -- the TANGENT sits below the
# observations 89-93% of the time. Signed offset computed from those columns:
# observations lie **+0.12 max_light above the tangent** through 88-96 deg, against
# **+0.055** for the logistic. The response the engine is actually handed
# systematically under-predicts twilight light, by about twice as much.
#
# That has the right sign for the observed bias: if the model expects less light at a
# given twilight zenith than really occurs, the fit resolves it by moving to where
# the sun sits higher -- the longer-day hemisphere, equatorward in northern winter.
#
# WHY RE-TEST SOMETHING ALREADY DECIDED. The logistic was compared against the
# tangent once before and lost badly (833 vs 354 km), and the tangent was adopted.
# But that comparison predates three changes that all bear on it: the
# calibration-window fix (pooled haul-out geometry + per-tag 15-day scale), tying the
# tangent's amplitude to `max_light` rather than the envelope's `amp`, and the
# cos(latitude) area correction. The original decision also rested partly on
# `spike_normaliser` behaving best when the expected curve attains both ends, which
# is an argument about the normaliser rather than about fit quality.
#
# REAL LIGHT, not the battery: the battery's noisy arm is GENERATED from the
# logistic, so testing response forms on it would be circular. Here neither form is
# correct, which is the situation that matters.
#
# PRE-REGISTERED: if the twilight under-prediction is doing real damage, the logistic
# should show less equatorward latitude bias. If the tangent still wins, the residual
# offset is not what drives the bias and the response form is settled.
#
# Truth = battery offset-0 Argos positions. Light = decimated archives, baselined
# exactly as the pipeline does: pmax(0, light - q05).
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))
suppressMessages(devtools::load_all(".", quiet = TRUE))
BAT <- readRDS("scratch/nes_calibration/light_battery.rds")
REC <- Filter(function(z) z$offset == 0, BAT)
arch <- lapply(readRDS("analysis/cache/nes/archives_v3_30min.rds"), `[[`, "main")
fv <- try(lapply(readRDS("analysis/cache/nes/fvilches_archives_v1_30min.rds"),
                 `[[`, "main"), silent = TRUE)
if (!inherits(fv, "try-error")) arch <- c(arch, fv)

OUT <- "scratch/nes_calibration/response_form_real.csv"
STEP_H <- 12; CELL <- 1; PSLAB <- 0.10; DIFF <- 110

fit_one <- function(b, arm) {
  g <- arch[[b$id]]
  if (is.null(g) || !nrow(g)) return(NULL)
  ML <- b$max_light; r <- b$response
  gt <- as.numeric(g$time); bt <- as.numeric(b$time)
  base <- as.numeric(quantile(as.numeric(g$light), 0.05, na.rm = TRUE))
  y <- pmax(0, approx(gt, as.numeric(g$light), bt, rule = 2)$y - base)

  cal <- if (arm == "logistic") r else {
    slope <- ML / (4 * r[4]); zero <- r[3] + 2 * r[4]
    c(slope * zero, slope)          # engine computes cal[1] - cal[2]*z, clamped
  }

  gg <- rast(xmin = floor(min(b$lon)) - 12, xmax = ceiling(max(b$lon)) + 12,
             ymin = max(-88, floor(min(b$lat)) - 12),
             ymax = min(88, ceiling(max(b$lat)) + 12),
             resolution = CELL, crs = "EPSG:4326")
  values(gg) <- 1
  lam <- 1 / (ML * 0.5)
  n <- length(bt)
  invisible(capture.output(
    f <- TwilightFreeGrid(b$time, y, grid = gg,
      start_lon = b$lon[1], start_lat = b$lat[1],
      end_lon = b$lon[n], end_lat = b$lat[n],
      step_hours = STEP_H, diffusion = DIFF, calibration = cal,
      likelihood_params = c(lam, ML, PSLAB), area_correction = TRUE)))

  gp <- grid_posterior(f); K <- nrow(gp$P)
  mu <- sdl <- rep(NA_real_, K)
  for (k in seq_len(K)) {
    w <- gp$P[k, ]; s <- sum(w)
    if (!is.finite(s) || s <= 0) next
    w <- w / s
    mu[k] <- sum(w * gp$lat); sdl[k] <- sqrt(sum(w * (gp$lat - mu[k])^2))
  }
  tk <- as.numeric(f$fit$time)
  tlat <- approx(bt, b$lat, tk, rule = 2)$y
  tlon <- approx(bt, lon360(b$lon), tk, rule = 2)$y
  sup  <- approx(bt, as.numeric(b$supported), tk, rule = 2)$y > 0.999
  e <- mu - tlat
  ok <- sup & is.finite(e) & is.finite(sdl) & sdl > 0
  data.table(id = b$id, arm = arm, n = sum(ok), log_z = f$log_z,
             bias = mean(e[ok]), bias_mode = mean((f$fit$lat - tlat)[ok]),
             rmse = sqrt(mean(e[ok]^2)), sd_lat = mean(sdl[ok]),
             cover = mean(abs(e[ok]) <= 1.96 * sdl[ok]),
             km = median(gc_km(lon360(f$fit$lon)[ok], f$fit$lat[ok],
                               tlon[ok], tlat[ok])))
}

done <- if (file.exists(OUT)) {
  d0 <- fread(OUT, colClasses = list(character = "id"))
  paste(d0$id, d0$arm)
} else character(0)
cat(sprintf("%d of %d fits already done\n\n", length(done), 2 * length(REC)))

for (a in c("tangent", "logistic")) for (b in REC) {
  key <- paste(b$id, a)
  if (key %in% done) next
  r <- try(fit_one(b, a), silent = TRUE)
  if (inherits(r, "try-error") || is.null(r)) {
    message("skip/failed ", key); next
  }
  fwrite(r, OUT, append = file.exists(OUT))
  message(key, " done")
}

D <- fread(OUT, colClasses = list(character = "id"))
cat(sprintf("\n%d fits, %d tags\n\n=== REAL LIGHT, scored against Argos ===\n",
            nrow(D), uniqueN(D$id)))
print(as.data.frame(D[, .(n = .N, bias = round(mean(bias), 3),
                          bias_mode = round(mean(bias_mode), 3),
                          rmse = round(mean(rmse), 3),
                          sd_lat = round(mean(sd_lat), 2),
                          cover = round(mean(cover), 3),
                          km = round(median(km)),
                          log_z = round(mean(log_z))), by = arm]), row.names = FALSE)

w <- dcast(D, id ~ arm, value.var = c("bias", "km", "cover", "log_z"))
if (all(c("bias_tangent", "bias_logistic") %in% names(w))) {
  w <- w[complete.cases(w)]
  cat(sprintf("\npaired within tag (n = %d):\n", nrow(w)))
  cat(sprintf("  d|bias|  %+0.3f deg   logistic better on %d/%d\n",
              mean(abs(w$bias_logistic)) - mean(abs(w$bias_tangent)),
              sum(abs(w$bias_logistic) < abs(w$bias_tangent)), nrow(w)))
  cat(sprintf("  d km     %+0.0f km      logistic better on %d/%d\n",
              mean(w$km_logistic) - mean(w$km_tangent),
              sum(w$km_logistic < w$km_tangent), nrow(w)))
  cat(sprintf("  d cover  %+0.3f\n", mean(w$cover_logistic) - mean(w$cover_tangent)))
  cat(sprintf("  d log_z  %+0.0f  (evidence prefers %s)\n",
              mean(w$log_z_logistic) - mean(w$log_z_tangent),
              if (mean(w$log_z_logistic) > mean(w$log_z_tangent)) "LOGISTIC" else "tangent"))
}
cat("\n  for reference: the earlier comparison, before the calibration-window fix,\n")
cat("  the max_light amplitude fix and the area correction, gave 833 km for the\n")
cat("  logistic against 354 for the tangent.\n")
