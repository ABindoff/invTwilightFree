# Per-tag observation noise, estimated rather than assumed.
#
# The target: latitude scatter varies fourfold across this panel (1.47 deg to
# 6.16) while the model hands every animal nearly the same uncertainty. That is
# why coverage runs 0.13 to 0.91 between tags.
#
# WHY NOT A TEMPERING PARAMETER. Dividing tag i's log-likelihood by tau_i is not
# identifiable from the likelihood it tempers: log_z is essentially monotone in
# tau, so there is nothing to maximise. That is the deeper reason `overcount`
# failed as a knob.
#
# WHY LAMBDA WORKS. The spike rate is a parameter INSIDE the likelihood -- the
# scale of the exponential spike around the expected light. Too large and the
# model insists observations match the curve tightly; too small and it learns
# nothing from them. So log_z has an INTERIOR optimum in lambda, and the grid
# HMM already computes log_z exactly. Selecting lambda per tag by maximising it
# is empirical Bayes: no ground truth, no feedback loop (log_z integrates over
# all tracks rather than conditioning on a fitted one), and it is identifiable.
#
# Three things this has to show before it is worth putting in the package:
#   1. log_z really does have an interior optimum in lambda, per tag
#   2. the selected lambda VARIES between tags -- otherwise a single panel value
#      is the honest answer and nothing is gained
#   3. selecting it improves COVERAGE. Accuracy is a secondary matter; the
#      problem being fixed is that noisy tags are handed over-confident intervals
#
# The default is lambda = 1/(max_light * 0.5); the sweep is a multiplier on it.
# RESUMABLE: appends per fit, skips what is already done.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))
OUT <- "scratch/nes_calibration/lambda_results.csv"

STEP_HOURS <- 12; DIFFUSION <- 110; CELL <- 1; CAL_DAYS <- 15
MULT <- c(0.25, 0.5, 1, 2, 4)

arch <- lapply(readRDS("analysis/cache/nes/archives_v3_30min.rds"), `[[`, "main")
meta <- nes_meta(); ar <- nes_argos()

tags <- list()
for (id in names(arch)) {
  m <- meta[match(id, meta$id), ]
  if (is.na(m$ptt) || !nzchar(m$ptt)) next
  df <- ar$deploy[Ptt == m$ptt]; a <- ar$fixes[Ptt == m$ptt]
  if (!nrow(df) || nrow(a) < 50) next
  x <- as.data.table(arch[[id]])
  t0 <- max(min(x$time), as.POSIXct(max(m$deploy_date, df$user_date), tz = "UTC"))
  t1 <- min(max(x$time), max(a$time))
  d <- x[time >= t0 & time <= t1]
  if (nrow(d) < 1000) next
  tf <- a[time >= max(time) - 72*3600]
  tags[[id]] <- list(id = id, light = d, argos = a,
                     p0 = c(df$deploy_lon, df$deploy_lat),
                     p1 = c(median(tf$lon), median(tf$lat)))
}
grid <- rast(xmin = 150, xmax = 250, ymin = 20, ymax = 70,
             resolution = CELL, crs = "EPSG:4326")
values(grid) <- 1

score <- function(id, r, mult) {
  g <- tags[[id]]
  lam <- mult / (r$max_light * 0.5)
  t <- system.time(invisible(capture.output(
    f <- TwilightFreeGrid(g$light$time, pmax(0, g$light$light - r$baseline), grid = grid,
      start_lon = unname(g$p0[1]), start_lat = unname(g$p0[2]),
      end_lon = unname(g$p1[1]), end_lat = unname(g$p1[2]),
      step_hours = STEP_HOURS, diffusion = DIFFUSION,
      calibration = r$calibration,
      likelihood_params = c(lam, r$max_light, 0.10)))))[3]
  gp <- grid_posterior(f)
  sdp <- sdl <- numeric(nrow(gp$P))
  for (i in seq_len(nrow(gp$P))) {
    w <- gp$P[i, ]; s <- sum(w)
    if (!is.finite(s) || s <= 0) { sdp[i] <- sdl[i] <- NA; next }
    w <- w / s
    sdl[i] <- sqrt(sum(w * (gp$lon - sum(w * gp$lon))^2))
    sdp[i] <- sqrt(sum(w * (gp$lat - sum(w * gp$lat))^2))
  }
  tm <- as.POSIXct(f$fit$time, origin = "1970-01-01", tz = "UTC")
  tr <- argos_at(g$argos, tm)
  ep <- f$fit$lat - tr$lat; el <- dlon(lon360(f$fit$lon), tr$lon)
  e <- gc_km(lon360(f$fit$lon), f$fit$lat, tr$lon, tr$lat)
  ok <- is.finite(ep) & is.finite(sdp) & sdp > 0
  data.frame(id = id, mult = mult, lambda = signif(lam, 4),
    log_z = round(f$log_z, 2), median_km = round(median(e[ok])),
    bias_lat = round(mean(ep[ok]), 2),
    rmse_lat = round(sqrt(mean(ep[ok]^2)), 2),
    post_sd_lat = round(median(sdp[ok]), 2),
    cover_lat = round(mean(abs(ep[ok]) <= 1.96 * sdp[ok]), 2),
    cover_lon = round(mean(abs(el[ok]) <= 1.96 * sdl[ok]), 2),
    secs = round(t), row.names = NULL)
}

scale_fits <- lapply(tags, function(g) {
  k <- g$light$time <= min(g$light$time) + CAL_DAYS * 86400
  fit_light_response(g$light$time[k], g$light$light[k], g$p0[1], g$p0[2])
})
geom_fits <- lapply(tags, function(g) {
  dep <- departure_time(g$light)
  k <- is.finite(dep) & as.numeric(g$light$time) < dep
  if (sum(k) < 200) return(NULL)
  fit_light_response(g$light$time[k], g$light$light[k], g$p0[1], g$p0[2])
})
responses <- pool_light_responses(geom_fits, scale_fits)

done <- if (file.exists(OUT)) {
  d <- data.table::fread(OUT); paste(d$id, d$mult) } else character(0)
cat(sprintf("%d fits on disk, %d to do\n", length(done),
            length(tags) * length(MULT) - length(done)))

for (mu in MULT) for (id in names(tags)) {
  if (paste(id, mu) %in% done) next
  message("lambda x", mu, " : ", id)
  data.table::fwrite(score(id, responses[[id]], mu), OUT, append = file.exists(OUT))
}

o <- data.table::fread(OUT)
cat("\n=== pooled over tags ===\n")
print(as.data.frame(o[, .(tags = .N, median_km = round(mean(median_km)),
  rmse_lat = round(mean(rmse_lat), 2), post_sd_lat = round(mean(post_sd_lat), 2),
  cover_lat = round(mean(cover_lat), 2), cover_lon = round(mean(cover_lon), 2),
  mean_log_z = round(mean(log_z), 1)), by = mult][order(mult)]), row.names = FALSE)

cat("\n=== 1. does log_z have an interior optimum, per tag? ===\n")
best <- o[, .(best_mult = mult[which.max(log_z)],
              interior = which.max(log_z) > 1 & which.max(log_z) < length(MULT),
              d_log_z = round(max(log_z) - log_z[mult == 1], 1)), by = id]
print(as.data.frame(best), row.names = FALSE)
cat(sprintf("\ninterior optimum for %d of %d tags\n", sum(best$interior), nrow(best)))

cat("\n=== 2. does the selected lambda vary between tags? ===\n")
print(table(best$best_mult))

cat("\n=== 3. what would selecting it per tag buy? ===\n")
sel <- merge(o, best[, .(id, mult = best_mult)], by = c("id", "mult"))
fix <- o[mult == 1]
cat(sprintf("fixed lambda (current):  median %3.0f km, rmse_lat %.2f, cover_lat %.2f\n",
            mean(fix$median_km), mean(fix$rmse_lat), mean(fix$cover_lat)))
cat(sprintf("per-tag lambda by log_z: median %3.0f km, rmse_lat %.2f, cover_lat %.2f\n",
            mean(sel$median_km), mean(sel$rmse_lat), mean(sel$cover_lat)))
cat("\nIf the selected lambda tracks the tags that are genuinely noisy, coverage\n")
cat("should even out between tags -- that is the point, not the mean error.\n")
print(as.data.frame(sel[, .(id, mult, log_z, median_km, rmse_lat, cover_lat)][order(mult)]),
      row.names = FALSE)
