# Can a tag say how noisy it is WITHOUT being told where it is?
#
# The target: per-tag latitude scatter varies fourfold across this panel (1.47
# deg for 2021033 to 6.16 for 2021032) while the model reports nearly the same
# uncertainty for every animal. That is the coverage problem. A per-tag
# dispersion inside the likelihood would fix it -- but only if something the tag
# knows about ITSELF predicts which tags are noisy.
#
# It has to be position-free, or we are back in the feedback loop that diverged
# twice in this project. So every index below is computed from the light series
# alone, with no zenith, no track and no Argos:
#
#   env_resid   how far the light sits below its own daily upper envelope
#   day_sd      day-to-day variability of the daily maximum
#   rough       lag-1 roughness of the within-day series, scaled by its range
#   dark_frac   fraction of daylight-hours observations that are near-dark
#   dive_atten  correlation between depth and light within a day (attenuation)
#
# Validation: correlate each against the MEASURED per-tag latitude scatter. If
# nothing correlates, no amount of implementation will help, and the honest
# answer is that the tag cannot tell.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))

k <- as.data.table(read.csv("analysis/output/nes/per_knot_errors.csv"))
k <- k[engine == "hier" & batch == "light" & is.finite(err_lat)]
k[, id := as.character(id)]
truth <- k[, .(offset = mean(err_lat), scatter = sd(err_lat),
               rmse_lat = sqrt(mean(err_lat^2)),
               cover = mean(abs(err_lat) <= 1.96 * lat_sd, na.rm = TRUE)), by = id]

arch <- lapply(readRDS("analysis/cache/nes/archives_v3_30min.rds"), `[[`, "main")
meta <- nes_meta(); ar <- nes_argos()

idx <- list()
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
  # drop the haul-out: we want the at-sea noise regime
  dep <- departure_time(d)
  if (is.finite(dep)) d <- d[as.numeric(time) >= dep]
  d[, day := as.integer(floor(as.numeric(time) / 86400))]

  rng <- diff(quantile(d$light, c(0.05, 0.95)))
  per_day <- d[, {
    n <- .N
    if (n < 24) .(env = NA_real_, mx = NA_real_, rr = NA_real_, dk = NA_real_,
                  da = NA_real_)
    else {
      # upper envelope of the day: a running max over a few hours
      w <- max(3, floor(n / 8))
      env <- sapply(seq_len(n), function(i)
        max(light[max(1, i-w):min(n, i+w)], na.rm = TRUE))
      bright <- env > quantile(env, 0.5)      # the daylight part of the day
      .(env = mean((env - light)[bright], na.rm = TRUE) / rng,
        mx  = max(light, na.rm = TRUE),
        rr  = mean(abs(diff(light)), na.rm = TRUE) / rng,
        dk  = mean(light[bright] < quantile(d$light, 0.2), na.rm = TRUE),
        da  = suppressWarnings(cor(depth_max, light, use = "complete.obs")))
    }
  }, by = day]
  per_day <- per_day[is.finite(env)]
  idx[[id]] <- data.frame(id = id, days = nrow(per_day),
    env_resid  = round(mean(per_day$env, na.rm = TRUE), 4),
    day_sd     = round(sd(per_day$mx, na.rm = TRUE) / rng, 4),
    rough      = round(mean(per_day$rr, na.rm = TRUE), 4),
    dark_frac  = round(mean(per_day$dk, na.rm = TRUE), 4),
    dive_atten = round(mean(per_day$da, na.rm = TRUE), 4), row.names = NULL)
}
I <- rbindlist(idx)
M <- merge(truth, I, by = "id")
M[, `:=`(offset = round(offset, 2), scatter = round(scatter, 2),
         rmse_lat = round(rmse_lat, 2), cover = round(cover, 2))]
cat("=== position-free noise indices against measured latitude scatter ===\n")
print(as.data.frame(M[order(-scatter)]), row.names = FALSE)

cat("\n=== correlations with the thing we want to predict (n = 10) ===\n")
for (v in c("env_resid", "day_sd", "rough", "dark_frac", "dive_atten")) {
  ct <- cor.test(M[[v]], M$scatter)
  cat(sprintf("  %-11s vs latitude scatter: r = %+.3f  (p = %.3f)\n",
              v, ct$estimate, ct$p.value))
}
cat("\n  ...and against coverage (which is what we actually want to fix):\n")
for (v in c("env_resid", "day_sd", "rough", "dark_frac", "dive_atten")) {
  ct <- cor.test(M[[v]], M$cover)
  cat(sprintf("  %-11s vs cover95_lat:      r = %+.3f  (p = %.3f)\n",
              v, ct$estimate, ct$p.value))
}
saveRDS(M, file.path(SCRATCH, "noiseindex.rds"))
