# Observed light minus theoretical light, at the KNOWN positions.
#
# This is the quantity the simulation is built to isolate, and it is worth
# looking at directly: it is the only place the model and the world can disagree
# once position is given. Everything else -- movement, sampler, domain, summary
# statistic -- operates downstream of it.
#
# The engine's likelihood treats the residual as a spike-and-slab: an asymmetric
# exponential about the expected curve (shading pulls DOWN, so the lower arm is
# slower) plus a uniform slab for junk. If that description is right, the
# residual should be a one-sided exponential whose scale does not depend on
# zenith. Anything systematic in the residual AS A FUNCTION OF ZENITH is a
# misspecification the fit cannot absorb, and near twilight it converts directly
# into latitude error.
#
# Four questions:
#   1. is the residual centred where the model thinks (i.e. is the curve right)?
#   2. does its shape change with zenith -- especially through twilight?
#   3. is the asymmetry the model assumes (down-shading) actually what happens?
#   4. does the twilight residual predict which tags are biased?
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))
CAL_DAYS <- 15

old <- lapply(readRDS("analysis/cache/nes/archives_v3_30min.rds"), `[[`, "main")
new <- lapply(readRDS("analysis/cache/nes/fvilches_archives_v1_30min.rds"), `[[`, "main")
new_argos <- readRDS("analysis/cache/nes/fvilches_argos_v1.rds")
new_man <- fread("scratch/nes_calibration/fvilches_manifest.csv",
                 colClasses = list(character = "topp"))
meta <- nes_meta(); ar <- nes_argos()
expected <- function(z, cal, maxl) pmin(pmax(cal[1] - cal[2]*z, 0), maxl)

tags <- list()
for (tg in names(old)) {
  m <- meta[match(tg, meta$id), ]
  if (is.na(m$ptt) || !nzchar(m$ptt)) next
  df <- ar$deploy[Ptt == m$ptt]; a <- ar$fixes[Ptt == m$ptt]
  if (!nrow(df) || nrow(a) < 50) next
  x <- as.data.table(old[[tg]])
  t0 <- max(min(x$time), as.POSIXct(max(m$deploy_date, df$user_date), tz = "UTC"))
  t1 <- min(max(x$time), max(a$time))
  d <- x[time >= t0 & time <= t1]; if (nrow(d) < 1000) next
  tags[[tg]] <- list(id = tg, light = d, argos = a, p0 = c(df$deploy_lon, df$deploy_lat))
}
for (tg in names(new)) {
  mm <- new_man[topp == tg]; if (!nrow(mm)) next
  a <- as.data.table(new_argos[[tg]]); d <- as.data.table(new[[tg]])
  d <- d[time >= min(a$time) & time <= max(a$time)]; if (nrow(d) < 1000) next
  tags[[tg]] <- list(id = tg, light = d, argos = a,
                     p0 = c(mm$deploy_lon[1], mm$deploy_lat[1]))
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

R <- rbindlist(lapply(names(tags), function(tg) {
  g <- tags[[tg]]; r <- responses[[tg]]; if (is.null(r)) return(NULL)
  a <- as.data.table(g$argos); setorder(a, time)
  lon <- approx(as.numeric(a$time), a$lon, as.numeric(g$light$time), rule = 2)$y
  lat <- approx(as.numeric(a$time), a$lat, as.numeric(g$light$time), rule = 2)$y
  z <- solar_zenith(as.numeric(g$light$time), lon, lat)
  obs <- pmax(0, g$light$light - r$baseline)
  mu <- expected(z, r$calibration, r$max_light)
  data.table(id = tg, z = z, obs = obs, mu = mu, res = obs - mu,
             rel = (obs - mu) / r$max_light, depth = g$light$depth_min)
}))
cat(sprintf("%d observations across %d deployments\n\n", nrow(R), uniqueN(R$id)))

cat("=== 1-2. the residual as a function of zenith ===\n")
R[, zb := cut(z, breaks = c(0, 60, 80, 86, 90, 94, 98, 102, 110, 180))]
print(as.data.frame(R[!is.na(zb), .(n = .N,
  mean_mu = round(mean(mu)), mean_obs = round(mean(obs)),
  mean_res = round(mean(res), 1), median_res = round(median(res), 1),
  sd_res = round(sd(res), 1),
  frac_above = round(mean(res > 0), 3)), by = zb][order(zb)]), row.names = FALSE)
cat("\nfrac_above is the share of observations BRIGHTER than the model expects.\n")
cat("The spike-and-slab assumes shading, so it expects this well below 0.5 and\n")
cat("roughly constant. A swing with zenith is a misspecified CURVE, not noise.\n")

cat("\n=== 3. is the asymmetry the model assumes what actually happens? ===\n")
tw <- R[z > 86 & z < 102]
cat(sprintf("through twilight (86-102 deg): mean %+.1f, median %+.1f, skew %.2f\n",
            mean(tw$res), median(tw$res),
            mean((tw$res - mean(tw$res))^3) / sd(tw$res)^3))
day <- R[z < 80]; night <- R[z > 110]
cat(sprintf("daytime  (<80):  mean %+.1f, frac brighter than expected %.3f\n",
            mean(day$res), mean(day$res > 0)))
cat(sprintf("night   (>110):  mean %+.1f, frac brighter than expected %.3f\n",
            mean(night$res), mean(night$res > 0)))

cat("\n=== 4. does the twilight residual predict which tags are biased? ===\n")
k <- fread("scratch/nes_calibration/all29_knots.csv")
b <- k[is.finite(err_lat), .(bias = mean(err_lat), scatter = sd(err_lat)), by = id]
b[, id := as.character(id)]
tt <- R[z > 86 & z < 102, .(tw_res = mean(res), tw_sd = sd(res),
                            tw_frac_above = mean(res > 0)), by = id]
M <- merge(b, tt, by = "id")
for (v in c("tw_res", "tw_sd", "tw_frac_above")) {
  c1 <- cor.test(M[[v]], M$bias); c2 <- cor.test(M[[v]], M$scatter)
  cat(sprintf("  %-13s vs bias r = %+.3f (p=%.3f)   vs scatter r = %+.3f (p=%.3f)\n",
              v, c1$estimate, c1$p.value, c2$estimate, c2$p.value))
}
cat("\nA twilight residual that predicts per-tag BIAS would be the first thing\n")
cat("found all session that does. Scatter has been predicted weakly before.\n")
saveRDS(R[sample(.N, min(.N, 2e5))], file.path(SCRATCH, "residuals.rds"))
