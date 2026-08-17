# How much light is being penalised that should not be penalised at all?
#
# The spike began as a Kuo-Mallick-style penalty on light appearing where the
# model expects darkness. That is a well-posed thing to penalise: a genuinely
# anomalous bright reading at night. The current likelihood instead applies the
# upper arm to EVERY observation above `expected`, at rate lam_hi = 2*lam_lo,
# wherever it occurs.
#
# The trouble is that `expected` is the clamped-linear tangent, which is zero
# above ~106 degrees. The tag is not zero there: its dark reading is 39-44 on
# the baselined scale. So every night observation is "light where darkness was
# expected" and is charged at the harsh rate, when nothing anomalous has
# happened at all. It is the sensor's own floor.
#
# Quantify it. For each calibration, per zenith band:
#   - what share of observations land on the UPPER arm
#   - the mean and total penalty in nats
#   - what share of the whole penalty budget that band consumes
# Then the spurious part: the penalty the tangent charges at night MINUS what a
# correctly-floored response charges for the same observations.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))
CAL_DAYS <- 15

old <- lapply(readRDS("analysis/cache/nes/archives_v3_30min.rds"), `[[`, "main")
new <- lapply(readRDS("analysis/cache/nes/fvilches_archives_v1_30min.rds"), `[[`, "main")
new_argos <- readRDS("analysis/cache/nes/fvilches_argos_v1.rds")
new_man <- fread("scratch/nes_calibration/fvilches_manifest.csv",
                 colClasses = list(character = "topp"))
meta <- nes_meta(); ar <- nes_argos()

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
src <- Filter(Negate(is.null), geom_fits)

R <- rbindlist(lapply(names(tags), function(tg) {
  g <- tags[[tg]]; r <- responses[[tg]]; if (is.null(r)) return(NULL)
  a <- as.data.table(g$argos); setorder(a, time)
  lon <- approx(as.numeric(a$time), a$lon, as.numeric(g$light$time), rule = 2)$y
  lat <- approx(as.numeric(a$time), a$lat, as.numeric(g$light$time), rule = 2)$y
  z <- solar_zenith(as.numeric(g$light$time), lon, lat)
  obs <- pmax(0, g$light$light - r$baseline)
  tab <- light_response_table(src, from = 30, to = 140, by = 1, max_light = r$max_light)
  data.table(id = tg, z = z, obs = obs, maxl = r$max_light,
             mu_tan = invTwilightFree:::.tf_expected_light(z, r$calibration, r$max_light),
             mu_tab = invTwilightFree:::.tf_expected_light(z, tab, r$max_light))
}))
cat(sprintf("%d observations, %d deployments\n", nrow(R), uniqueN(R$id)))

# penalty in nats on the upper arm: lam_hi * (obs - mu), lam_hi = 2*lam_lo,
# lam_lo = 1/(0.5*max_light), so lam_hi = 4/max_light
R[, lam_hi := 4 / maxl]
R[, pen_tan := pmax(obs - mu_tan, 0) * lam_hi]
R[, pen_tab := pmax(obs - mu_tab, 0) * lam_hi]
R[, band := cut(z, c(0, 80, 86, 94, 100, 106, 180),
                labels = c("day <80","80-86","86-94","94-100","100-106","night >106"))]

cat("\n=== the TANGENT's penalty budget, by zenith band ===\n")
tot <- sum(R$pen_tan)
print(as.data.frame(R[!is.na(band), .(
  n = .N, pct_obs = round(100*.N/nrow(R), 1),
  mean_mu = round(mean(mu_tan), 1), mean_obs = round(mean(obs), 1),
  pct_upper_arm = round(100*mean(obs > mu_tan), 1),
  mean_nats = round(mean(pen_tan), 3),
  pct_of_budget = round(100*sum(pen_tan)/tot, 1)), by = band][order(band)]),
  row.names = FALSE)

cat("\n=== the same observations under a correctly-floored (table) response ===\n")
tot2 <- sum(R$pen_tab)
print(as.data.frame(R[!is.na(band), .(
  n = .N, mean_mu = round(mean(mu_tab), 1),
  pct_upper_arm = round(100*mean(obs > mu_tab), 1),
  mean_nats = round(mean(pen_tab), 3),
  pct_of_budget = round(100*sum(pen_tab)/tot2, 1)), by = band][order(band)]),
  row.names = FALSE)

cat(sprintf("\ntotal penalty: tangent %.0f nats, table %.0f nats (%.0f%% less)\n",
            tot, tot2, 100*(1 - tot2/tot)))
nt <- R[z > 106]
cat(sprintf("\nNIGHT ONLY (%.0f%% of all observations):\n", 100*nrow(nt)/nrow(R)))
cat(sprintf("  tangent expects %.1f, tag reads %.1f, charges %.3f nats each\n",
            mean(nt$mu_tan), mean(nt$obs), mean(nt$pen_tan)))
cat(sprintf("  table   expects %.1f, same obs,      charges %.3f nats each\n",
            mean(nt$mu_tab), mean(nt$pen_tab)))
cat(sprintf("  SPURIOUS night penalty = %.0f%% of the tangent's ENTIRE budget\n",
            100*(sum(nt$pen_tan) - sum(nt$pen_tab))/tot))
cat(sprintf("\nper 12-hour knot (%.1f observations of which %.1f at night):\n",
            24, 24*nrow(nt)/nrow(R)))
cat(sprintf("  tangent charges %.2f nats per knot, %.2f of it at night\n",
            24*mean(R$pen_tan), 24*(nrow(nt)/nrow(R))*mean(nt$pen_tan)))
cat("\nA Kuo-Mallick spike is meant to fire on ANOMALOUS light. Whatever share of\n")
cat("the budget is night dark-current is the spike firing on nothing.\n")
