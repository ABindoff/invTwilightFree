# Can a Kuo-Mallick darkness indicator carry the day-length signal?
#
# The tangent gets its latitude information by accident: it expects ZERO at
# night, so every non-zero reading is charged, and the charge differs between a
# position that predicts night and one that predicts twilight. Give the response
# its true floor and that lever vanishes (451 km against 242).
#
# The principled replacement is an indicator: does this observation EXCEED what
# darkness can produce? Genuine dark readings run from 0 (diving) up to the dark
# envelope; anything above that is either not-dark or contaminated (ALAN, moon).
# So the discriminant is P(obs > dark envelope), and it is only useful if that
# probability differs sharply between night and twilight.
#
# Measure it, and measure the contamination the indicator would have to absorb.
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
  tab <- light_response_table(src, from = 30, to = 140, by = 1, max_light = r$max_light)
  data.table(id = tg, z = z, maxl = r$max_light,
             obs = pmax(0, g$light$light - r$baseline),
             dark_env = tail(tab, 1),          # the table's night value
             mu = invTwilightFree:::.tf_expected_light(z, tab, r$max_light))
}))
cat(sprintf("%d observations, %d deployments\n", nrow(R), uniqueN(R$id)))

cat("\n=== does 'exceeds the dark envelope' separate night from twilight? ===\n")
R[, band := cut(z, c(0,94,98,102,106,110,120,180))]
print(as.data.frame(R[!is.na(band), .(
  n = .N,
  mean_obs = round(mean(obs), 1),
  dark_env = round(mean(dark_env), 1),
  pct_over_env = round(100*mean(obs > dark_env), 1),
  pct_over_1.5x = round(100*mean(obs > 1.5*dark_env), 1),
  q99 = round(quantile(obs, 0.99), 1)), by = band][order(band)]), row.names = FALSE)
cat("\npct_over_env is the indicator. It is a usable discriminant only if it\n")
cat("falls sharply from twilight into full night.\n")

cat("\n=== how much contamination must the slab absorb? ===\n")
nt <- R[z > 110]
cat(sprintf("deep night (z>110), %d obs across %d tags:\n", nrow(nt), uniqueN(nt$id)))
for (m in c(1, 1.25, 1.5, 2)) {
  cat(sprintf("  obs > %.2f x dark envelope: %.2f%% of observations\n",
              m, 100*mean(nt$obs > m*nt$dark_env)))
}
per <- nt[, .(pct = 100*mean(obs > 1.5*dark_env)), by = id][order(-pct)]
cat(sprintf("\nper-tag rate of bright-when-dark (>1.5x env): median %.2f%%, max %.2f%% (%s)\n",
            median(per$pct), per$pct[1], per$id[1]))
cat("A slab weight has to cover the worst tag, not the median, or that tag breaks.\n")
cat("This is the ALAN / moonlight budget the indicator has to absorb.\n")
