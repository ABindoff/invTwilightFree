# Is 2021032 contaminated by artificial light at night?
#
# That animal has resisted everything: coverage 0.13 at the default and 0.27-0.36
# across the lambda sweep, unmoved by recalibration, by the movement model, by
# state-switching. Four mechanisms have missed it, which suggests the cause is
# not any of the things being varied.
#
# ALAN would do exactly this. Light where the model expects darkness corrupts the
# twilight edge, and latitude is read from that edge. It would be sporadic
# (whenever the animal is near a vessel), so it would inflate SCATTER rather than
# impose a constant offset -- which matches: 2021032's error is scatter 6.16
# against an offset of 2.53.
#
# And the location is right. These animals forage in the North Pacific Transition
# Zone, which carries one of the world's largest squid-jigging fleets; those
# vessels run banks of metal-halide lamps bright enough to be a standard feature
# of night-time satellite imagery.
#
# The test: solar zenith at the KNOWN Argos position, so "night" is defined by
# geometry rather than by the model. Astronomical night is zenith > 108 degrees.
# Any light above the tag's own dark level there has no natural explanation
# beyond moonlight and bioluminescence.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))
CAL_DAYS <- 15

arch <- lapply(readRDS("analysis/cache/nes/archives_v3_30min.rds"), `[[`, "main")
meta <- nes_meta(); ar <- nes_argos()

dat <- list()
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
  tr <- argos_at(a, d$time)
  d[, `:=`(tlon = tr$lon, tlat = tr$lat)]
  d <- d[is.finite(tlon)]
  d[, zen := solar_zenith(as.numeric(time), tlon, tlat)]
  k <- d$time <= min(d$time) + CAL_DAYS * 86400
  r <- fit_light_response(d$time[k], d$light[k], df$deploy_lon, df$deploy_lat)
  d[, dark := if (is.null(r)) as.numeric(quantile(light, .05)) else r$dark_level]
  d[, id := id]
  dat[[id]] <- d
}
D <- rbindlist(dat, fill = TRUE)

# Night light, defined geometrically and measured against the tag's own dark level
night <- D[zen > 108]
tab <- night[, .(night_obs = .N,
  median_night = round(median(light)),
  dark_level = round(median(dark)),
  excess_med = round(median(light - dark), 1),
  frac_above_10 = round(mean(light - dark > 10), 4),
  frac_above_30 = round(mean(light - dark > 30), 4),
  max_excess = round(max(light - dark))), by = id]

k <- as.data.table(read.csv("analysis/output/nes/per_knot_errors.csv"))
k <- k[engine == "hier" & batch == "light" & is.finite(err_lat)]
k[, id := as.character(id)]
sc <- k[, .(scatter = round(sd(err_lat), 2),
            cover = round(mean(abs(err_lat) <= 1.96 * lat_sd, na.rm = TRUE), 2)), by = id]
M <- merge(tab, sc, by = "id")
cat("=== light during astronomical night (solar zenith > 108 deg, at TRUE positions) ===\n")
print(as.data.frame(M[order(-frac_above_10)]), row.names = FALSE)
cat(sprintf("\ncorrelation of frac_above_10 with latitude scatter: %+.3f (p = %.3f)\n",
            cor(M$frac_above_10, M$scatter),
            cor.test(M$frac_above_10, M$scatter)$p.value))
cat(sprintf("correlation of frac_above_10 with coverage:         %+.3f (p = %.3f)\n",
            cor(M$frac_above_10, M$cover),
            cor.test(M$frac_above_10, M$cover)$p.value))

# ---- figure -----------------------------------------------------------------
FOCUS <- "2021032"; CLEAN <- "2021028"
png(file.path("scratch/nes_calibration", "alan_2021032.png"),
    width = 1500, height = 1150, res = 110)
layout(matrix(1:4, 2, 2, byrow = TRUE))
op <- par(mar = c(4.2, 4.2, 2.6, 1))

for (tg in c(FOCUS, CLEAN)) {
  d <- D[id == tg]
  plot(d$zen, d$light, pch = 16, cex = 0.18,
       col = grDevices::adjustcolor(ifelse(d$depth_max > 50, "#2980b9", "#e67e22"), 0.25),
       xlab = "solar zenith at the true position (deg)", ylab = "light (tag units)",
       main = sprintf("%s%s: light against zenith", tg,
                      if (tg == FOCUS) "  (the problem tag)" else "  (a clean tag)"))
  abline(v = c(90, 108), lty = c(2, 3), col = c("grey40", "red"))
  abline(h = median(d$dark), col = "darkgreen", lwd = 1.5)
  legend("topright", c("shallow (<50 m)", "diving", "dark level", "horizon", "astro night"),
         col = c("#e67e22", "#2980b9", "darkgreen", "grey40", "red"),
         pch = c(16, 16, NA, NA, NA), lty = c(NA, NA, 1, 2, 3), bty = "n", cex = 0.75)
}
# night-light excess through the deployment
for (tg in c(FOCUS, CLEAN)) {
  d <- D[id == tg & zen > 108]
  plot(d$time, d$light - d$dark, pch = 16, cex = 0.3,
       col = grDevices::adjustcolor("#8e44ad", 0.45),
       xlab = "", ylab = "night light above dark level",
       main = sprintf("%s: night-time excess through the trip", tg))
  abline(h = 0, col = "darkgreen")
  abline(h = 10, lty = 3, col = "red")
}
par(op); dev.off()
cat("\nwrote scratch/nes_calibration/alan_2021032.png\n")
