# What shape IS the Mk9 light response, measured rather than assumed?
#
# The Mk9 is understood to carry an optical filter that maximises sensitivity
# around twilight, for template-fitting methods. Twilight light is blue-shifted
# (Chappuis absorption), so a blue-weighted channel reads relatively MORE light
# through the twilight band than a broadband sensor would. That is a genuine
# enhancement over a zenith band, and NO MONOTONE SIGMOID CAN REPRESENT IT --
# neither the clamped-linear tangent the engine is given nor the four-parameter
# logistic it is fitted from.
#
# The residual analysis is consistent: observed minus theoretical runs -22, -34,
# -14 through the day (shading, as assumed), then +3.7, +9.9, +3.4 across
# 86-98 degrees, then +12.9 at night. That is a BUMP through twilight plus a
# separate floor offset, not a shifted curve.
#
# Measure it. Pool all 29 deployments at their TRUE Argos positions, take the
# upper envelope per zenith bin (shading only ever subtracts, so the clear-sky
# response is the top of the cloud), and compare with the two parametric forms
# the package can currently express. The gap is the shape the model cannot fit.
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
  tags[[tg]] <- list(id = tg, light = d, argos = a, p0 = c(df$deploy_lon, df$deploy_lat),
                     fam = "Mk9_219")
}
for (tg in names(new)) {
  mm <- new_man[topp == tg]; if (!nrow(mm)) next
  a <- as.data.table(new_argos[[tg]]); d <- as.data.table(new[[tg]])
  d <- d[time >= min(a$time) & time <= max(a$time)]; if (nrow(d) < 1000) next
  ser <- sub("^[0-9]+_(.+)-out-Archive[.]csv$", "\\1",
             list.files("fvilches/extracted/TDR raw",
                        pattern = paste0("^", tg, "_.*Archive[.]csv$")))[1]
  tags[[tg]] <- list(id = tg, light = d, argos = a,
                     p0 = c(mm$deploy_lon[1], mm$deploy_lat[1]),
                     fam = if (grepl("^18A", ser)) "F18A" else "Mk9_219")
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
P <- attr(responses, "pooled")

# observed light against TRUE zenith, all deployments
R <- rbindlist(lapply(names(tags), function(tg) {
  g <- tags[[tg]]; r <- responses[[tg]]; if (is.null(r)) return(NULL)
  a <- as.data.table(g$argos); setorder(a, time)
  lon <- approx(as.numeric(a$time), a$lon, as.numeric(g$light$time), rule = 2)$y
  lat <- approx(as.numeric(a$time), a$lat, as.numeric(g$light$time), rule = 2)$y
  data.table(id = tg, fam = g$fam,
             z = solar_zenith(as.numeric(g$light$time), lon, lat),
             # on the scale the engine works in: baselined, so 0 = the 5th pct
             y = pmax(0, g$light$light - r$baseline),
             maxl = r$max_light,
             cal1 = r$calibration[1], cal2 = r$calibration[2],
             z50 = r$z50, sc = r$scale, floor = r$floor, amp = r$amp)
}))
cat(sprintf("%d observations, %d deployments\n", nrow(R), uniqueN(R$id)))

# upper envelope per 1-degree zenith bin, per family
R[, zb := round(z)]
env <- R[zb >= 40 & zb <= 130, .(n = .N,
  env95 = quantile(y, 0.95), env90 = quantile(y, 0.90),
  med = median(y),
  # what the engine is currently given: the clamped-linear tangent
  tangent = pmin(pmax(mean(cal1) - mean(cal2)*mean(zb), 0), mean(maxl)),
  # what the logistic form would give, on the same baselined scale
  logistic = mean(floor) + mean(amp)/(1 + exp((mean(zb) - mean(z50))/mean(sc))) -
             mean(floor)),
  by = .(fam, zb)][n >= 200]
env[, gap_tangent := round(env95 - tangent, 1)]
env[, gap_logistic := round(env95 - logistic, 1)]

cat("\n=== measured clear-sky response vs what the package can express (Mk9) ===\n")
e <- env[fam == "Mk9_219" & zb %in% seq(40, 130, by = 5)]
print(as.data.frame(e[, .(zenith = zb, n, measured = round(env95),
                          tangent = round(tangent), gap_tangent,
                          logistic = round(logistic), gap_logistic)]),
      row.names = FALSE)

cat("\n=== where is the gap, and is it a bump? ===\n")
tw <- env[fam == "Mk9_219" & zb >= 84 & zb <= 100]
dy <- env[fam == "Mk9_219" & zb >= 50 & zb <= 75]
nt <- env[fam == "Mk9_219" & zb >= 110]
cat(sprintf("daytime  (50-75):  mean gap to tangent %+.1f\n", mean(dy$gap_tangent)))
cat(sprintf("twilight (84-100): mean gap to tangent %+.1f   to logistic %+.1f\n",
            mean(tw$gap_tangent), mean(tw$gap_logistic)))
cat(sprintf("night    (>110):   mean gap to tangent %+.1f   to logistic %+.1f\n",
            mean(nt$gap_tangent), mean(nt$gap_logistic)))
cat(sprintf("\npeak twilight gap: %+.1f at zenith %d\n",
            max(tw$gap_tangent), tw$zb[which.max(tw$gap_tangent)]))
cat(sprintf("as a fraction of the response range (%.0f): %.1f%%\n",
            mean(R$maxl), 100*max(tw$gap_tangent)/mean(R$maxl)))

cat("\nA gap that is NEGATIVE in daylight (shading) and POSITIVE through twilight\n")
cat("is a bump the sigmoid cannot make. Both forms are monotone in zenith; a\n")
cat("filtered channel is not.\n")

cat("\n=== does the 18A family show the same shape? ===\n")
if (nrow(env[fam == "F18A"])) {
  tw2 <- env[fam == "F18A" & zb >= 84 & zb <= 100]
  cat(sprintf("18A twilight gap to tangent: %+.1f (Mk9: %+.1f)\n",
              mean(tw2$gap_tangent), mean(tw$gap_tangent)))
  cat("A different tag family sharing the bump argues for physics (twilight\n")
  cat("spectrum) over one instrument's filter.\n")
} else cat("no 18A bins with enough data\n")
fwrite(env, "scratch/nes_calibration/empirical_response.csv")
