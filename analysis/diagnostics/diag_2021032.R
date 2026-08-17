# Six mechanisms have now failed to explain 2021032: ALAN, response calibration
# (its pooled curve is the BEST-matched of the ten), the movement model,
# state-switching, the observation-noise scale, and the search domain.
#
# The track map suggests why they all missed: the failure is not diffuse. The
# estimate makes one large excursion to about 58 N while truth is near 47. If
# the error is concentrated in a few weeks rather than spread over 247 days,
# then it is an EPISODE, and no global parameter could ever have fixed it.
#
# So: locate it in time, then ask what the light was doing while it happened.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))
FOCUS <- "2021032"; CLEAN <- "2021028"

k <- as.data.table(read.csv("analysis/output/nes/per_knot_errors.csv"))
k <- k[engine == "hier" & batch == "light" & is.finite(err_lat)]
k[, `:=`(id = as.character(id), time = as.POSIXct(time, tz = "UTC"))]

cat("=== is the error concentrated in time? ===\n")
conc <- k[, {
  s <- sort(err_lat^2, decreasing = TRUE)
  .(knots = .N, rmse = round(sqrt(mean(err_lat^2)), 2),
    pct_of_error_in_worst_10pct = round(100 * sum(s[1:ceiling(.N*0.1)]) / sum(s)),
    pct_in_worst_25pct = round(100 * sum(s[1:ceiling(.N*0.25)]) / sum(s)))
}, by = id]
print(as.data.frame(conc[order(-pct_of_error_in_worst_10pct)]), row.names = FALSE)
cat("\nIf a tag carries most of its squared error in a tenth of its knots, the\n")
cat("problem is an episode; if it is spread evenly, it is a property.\n")

f <- k[id == FOCUS][order(time)]
f[, mon := format(time, "%Y-%m")]
cat(sprintf("\n=== %s, latitude error by month ===\n", FOCUS))
print(as.data.frame(f[, .(knots = .N, mean_true_lat = round(mean(true_lat), 1),
  bias = round(mean(err_lat), 2), rmse = round(sqrt(mean(err_lat^2)), 2),
  worst = round(max(abs(err_lat)), 1)), by = mon]), row.names = FALSE)

# what was the light doing during the worst window?
arch <- lapply(readRDS("analysis/cache/nes/archives_v3_30min.rds"), `[[`, "main")
meta <- nes_meta(); ar <- nes_argos()
m <- meta[match(FOCUS, meta$id), ]
a <- ar$fixes[Ptt == m$ptt]
x <- as.data.table(arch[[FOCUS]])
tr <- argos_at(a, x$time)
x[, `:=`(tlon = tr$lon, tlat = tr$lat)]
x <- x[is.finite(tlon)]
x[, zen := solar_zenith(as.numeric(time), tlon, tlat)]
x[, mon := format(time, "%Y-%m")]

bad_mon <- f[, .(r = sqrt(mean(err_lat^2))), by = mon][order(-r)][1:2]$mon
cat(sprintf("\n=== light behaviour by month (worst months: %s) ===\n",
            paste(bad_mon, collapse = ", ")))
rng <- diff(quantile(x$light, c(.05, .95)))
print(as.data.frame(x[, .(obs = .N,
  day_max   = round(quantile(light, 0.99)),
  night_med = round(median(light[zen > 108])),
  twilight_n = sum(zen > 88 & zen < 100),
  # how cleanly does the light track zenith through twilight?
  tw_spread = round(sd(light[zen > 88 & zen < 100]) / rng, 3),
  med_depth = round(median(depth_max)),
  frac_deep = round(mean(depth_max > 400), 2),
  worst = mon %in% bad_mon), by = mon][order(mon)]), row.names = FALSE)

cat("\n=== the same, for a clean tag, as a control ===\n")
m2 <- meta[match(CLEAN, meta$id), ]; a2 <- ar$fixes[Ptt == m2$ptt]
x2 <- as.data.table(arch[[CLEAN]]); tr2 <- argos_at(a2, x2$time)
x2[, `:=`(tlon = tr2$lon, tlat = tr2$lat)]; x2 <- x2[is.finite(tlon)]
x2[, zen := solar_zenith(as.numeric(time), tlon, tlat)][, mon := format(time, "%Y-%m")]
rng2 <- diff(quantile(x2$light, c(.05, .95)))
print(as.data.frame(x2[, .(obs = .N, day_max = round(quantile(light, 0.99)),
  night_med = round(median(light[zen > 108])),
  tw_spread = round(sd(light[zen > 88 & zen < 100]) / rng2, 3),
  med_depth = round(median(depth_max)),
  frac_deep = round(mean(depth_max > 400), 2)), by = mon][order(mon)]), row.names = FALSE)
