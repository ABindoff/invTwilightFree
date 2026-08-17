# The declination mechanism, tested properly this time.
#
# Earlier I proposed that a fixed error in the assumed threshold zenith produces
# opposite-signed latitude bias either side of the equinox, and then "falsified"
# it with a test pooled across all ten animals (r = -0.001). That test was
# mis-specified: each tag has its OWN width error, so its own swing magnitude and
# SIGN, and pooling averages them to nothing.
#
# 2021032 shows it unmistakably on its own: bias -3.8, -4.8, -3.2 through
# July-September, then +8.5, +10.8, +8.2 through October-December.
#
# The prediction, per tag. The engine is handed the pooled transition width;
# through the tangent construction, the assumed zero-crossing is
#     zero_at = z50 + 2 * scale,  scale = width / 4.394
# so a pooled width WIDER than the tag's true width puts the assumed threshold
# too deep. Too deep means the model expects light to persist past the true
# horizon, which reads as too much daylight in summer (biased south) and too
# little in winter (biased north).
#
# So: swing = (winter bias) - (summer bias) should be POSITIVE when the pooled
# zero_at exceeds the tag's true zero_at, and should scale with the size of that
# error. That is a signed, quantitative prediction across ten animals.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))

declination <- function(t) {
  doy <- as.numeric(format(t, "%j"))
  hr  <- as.numeric(format(t, "%H")) + as.numeric(format(t, "%M")) / 60
  g <- 2*pi/365 * (doy - 1 + (hr - 12)/24)
  (0.006918 - 0.399912*cos(g) + 0.070257*sin(g) - 0.006758*cos(2*g) +
     0.000907*sin(2*g) - 0.002697*cos(3*g) + 0.00148*sin(3*g)) * 180/pi
}

k <- as.data.table(read.csv("analysis/output/nes/per_knot_errors.csv"))
k <- k[engine == "hier" & batch == "light" & is.finite(err_lat)]
k[, `:=`(id = as.character(id), time = as.POSIXct(time, tz = "UTC"))]
k[, decl := declination(time)]

sw <- k[, .(n_sum = sum(decl > 5), n_win = sum(decl < -5),
            bias_summer = round(mean(err_lat[decl > 5]), 2),
            bias_winter = round(mean(err_lat[decl < -5]), 2)), by = id]
sw[, swing := round(bias_winter - bias_summer, 2)]

# the tag's true at-sea response, and the zero_at the engine was actually given
arch <- lapply(readRDS("analysis/cache/nes/archives_v3_30min.rds"), `[[`, "main")
meta <- nes_meta(); ar <- nes_argos()
rows <- list()
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
  dep <- departure_time(d)
  at_sea <- if (is.finite(dep)) d[as.numeric(time) >= dep] else d
  tr <- argos_at(a, at_sea$time); ok <- is.finite(tr$lon)
  r <- if (sum(ok) >= 200)
    fit_light_response(at_sea$time[ok], at_sea$light[ok], tr$lon[ok], tr$lat[ok]) else NULL
  if (is.null(r)) next
  rows[[length(rows)+1]] <- data.frame(id = id,
    true_z50 = r$z50, true_scale = r$scale,
    true_zero = r$z50 + 2 * r$scale, row.names = NULL)
}
G <- rbindlist(rows)
Z_POOL <- 92.09; S_POOL <- 27.0 / 4.394
pooled_zero <- Z_POOL + 2 * S_POOL
G[, zero_err := round(pooled_zero - true_zero, 2)]   # positive = assumed too deep

M <- merge(sw, G, by = "id")
cat(sprintf("pooled zero-crossing handed to the engine: %.2f deg\n\n", pooled_zero))
cat("zero_err > 0 means the assumed threshold is DEEPER than the tag's true one,\n")
cat("which predicts summer bias south, winter bias north, so swing > 0.\n\n")
print(as.data.frame(M[order(-zero_err),
  .(id, true_zero = round(true_zero, 2), zero_err, bias_summer, bias_winter, swing)]),
  row.names = FALSE)

ct <- cor.test(M$zero_err, M$swing)
cat(sprintf("\ncorrelation of seasonal swing with the zero-crossing error: r = %+.3f (p = %.4f)\n",
            ct$estimate, ct$p.value))
f <- lm(swing ~ zero_err, data = M)
cat(sprintf("slope %.2f deg of latitude swing per deg of zenith error (R2 = %.2f)\n",
            coef(f)[2], summary(f)$r.squared))
cat(sprintf("\ntags with swing of the predicted sign: %d of %d\n",
            sum(sign(M$swing) == sign(M$zero_err)), nrow(M)))
cat(sprintf("mean |swing|: %.2f deg -- for comparison, mean per-tag offset is %.2f\n",
            mean(abs(M$swing)), mean(abs(k[, mean(err_lat), by = id]$V1))))
