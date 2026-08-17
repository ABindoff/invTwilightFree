# Does the LIGHT RESPONSE drift over a deployment?
#
# The model fits one response for 240 days. Everything else has been excluded:
# the movement model, the samplers, the search domain, the observation-noise
# scale, ALAN, and now depth attenuation (96% of windows reach within 10 m of
# the surface, so the decimated maximum IS surface light and there is nothing
# left to attenuate).
#
# What remains has to vary on a multi-week timescale, differ between animals,
# and change sign. A drifting response does all three: a static curve fitted
# once is right at one end of the trip and wrong at the other, and which end it
# is right at depends on where in the record the calibration window sat.
#
# The hint: 2021032's twilight scatter rose from 0.161 (Jul-Sep) to 0.224
# (Oct-Dec) while its latitude bias flipped from -3.8 to +8.0.
#
# TEST. Fit the response separately in each third of the deployment, at the TRUE
# Argos positions so the fit is honest, and ask two things:
#   1. does the response actually move over the deployment?
#   2. does the movement predict the bias change over the same period?
# The second matters more. A response that drifts but does not track the bias is
# a curiosity; one that tracks it is the mechanism.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))

arch <- lapply(readRDS("analysis/cache/nes/archives_v3_30min.rds"), `[[`, "main")
meta <- nes_meta(); ar <- nes_argos()
k <- as.data.table(read.csv("analysis/output/nes/per_knot_errors.csv"))
k <- k[engine == "hier" & batch == "light" & is.finite(err_lat)]
k[, `:=`(id = as.character(id), time = as.POSIXct(time, tz = "UTC"))]

rows <- list()
for (tg in names(arch)) {
  m <- meta[match(tg, meta$id), ]
  if (is.na(m$ptt) || !nzchar(m$ptt)) next
  df <- ar$deploy[Ptt == m$ptt]; a <- ar$fixes[Ptt == m$ptt]
  if (!nrow(df) || nrow(a) < 50) next
  x <- as.data.table(arch[[tg]])
  t0 <- max(min(x$time), as.POSIXct(max(m$deploy_date, df$user_date), tz = "UTC"))
  t1 <- min(max(x$time), max(a$time))
  d <- x[time >= t0 & time <= t1]
  if (nrow(d) < 1000) next
  dep <- departure_time(d)
  d <- if (is.finite(dep)) d[as.numeric(time) >= dep] else d      # at sea only
  tr <- argos_at(a, d$time); ok <- is.finite(tr$lon)
  d <- d[ok]; tlon <- tr$lon[ok]; tlat <- tr$lat[ok]
  n <- nrow(d); if (n < 900) next
  cuts <- floor(seq(1, n + 1, length.out = 4))
  kk <- k[id == tg]   # `id` alone would resolve to the COLUMN, not the loop variable
  for (p in 1:3) {
    idx <- cuts[p]:(cuts[p+1] - 1)
    r <- fit_light_response(d$time[idx], d$light[idx], tlon[idx], tlat[idx])
    if (is.null(r)) next
    win <- range(d$time[idx])
    b <- kk[time >= win[1] & time <= win[2], mean(err_lat)]
    rows[[length(rows)+1]] <- data.frame(id = tg, third = p,
      z50 = round(r$z50, 2), width = round(r$width_deg, 1),
      zero_at = round(r$zero_at, 2),
      bias = if (length(b) && is.finite(b)) round(b, 2) else NA_real_,
      row.names = NULL)
  }
}
R <- rbindlist(rows)
cat("=== response fitted separately in each third, at TRUE positions ===\n")
print(as.data.frame(R), row.names = FALSE)

cat("\n=== 1. does the response move over a deployment? ===\n")
w <- dcast(R, id ~ third, value.var = "zero_at")
setnames(w, c("id", "t1", "t2", "t3"))
w[, drift := round(t3 - t1, 2)]
print(as.data.frame(w), row.names = FALSE)
cat(sprintf("\nmean |drift| in the zero-crossing over a deployment: %.2f deg\n",
            mean(abs(w$drift), na.rm = TRUE)))
cat(sprintf("for comparison, the whole between-animal spread is 1.85 deg\n"))

cat("\n=== 2. does the drift track the bias? (within tag, third to third) ===\n")
R[, `:=`(zero_c = zero_at - mean(zero_at), bias_c = bias - mean(bias)), by = id]
ok <- R[is.finite(bias_c) & is.finite(zero_c)]
if (nrow(ok) > 6) {
  ct <- cor.test(ok$zero_c, ok$bias_c)
  cat(sprintf("within-tag correlation of zero-crossing with latitude bias: %+.3f (p = %.3f, n = %d)\n",
              ct$estimate, ct$p.value, nrow(ok)))
  f <- lm(bias_c ~ zero_c, data = ok)
  cat(sprintf("slope %.2f deg of latitude per deg of zero-crossing (R2 = %.2f)\n",
              coef(f)[2], summary(f)$r.squared))
  cat("\nA deeper assumed threshold than the truth shortens the apparent day, so a\n")
  cat("POSITIVE slope means the drift is pushing latitude the way the physics says.\n")
} else cat("too few complete thirds\n")
