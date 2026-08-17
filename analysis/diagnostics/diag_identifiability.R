# Do the shading arm and the false-light slab absorb each other?
#
# The shading spike explains observations BELOW the expected curve; the slab
# explains observations ABOVE it. Between them they can explain anything. Make
# both generous enough and the likelihood goes flat: position stops being
# identified, the posterior widens toward the prior, and nothing looks broken,
# because nothing IS broken -- the model has simply stopped saying anything.
#
# That corner is not far away. The shading rate is already 3.1x looser than the
# measured residual scale, so the lower arm is over-generous before the night
# slab is touched. Map it.
#
# TWO QUANTITIES, both per 12-hour knot, both in nats:
#
#   SIGNAL     logL(true position) - logL(true position + 3 degrees latitude).
#              This IS the day-length information, measured. It is what the
#              clamped-linear response currently obtains by accident, through a
#              floor that is wrong by 40 units.
#
#   ROBUSTNESS the cost of one injected false-light reading at 0.9*max_light in
#              a night window. Low is good: the model absorbs contamination
#              instead of moving the track to accommodate it.
#
# The two trade off. The useful output is the frontier, and the region where
# signal collapses -- which is the failure mode, drawn.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))
set.seed(5)
CAL_DAYS <- 15; STEP_HOURS <- 12; DLAT <- 3

old <- lapply(readRDS("analysis/cache/nes/archives_v3_30min.rds"), `[[`, "main")
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
tags <- tags[seq_len(min(6, length(tags)))]
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

# per-tag windows at the engine's own knot spacing, positions from Argos
W <- list()
for (tg in names(tags)) {
  g <- tags[[tg]]; r <- responses[[tg]]; if (is.null(r)) next
  a <- as.data.table(g$argos); setorder(a, time)
  tt <- as.numeric(g$light$time)
  obs <- pmax(0, g$light$light - r$baseline)
  brk <- seq(min(tt), max(tt), by = STEP_HOURS * 3600)
  k <- findInterval(tt, brk)
  lonk <- approx(as.numeric(a$time), a$lon, brk, rule = 2)$y
  latk <- approx(as.numeric(a$time), a$lat, brk, rule = 2)$y
  tab <- light_response_table(src, from = 30, to = 140, by = 1, max_light = r$max_light)
  keep <- which(tabulate(k, length(brk)) >= 12)
  keep <- keep[seq(1, length(keep), length.out = min(60, length(keep)))]
  W[[tg]] <- list(tt = tt, obs = obs, k = k, brk = brk, lon = lonk, lat = latk,
                  keep = keep, tab = tab, cal = r$calibration, ml = r$max_light)
}
cat(sprintf("%d deployments, %d windows\n\n", length(W), sum(vapply(W, function(w) length(w$keep), 0L))))

score <- function(cal, lp) {
  sig <- rob <- numeric(0)
  for (w in W) {
    for (i in w$keep) {
      j <- which(w$k == i); if (length(j) < 12) next
      v <- eval_logpk_grid(c(w$lon[i], w$lon[i]), c(w$lat[i], w$lat[i] + DLAT),
                           w$tt[j], w$obs[j], cal, lp, 2)
      sig <- c(sig, v[1] - v[2])
      # inject one false-light reading into the darkest observation of the window
      y <- w$obs[j]; y[which.min(y)] <- 0.9 * w$ml
      v2 <- eval_logpk_grid(w$lon[i], w$lat[i], w$tt[j], y, cal, lp, 2)
      rob <- c(rob, v[1] - v2)
    }
  }
  # Mean signal is the wrong summary on its own. A response that makes larger
  # claims AND is wrong more confidently can have a bigger mean and be worse to
  # fit with. What the HMM actually accumulates is signal against its own
  # spread, so report that too.
  c(signal = mean(sig, na.rm = TRUE), robustness = mean(rob, na.rm = TRUE),
    frac_neg = mean(sig < 0, na.rm = TRUE),
    snr = mean(sig, na.rm = TRUE) / stats::sd(sig, na.rm = TRUE))
}

w1 <- W[[1]]
cat("=== reference: the clamped-linear tangent, whose signal is accidental ===\n")
ref <- score(w1$cal, c(1/(w1$ml*0.5), w1$ml, 0.10))
cat(sprintf("signal %.2f nats/knot   false-light cost %.2f nats   sign wrong on %.0f%% of windows\n\n",
            ref["signal"], ref["robustness"], 100*ref["frac_neg"]))

# CONTROL. The map is only worth reading if this metric orders the arms the way
# the real fits did. Two are already scored against Argos on all 29:
#   tangent     242 km  (no darkness regime)
#   table_flat  451 km  (no darkness regime, corrected floor)
# If snr ranks those the same way round, it is measuring something real, and its
# verdict on the darkness regime is a prediction rather than a curiosity.
cat("=== control: the two arms already scored against Argos ===\n")
ctl <- score(w1$tab, c(1/(w1$ml*0.5), w1$ml, 0.10))
cat(sprintf("tangent    (242 km on 29 tags): snr %.3f, signal %.2f, wrong sign %.0f%%\n",
            ref["snr"], ref["signal"], 100*ref["frac_neg"]))
cat(sprintf("table_flat (451 km on 29 tags): snr %.3f, signal %.2f, wrong sign %.0f%%\n\n",
            ctl["snr"], ctl["signal"], 100*ctl["frac_neg"]))

cat("=== the (shading rate, night slab) plane, on the corrected-floor table ===\n")
cat("lam_mult 1 is the current setting; 3.1 matches the measured residual scale.\n\n")
G <- expand.grid(lam_mult = c(0.5, 1, 2, 3, 4), pslab = c(0.05, 0.15, 0.35, 0.60))
res <- rbindlist(lapply(seq_len(nrow(G)), function(i) {
  lp <- c(G$lam_mult[i]/(w1$ml*0.5), w1$ml, 0.10, 0.35, 20, G$pslab[i])
  s <- score(w1$tab, lp)
  data.table(lam_mult = G$lam_mult[i], pslab = G$pslab[i],
             signal = round(s["signal"], 2), robust = round(s["robustness"], 2),
             pct_wrong_sign = round(100*s["frac_neg"]),
             snr = round(s["snr"], 3))
}))
print(as.data.frame(dcast(res, lam_mult ~ pslab, value.var = "signal")), row.names = FALSE)
cat("  ^ SIGNAL (nats per knot for a 3 degree latitude displacement)\n\n")
print(as.data.frame(dcast(res, lam_mult ~ pslab, value.var = "robust")), row.names = FALSE)
cat("  ^ COST OF ONE FALSE-LIGHT READING (nats; lower = better absorbed)\n\n")
print(as.data.frame(dcast(res, lam_mult ~ pslab, value.var = "pct_wrong_sign")), row.names = FALSE)
cat("  ^ % OF WINDOWS PREFERRING THE WRONG POSITION\n\n")
print(as.data.frame(dcast(res, lam_mult ~ pslab, value.var = "snr")), row.names = FALSE)
cat("  ^ SIGNAL / ITS OWN SD ACROSS WINDOWS -- what the HMM actually accumulates\n")

cat(sprintf("\ntangent reference: signal %.2f, false-light cost %.2f, snr %.3f, wrong sign %.0f%%\n",
            ref["signal"], ref["robustness"], ref["snr"], 100*ref["frac_neg"]))
cat("\nThe failure mode is signal -> 0 with the slab large and the shading arm\n")
cat("loose. If signal holds up across the row, they are not eating each other.\n")
fwrite(res, "scratch/nes_calibration/identifiability_map.csv")
