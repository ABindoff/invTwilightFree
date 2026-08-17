# Is the latitude likelihood ASYMMETRIC, and does lambda set the asymmetry?
#
# THE PICTURE UNDER TEST (user's, and it makes a falsifiable prediction). Longitude
# is well identified (rmse ~1.15 deg, coverage 0.96-1.00), so treat the problem as
# one-dimensional along a meridian, where latitude IS day length. In northern summer:
#
#   POLEWARD of truth  -> model predicts a LONGER day -> more light than observed
#                      -> observations sit BELOW the expected curve
#                      -> absorbed by the SHADING arm, rate lam_lo. CHEAP.
#   EQUATORWARD        -> model predicts a SHORTER day -> less light than observed
#                      -> observations sit ABOVE the curve
#                      -> upper arm at shade_ratio * lam_lo (2x), or the slab. DEAR.
#
# So the profile should have a long cheap poleward tail and a short expensive
# equatorward one, and since the per-observation cost scales with lambda, the
# asymmetry in nats should GROW with lambda -- pulling the posterior mean poleward.
#
# PREDICTIONS, stated before running:
#   P1  logL(truth + d) > logL(truth - d): cheaper to be wrong poleward.
#   P2  the gap grows with lambda, roughly proportionally.
#   P3  the profile is right-skewed, so mean sits POLEWARD of mode, and that gap
#       also grows with lambda.
#   P4  at very loose lambda all of this collapses: the spike flattens, the slab
#       dominates, latitude stops being identified and the profile goes flat. This
#       is the second regime that makes the bias-vs-lambda curve non-monotone.
#
# Sign convention: everything is expressed as POLEWARD-POSITIVE, so the southern
# hemisphere would flip. All 29 deployments here are northern, so poleward = north.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))
set.seed(7)
CAL_DAYS <- 15; STEP_HOURS <- 12
LATS <- seq(20, 70, by = 0.25)
LAMS <- 2^seq(-3, 4, by = 1)
SHADE_RATIO <- 2
DPROBE <- 3          # degrees either side of truth for the P1/P2 gap

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
resp <- lapply(tags, function(g) {
  k <- g$light$time <= min(g$light$time) + CAL_DAYS * 86400
  fit_light_response(g$light$time[k], g$light$light[k], g$p0[1], g$p0[2])
})
cat(sprintf("%d deployments\n\n", length(tags)))

R <- list()
for (tg in names(tags)) {
  g <- tags[[tg]]; r <- resp[[tg]]; if (is.null(r)) next
  a <- as.data.table(g$argos); setorder(a, time)
  tt <- as.numeric(g$light$time); obs <- pmax(0, g$light$light - r$baseline)
  brk <- seq(min(tt), max(tt), by = STEP_HOURS * 3600)
  kk <- findInterval(tt, brk)
  tlon <- approx(as.numeric(a$time), a$lon, brk, rule = 2)$y
  tlat <- approx(as.numeric(a$time), a$lat, brk, rule = 2)$y
  lam0 <- 1 / (r$max_light * 0.5)
  use <- which(tabulate(kk, length(brk)) >= 16)
  if (!length(use)) next
  use <- use[seq(1, length(use), length.out = min(20, length(use)))]
  for (i in use) {
    j <- which(kk == i); if (length(j) < 16) next
    for (lm in LAMS) {
      lp <- c(lam0 * lm, r$max_light, 0.10)
      ll <- eval_logpk_grid(rep(tlon[i], length(LATS)), LATS, tt[j], obs[j],
                            r$calibration, lp, SHADE_RATIO)
      if (!any(is.finite(ll))) next
      f <- function(x) ll[which.min(abs(LATS - x))]
      pole <- f(tlat[i] + DPROBE); equa <- f(tlat[i] - DPROBE)   # northern: pole = +
      w <- exp(ll - max(ll)); s <- sum(w); if (!is.finite(s) || s <= 0) next
      w <- w / s
      mu <- sum(w * LATS); md <- LATS[which.max(ll)]
      R[[length(R) + 1]] <- data.table(
        id = tg, win = i, lam_mult = lm, true_lat = tlat[i],
        gap = pole - equa,            # P1/P2: >0 means cheaper poleward
        skew = mu - md,               # P3: >0 means mean poleward of mode
        prof_sd = sqrt(sum(w * (LATS - mu)^2)),
        bias = mu - tlat[i])
    }
  }
}
D <- rbindlist(R)
fwrite(D, "scratch/nes_calibration/profile_asymmetry.csv")

S <- D[, .(n = .N,
           gap = round(median(gap, na.rm = TRUE), 2),
           frac_pole_cheaper = round(mean(gap > 0, na.rm = TRUE), 3),
           skew = round(median(skew, na.rm = TRUE), 3),
           prof_sd = round(median(prof_sd, na.rm = TRUE), 2),
           bias = round(mean(bias, na.rm = TRUE), 3)), by = lam_mult][order(lam_mult)]
cat("gap  = logL(truth+3) - logL(truth-3), nats. >0 means being wrong POLEWARD is cheaper.\n")
cat("skew = profile mean - mode. >0 means the poleward tail drags the mean.\n\n")
cat(sprintf("%9s %6s %9s %14s %9s %9s %9s\n",
            "lam_mult", "n", "gap", "frac_pole_ch", "skew", "prof_sd", "mean_bias"))
for (i in seq_len(nrow(S)))
  cat(sprintf("%9.3f %6d %+9.2f %14.3f %+9.3f %9.2f %+9.3f\n", S$lam_mult[i], S$n[i],
              S$gap[i], S$frac_pole_cheaper[i], S$skew[i], S$prof_sd[i], S$bias[i]))

cat("\n=== verdicts ===\n")
cat(sprintf("P1 poleward is cheaper at the default lambda: %s (median gap %+.2f nats, %.0f%% of windows)\n",
            if (S[lam_mult == 1]$gap > 0) "YES" else "NO",
            S[lam_mult == 1]$gap, 100 * S[lam_mult == 1]$frac_pole_cheaper))
ct <- suppressWarnings(cor.test(log2(D$lam_mult), D$gap))
cat(sprintf("P2 gap grows with lambda: r = %+.3f, p = %.3g\n", ct$estimate, ct$p.value))
sl <- coef(lm(gap ~ lam_mult, data = D))[2]
cat(sprintf("   slope %.2f nats per unit lambda multiplier (proportionality would be linear)\n", sl))
ct2 <- suppressWarnings(cor.test(log2(D$lam_mult), D$skew))
cat(sprintf("P3 skew grows with lambda: r = %+.3f, p = %.3g\n", ct2$estimate, ct2$p.value))
cat(sprintf("P4 profile flattens at loose lambda: sd %.1f at %.3fx vs %.1f at %.0fx\n",
            S[1]$prof_sd, S[1]$lam_mult, S[nrow(S)]$prof_sd, S[nrow(S)]$lam_mult))
cat("\nIf P1-P4 all hold, the asymmetric-shading picture is the mechanism and lambda is\n")
cat("its gain control: tight lambda buys latitude information at the price of a\n")
cat("poleward pull, loose lambda hands the answer to the prior.\n")
