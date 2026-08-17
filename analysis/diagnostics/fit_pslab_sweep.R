# Does prob_slab move the lambda at which the latitude bias flips sign?
#
# THE MECHANISM UNDER TEST. Shading puts observations BELOW clear sky. The model
# reduces that mismatch by predicting less light, which in northern summer means a
# shorter day, i.e. EQUATORWARD -- the -1.5 deg we have been chasing. But a shaded
# observation cannot cost more than the slab: its price is
#       min( lambda * mismatch , -log(prob_slab / max_light) )
# so once lambda is large enough that shading exceeds the slab ceiling, shaded
# observations are absorbed as contamination and stop pulling position at all. The
# unshaded remainder then takes over and the bias flips poleward. Measured: the flip
# sits near lambda ~5 at prob_slab = 0.10.
#
# THE PREDICTION, which is quantitative rather than directional. The flip is where
#       lambda_flip * Delta  =  -log(prob_slab / max_light)
# so lambda_flip should FALL as prob_slab RISES, and should be roughly LINEAR in
#       -log(prob_slab)
# A cheaper slab (large prob_slab) saturates sooner, so the flip arrives at smaller
# lambda. If instead the flip does not move, the slab is not the mechanism and the
# sign change has some other cause.
#
# CONTROL: prob_slab = 0.10, lambda = 1 is the shipped setting and must reproduce the
# known all-band emission bias of about -1.2 to -1.3 deg.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))
OUT <- "scratch/nes_calibration/pslab_sweep.csv"

CAL_DAYS <- 15; STEP_HOURS <- 12
LATS   <- seq(20, 70, by = 0.5)
LAMS   <- 2^seq(-2, 5, by = 0.5)
PSLABS <- c(0.02, 0.05, 0.10, 0.20, 0.35)
MAXWIN <- 20
SHADE_RATIO <- 2

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
resp <- lapply(tags, function(g) {
  k <- g$light$time <= min(g$light$time) + CAL_DAYS * 86400
  fit_light_response(g$light$time[k], g$light$light[k], g$p0[1], g$p0[2])
})
cat(sprintf("%d deployments | %d lambdas x %d prob_slab = %d configs\n\n",
            length(tags), length(LAMS), length(PSLABS), length(LAMS) * length(PSLABS)))

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
  use <- use[seq(1, length(use), length.out = min(MAXWIN, length(use)))]
  message(tg, " (", length(use), " windows)")
  for (i in use) {
    j <- which(kk == i); if (length(j) < 16) next
    for (ps in PSLABS) for (lm in LAMS) {
      ll <- eval_logpk_grid(rep(tlon[i], length(LATS)), LATS, tt[j], obs[j],
                            r$calibration, c(lam0 * lm, r$max_light, ps), SHADE_RATIO)
      w <- exp(ll - max(ll)); s <- sum(w)
      if (!is.finite(s) || s <= 0) next
      w <- w / s; mu <- sum(w * LATS)
      R[[length(R) + 1]] <- data.table(id = tg, win = i, lam_mult = lm, pslab = ps,
                                       bias = mu - tlat[i],
                                       prof_sd = sqrt(sum(w * (LATS - mu)^2)),
                                       slab_cost = -log(ps / r$max_light))
    }
  }
}
D <- rbindlist(R)
fwrite(D, OUT)

S <- D[, .(bias = mean(bias, na.rm = TRUE), sd = mean(prof_sd, na.rm = TRUE),
           slab_cost = mean(slab_cost)), by = .(pslab, lam_mult)]
cat("\n=== mean latitude bias (rows lambda, cols prob_slab) ===\n")
print(as.data.frame(dcast(S, lam_mult ~ pslab, value.var = "bias")), row.names = FALSE, digits = 3)

cat(sprintf("\nCONTROL pslab 0.10, lambda 1: bias %+0.3f (known about -1.2 to -1.3) -- %s\n",
            S[pslab == 0.10 & lam_mult == 1]$bias,
            if (abs(S[pslab == 0.10 & lam_mult == 1]$bias + 1.25) < 0.35) "reproduces" else "*** DIFFERS"))

xing <- function(lam, y) {
  x <- log2(lam); o <- order(x); x <- x[o]; y <- y[o]
  k <- which(y[-1] * y[-length(y)] < 0); if (!length(k)) return(NA_real_)
  i <- tail(k, 1)      # the UPPER root: the equatorward -> poleward flip
  2^(x[i] - y[i] * (x[i+1] - x[i]) / (y[i+1] - y[i]))
}
FL <- S[, .(lam_flip = xing(lam_mult, bias), slab_cost = mean(slab_cost)), by = pslab][order(pslab)]
cat("\n=== the flip point vs prob_slab ===\n")
cat(sprintf("%8s %12s %12s\n", "pslab", "lam_flip", "slab_cost"))
for (i in seq_len(nrow(FL)))
  cat(sprintf("%8.2f %12s %12.2f\n", FL$pslab[i],
              if (is.na(FL$lam_flip[i])) "--" else sprintf("%.2f", FL$lam_flip[i]),
              FL$slab_cost[i]))
G <- FL[is.finite(lam_flip)]
if (nrow(G) >= 3) {
  ct <- suppressWarnings(cor.test(log(G$pslab), log(G$lam_flip)))
  cat(sprintf("\n  cor(log pslab, log lam_flip) = %+.3f, p = %.4f, n = %d\n",
              ct$estimate, ct$p.value, nrow(G)))
  ct2 <- suppressWarnings(cor.test(G$slab_cost, G$lam_flip))
  cat(sprintf("  cor(slab cost, lam_flip)     = %+.3f, p = %.4f  <- the mechanism's own variable\n",
              ct2$estimate, ct2$p.value))
  cat(sprintf("  flip moves %.1fx across a %.1fx range of prob_slab\n",
              max(G$lam_flip)/min(G$lam_flip), max(G$pslab)/min(G$pslab)))
  cat("\n  PREDICTED: negative cor(log pslab, log lam_flip) -- a cheaper slab saturates\n")
  cat("  sooner, so the flip arrives at smaller lambda. A flat flip point would mean\n")
  cat("  the slab is NOT the mechanism.\n")
}
