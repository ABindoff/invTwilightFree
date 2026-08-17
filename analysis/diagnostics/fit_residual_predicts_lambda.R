# Does per-window RESIDUAL SPREAD (truth-free) predict that window's lambda_flip
# (which needs Argos)? If so, lambda can be set per window from the tag alone.
#
# WHY THIS IS THE TEST THAT MATTERS. lambda is not a per-tag constant: within a
# deployment the bias-zeroing lambda varies by a median factor of 34.6, and only 29%
# of the variance is between deployments. So any usable estimator has to be
# per-window, and it has to work without ground truth. The prob_slab sweep showed the
# bias sign-flip sits at lambda ~5 almost regardless of the slab price (flip moves
# 1.2x across a 17.5x range of prob_slab), which says the flip is pinned by the DATA
# -- by how much shading is actually present -- rather than by the model's own
# constants.
#
# QUANTITATIVE PREDICTION, not just directional: if the flip is where lambda matches
# the residual scale, then lambda_flip is proportional to 1/spread, so
#         log(lambda_flip) vs log(spread) has slope about -1.
# A slope near 0 means residual spread carries no information about lambda; a slope
# near -1 means it carries essentially all of it.
#
# TRUTH-FREE CONSTRUCTION. The position for each window is estimated by maximum
# likelihood from that window's own light (coarse 2-D pass at lambda = 1), never from
# Argos. Residuals are observed minus expected at that estimated position. Argos is
# used ONLY to compute the target, lambda_flip, which is what we are trying to
# predict.
#
# CONTROLS.
#  - `mean_light` is carried as a NEGATIVE control: window brightness is a
#    season/geometry property, not a shading property, so it should NOT predict
#    lambda_flip. If it predicts as well as the residual statistics, they are all
#    proxying something else (time of year, latitude) and the result is void.
#  - several spread statistics are computed rather than one, because "residual
#    spread" is ambiguous: shading is one-sided, so the spread of the NEGATIVE
#    residuals is the physically meaningful one and the symmetric spread is the
#    naive one. They are reported side by side.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))
OUT <- "scratch/nes_calibration/residual_vs_lambda.csv"

CAL_DAYS <- 15; STEP_HOURS <- 12
LATS <- seq(20, 70, by = 0.5); LONS <- seq(150, 250, by = 2)
MAXWIN <- 20; SHADE_RATIO <- 2

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

# expected light from the clamped-linear response, mirroring the engine
expect <- function(z, cal, ml) pmin(pmax(cal[1] - cal[2] * z, 0), ml)

R <- list()
for (tg in names(tags)) {
  g <- tags[[tg]]; r <- resp[[tg]]; if (is.null(r)) next
  tt <- as.numeric(g$light$time); obs <- pmax(0, g$light$light - r$baseline)
  brk <- seq(min(tt), max(tt), by = STEP_HOURS * 3600)
  kk <- findInterval(tt, brk)
  lam0 <- 1 / (r$max_light * 0.5)
  use <- which(tabulate(kk, length(brk)) >= 16)     # identical selection to the sweep
  if (!length(use)) next
  use <- use[seq(1, length(use), length.out = min(MAXWIN, length(use)))]
  message(tg, " (", length(use), " windows)")
  for (i in use) {
    j <- which(kk == i); if (length(j) < 16) next
    # --- position from THIS WINDOW'S OWN LIGHT, no Argos ---------------------
    gr <- expand.grid(lon = LONS, lat = LATS)
    ll <- eval_logpk_grid(gr$lon, gr$lat, tt[j], obs[j], r$calibration,
                          c(lam0, r$max_light, 0.10), SHADE_RATIO)
    bi <- which.max(ll)
    z  <- solar_zenith(tt[j], rep(gr$lon[bi], length(j)), rep(gr$lat[bi], length(j)))
    mu <- expect(z, r$calibration, r$max_light)
    res <- obs[j] - mu
    lit <- mu > 0.05 * r$max_light        # only where light was expected at all;
                                          # a residual in the dark is not shading
    if (sum(lit) < 8) next
    rl <- res[lit]
    R[[length(R) + 1]] <- data.table(
      id = tg, win = i,
      # symmetric (naive) spread
      mad_all  = stats::mad(rl),
      # one-sided: shading only pushes light DOWN, so this is the physical one
      mad_neg  = stats::mad(rl[rl < 0]),
      mean_neg = -mean(pmin(rl, 0)),
      # fraction of lit observations well below expectation
      frac_sh  = mean(rl < -0.1 * r$max_light),
      # NEGATIVE CONTROL: brightness is geometry/season, not shading
      mean_light = mean(mu),
      max_light = r$max_light)
  }
}
X <- rbindlist(R)

# --- target: per-window lambda_flip from the pslab sweep (pslab = 0.10) -------
S <- fread("scratch/nes_calibration/pslab_sweep.csv")[pslab == 0.10]
S[, id := as.character(id)]; X[, id := as.character(id)]
xing <- function(lam, y) {
  x <- log2(lam); o <- order(x); x <- x[o]; y <- y[o]
  k <- which(y[-1] * y[-length(y)] < 0); if (!length(k)) return(NA_real_)
  i <- tail(k, 1)
  2^(x[i] - y[i] * (x[i+1] - x[i]) / (y[i+1] - y[i]))
}
FL <- S[, .(lam_flip = xing(lam_mult, bias)), by = .(id, win)][is.finite(lam_flip)]
D <- merge(X, FL, by = c("id", "win"))
fwrite(D, OUT)
cat(sprintf("\n%d windows with both a residual statistic and a lambda_flip (%d tags)\n\n",
            nrow(D), uniqueN(D$id)))

# scale-free: express spreads relative to the tag's own max_light
D[, `:=`(s_all = mad_all / max_light, s_neg = mad_neg / max_light,
         s_mean = mean_neg / max_light, m_light = mean_light / max_light)]
D <- D[is.finite(s_neg) & s_neg > 0 & is.finite(lam_flip)]

cat(sprintf("%-24s %9s %9s %9s %11s\n", "predictor", "r", "p", "slope", "verdict"))
tst <- function(v, lab, ctrl = FALSE) {
  ok <- is.finite(v) & v > 0
  if (sum(ok) < 10) { cat(sprintf("%-24s   too few\n", lab)); return(invisible()) }
  ct <- suppressWarnings(cor.test(log(v[ok]), log(D$lam_flip[ok])))
  sl <- coef(lm(log(D$lam_flip[ok]) ~ log(v[ok])))[2]
  cat(sprintf("%-24s %+9.3f %9.2g %+9.2f %11s\n", lab, ct$estimate, ct$p.value, sl,
              if (ctrl) (if (ct$p.value > 0.05) "clean" else "*** CONFOUND")
              else (if (ct$p.value < 0.05) "predicts" else "no signal")))
}
tst(D$s_neg,   "MAD of neg residuals")
tst(D$s_all,   "MAD of all residuals")
tst(D$s_mean,  "mean shortfall")
tst(D$frac_sh, "frac heavily shaded")
tst(D$m_light, "mean expected light [CTRL]", ctrl = TRUE)

cat("\n  PREDICTION: slope near -1 if lambda_flip tracks 1/residual-scale.\n")
cat("  Slope near 0 means residual spread carries no information about lambda.\n")

# within-deployment: removes any between-tag confound
W <- D[, .(x = log(s_neg) - mean(log(s_neg)), y = log(lam_flip) - mean(log(lam_flip))), by = id]
W <- W[is.finite(x) & is.finite(y)]
if (nrow(W) > 20) {
  ct <- suppressWarnings(cor.test(W$x, W$y))
  cat(sprintf("\n  WITHIN-deployment (both centred): r = %+.3f, p = %.3g, slope %+.2f, n = %d\n",
              ct$estimate, ct$p.value, coef(lm(y ~ x, data = W))[2], nrow(W)))
  cat("  this is the one that matters: lambda must vary WITHIN a deployment, so a\n")
  cat("  between-tag correlation would not give a usable per-window estimator.\n")
}
