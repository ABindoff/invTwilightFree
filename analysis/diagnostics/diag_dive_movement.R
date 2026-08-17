# Does the DEPTH channel inform the movement model rather than the bathymetry?
#
# The argument is a time budget: hours spent diving are hours not spent
# travelling, so dive effort over a knot should predict how far the animal
# actually got. If it does, the depth record supplies a per-knot movement scale
# driven by an OBSERVED covariate -- which is a much better position than a
# latent ARS/transit state, because there is nothing to identify. It would also
# go straight at the scale problem: the one-step scale here is 31 km/12 h and
# the displacement-equivalent scale 316, and a memoryless model splits the
# difference at ~90, which is neither quantity.
#
# Tested against Argos truth, so this is the real step length, not a modelled one.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))

STEP_H <- 12
arch <- lapply(readRDS("analysis/cache/nes/archives_v3_30min.rds"), `[[`, "main")
meta <- nes_meta(); ar <- nes_argos()

rows <- list()
for (id in names(arch)) {
  m <- meta[match(id, meta$id), ]
  if (is.na(m$ptt) || !nzchar(m$ptt)) next
  a <- ar$fixes[Ptt == m$ptt]; df <- ar$deploy[Ptt == m$ptt]
  if (!nrow(df) || nrow(a) < 50) next
  x <- as.data.table(arch[[id]])
  t0 <- max(min(x$time), as.POSIXct(max(m$deploy_date, df$user_date), tz = "UTC"))
  t1 <- min(max(x$time), max(a$time))
  x <- x[time >= t0 & time <= t1]
  if (nrow(x) < 1000) next

  # knot edges, and dive effort summarised within each knot
  brk <- seq(as.numeric(t0), as.numeric(t1), by = STEP_H * 3600)
  if (length(brk) < 6) next
  x[, kn := findInterval(as.numeric(time), brk)]
  eff <- x[kn >= 1 & kn < length(brk), .(
    n_win      = .N,
    mean_dmax  = mean(depth_max, na.rm = TRUE),   # mean of the 30-min maxima
    med_dmax   = median(depth_max, na.rm = TRUE),
    frac_deep  = mean(depth_max > 200, na.rm = TRUE),
    frac_surf  = mean(depth_min < 5, na.rm = TRUE),
    max_dmax   = max(depth_max, na.rm = TRUE)), by = kn]
  eff <- eff[n_win >= STEP_H]                      # a reasonably complete knot

  # true displacement over the same knot, from Argos
  mid <- as.POSIXct(brk[-length(brk)] + STEP_H*1800, origin = "1970-01-01", tz = "UTC")
  p0 <- argos_at(a, as.POSIXct(brk[-length(brk)], origin = "1970-01-01", tz = "UTC"))
  p1 <- argos_at(a, as.POSIXct(brk[-1], origin = "1970-01-01", tz = "UTC"))
  step <- gc_km(p0$lon, p0$lat, p1$lon, p1$lat)
  d <- data.table(id = id, kn = seq_along(step), step_km = step,
                  time = mid, lat = p0$lat)
  rows[[length(rows)+1]] <- merge(d, eff, by = "kn")
}
D <- rbindlist(rows)
D <- D[is.finite(step_km) & is.finite(mean_dmax) & step_km < 200]
cat(sprintf("%d knots across %d animals\n", nrow(D), length(unique(D$id))))
cat(sprintf("step length: mean %.1f km / %d h, sd %.1f\n",
            mean(D$step_km), STEP_H, sd(D$step_km)))

cat("\n=== step length by dive effort (mean of the 30-min depth maxima) ===\n")
D[, ebin := cut(mean_dmax, breaks = quantile(mean_dmax, seq(0, 1, 0.2)),
                include.lowest = TRUE)]
print(as.data.frame(D[, .(knots = .N, mean_depth_m = round(mean(mean_dmax)),
  frac_deep = round(mean(frac_deep), 2), step_km = round(mean(step_km), 1),
  sd_step = round(sd(step_km), 1)), by = ebin][order(ebin)]), row.names = FALSE)

cat("\n=== and by the fraction of the knot spent below 200 m ===\n")
D[, fbin := cut(frac_deep, breaks = c(-0.01, 0.2, 0.4, 0.6, 0.8, 1))]
print(as.data.frame(D[, .(knots = .N, step_km = round(mean(step_km), 1),
  sd_step = round(sd(step_km), 1)), by = fbin][order(fbin)]), row.names = FALSE)

cat("\n=== correlations ===\n")
cc <- function(v, nm) cat(sprintf("  %-12s pooled %+.3f   within-animal %+.3f\n", nm,
  cor(D[[v]], D$step_km, use = "complete.obs"),
  mean(D[, .(r = if (.N > 30) cor(get(v), step_km, use = "complete.obs") else NA_real_),
          by = id]$r, na.rm = TRUE)))
for (v in c("mean_dmax", "med_dmax", "frac_deep", "frac_surf", "max_dmax")) cc(v, v)

cat("\n=== how much of the step-length variance does dive effort explain? ===\n")
m0 <- lm(log(pmax(step_km, 1)) ~ factor(id), data = D)
m1 <- lm(log(pmax(step_km, 1)) ~ factor(id) + mean_dmax + frac_deep + frac_surf, data = D)
cat(sprintf("animal only          R2 = %.3f\n", summary(m0)$r.squared))
cat(sprintf("animal + dive effort R2 = %.3f  (increment %.3f, p = %.3g)\n",
            summary(m1)$r.squared, summary(m1)$r.squared - summary(m0)$r.squared,
            anova(m0, m1)$`Pr(>F)`[2]))
print(round(summary(m1)$coefficients[c("mean_dmax", "frac_deep", "frac_surf"), ], 5))

cat("\n=== what a covariate-driven sigma would buy ===\n")
# residual sd of log step, with and without the covariate: this is the factor by
# which a per-knot movement prior could be tightened
s0 <- sd(residuals(m0)); s1 <- sd(residuals(m1))
cat(sprintf("residual sd of log step: %.3f -> %.3f, a %.1f%% reduction\n",
            s0, s1, 100 * (1 - s1/s0)))
cat(sprintf("implied sigma ratio between the least and most active quintiles: %.2fx\n",
            exp(diff(range(D[, .(m = mean(log(pmax(step_km,1)))), by = ebin]$m)))))
saveRDS(D, file.path(SCRATCH, "dive_movement.rds"))
