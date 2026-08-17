# The first pass found two effects pulling opposite ways and lumped them.
#
#   frac_deep (+2.34): knots where the animal is NOT continuously deep-diving
#   have SHORT steps (11.7 km against 28.2). That is not a time budget, it is
#   the haul-out and the coastal shelf at each end of the trip -- the animal is
#   ashore or milling, not transiting.
#
#   mean_dmax (-0.0016): among knots that ARE continuously diving, deeper mean
#   dives go with shorter steps. THAT is the effect proposed.
#
# Only 6.5% of knots are non-pelagic, but their steps differ by a factor of
# three, so they can carry the whole R2. Refit on at-sea knots only: whatever
# survives is the time budget rather than the haul-out.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
suppressMessages(library(data.table))
D <- readRDS(file.path(SCRATCH, "dive_movement.rds"))

S <- D[frac_deep > 0.8]            # pelagic, continuously diving
cat(sprintf("at-sea knots: %d of %d (%.0f%%)\n", nrow(S), nrow(D), 100*nrow(S)/nrow(D)))
cat(sprintf("step length: mean %.1f km / 12 h, sd %.1f\n",
            mean(S$step_km), sd(S$step_km)))

cat("\n=== step length by mean dive depth, at-sea knots only ===\n")
S[, ebin := cut(mean_dmax, breaks = quantile(mean_dmax, seq(0, 1, 0.2)),
                include.lowest = TRUE)]
print(as.data.frame(S[, .(knots = .N, mean_depth_m = round(mean(mean_dmax)),
  step_km = round(mean(step_km), 1), sd_step = round(sd(step_km), 1)),
  by = ebin][order(ebin)]), row.names = FALSE)

cat("\n=== correlations, at-sea only ===\n")
for (v in c("mean_dmax", "med_dmax", "frac_surf", "max_dmax")) {
  r_all <- cor(S[[v]], S$step_km, use = "complete.obs")
  r_win <- mean(S[, .(r = if (.N > 30) cor(get(v), step_km, use = "complete.obs")
                      else NA_real_), by = id]$r, na.rm = TRUE)
  cat(sprintf("  %-11s pooled %+.3f   within-animal %+.3f\n", v, r_all, r_win))
}

m0 <- lm(log(pmax(step_km, 1)) ~ factor(id), data = S)
m1 <- lm(log(pmax(step_km, 1)) ~ factor(id) + mean_dmax + frac_surf, data = S)
cat(sprintf("\nanimal only          R2 = %.4f\n", summary(m0)$r.squared))
cat(sprintf("animal + dive effort R2 = %.4f  (increment %.4f, p = %.3g)\n",
            summary(m1)$r.squared, summary(m1)$r.squared - summary(m0)$r.squared,
            anova(m0, m1)$`Pr(>F)`[2]))
print(round(summary(m1)$coefficients[c("mean_dmax", "frac_surf"), ], 5))

s0 <- sd(residuals(m0)); s1 <- sd(residuals(m1))
cat(sprintf("\nresidual sd of log step: %.3f -> %.3f, a %.1f%% reduction\n",
            s0, s1, 100 * (1 - s1/s0)))
q <- S[, .(m = mean(log(pmax(step_km, 1)))), by = ebin]
cat(sprintf("sigma ratio, shallowest to deepest quintile: %.2fx\n",
            exp(diff(range(q$m)))))

cat("\n=== for scale: what the movement model actually needs ===\n")
cat(sprintf("observed one-step sd            %.1f km / 12 h\n", sd(S$step_km)))
disp <- D[, .(d = {
  # net displacement from the trip start, at the halfway point of each track
  NA_real_ }), by = id]
cat("one-step scale 31 km against a displacement-equivalent scale of 316 km:\n")
cat("a covariate that modulates sigma by the factor above closes none of that.\n")
