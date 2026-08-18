# HOW MUCH WOULD MARGINALISING OVER THE CALIBRATION POSTERIOR ACTUALLY WIDEN THINGS?
#
# THE QUESTION. The pipeline fits the light response once and plugs the point
# estimate into the engine, so the position posterior is CONDITIONAL on a calibration
# taken as known. Uncertainty in the calibration is expressed nowhere. The
# statistically clean alternative is to treat the calibration as uncertain and
# marginalise: p(track | y) = integral p(track | y, cal) p(cal | y) d cal.
#
# That is correct, and it is the right thing to do if the paper claims calibrated
# uncertainty. The question worth settling BEFORE building it is how big the effect
# is -- specifically whether it can touch the coverage failure (latitude 0.64 against
# a nominal 0.95, needing intervals roughly 2.6x wider).
#
# We can size it without any new fitting. `diag_profile_recovery.R` already evaluated
# both `log_z` and the latitude bias on a grid of z50 offsets, per tag. Normalised
# evidence weights on that grid ARE a calibration posterior (a flat prior on the
# grid), and the spread of the position estimate under those weights IS the
# marginalisation effect, to first order.
#
#   sd_marg(position) = sd of the fitted latitude across the calibration posterior
#
# Combined in quadrature with the within-calibration posterior sd, that gives the
# marginalised interval. If it moves the total by a few percent, marginalising is a
# correctness improvement and not a fix for coverage; if it moves it by tens of
# percent, it is both.
#
# CAVEAT, stated up front: this varies z50 only. The full calibration posterior also
# covers floor, amplitude, scale, lambda and prob_slab -- and, more importantly, the
# response FORM itself (clamped line vs logistic vs lookup table), which is model
# uncertainty rather than parameter uncertainty and is not represented on any
# parametric grid. So this is a LOWER BOUND on what full marginalisation would do.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
suppressMessages(library(data.table))

D <- fread("scratch/nes_calibration/profile_recovery.csv",
           colClasses = list(character = "id"))
setorder(D, id, dz)

cat("=== per tag: the calibration posterior implied by log_z, and what it costs ===\n")
res <- D[, {
  w <- exp(log_z - max(log_z)); w <- w / sum(w)
  z_mean <- sum(w * dz)
  z_sd   <- sqrt(sum(w * (dz - z_mean)^2))
  b_mean <- sum(w * bias)
  b_sd   <- sqrt(sum(w * (bias - b_mean)^2))     # position spread from calibration
  .(z_hat = z_mean, z_sd = z_sd, bias_marg = b_mean, sd_from_cal = b_sd)
}, by = id]
print(as.data.frame(res[, lapply(.SD, function(x)
  if (is.numeric(x)) round(x, 4) else x)]), row.names = FALSE)

cat(sprintf("\n  median calibration posterior sd : %.3f deg of z50\n", median(res$z_sd)))
cat(sprintf("  median induced position sd      : %.4f deg of latitude\n",
            median(res$sd_from_cal)))

# within-calibration posterior sd, from the same fits at the evidence peak
sd_within <- 2.35   # mean per-knot posterior sd on this arm (diag_axis_prior_fits)
tot <- sqrt(sd_within^2 + median(res$sd_from_cal)^2)
cat(sprintf("\n=== combined in quadrature ===\n"))
cat(sprintf("  within-calibration posterior sd : %.3f deg\n", sd_within))
cat(sprintf("  calibration contribution        : %.4f deg\n", median(res$sd_from_cal)))
cat(sprintf("  marginalised total              : %.4f deg\n", tot))
cat(sprintf("  widening                        : %.3f%%\n", 100 * (tot / sd_within - 1)))

cat("\n=== what coverage actually needs ===\n")
cat("  real-data latitude coverage is 0.64 against a nominal 0.95; error sd 3.74 deg\n")
cat("  against a reported 1.38, so intervals need to be about 2.6x wider.\n")
cat(sprintf("  marginalising over z50 delivers %.3fx.\n", tot / sd_within))

cat("\n=== how sharp is the evidence, really? ===\n")
cat("  (a 2-unit drop in log_z is roughly a 1-sd interval on z50)\n")
P <- D[, .(log_z = sum(log_z)), by = dz][order(dz)]
P[, rel := log_z - max(log_z)]
print(as.data.frame(P[, .(dz, rel_logz = round(rel, 1))]), row.names = FALSE)
cat(sprintf("\n  pooled across 6 tags the evidence falls %.0f units by dz = +-0.5,\n",
            -P[dz == 0.5]$rel))
cat("  so the pooled calibration posterior is very tight indeed -- which is the\n")
cat("  same result the Fisher analysis gave (se(z50) = 0.045 deg over six tags).\n")
cat("  A tight posterior is exactly what makes marginalising over it cheap.\n")
