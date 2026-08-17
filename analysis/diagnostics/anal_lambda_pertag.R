# Is there a good DEFAULT lambda, or does it have to be per-tag?
#
# Stage 1 already showed a global lambda removes the MEAN bias (zero near 3.7x
# tighter) while per-logger bias stays scattered -2.6 to +3.5. That is suggestive but
# not the answer: the scatter could be irreducible noise in each logger's bias, or it
# could be that each logger wants a DIFFERENT lambda. Those have opposite
# implications -- the first says ship a default and quote a residual spread, the
# second says lambda must be estimated per tag and a default cannot work.
#
# Distinguish them: for each logger, find the lambda that zeroes ITS bias, then ask
# how tightly those agree. Tight agreement => a default exists. Wide disagreement
# => it does not, and the question becomes how to estimate lambda without truth.
#
# Loggers, not deployments (15 vs 29): lambda is a property of the sensor's residual
# scale, and two deployments of one logger share it.
suppressMessages(library(data.table))
# run from the package root (setwd removed for portability)
D <- fread("scratch/nes_calibration/bias_gridsearch.csv")
D <- D[band == "all" & delta == 0 & !is.na(serial)]

per <- D[, .(bias = mean(bias, na.rm = TRUE), sd = mean(sd, na.rm = TRUE)),
         by = .(serial, lam_mult)][order(serial, lam_mult)]

# zero crossing in log2(lambda), linearly interpolated; NA if it never crosses
cross <- per[, {
  x <- log2(lam_mult); y <- bias
  o <- order(x); x <- x[o]; y <- y[o]
  k <- which(y[-1] * y[-length(y)] < 0)
  if (!length(k)) {
    .(lam_opt = NA_real_, note = if (all(y < 0)) ">4 (never reached)" else "<0.25")
  } else {
    i <- k[1]
    .(lam_opt = 2^(x[i] + (0 - y[i]) * (x[i+1] - x[i]) / (y[i+1] - y[i])), note = "interior")
  }
}, by = serial]

cat("=== lambda that zeroes each logger's own latitude bias ===\n")
cat("(multiplier on the shipped default; >1 means the default is too loose)\n\n")
cat(sprintf("%-12s %10s  %s\n", "logger", "lam_opt", "note"))
for (i in seq_len(nrow(cross)))
  cat(sprintf("%-12s %10s  %s\n", cross$serial[i],
              if (is.na(cross$lam_opt[i])) "--" else sprintf("%.2f", cross$lam_opt[i]),
              cross$note[i]))

ok <- cross[!is.na(lam_opt)]
cat(sprintf("\ninterior optima: %d of %d loggers\n", nrow(ok), nrow(cross)))
if (nrow(ok) >= 3) {
  cat(sprintf("  median %.2f | IQR %.2f-%.2f | range %.2f-%.2f | fold-spread %.1fx\n",
              median(ok$lam_opt), quantile(ok$lam_opt, .25), quantile(ok$lam_opt, .75),
              min(ok$lam_opt), max(ok$lam_opt), max(ok$lam_opt)/min(ok$lam_opt)))
  cat(sprintf("  sd of log2(lam_opt) = %.2f  (1.0 means a factor-of-2 spread)\n",
              sd(log2(ok$lam_opt))))
}

# What does a single default cost each logger, versus its own optimum?
cat("\n=== cost of a global default ===\n")
grid_lams <- sort(unique(per$lam_mult))
for (L in grid_lams) {
  b <- per[lam_mult == L]
  cat(sprintf("  lam = %4.2f : mean |bias| %.3f | worst logger %+0.2f | %d/%d loggers within 1 deg\n",
              L, mean(abs(b$bias)), b$bias[which.max(abs(b$bias))],
              sum(abs(b$bias) <= 1), nrow(b)))
}
best <- grid_lams[which.min(vapply(grid_lams, function(L) mean(abs(per[lam_mult == L]$bias)), 0))]
cat(sprintf("\n  best single default on this grid: %.2f\n", best))

# The decisive comparison: how much of the bias does a per-logger lambda remove that
# the best global one cannot?
if (nrow(ok) >= 3) {
  glob <- mean(abs(per[lam_mult == best]$bias))
  cat(sprintf("\n  mean |bias| at the best GLOBAL lambda   : %.3f\n", glob))
  cat("  mean |bias| at each logger's OWN lambda : 0.000 by construction\n")
  cat(sprintf("  loggers whose own optimum lies outside the swept range: %d\n",
              nrow(cross) - nrow(ok)))
  cat("\n  If the interior optima cluster tightly, a default works and the residual\n")
  cat("  scatter is noise. If they span a large factor, lambda is genuinely per-tag\n")
  cat("  and the problem becomes estimating it WITHOUT ground truth.\n")
}

# Does a logger's optimal lambda relate to anything measurable WITHOUT truth?
# If it correlates with something observable from the light record alone, that is the
# beginning of a truth-free estimator.
cat("\n=== is lam_opt predictable from the profile width (no truth needed)? ===\n")
w <- merge(ok, per[lam_mult == 1, .(serial, sd_at_default = sd)], by = "serial")
if (nrow(w) >= 4) {
  ct <- suppressWarnings(cor.test(log2(w$lam_opt), w$sd_at_default))
  cat(sprintf("  cor(log2 lam_opt, profile sd at default) = %+.3f, p = %.3f, n = %d\n",
              ct$estimate, ct$p.value, nrow(w)))
  cat("  profile sd is computable from the tag's own light with no Argos at all,\n")
  cat("  so a real correlation here would be a truth-free handle on lambda.\n")
}
