# Does the dawn/dusk criterion pick the same lambda that Argos does?
#
# The criterion is truth-free by construction. This script is the VALIDATION: it
# compares the lambda that makes sd(z) = 1 against the lambda that zeroes the tag's
# latitude bias measured against Argos. If they agree, we have an estimator usable on
# tags with no ground truth. If they disagree, the criterion calibrates DISPERSION but
# not BIAS -- which would be a real finding, since those are different failures and
# nothing guarantees one lambda fixes both.
suppressMessages(library(data.table))
# run from the package root (setwd removed for portability)
D <- fread("scratch/nes_calibration/dawndusk.csv")

cat("=== pooled: sd(z) and mean(z) against lambda ===\n")
cat("sd(z) = 1 is calibrated; > 1 over-confident (lambda too tight); < 1 under-confident\n")
cat("mean(z) far from 0 means a dawn/dusk ASYMMETRY, which lambda cannot fix\n\n")
P <- D[, .(n = .N, sd_z = sd(z, na.rm = TRUE), mean_z = mean(z, na.rm = TRUE),
           med_sd = median(sqrt(sd_dawn^2 + sd_dusk^2), na.rm = TRUE)), by = lam_mult][order(lam_mult)]
cat(sprintf("%9s %7s %9s %9s %10s\n", "lam_mult", "n", "sd(z)", "mean(z)", "claimed_sd"))
for (i in seq_len(nrow(P)))
  cat(sprintf("%9.3f %7d %9.3f %+9.3f %10.2f\n", P$lam_mult[i], P$n[i], P$sd_z[i],
              P$mean_z[i], P$med_sd[i]))

xing <- function(lam, y, target = 1) {
  x <- log2(lam); o <- order(x); x <- x[o]; y <- y[o]
  k <- which((y[-1] - target) * (y[-length(y)] - target) < 0)
  if (!length(k)) return(NA_real_)
  i <- k[1]
  2^(x[i] + (target - y[i]) * (x[i+1] - x[i]) / (y[i+1] - y[i]))
}
lam_pool <- xing(P$lam_mult, P$sd_z)
cat(sprintf("\npooled lambda with sd(z) = 1: %s\n",
            if (is.na(lam_pool)) "no crossing in range" else sprintf("%.2f", lam_pool)))

cat("\n=== per logger ===\n")
L <- D[!is.na(serial), .(sd_z = sd(z, na.rm = TRUE), n = .N), by = .(serial, lam_mult)]
lam_dd <- L[, .(lam_dd = xing(lam_mult, sd_z), ndays = max(n)), by = serial][order(serial)]
cat(sprintf("%-12s %10s %8s\n", "logger", "lam_dd", "n"))
for (i in seq_len(nrow(lam_dd)))
  cat(sprintf("%-12s %10s %8d\n", lam_dd$serial[i],
              if (is.na(lam_dd$lam_dd[i])) "--" else sprintf("%.2f", lam_dd$lam_dd[i]),
              lam_dd$ndays[i]))
gd <- lam_dd[!is.na(lam_dd)]
cat(sprintf("\ninterior crossings: %d of %d loggers", nrow(gd), nrow(lam_dd)))
if (nrow(gd) >= 3)
  cat(sprintf(" | median %.2f | range %.2f-%.2f | sd log2 %.2f\n",
              median(gd$lam_dd), min(gd$lam_dd), max(gd$lam_dd), sd(log2(gd$lam_dd))))
cat("\n")

# ---- VALIDATION against the Argos-optimal lambda ---------------------------
f <- "scratch/nes_calibration/bias_lambda_wide.csv"
if (!file.exists(f)) {
  cat("=== VALIDATION SKIPPED: the wide Argos lambda sweep has not finished ===\n")
} else {
  W <- fread(f)[band == "all" & !is.na(serial)]
  per <- W[, .(bias = mean(bias, na.rm = TRUE)), by = .(serial, lam_mult)]
  lam_ar <- per[, .(lam_argos = xing(lam_mult, bias, target = 0)), by = serial]
  M <- merge(lam_dd[, .(serial, lam_dd)], lam_ar, by = "serial")
  M <- M[!is.na(lam_dd) & !is.na(lam_argos)]
  cat("=== VALIDATION: truth-free lambda vs Argos-optimal lambda ===\n")
  cat(sprintf("%-12s %10s %12s %9s\n", "logger", "lam_dd", "lam_argos", "log2 gap"))
  for (i in seq_len(nrow(M)))
    cat(sprintf("%-12s %10.2f %12.2f %+9.2f\n", M$serial[i], M$lam_dd[i], M$lam_argos[i],
                log2(M$lam_dd[i]) - log2(M$lam_argos[i])))
  if (nrow(M) >= 4) {
    ct <- suppressWarnings(cor.test(log2(M$lam_dd), log2(M$lam_argos)))
    cat(sprintf("\n  n = %d loggers with both\n", nrow(M)))
    cat(sprintf("  cor(log2 lam_dd, log2 lam_argos) = %+.3f, p = %.4f\n",
                ct$estimate, ct$p.value))
    cat(sprintf("  median log2 gap %+.2f (a constant offset is correctable; scatter is not)\n",
                median(log2(M$lam_dd) - log2(M$lam_argos))))
    cat(sprintf("  sd of log2 gap  %.2f  (0.5 = within a factor of 1.4; 1.0 = a factor of 2)\n",
                sd(log2(M$lam_dd) - log2(M$lam_argos))))
    cat("\n  A strong positive correlation means the criterion tracks the right\n")
    cat("  answer and can be used without truth. A null correlation means it\n")
    cat("  calibrates dispersion but says nothing about bias.\n")
  } else cat("\n  too few loggers with both to correlate\n")
}
