# Read the 2-D (delta, lambda) surface from fit_bias_gridsearch.R.
# Separate from the sweep so the surface can be re-analysed without re-running it.
suppressMessages({ library(data.table) })
# run from the package root (setwd removed for portability)
D <- fread("scratch/nes_calibration/bias_gridsearch.csv")

# Logger-level first (15 units, not 29): average within logger, then across loggers.
per_log <- D[, .(bias = mean(bias, na.rm = TRUE), sd = mean(sd, na.rm = TRUE)),
             by = .(delta, lam_mult, band, serial)]
S <- per_log[, .(bias = mean(bias, na.rm = TRUE), sd = mean(sd, na.rm = TRUE),
                 loggers = .N), by = .(delta, lam_mult, band)]

cat("=== CONTROL 1: delta = 0, lam_mult = 1 must BE the shipped tangent ===\n")
cat("known from the band decomposition: day ~ +0.36, twilight ~ -1.25, night ~ -1.68\n")
ref <- S[delta == 0 & lam_mult == 1]
for (b in c("day", "twilight", "night", "all"))
  cat(sprintf("  %-9s bias %+0.3f (sd %5.2f)\n", b, ref[band == b]$bias, ref[band == b]$sd))
ok <- abs(ref[band == "twilight"]$bias - (-1.25)) < 0.25 &&
      abs(ref[band == "day"]$bias - 0.36) < 0.25
cat(sprintf("  verdict: %s\n\n", if (ok) "reproduces -- surface is readable"
            else "*** DOES NOT REPRODUCE: harness differs, stop here"))

pr <- function(bnd, what) {
  M <- dcast(S[band == bnd], delta ~ lam_mult, value.var = what)
  cat(sprintf("--- %s : %s (rows delta, cols lambda multiplier) ---\n", bnd, what))
  print(as.data.frame(M), row.names = FALSE, digits = 3)
  cat("\n")
}
pr("all", "bias"); pr("all", "sd")
pr("twilight", "bias"); pr("twilight", "sd")

cat("=== CONTROL 2: does any config buy low bias by flattening the profile? ===\n")
ref_sd <- ref[band == "all"]$sd
A <- S[band == "all"][order(abs(bias))]
cat(sprintf("reference profile sd at the shipped setting: %.2f\n", ref_sd))
cat(sprintf("%8s %9s %9s %8s %10s\n", "delta", "lam_mult", "bias", "sd", "verdict"))
for (i in seq_len(min(10, nrow(A))))
  cat(sprintf("%8.1f %9.2f %+9.3f %8.2f %10s\n", A$delta[i], A$lam_mult[i], A$bias[i],
              A$sd[i], if (A$sd[i] <= ref_sd * 1.05) "admissible" else "FLATTER"))

cat("\n=== in-sample optimum (all loggers), admissible only ===\n")
adm <- A[sd <= ref_sd * 1.05]
best <- adm[1]
cat(sprintf("  delta = %+0.1f deg, lam_mult = %.2f  ->  bias %+0.3f (sd %.2f)\n",
            best$delta, best$lam_mult, best$bias, best$sd))
cat(sprintf("  shipped setting                     ->  bias %+0.3f (sd %.2f)\n",
            ref[band == "all"]$bias, ref_sd))

cat("\n=== LEAVE-ONE-LOGGER-OUT (the honest number) ===\n")
cat("29 deployments come from 15 loggers, so leaving out a DEPLOYMENT would leak\n")
cat("through the shared sensor. Tune on 14 loggers, score on the 15th.\n\n")
pl <- per_log[band == "all"]
logs <- unique(pl$serial); logs <- logs[!is.na(logs)]
rows <- rbindlist(lapply(logs, function(L) {
  tr <- pl[serial != L, .(bias = mean(bias, na.rm = TRUE), sd = mean(sd, na.rm = TRUE)),
           by = .(delta, lam_mult)]
  tr <- tr[sd <= ref_sd * 1.05][order(abs(bias))]
  if (!nrow(tr)) return(NULL)
  ch <- tr[1]
  te <- pl[serial == L & delta == ch$delta & lam_mult == ch$lam_mult]
  base <- pl[serial == L & delta == 0 & lam_mult == 1]
  data.table(serial = L, delta = ch$delta, lam_mult = ch$lam_mult,
             held_out = te$bias[1], shipped = base$bias[1])
}))
print(as.data.frame(rows), row.names = FALSE, digits = 3)
cat(sprintf("\n  held-out mean |bias| tuned   : %.3f\n", mean(abs(rows$held_out), na.rm = TRUE)))
cat(sprintf("  held-out mean |bias| shipped : %.3f\n", mean(abs(rows$shipped), na.rm = TRUE)))
cat(sprintf("  in-sample  mean |bias| tuned : %.3f  (gap to held-out = overfitting)\n",
            abs(best$bias)))
cat(sprintf("  chosen delta: %s\n", paste(sort(unique(rows$delta)), collapse = ", ")))
cat(sprintf("  chosen lam_mult: %s\n", paste(sort(unique(rows$lam_mult)), collapse = ", ")))
cat("  a STABLE choice across folds means the surface has a real optimum;\n")
cat("  a scattered one means it is flat and the tuning is fitting noise.\n")
if (nrow(rows) > 2) {
  p <- suppressWarnings(wilcox.test(abs(rows$held_out), abs(rows$shipped), paired = TRUE)$p.value)
  cat(sprintf("  paired Wilcoxon, held-out tuned vs shipped: p = %.4f\n", p))
}

cat("\n=== which knob does the work? ===\n")
b_delta <- S[band == "all" & lam_mult == 1][order(delta)]
b_lam   <- S[band == "all" & delta == 0][order(lam_mult)]
cat("  delta alone (lam_mult = 1):\n")
cat(sprintf("    %s\n", paste(sprintf("%+0.1f:%+0.2f", b_delta$delta, b_delta$bias), collapse = "  ")))
cat("  lambda alone (delta = 0):\n")
cat(sprintf("    %s\n", paste(sprintf("%.2fx:%+0.2f", b_lam$lam_mult, b_lam$bias), collapse = "  ")))
cat(sprintf("\n  bias range from delta: %.2f deg | from lambda: %.2f deg\n",
            diff(range(b_delta$bias, na.rm = TRUE)), diff(range(b_lam$bias, na.rm = TRUE))))
cat("  the larger range is the knob that controls the bias.\n")
