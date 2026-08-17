# Is the per-logger Argos-optimal lambda actually ESTIMABLE?
#
# The dawn/dusk criterion gives a stable per-logger lambda (15/15 crossings, spread
# only 2.7x, sd log2 0.51) but correlates -0.42 with the Argos-optimal lambda, whose
# spread is 256x. Two readings, with opposite consequences:
#
#   A. The criterion genuinely fails per tag. Then no constant offset can rescue it.
#   B. The TARGET is noise. If each logger's bias(lambda) curve is shallow near its
#      zero crossing, the crossing is barely identified, the 256x spread is sampling
#      noise, and the low correlation is attenuation from measurement error in the
#      target rather than a failure of the criterion. Then a global offset applied to
#      the dawn/dusk lambda may be exactly right.
#
# Distinguish them by bootstrapping the crossing. Each logger's bias at a given
# lambda is a mean over knot windows, so resample WINDOWS within logger and recompute
# the crossing. A tight bootstrap interval means the target is real (reading A); a
# wide one means it is not (reading B).
#
# Also report the slope d(bias)/d(log2 lambda) at the crossing: a shallow slope IS
# the mechanism for a poorly determined crossing, and it says how much bias a
# factor-of-2 error in lambda actually costs -- which is the practically important
# number regardless of which reading wins.
suppressMessages(library(data.table))
# run from the package root (setwd removed for portability)
set.seed(1)
W <- fread("scratch/nes_calibration/bias_lambda_wide.csv")[band == "all" & !is.na(serial)]
W[, win := rleid(id, bias)]                       # a row per (window, lambda)

xing <- function(lam, y) {
  x <- log2(lam); o <- order(x); x <- x[o]; y <- y[o]
  k <- which(y[-1] * y[-length(y)] < 0)
  if (!length(k)) return(c(NA_real_, NA_real_))
  i <- k[1]
  sl <- (y[i+1] - y[i]) / (x[i+1] - x[i])          # deg of bias per doubling
  c(2^(x[i] - y[i] / sl), sl)
}

B <- 300
res <- rbindlist(lapply(unique(W$serial), function(L) {
  d <- W[serial == L]
  keys <- unique(d[, .(id, win)])
  pt <- d[, .(bias = mean(bias, na.rm = TRUE)), by = lam_mult]
  o  <- xing(pt$lam_mult, pt$bias)
  bs <- replicate(B, {
    s <- keys[sample.int(nrow(keys), nrow(keys), replace = TRUE)]
    dd <- merge(d, s[, .N, by = .(id, win)], by = c("id", "win"))
    m <- dd[, .(bias = weighted.mean(bias, N, na.rm = TRUE)), by = lam_mult]
    xing(m$lam_mult, m$bias)[1]
  })
  bs <- bs[is.finite(bs)]
  data.table(serial = L, n_win = nrow(keys), lam_hat = o[1], slope = o[2],
             lo = if (length(bs) > 20) quantile(bs, .1) else NA_real_,
             hi = if (length(bs) > 20) quantile(bs, .9) else NA_real_,
             frac_crossing = length(bs) / B)
}))

cat("=== is each logger's Argos lambda identified? (bootstrap over knot windows) ===\n")
cat(sprintf("%-12s %6s %9s %9s %20s %8s\n",
            "logger", "n_win", "lam_hat", "slope", "80% interval", "cross%"))
for (i in seq_len(nrow(res)))
  cat(sprintf("%-12s %6d %9s %9s %20s %8.0f\n", res$serial[i], res$n_win[i],
              if (is.na(res$lam_hat[i])) "--" else sprintf("%.2f", res$lam_hat[i]),
              if (is.na(res$slope[i])) "--" else sprintf("%+.2f", res$slope[i]),
              if (is.na(res$lo[i])) "--" else sprintf("%.2f - %.2f", res$lo[i], res$hi[i]),
              100 * res$frac_crossing[i]))

g <- res[is.finite(lam_hat) & is.finite(lo)]
if (nrow(g) >= 3) {
  wid <- log2(g$hi) - log2(g$lo)
  cat(sprintf("\n  median 80%% interval width: %.2f in log2, i.e. a factor of %.1f\n",
              median(wid), 2^median(wid)))
  cat(sprintf("  observed between-logger spread: %.2f in log2 (factor %.0f)\n",
              sd(log2(g$lam_hat)), 2^(2*sd(log2(g$lam_hat)))))
  cat(sprintf("  median |slope| at the crossing: %.2f deg of bias per doubling of lambda\n",
              median(abs(g$slope), na.rm = TRUE)))
  cat("\n  VERDICT:\n")
  if (median(wid) > 2 * sd(log2(g$lam_hat)))
    cat("  the per-logger crossing is NOT identified -- its uncertainty swamps the\n  between-logger spread, so the 256x range is mostly noise and the low\n  correlation with dawn/dusk is attenuation, not failure.\n")
  else if (median(wid) < sd(log2(g$lam_hat)))
    cat("  the per-logger crossing IS identified -- real between-logger variation\n  exceeds its uncertainty, so the dawn/dusk criterion genuinely does not\n  track it and no constant offset will fix that.\n")
  else
    cat("  borderline: uncertainty and between-logger spread are comparable, so the\n  target is only partly estimable and any offset will be weakly determined.\n")
}

cat("\n=== practical cost of getting lambda wrong ===\n")
pt <- W[, .(bias = mean(bias, na.rm = TRUE)), by = lam_mult][order(lam_mult)]
cat("pooled bias vs lambda (deg):\n")
cat(sprintf("  %s\n", paste(sprintf("%.2fx:%+0.2f", pt$lam_mult, pt$bias), collapse = "  ")))
cat(sprintf("\n  a factor-of-2 error in lambda costs roughly %.2f deg of latitude bias\n",
            median(abs(res$slope), na.rm = TRUE)))
