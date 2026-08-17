# Is lambda even a per-TAG quantity? Test stationarity within deployments.
#
# Everything so far has estimated ONE lambda per tag from the whole record. That
# assumes lambda is constant within a deployment -- which is assumed, never tested,
# and the start-month correlation (r = -0.49, p = 0.024) is exactly what you would
# see if it were FALSE. All 29 trips are the same post-moult migration departing
# May/June and returning Jan/Feb, so two deployments a month apart sample the same
# seasonal sequence at different phases. If lambda varies through the trip, their
# whole-record averages differ for that reason alone, and "per-tag lambda" is a
# smeared average of lambda(t) rather than a property of the tag.
#
# THE TEST. Split each deployment into calendar months and find the bias-zeroing
# lambda within each month. Then compare:
#     within-deployment variance of log2(lambda_opt)  vs  between-deployment variance
# If within >= between, lambda is not a per-tag constant and the entire per-tag
# framing is wrong. If within << between, the per-tag view survives and the start
# month effect is something else.
#
# SECOND QUESTION, which is the useful one if the first says non-stationary: does
# lambda_opt track a KNOWN function of time -- solar declination, day length, the
# animal's latitude? Those are all computable without Argos (declination exactly,
# the others from the fit itself), so a time-varying lambda would still be usable.
suppressMessages(library(data.table))
# run from the package root (setwd removed for portability)
D <- fread("scratch/nes_calibration/bias_lambda_time.csv")[band == "all"]
D[, wtime := as.POSIXct(wtime, origin = "1970-01-01", tz = "UTC")]
D[, mon := format(wtime, "%Y-%m")]
D[, doy := as.integer(format(wtime, "%j"))]

xing <- function(lam, y) {
  x <- log2(lam); o <- order(x); x <- x[o]; y <- y[o]
  k <- which(y[-1] * y[-length(y)] < 0); if (!length(k)) return(NA_real_)
  i <- k[1]; 2^(x[i] - y[i] * (x[i+1] - x[i]) / (y[i+1] - y[i]))
}

# per deployment-month, needing enough windows for the mean to mean anything
cnt <- D[lam_mult == 1, .N, by = .(id, mon)][N >= 4]
DM <- merge(D, cnt[, .(id, mon)], by = c("id", "mon"))
pm <- DM[, .(bias = mean(bias, na.rm = TRUE), nwin = .N), by = .(id, mon, lam_mult)]
lm_month <- pm[, .(lam_opt = xing(lam_mult, bias), nwin = max(nwin)), by = .(id, mon)]
lm_month <- lm_month[is.finite(lam_opt)]
lm_month[, `:=`(l2 = log2(lam_opt), doy = as.integer(format(as.Date(paste0(mon, "-15")), "%j")))]

cat(sprintf("deployment-months with an identifiable lambda: %d across %d deployments\n\n",
            nrow(lm_month), uniqueN(lm_month$id)))

keep <- lm_month[, .N, by = id][N >= 3]$id
G <- lm_month[id %in% keep]
cat(sprintf("=== stationarity: %d deployments with >=3 monthly estimates ===\n", uniqueN(G$id)))
if (uniqueN(G$id) >= 3) {
  a <- aov(l2 ~ factor(id), data = G)
  ss <- summary(a)[[1]][["Sum Sq"]]
  within_sd  <- sd(G[, l2 - mean(l2), by = id]$V1)
  between_sd <- sd(G[, mean(l2), by = id]$V1)
  cat(sprintf("  within-deployment  sd of log2(lambda): %.2f  (factor %.1f)\n",
              within_sd, 2^within_sd))
  cat(sprintf("  between-deployment sd of log2(lambda): %.2f  (factor %.1f)\n",
              between_sd, 2^between_sd))
  cat(sprintf("  between-deployment share of variance: %.0f%%\n", 100*ss[1]/sum(ss)))
  cat("\n  VERDICT: ")
  if (within_sd >= between_sd)
    cat("lambda is NOT a per-tag constant. Within-deployment variation matches or\n  exceeds between-deployment variation, so a single per-tag lambda is a smeared\n  average of lambda(t) and the per-tag framing is wrong.\n")
  else
    cat("the per-tag view survives: between-deployment variation dominates.\n")
  cat(sprintf("\n  per-deployment range of lambda_opt (median across deployments): factor %.1f\n",
              2^median(G[, diff(range(l2)), by = id]$V1)))
}

cat("\n=== does lambda track time of year rather than the tag? ===\n")
if (nrow(lm_month) >= 8) {
  decl <- function(doy) {
    g <- 2*pi/365 * (doy - 1)
    (0.006918 - 0.399912*cos(g) + 0.070257*sin(g) - 0.006758*cos(2*g) +
       0.000907*sin(2*g) - 0.002697*cos(3*g) + 0.00148*sin(3*g)) * 180/pi
  }
  lm_month[, dec := decl(doy)]
  for (v in c("doy", "dec")) {
    ct <- suppressWarnings(cor.test(lm_month[[v]], lm_month$l2))
    cat(sprintf("  cor(%-4s, log2 lambda_opt) = %+.3f, p = %.4f, n = %d\n",
                v, ct$estimate, ct$p.value, nrow(lm_month)))
  }
  # within-deployment only: removes any between-tag confound entirely
  W <- lm_month[id %in% keep][, .(dec = dec - mean(dec), l2 = l2 - mean(l2)), by = id]
  ct <- suppressWarnings(cor.test(W$dec, W$l2))
  cat(sprintf("  WITHIN-deployment (both centred): cor = %+.3f, p = %.4f, n = %d\n",
              ct$estimate, ct$p.value, nrow(W)))
  cat("  the within-deployment correlation is the clean one: it cannot be produced\n")
  cat("  by tags differing from each other, only by lambda moving through the year.\n")
  cat("\n  monthly medians of lambda_opt:\n")
  mm <- lm_month[, .(n = .N, lam = round(2^median(l2), 2)), by = .(m = substr(mon, 6, 7))][order(m)]
  cat(sprintf("    %s\n", paste(sprintf("%s:%.2f(n=%d)", month.abb[as.integer(mm$m)], mm$lam, mm$n),
                                collapse = "  ")))
}
