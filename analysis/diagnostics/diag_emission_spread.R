# EMISSION SPREAD: where the assumed spike-and-slab departs from real light.
#
# WHY THIS AND NOT ANOTHER RESPONSE. Coverage is the last unexplained failure and no
# response shape touches it: across the three families tested at fit level it runs
# 0.593 (tangent) / 0.338 (dark regime) / 0.077 (gompertz) and gets WORSE as bias
# improves. On the battery's noisy arm, where the emission is correct by
# construction, coverage is 0.996. So the coverage failure is the emission's SPREAD,
# not its centre and not its shape.
#
# WHAT THE MODEL ASSUMES (lib.rs `spike_density` / `LightMix`):
#   y | mu ~ (1-p) * AsymLaplace(mu; lam below, shade_ratio*lam above,
#            truncated to [0, max_light])  +  p * Uniform[0, max_light]
# with lam = 1/(0.5*max_light), shade_ratio = 2, p = 0.10, CONSTANT everywhere.
# Consequences: mean |residual| is 1/lam below and 1/(shade_ratio*lam) above, so
# their ratio recovers shade_ratio; and tail mass beyond the spike recovers p.
#
# MEASURED AGAINST THE PIPELINE'S OWN CALIBRATION, not the battery's. Uses
# `panel_v1.rds`: pooled-within-family haul-out geometry, per-tag 15-day scale,
# tangent tied to max_light, light baselined as `pmax(0, light - r$baseline)`,
# truth from `argos_at` (the strict scorer -- NA across Argos gaps over 24 h).
#
# CLAMPING IS FLAGGED, NOT IGNORED. The tangent clamps to 0 beyond zenith ~106, so
# residuals there CANNOT be negative and any asymmetry statistic is an artefact. An
# earlier version of this measurement reported a ratio of 0.27 at night without
# noting that. Bands where the expected curve is clamped are marked and excluded
# from the recommendation.
#
# THE POINT is to state what the engine would have to express, and whether it can.
# The engine has ONE global `shade_ratio`, ONE global `prob_slab`, and ONE darkness
# regime (a single threshold `dark_frac`, below which `shade_ratio_dark` and
# `prob_slab_dark` take over). If the requirement is monotone in zenith, that is
# expressible; if it needs two thresholds in opposite directions, it is not.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))
suppressMessages(devtools::load_all(".", quiet = TRUE))

P <- readRDS("analysis/cache/nes/panel_v1.rds")
tags <- P$tags; responses <- P$responses
SHADE <- 2; PSLAB <- 0.10

R <- list()
for (tg in names(tags)) {
  g <- tags[[tg]]; r <- responses[[tg]]
  if (is.null(r)) next
  ML <- r$max_light
  y <- pmax(0, as.numeric(g$light$light) - r$baseline)
  tm <- g$light$time
  tr <- argos_at(as.data.table(g$argos), tm)     # strict: NA across >24 h gaps
  ok <- is.finite(tr$lat) & is.finite(y)
  if (sum(ok) < 2000) next
  z  <- solar_zenith(as.numeric(tm)[ok], tr$lon[ok], tr$lat[ok])
  mu <- pmin(pmax(r$calibration[1] - r$calibration[2] * z, 0), ML)   # the tangent
  R[[tg]] <- data.table(id = tg, z = z, y = y[ok], mu = mu, ML = ML,
                        clamped = mu <= 1e-9 | mu >= ML - 1e-9)
}
D <- rbindlist(R)
D[, res := (y - mu) / ML]                        # residual in units of max_light
cat(sprintf("%d tags, %d observations with supported truth\n\n", uniqueN(D$id), nrow(D)))

emp <- function(x) {
  lo <- x[x < 0]; hi <- x[x > 0]
  list(p_below = mean(x < 0),
       mad_lo = if (length(lo)) mean(abs(lo)) else NA_real_,
       mad_hi = if (length(hi)) mean(hi) else NA_real_,
       ratio  = if (length(lo) && length(hi)) mean(abs(lo)) / mean(hi) else NA_real_,
       tail   = mean(abs(x) > 0.5))
}

cat("MODEL ASSUMES, in every band:\n")
cat(sprintf("  P(below) %.3f | mad_lo %.3f | mad_hi %.3f | ratio %.2f | tail(slab) %.2f\n\n",
            SHADE/(1+SHADE), 1/2, 1/(2*SHADE), SHADE, PSLAB))

D[, zb := cut(z, c(0, 80, 86, 90, 94, 98, 102, 180),
              labels = c("day <80","80-86","86-90","90-94","94-98","98-102","night >102"))]
T <- D[!is.na(zb), c(list(n = .N, clamped = round(mean(clamped), 2)), emp(res)), by = zb][order(zb)]
T[, `:=`(p_below = round(p_below, 3), mad_lo = round(mad_lo, 3), mad_hi = round(mad_hi, 3),
         ratio = round(ratio, 2), tail = round(tail, 3),
         lam_lo_impl = round(1 / mad_lo, 2), lam_hi_impl = round(1 / mad_hi, 2))]
cat("=== measured, by zenith (residual in units of max_light) ===\n")
print(as.data.frame(T[, .(zb, n, clamped, p_below, mad_lo, mad_hi, ratio, tail)]),
      row.names = FALSE)
cat("\n  `clamped` = fraction of the band where the tangent is at 0 or max_light.\n")
cat("  Where that is high the asymmetry statistics are ARTEFACTS -- residuals cannot\n")
cat("  take one sign -- and the row must not be used for tuning.\n")

use <- T[clamped < 0.20]
cat("\n=== usable bands only (clamped < 0.20) ===\n")
print(as.data.frame(use[, .(zb, n, ratio, implied_shade_ratio = ratio,
                            implied_lam_x = round(0.5 / ((mad_lo + mad_hi) / 2), 2),
                            tail_vs_pslab = round(tail / PSLAB, 2))]), row.names = FALSE)
cat(sprintf("\n  shipped shade_ratio %.1f ; measured spans %.2f to %.2f across usable bands\n",
            SHADE, min(use$ratio, na.rm = TRUE), max(use$ratio, na.rm = TRUE)))
cat(sprintf("  shipped lambda is 1/(0.5*ML); measured implies %.2f to %.2f times that\n",
            min(0.5 / ((use$mad_lo + use$mad_hi) / 2)), max(0.5 / ((use$mad_lo + use$mad_hi) / 2))))
cat(sprintf("  shipped prob_slab %.2f ; measured tail mass %.3f to %.3f\n",
            PSLAB, min(use$tail), max(use$tail)))

cat("\n=== can the engine express what is required? ===\n")
cat("  It has ONE global shade_ratio, ONE global prob_slab, and ONE darkness\n")
cat("  threshold below which both switch to a second value.\n")
mono <- all(diff(use$ratio[!is.na(use$ratio)]) <= 0) || all(diff(use$ratio[!is.na(use$ratio)]) >= 0)
cat(sprintf("  required shade_ratio monotone in zenith across usable bands: %s\n", mono))
cat(sprintf("  => %s\n", if (mono)
  "expressible by the existing darkness regime (one threshold, one alternative)." else
  "NOT expressible: it would need the ratio to move in both directions, i.e. two\n     thresholds. Any single-threshold tuning can only fit part of it."))
