# HOW DOES REAL LIGHT DEPART FROM THE SPIKE-AND-SLAB THE MODEL ASSUMES?
#
# WHY THIS IS THE TARGET. On the battery's `noisy` arm -- where the emission is
# correct by construction -- the engine gives latitude bias -0.51 deg, 189 km and
# coverage 0.996. On real data it gives -1.54 deg, ~244 km and coverage 0.64. The
# coverage failure does not reproduce AT ALL when the emission is right. So roughly a
# degree of bias and essentially the whole coverage problem is emission fidelity to
# real light. Everything else -- inference, geometry, movement prior, readout,
# calibration identifiability -- has been tested and is not it.
#
# WHAT THE MODEL ASSUMES, precisely (lib.rs `spike_density`, `LightMix::logden`):
#   y | mu ~ (1 - pslab) * AsymLaplace(mu; lam_lo below, lam_hi = shade_ratio*lam_lo
#            above, truncated to [0, max_light])  +  pslab * Uniform[0, max_light]
# with lam_lo = 1/(0.5*max_light) and shade_ratio = 2, CONSTANT over the whole record.
#
# Three consequences that are directly checkable against real light:
#   (a) P(y < mu) = lam_hi/(lam_lo + lam_hi) = 2/3 exactly, for shade_ratio = 2.
#   (b) Each side is EXPONENTIAL, so the mean absolute residual on the low side is
#       1/lam_lo and on the high side 1/lam_hi; their ratio recovers shade_ratio.
#   (c) None of this depends on zenith, depth or season -- the emission is assumed
#       identical everywhere.
#
# (c) IS THE ONE THAT MATTERS. Latitude is read from the shape of the light against
# zenith. If the residual distribution itself varies with zenith, the model is
# misspecified in a zenith-dependent way, and that maps straight into latitude. A
# diving seal shades the sensor in DAYLIGHT; at night there is no sky to shade. So
# there is every reason to expect (c) to fail, and the size and shape of the failure
# is what tells us what to replace the emission with.
#
# Residuals are reported against BOTH response forms, because they differ and the
# gap is itself part of the misspecification:
#   * the fitted LOGISTIC, the better description of the light;
#   * its TANGENT at half amplitude (slope = max_light/(4*scale),
#     zero = z50 + 2*scale), which is what the engine is actually handed.
#
# Truth is the battery's offset-0 positions (Argos, clipped per deployment). Real
# light comes from the decimated archives at the same timestamps. Aggregates only.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))
suppressMessages(devtools::load_all(".", quiet = TRUE))
BAT <- readRDS("scratch/nes_calibration/light_battery.rds")
REC <- Filter(function(z) z$offset == 0, BAT)

arch <- lapply(readRDS("analysis/cache/nes/archives_v3_30min.rds"), `[[`, "main")
fv <- try(lapply(readRDS("analysis/cache/nes/fvilches_archives_v1_30min.rds"),
                 `[[`, "main"), silent = TRUE)
if (!inherits(fv, "try-error")) arch <- c(arch, fv)
cat(sprintf("archives available: %d tags\n", length(arch)))

PSLAB <- 0.10; RATIO <- 2

R <- list()
for (b in REC) {
  g <- arch[[b$id]]
  if (is.null(g) || !nrow(g)) next
  ML <- b$max_light; r <- b$response
  gt <- as.numeric(g$time); bt <- as.numeric(b$time)
  y <- approx(gt, as.numeric(g$light), bt, rule = 2)$y
  dep <- if (!is.null(g$depth)) approx(gt, as.numeric(g$depth), bt, rule = 2)$y else NA_real_
  # THE ENGINE IS FED `pmax(0, light - baseline)`, with baseline the 5th percentile
  # of the tag's own light (fit_light_response returns q05). Comparing RAW light to
  # the response puts a constant positive offset in every residual -- an earlier
  # version of this script did exactly that and reported P(y < mu) = 0.057 against a
  # model value of 0.667. Reproduce what the engine actually sees.
  base <- as.numeric(quantile(as.numeric(g$light), 0.05, na.rm = TRUE))
  y <- pmax(0, y - base)
  ok <- b$supported & is.finite(y)
  if (sum(ok) < 2000) next
  tt <- bt[ok]
  z <- solar_zenith(tt, lon360(b$lon)[ok], b$lat[ok])
  mu_log <- pmin(pmax(r[1] + r[2] / (1 + exp((z - r[3]) / r[4])), 0), ML)
  slope <- ML / (4 * r[4]); zero <- r[3] + 2 * r[4]
  mu_tan <- pmin(pmax(slope * (zero - z), 0), ML)
  R[[b$id]] <- data.table(id = b$id, z = z, y = y[ok], ML = ML,
                          depth = if (length(dep) > 1) dep[ok] else NA_real_,
                          decl = solar_declination(as.POSIXct(tt, origin = "1970-01-01",
                                                              tz = "UTC")),
                          mu_log = mu_log, mu_tan = mu_tan)
}
D <- rbindlist(R)
cat(sprintf("%d tags, %d supported observations\n\n", uniqueN(D$id), nrow(D)))

D[, `:=`(r_log = (y - mu_log) / ML, r_tan = (y - mu_tan) / ML)]   # in units of max_light
LAM <- 2                       # lam_lo * max_light = 1/0.5 = 2, so scale = 1/2 of ML
cat(sprintf("MODEL SAYS: P(y < mu) = %.3f ; mean|resid| low = %.3f ML, high = %.3f ML ; ratio %.1f\n\n",
            RATIO / (1 + RATIO), 1 / LAM, 1 / (LAM * RATIO), RATIO))

emp <- function(x) {
  lo <- x[x < 0]; hi <- x[x > 0]
  data.table(n = length(x),
             p_below = mean(x < 0),
             mad_lo = if (length(lo)) mean(abs(lo)) else NA_real_,
             mad_hi = if (length(hi)) mean(hi) else NA_real_,
             ratio = if (length(lo) && length(hi)) mean(abs(lo)) / mean(hi) else NA_real_,
             far = mean(abs(x) > 0.5))          # tail mass the slab is meant to absorb
}

cat("=== (a,b) OVERALL, against each response form ===\n")
for (v in c("r_log", "r_tan")) {
  e <- emp(D[[v]])
  cat(sprintf("  %s : P(below) %.3f | mean|lo| %.3f | mean|hi| %.3f | implied ratio %.2f | P(|r|>0.5) %.3f\n",
              v, e$p_below, e$mad_lo, e$mad_hi, e$ratio, e$far))
}
cat(sprintf("  model : P(below) %.3f | mean|lo| %.3f | mean|hi| %.3f | ratio %.2f\n",
            RATIO / (1 + RATIO), 1 / LAM, 1 / (LAM * RATIO), RATIO))

cat("\n=== (c) THE ONE THAT MATTERS: does the residual depend on ZENITH? ===\n")
D[, zb := cut(z, c(0, 80, 88, 92, 96, 102, 180),
              labels = c("day <80", "80-88", "88-92", "92-96", "96-102", "night >102"))]
for (v in c("r_log", "r_tan")) {
  cat(sprintf("\n  -- %s --\n", v))
  t <- D[!is.na(zb), emp(get(v)), by = zb][order(zb)]
  t[, `:=`(p_below = round(p_below, 3), mad_lo = round(mad_lo, 3),
           mad_hi = round(mad_hi, 3), ratio = round(ratio, 2), far = round(far, 3))]
  print(as.data.frame(t), row.names = FALSE)
}
cat(sprintf("\n  model is CONSTANT across every row: P(below) %.3f, ratio %.1f\n",
            RATIO / (1 + RATIO), RATIO))

cat("\n=== does it depend on DEPTH? (the shading mechanism) ===\n")
if (any(is.finite(D$depth))) {
  D[, db := cut(depth, c(-Inf, 5, 50, 200, 500, Inf),
                labels = c("<5m", "5-50", "50-200", "200-500", ">500"))]
  t <- D[!is.na(db) & z < 92, emp(r_log), by = db][order(db)]   # daylight only
  t[, `:=`(p_below = round(p_below, 3), mad_lo = round(mad_lo, 3),
           mad_hi = round(mad_hi, 3), ratio = round(ratio, 2), far = round(far, 3))]
  cat("  (daylight only, z < 92)\n")
  print(as.data.frame(t), row.names = FALSE)
} else cat("  depth not available in the cached archives\n")

cat("\n=== and on SEASON? ===\n")
D[, sb := cut(decl, c(-24, -10, 10, 24), labels = c("NH winter", "equinox", "NH summer"))]
t <- D[!is.na(sb), emp(r_log), by = sb][order(sb)]
t[, `:=`(p_below = round(p_below, 3), mad_lo = round(mad_lo, 3),
         mad_hi = round(mad_hi, 3), ratio = round(ratio, 2), far = round(far, 3))]
print(as.data.frame(t), row.names = FALSE)

cat("\nREADING\n")
cat("  P(below) far from 0.667, or an implied ratio far from 2, means the assumed\n")
cat("  asymmetry is wrong. Those varying strongly with ZENITH is the serious finding:\n")
cat("  it means one constant emission cannot describe the record, and the error is\n")
cat("  organised along exactly the axis latitude is read from.\n")
