# DOES SOLAR GEOLOCATION RECOVER HABITAT OCCUPANCY, EVEN WHEN POSITIONS ARE POOR?
#
# THE QUESTION THIS ACTUALLY ANSWERS. Per-knot position error is the wrong metric
# for the science these tags are deployed to do. Nobody tracks an elephant seal to
# learn where it was on a Tuesday; they track it to learn which water it used and
# when. So the assessment that belongs in the manuscript is whether the OCCUPANCY
# DISTRIBUTION -- proportion of time in each latitude band, by month -- is
# recovered, not whether each fix is close.
#
# Latitude is the right axis to test on. It carries 230 of the engine's 243 km of
# error, and North Pacific habitat structure is largely zonal (the transition zone
# chlorophyll front, the subarctic boundary), so it is simultaneously the worst
# coordinate and the one the ecology is organised along.
#
# THE PART OF THE ARGUMENT THAT DOES NOT HOLD, and which this measures directly:
# random error averages out over animals, bias does not. Worse, the calibration
# here is POOLED across tags by construction, so every Mk9 shares one z50 and their
# errors are CORRELATED rather than independent. Pooling 100 tracks divides the
# random part by 10 and the shared part by 1. If occupancy is recovered it will be
# because the bias is small against the habitat scale, not because n was large --
# and those two claims need to be distinguished in print.
#
# Metrics:
#   overlap  sum_b min(p_true[b], p_est[b]), the standard occupancy overlap, 1 = identical
#   per-tag  computed within tag, then summarised
#   pooled   computed on the population distribution (all knots from all tags)
# If pooling helps, pooled overlap exceeds the per-tag median. If the error is a
# shared bias, it will not.
suppressMessages({ library(data.table) })

K <- fread("scratch/nes_calibration/final29_knots.csv", colClasses = list(character = "id"))
K <- K[arm == "A_shipped" & is.finite(true_lat) & is.finite(err_lat)]
K[, est_lat := true_lat + err_lat]
K[, month := as.integer(format(as.POSIXct(time, tz = "UTC"), "%m"))]
cat(sprintf("%d knots, %d tags\n", nrow(K), uniqueN(K$id)))
cat(sprintf("latitude range: Argos %.1f to %.1f, light %.1f to %.1f\n\n",
            min(K$true_lat), max(K$true_lat), min(K$est_lat), max(K$est_lat)))

overlap <- function(a, b, brk) {
  pa <- prop.table(table(cut(a, brk))); pb <- prop.table(table(cut(b, brk)))
  sum(pmin(pa, pb))
}

for (W in c(5, 2)) {
  brk <- seq(10, 80, by = W)
  cat(sprintf("=== %d-degree latitude bands ===\n", W))
  per <- K[, .(ov = overlap(true_lat, est_lat, brk), n = .N), by = id]
  pool <- overlap(K$true_lat, K$est_lat, brk)
  cat(sprintf("  per-tag overlap  (q25|med|q75): %s\n",
              paste(round(quantile(per$ov, c(.25, .5, .75)), 3), collapse = " | ")))
  cat(sprintf("  POOLED overlap                : %.3f\n", pool))
  cat(sprintf("  pooling gain over per-tag median: %+.3f  -> %s\n",
              pool - median(per$ov),
              if (pool > median(per$ov) + 0.02) "pooling helps (error partly independent)"
              else "pooling does NOT help (error is a shared bias)"))
  cat("\n")
}

cat("=== population occupancy by 5-degree band, Argos vs light ===\n")
brk <- seq(10, 80, by = 5)
K[, b_true := cut(true_lat, brk)][, b_est := cut(est_lat, brk)]
tb <- merge(K[, .(argos = .N / nrow(K)), by = .(band = b_true)],
            K[, .(light = .N / nrow(K)), by = .(band = b_est)], by = "band", all = TRUE)
tb[is.na(tb)] <- 0
setorder(tb, band)
tb <- tb[argos > 0.005 | light > 0.005]
print(as.data.frame(tb[, .(band, argos = round(argos, 3), light = round(light, 3),
                           diff = round(light - argos, 3))]), row.names = FALSE)
cat(sprintf("\n  population mean latitude: Argos %.2f, light %.2f, shift %+.2f deg (%+.0f km)\n",
            mean(K$true_lat), mean(K$est_lat), mean(K$est_lat) - mean(K$true_lat),
            (mean(K$est_lat) - mean(K$true_lat)) * 111.19))

cat("\n=== seasonal: overlap of the month x band occupancy surface ===\n")
mb <- function(lat, mon) prop.table(table(cut(lat, brk), mon))
pa <- mb(K$true_lat, K$month); pb <- mb(K$est_lat, K$month)
cat(sprintf("  month x 5-deg-band overlap: %.3f\n", sum(pmin(pa, pb))))
ms <- K[, .(argos = mean(true_lat), light = mean(est_lat), n = .N), by = month][order(month)]
ms[, shift := light - argos]
cat("\n  monthly mean latitude:\n")
print(as.data.frame(ms[, .(month, n, argos = round(argos, 1), light = round(light, 1),
                           shift = round(shift, 2))]), row.names = FALSE)
cat("\n  a CONSTANT shift across months is a calibration bias and is correctable;\n")
cat("  one that varies with season is not, and would distort seasonal ecology.\n")
