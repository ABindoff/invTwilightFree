# The bias is structured ALONG the track, not accumulating. What is it tracking?
#
# The filtered/smoothed profiles are near zero for the first 40% of the track, dive
# to about -1.9 deg in the middle, and partly recover. Two candidates, and they are
# confounded in these records because the seals go north and come back over the
# same months:
#
#   LATITUDE. A per-knot TILT of the marginal by cos(latitude) shifts the mean by
#     about -sigma_post^2 * tan(lat), so it grows with latitude. The animals are at
#     their northernmost mid-trip, which is exactly where the bias peaks. This
#     predicts bias ~ tan(lat), and more sharply bias ~ sigma_post^2 * tan(lat).
#
#   SEASON. Offset-0 records run May to December, so mid-track is near the
#     September equinox, where day length carries least latitude information.
#
# These are separable HERE, because the battery holds four date offsets per track:
# the same movement placed in four seasons. If the bias follows latitude it will
# look the same at every offset; if it follows season it will move with the offset.
# This script uses the offset-0 fits already on disk and states what the offset
# comparison would add.
suppressMessages(library(data.table))
D <- fread("scratch/nes_calibration/filtered_vs_smoothed.csv",
           colClasses = list(character = "id"))
D <- D[sup == TRUE & is.finite(smoothed) & is.finite(truth)]
D[, e := smoothed - truth]
D[, tanlat := tan(truth * pi / 180)]

cat(sprintf("%d knots, %d tracks, drift arms: %s\n\n",
            nrow(D[drift == FALSE]), length(unique(D$id)),
            paste(unique(D$drift), collapse = "/")))

X <- D[drift == FALSE]

cat("=== bias by TRUE latitude band ===\n")
X[, band := cut(truth, c(-Inf, 40, 45, 50, 55, Inf),
                labels = c("<40", "40-45", "45-50", "50-55", ">55"))]
b <- X[, .(n = .N, lat = round(mean(truth), 1), bias = round(mean(e), 3),
           tanlat = round(mean(tanlat), 3)), by = band][order(band)]
print(as.data.frame(b), row.names = FALSE)

cat("\n=== bias by position along track (for comparison) ===\n")
X[, q := cut(frac, seq(0, 1, 0.2), include.lowest = TRUE, labels = FALSE)]
p <- X[, .(n = .N, lat = round(mean(truth), 1), bias = round(mean(e), 3)), by = q][order(q)]
print(as.data.frame(p), row.names = FALSE)

cat("\n=== which predicts the bias: latitude, or position along the track? ===\n")
m1 <- lm(e ~ tanlat, data = X)
m2 <- lm(e ~ frac,   data = X)
m3 <- lm(e ~ tanlat + frac, data = X)
cat(sprintf("  e ~ tan(lat)        : R2 %.4f, slope %+0.3f\n",
            summary(m1)$r.squared, coef(m1)[2]))
cat(sprintf("  e ~ frac            : R2 %.4f, slope %+0.3f\n",
            summary(m2)$r.squared, coef(m2)[2]))
cat(sprintf("  e ~ tan(lat) + frac : R2 %.4f | tan(lat) %+0.3f (p=%.2g), frac %+0.3f (p=%.2g)\n",
            summary(m3)$r.squared,
            coef(m3)[2], summary(m3)$coefficients[2, 4],
            coef(m3)[3], summary(m3)$coefficients[3, 4]))

cat("\n=== within-track, to remove between-animal differences ===\n")
for (tg in sort(unique(X$id))) {
  x <- X[id == tg]
  r <- suppressWarnings(cor(x$e, x$tanlat, method = "spearman"))
  cat(sprintf("  %-9s n=%4d  lat %.1f-%.1f  bias %+0.2f  rho(e, tan lat) %+0.2f\n",
              tg, nrow(x), min(x$truth), max(x$truth), mean(x$e), r))
}

cat("\nNOTE the confound: in these offset-0 records the animals are at their\n")
cat("northernmost mid-trip, which is also near the September equinox. Latitude and\n")
cat("season cannot be separated from offset 0 alone. The battery's other three date\n")
cat("offsets place the SAME movement in different seasons and would separate them;\n")
cat("that is 18 further fits and is the obvious next run if tan(lat) looks strong.\n")
