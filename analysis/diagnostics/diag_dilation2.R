# The declination sign-flip idea is dead: at a self-crossing the latitude error
# does not track the declination difference at all (r = -0.001, p = 0.88). So
# the dilation is not a signed calibration error flipping either side of the
# equinox. Three better candidates, all measurable:
#
#   CENTRE PULL. At the equinox day length is ~12 h everywhere, so latitude is
#   nearly unidentifiable and the posterior mean relaxes toward the middle of
#   whatever search domain it was given. The grid here spans 29-65 N, centre 47.
#   Mean truth at the equinox is 45.0 and the mean estimate 47.5, which is
#   suspiciously close to the centre. PREDICTION: |estimate - 47| < |truth - 47|
#   near the equinox and not at the solstices, and the effect scales with how
#   far truth sits from the centre.
#
#   NOISE INFLATION. Estimated latitude is truth plus error of sd ~3 deg. The
#   RANGE of a noisy series exceeds the range of a clean one, so any loop looks
#   fatter without any bias at all. PREDICTION: latitude span of the estimate
#   exceeds truth's at matched longitude, by roughly the error sd.
#
#   LIMB ASYMMETRY. The thing actually reported: outbound displaced one way and
#   return the other. PREDICTION: within a tag, the two limbs carry
#   opposite-signed latitude bias.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
suppressMessages({ library(data.table) })

declination <- function(t) {
  doy <- as.numeric(format(t, "%j"))
  hr  <- as.numeric(format(t, "%H")) + as.numeric(format(t, "%M")) / 60
  g <- 2*pi/365 * (doy - 1 + (hr - 12)/24)
  (0.006918 - 0.399912*cos(g) + 0.070257*sin(g) - 0.006758*cos(2*g) +
     0.000907*sin(2*g) - 0.002697*cos(3*g) + 0.00148*sin(3*g)) * 180/pi
}

k <- as.data.table(read.csv("analysis/output/nes/per_knot_errors.csv"))
k[, time := as.POSIXct(time, tz = "UTC")]
k <- k[engine == "hier" & batch == "light" & is.finite(err_lat) & is.finite(true_lat)]
k[, decl := declination(time)]
setorder(k, id, time)
GRID_MID <- (29 + 65) / 2

cat("=== 1. centre pull: does the estimate move toward the middle of the box? ===\n")
k[, `:=`(d_true = abs(true_lat - GRID_MID), d_est = abs(est_lat - GRID_MID))]
k[, abin := cut(abs(decl), breaks = c(0, 5, 10, 15, 20, 24))]
print(as.data.frame(k[!is.na(abin), .(knots = .N,
  mean_true_lat = round(mean(true_lat), 1), mean_est_lat = round(mean(est_lat), 1),
  dist_true = round(mean(d_true), 2), dist_est = round(mean(d_est), 2),
  pulled_in = round(mean(d_est < d_true), 2)),
  by = abin][order(abin)]), row.names = FALSE)
cat("\n`pulled_in` is the fraction of knots where the ESTIMATE sits closer to the\n")
cat("box centre than the TRUTH does. 0.5 is no pull.\n")

# regression of the estimate on truth: slope < 1 is shrinkage toward a centre,
# slope > 1 is dilation away from it
cat("\n=== 2. est_lat regressed on true_lat, by |declination| ===\n")
print(as.data.frame(k[!is.na(abin), {
  f <- lm(est_lat ~ true_lat)
  .(knots = .N, slope = round(coef(f)[2], 3),
    intercept = round(coef(f)[1], 2),
    fixed_point = round(coef(f)[1] / (1 - coef(f)[2]), 1))
}, by = abin][order(abin)]), row.names = FALSE)
cat("slope < 1 means the estimate is COMPRESSED toward `fixed_point`;\n")
cat("slope > 1 means it is DILATED away from it.\n")

cat("\n=== 3. noise inflation: latitude span at matched longitude ===\n")
sp <- k[, {
  b <- cut(true_lon, breaks = seq(160, 265, by = 5))
  d <- data.table(b = b, tl = true_lat, el = est_lat)[!is.na(b)]
  s <- d[, .(span_true = diff(range(tl)), span_est = diff(range(el)), n = .N),
         by = b][n >= 5]
  .(span_true = mean(s$span_true), span_est = mean(s$span_est), bins = nrow(s))
}, by = id]
sp[, ratio := round(span_est / span_true, 2)]
print(as.data.frame(sp[, .(id, bins, span_true = round(span_true, 2),
                           span_est = round(span_est, 2), ratio)]), row.names = FALSE)
cat(sprintf("\nmean span ratio %.2f -- above 1 means the loop is drawn fatter than it is\n",
            mean(sp$ratio)))

cat("\n=== 4. limb asymmetry: outbound against return, within tag ===\n")
lb <- k[, {
  # the turning point is the westernmost knot of the TRUE track
  turn <- which.min(true_lon)
  limb <- ifelse(seq_len(.N) <= turn, "outbound", "return")
  .(limb = limb, err_lat = err_lat, true_lat = true_lat, decl = decl)
}, by = id]
tab <- lb[, .(knots = .N, mean_true_lat = round(mean(true_lat), 1),
              bias = round(mean(err_lat), 2)), by = .(id, limb)]
w <- dcast(tab, id ~ limb, value.var = c("mean_true_lat", "bias"))
w[, opposite_sign := sign(bias_outbound) != sign(bias_return)]
w[, dilation := round(bias_outbound - bias_return, 2)]
print(as.data.frame(w), row.names = FALSE)
cat(sprintf("\n%d of %d tags carry opposite-signed bias on the two limbs\n",
            sum(w$opposite_sign), nrow(w)))
cat(sprintf("mean outbound-minus-return bias: %.2f deg\n", mean(w$dilation)))
cat("Positive means the outbound limb is displaced north relative to the return,\n")
cat("i.e. the loop is opened out rather than shifted bodily.\n")
