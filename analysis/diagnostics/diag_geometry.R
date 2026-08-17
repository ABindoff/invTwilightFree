# Is the "flat" relationship between the movement scale and accuracy a real
# regime, or was it an accident of one configuration?
#
# ABSTRACT SETUP. Near a knot the light likelihood is approximately Gaussian with
# covariance diag(tau_lon^2, tau_lat^2); the prior is an isotropic random walk
# with per-step sd sigma. That is a linear-Gaussian smoother, whose effective
# bandwidth is b ~ tau/sigma knots. Three predictions follow:
#
#   1. ANISOTROPY. One isotropic sigma gives DIFFERENT bandwidths per axis,
#      because tau does. b_lon = tau_lon/sigma, b_lat = tau_lat/sigma.
#      With tau_lat >> tau_lon the prior smooths latitude and barely touches
#      longitude -- which is why longitude intervals are conservative and
#      latitude ones are not.
#   2. WIDTH. Averaging b knots gives sd ~ tau/sqrt(b) = sqrt(tau*sigma), so the
#      reported posterior sd should scale as sqrt(sigma).
#   3. ACCURACY. Write the error as bias + independent noise. Smoothing divides
#      the noise by sqrt(b) and leaves the bias alone, while the model reports a
#      width reduced by sqrt(b) on both. So accuracy is flat in sigma IFF the
#      bias term dominates -- i.e. iff the error is correlated over lags shorter
#      than the bandwidth.
#
# Prediction 3 is the one in question and it is decidable: measure the
# autocorrelation of the error. If it decays over more knots than the bandwidth,
# smoothing cannot help and "flat" is real. If it decays faster, smoothing should
# help and my claim was wrong.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))

k <- as.data.table(read.csv("analysis/output/nes/per_knot_errors.csv"))
k <- k[engine == "hier" & batch == "light" & is.finite(err_lat) & is.finite(err_lon)]
k[, `:=`(id = as.character(id), time = as.POSIXct(time, tz = "UTC"))]
setorder(k, id, time)
KM_PER_DEG <- 111.32

cat("=== observed error scales (deg, and km at 45 N) ===\n")
sc <- k[, .(sd_lat = round(sd(err_lat), 2), sd_lon = round(sd(err_lon), 2),
            km_lat = round(sd(err_lat) * KM_PER_DEG),
            km_lon = round(sd(err_lon) * KM_PER_DEG * cos(45*pi/180))), by = id]
print(as.data.frame(sc), row.names = FALSE)
cat(sprintf("\npooled: latitude %.0f km, longitude %.0f km -> anisotropy %.1fx\n",
            sd(k$err_lat) * KM_PER_DEG,
            sd(k$err_lon) * KM_PER_DEG * cos(45*pi/180),
            (sd(k$err_lat)) / (sd(k$err_lon) * cos(45*pi/180))))

cat("\n=== autocorrelation of the error, per tag (lag in 12-h knots) ===\n")
acf_tab <- k[, {
  al <- acf(err_lat, lag.max = 20, plot = FALSE)$acf[, , 1]
  ao <- acf(err_lon, lag.max = 20, plot = FALSE)$acf[, , 1]
  # correlation length: first lag at which the acf drops below 1/e
  cl <- which(al < exp(-1))[1] - 1
  co <- which(ao < exp(-1))[1] - 1
  .(lat_lag1 = round(al[2], 2), lat_lag4 = round(al[5], 2), lat_lag8 = round(al[9], 2),
    lat_corr_len = cl,
    lon_lag1 = round(ao[2], 2), lon_lag4 = round(ao[5], 2), lon_corr_len = co)
}, by = id]
print(as.data.frame(acf_tab), row.names = FALSE)

cat("\n=== the test ===\n")
sigma_km <- 110 * sqrt(0.5)           # D = 110 km/sqrt(day), 12-h step
tau_lat  <- sd(k$err_lat) * KM_PER_DEG
tau_lon  <- sd(k$err_lon) * KM_PER_DEG * cos(45*pi/180)
cat(sprintf("sigma (D = 110, 12 h)      : %.0f km\n", sigma_km))
cat(sprintf("bandwidth b_lat = tau/sigma: %.1f knots\n", tau_lat / sigma_km))
cat(sprintf("bandwidth b_lon = tau/sigma: %.1f knots\n", tau_lon / sigma_km))
cat(sprintf("median latitude correlation length: %.0f knots\n",
            median(acf_tab$lat_corr_len, na.rm = TRUE)))
cat(sprintf("median longitude correlation length: %.0f knots\n",
            median(acf_tab$lon_corr_len, na.rm = TRUE)))
cat("\nIf the latitude correlation length EXCEEDS b_lat, the smoother is averaging\n")
cat("over a window inside which the error barely changes: it cannot cancel, so\n")
cat("accuracy is flat in sigma while the reported width still shrinks as\n")
cat("sqrt(sigma). That is the regime, and it is what makes coverage degrade as\n")
cat("the prior tightens.\n")

cat("\n=== how much of the latitude error could ANY amount of smoothing remove? ===\n")
# variance of the knot-level error, split into a slowly varying part (a moving
# average over the bandwidth) and the residual the smoother could cancel
for (b in c(2, 4, 8, 16)) {
  v <- k[, {
    s <- stats::filter(err_lat, rep(1/b, b), sides = 2)
    .(slow = var(s, na.rm = TRUE), tot = var(err_lat))
  }, by = id]
  cat(sprintf("  bandwidth %2d knots: %.0f%% of latitude error variance survives averaging\n",
              b, 100 * mean(v$slow / v$tot)))
}
