# Is the SEASONAL swing the max_light CLAMP on the expected curve?
#
# Censoring nearly zeroes the mean tilt (-0.0548 -> +0.0052) but leaves the seasonal
# swing intact (0.2069 -> 0.2360). A seasonally varying bias is worse than a constant
# one for habitat and climate inference, because it manufactures apparent seasonal
# movement, so this is the part that has to go.
#
# THE INCOHERENCE. `expected_light` clamps mu to [0, max_light] (lib.rs:221). Under
# TRUNCATION it has to: the normaliser Z(mu) is the mass over [0, ml] and is not
# defined for mu outside it. But the clamp is physically wrong -- true daylight
# routinely EXCEEDS a sensor's ceiling; the sensor just reports the ceiling. And the
# clamped FRACTION of the daily curve varies with maximum sun elevation, which is a
# seasonal quantity. That is a seasonal, latitude-dependent distortion of the
# expected curve, sitting exactly where the tilt is made (daylight, not twilight).
#
# Under CENSORING the support is unbounded, so mu may exceed max_light freely and the
# atom at the ceiling absorbs it correctly. Censoring therefore MAKES UNCLAMPING
# AVAILABLE, which truncation does not.
#
# PREDICTION, WRITTEN BEFORE RUNNING: censored + unclamped should reduce the SEASONAL
# SWING substantially (target: below 0.08 nats/deg, i.e. under 40% of the current
# 0.207), while keeping the mean tilt near zero. If the swing does not fall, the
# clamp is not the seasonal driver and this route is dead.
#
# Variants, all on the validated R twin:
#   trunc_clamp   current engine
#   cens_clamp    censored, mu still clamped
#   cens_free     censored, mu UNCLAMPED            <- the physically coherent model
#   cens_free_sym censored, unclamped, symmetric arms (what is left after everything)
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))
suppressMessages(devtools::load_all(".", quiet = TRUE))
BAT <- readRDS("scratch/nes_calibration/light_battery.rds")
PSLAB <- 0.10

variants <- function(b, ratio) {
  ML <- b$max_light; LAM <- 1 / (ML * 0.5); LH <- LAM * ratio
  CC <- LAM * LH / (LAM + LH)
  r <- b$response
  mu_clamp <- function(z) pmin(pmax(r[1] + r[2] / (1 + exp((z - r[3]) / r[4])), 0), ML)
  mu_free  <- function(z) r[1] + r[2] / (1 + exp((z - r[3]) / r[4]))

  ltrunc <- function(y, mu) {
    raw <- ifelse(y <= mu, LAM * exp(-LAM * (mu - y)), LAM * exp(-LH * (y - mu)))
    lo <- (1 - exp(-LAM * mu)) / LAM
    hi <- (1 - exp(-LH * (ML - mu))) / LH
    log((1 - PSLAB) * (raw / pmax((lo + hi) * LAM, 1e-12)) + PSLAB / ML)
  }
  lcens <- function(y, mu) {
    dens <- ifelse(y <= mu, CC * exp(-LAM * (mu - y)), CC * exp(-LH * (y - mu)))
    a_lo <- (LH  / (LAM + LH)) * exp(-LAM * mu)
    a_hi <- (LAM / (LAM + LH)) * exp(-LH * (ML - mu))
    sp <- ifelse(y >= ML - 1e-9, a_hi, ifelse(y <= 1e-9, a_lo, dens))
    log((1 - PSLAB) * sp + PSLAB / ML)
  }
  list(trunc_clamp = list(f = ltrunc, mu = mu_clamp),
       cens_clamp  = list(f = lcens,  mu = mu_clamp),
       cens_free   = list(f = lcens,  mu = mu_free))
}

tilt_record <- function(b, ratio = 2, want = NULL) {
  V <- variants(b, ratio)
  if (!is.null(want)) V <- V[want]
  tn <- as.numeric(b$time)
  kb <- seq(min(tn), max(tn), by = 12 * 3600)
  dl <- 0.5; out <- list()
  for (k in seq_len(length(kb) - 1)) {
    sel <- tn > kb[k] & tn <= kb[k + 1]
    if (sum(sel) < 6) next
    t_k <- tn[sel]; y_k <- b$noisy[sel]; m <- sum(sel)
    lat0 <- approx(tn, b$lat, mean(t_k), rule = 2)$y
    lon0 <- approx(tn, b$lon, mean(t_k), rule = 2)$y
    z_up <- solar_zenith(t_k, rep(lon0, m), rep(lat0 + dl, m))
    z_dn <- solar_zenith(t_k, rep(lon0, m), rep(lat0 - dl, m))
    row <- data.table(
      id = b$id, offset = b$offset,
      decl = solar_declination(as.POSIXct(mean(t_k), origin = "1970-01-01", tz = "UTC")))
    for (nm in names(V)) {
      v <- V[[nm]]
      row[[nm]] <- (sum(v$f(y_k, v$mu(z_up))) - sum(v$f(y_k, v$mu(z_dn)))) / (2 * dl)
    }
    out[[length(out) + 1]] <- row
  }
  rbindlist(out)
}

RECS <- Filter(function(z) z$id %in% c("2021033", "2023032"), BAT)
D  <- rbindlist(lapply(RECS, tilt_record))
DS <- rbindlist(lapply(RECS, tilt_record, ratio = 1, want = "cens_free"))
setnames(DS, "cens_free", "cens_free_sym")
D[, cens_free_sym := DS$cens_free_sym]

D[, sband := cut(decl, c(-24, -15, -7, 7, 15, 24), include.lowest = TRUE,
                 labels = c("NH winter", "-15..-7", "equinox", "+7..+15", "NH summer"))]
VN <- c("trunc_clamp", "cens_clamp", "cens_free", "cens_free_sym")

cat(sprintf("%d knots\n\n=== per-knot latitude tilt (nats/deg) ===\n", nrow(D)))
tab <- D[, lapply(.SD, function(x) round(mean(x), 4)), by = sband, .SDcols = VN][order(sband)]
print(as.data.frame(tab), row.names = FALSE)

cat("\n=== summary ===\n")
cat(sprintf("%-16s %10s %10s %14s\n", "variant", "mean", "swing", "vs current"))
sw0 <- NA
for (v in VN) {
  m <- mean(D[[v]]); s <- diff(range(tab[[v]]))
  if (v == "trunc_clamp") sw0 <- s
  cat(sprintf("%-16s %+10.4f %10.4f %13s\n", v, m, s,
              if (v == "trunc_clamp") "-" else sprintf("%.0f%% swing", 100 * s / sw0)))
}

cat("\n=== fraction of the expected curve that the clamp is binding on ===\n")
for (b in RECS[c(1, 3)]) {
  r <- b$response; ML <- b$max_light
  z <- solar_zenith(as.numeric(b$time), b$lon, b$lat)
  mf <- r[1] + r[2] / (1 + exp((z - r[3]) / r[4]))
  d <- solar_declination(b$time)
  cat(sprintf("  %s off %3d: clamped %.1f%% overall | winter %.1f%% | summer %.1f%%\n",
              b$id, b$offset, 100 * mean(mf > ML),
              100 * mean(mf[d < -15] > ML), 100 * mean(mf[d > 15] > ML)))
}

cat("\nVERDICT\n")
s_free <- diff(range(tab$cens_free))
if (s_free < 0.4 * sw0) {
  cat(sprintf("  GATE PASSED: censored+unclamped cuts the seasonal swing to %.4f (%.0f%% of\n",
              s_free, 100 * s_free / sw0))
  cat("  current). The clamp on the expected curve was the seasonal driver.\n")
} else {
  cat(sprintf("  GATE FAILED: swing %.4f is %.0f%% of current. The clamp is NOT the\n",
              s_free, 100 * s_free / sw0))
  cat("  seasonal driver; do not pursue unclamping on this basis.\n")
}
