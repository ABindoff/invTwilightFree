# Is our solar zenith right? An INDEPENDENT reference, not another transcription.
#
# The package's `solar_zenith()` (src/rust/src/lib.rs:11) is the USNO Astronomical
# Almanac LOW-PRECISION formula -- the same one SGAT uses, because it came from
# there. Comparing it against SGAT would therefore test only whether it was copied
# correctly, which is worth knowing but is not a test of the algorithm.
#
# So this implements NOAA's Solar Calculator (Meeus, Astronomical Algorithms ch.25)
# from the published coefficients: equation of centre, nutation-corrected apparent
# longitude, corrected obliquity, and the equation of time. It is a genuinely
# different algorithm with roughly two orders of magnitude better accuracy, so a
# disagreement is informative about OURS.
#
# It also reports the ATMOSPHERIC REFRACTION term separately. Ours computes the
# GEOMETRIC zenith with no refraction at all. Refraction is ~0.57 deg at the
# horizon and is the largest physical term in twilight geolocation, so whether it
# matters here is a question worth answering with a number rather than an opinion.
#
# NOTE ON WHAT THIS CAN AND CANNOT EXPLAIN. In the synthetic battery the light is
# GENERATED with the same `solar_zenith()` the engine fits with, so any zenith
# error is common to both sides and cancels exactly -- the truth stays the exact
# maximiser. This diagnostic therefore says nothing about the null gate's -0.64 deg.
# On REAL data there is no such cancellation: the sun is where it actually is, and
# a zenith error displaces the fitted position to compensate. That is the arm this
# is for.
suppressMessages(devtools::load_all(".", quiet = TRUE))

deg <- function(x) x * 180 / pi
rad <- function(x) x * pi / 180

# --- NOAA / Meeus reference -------------------------------------------------
noaa_zenith <- function(unix_time, lon, lat, refraction = FALSE) {
  jd <- unix_time / 86400 + 2440587.5
  T  <- (jd - 2451545) / 36525

  L0 <- (280.46646 + T * (36000.76983 + T * 0.0003032)) %% 360
  M  <- 357.52911 + T * (35999.05029 - 0.0001537 * T)
  e  <- 0.016708634 - T * (0.000042037 + 0.0000001267 * T)

  C <- sin(rad(M))     * (1.914602 - T * (0.004817 + 0.000014 * T)) +
       sin(rad(2 * M)) * (0.019993 - 0.000101 * T) +
       sin(rad(3 * M)) *  0.000289

  true_long <- L0 + C
  omega     <- 125.04 - 1934.136 * T
  app_long  <- true_long - 0.00569 - 0.00478 * sin(rad(omega))

  e0        <- 23 + (26 + (21.448 - T * (46.815 + T * (0.00059 - T * 0.001813))) / 60) / 60
  obliq     <- e0 + 0.00256 * cos(rad(omega))

  decl <- deg(asin(sin(rad(obliq)) * sin(rad(app_long))))

  # equation of time (minutes)
  y      <- tan(rad(obliq / 2))^2
  eqtime <- 4 * deg(y * sin(2 * rad(L0)) - 2 * e * sin(rad(M)) +
                    4 * e * y * sin(rad(M)) * cos(2 * rad(L0)) -
                    0.5 * y * y * sin(4 * rad(L0)) -
                    1.25 * e * e * sin(2 * rad(M)))

  mins <- (unix_time %% 86400) / 60          # minutes past UTC midnight
  tst  <- (mins + eqtime + 4 * lon) %% 1440  # true solar time
  ha   <- tst / 4 - 180                      # hour angle, degrees

  cosz <- sin(rad(lat)) * sin(rad(decl)) + cos(rad(lat)) * cos(rad(decl)) * cos(rad(ha))
  z    <- deg(acos(pmin(pmax(cosz, -1), 1)))

  if (!refraction) return(z)

  # NOAA refraction on the APPARENT elevation, in arc seconds
  el <- 90 - z
  te <- tan(rad(el))
  r  <- ifelse(el > 85, 0,
        ifelse(el > 5,  58.1 / te - 0.07 / te^3 + 0.000086 / te^5,
        ifelse(el > -0.575,
               1735 + el * (-518.2 + el * (103.4 + el * (-12.79 + el * 0.711))),
               -20.772 / te)))
  z - r / 3600
}

# --- comparison grid: the regime these deployments actually occupy ----------
set.seed(1)
n   <- 200000
t0  <- as.numeric(as.POSIXct("2021-01-01", tz = "UTC"))
tt  <- t0 + runif(n, 0, 3 * 365.25 * 86400)   # 2021-2023, the three seasons
lat <- runif(n, 30, 60)                        # NES range
lon <- runif(n, -180, 180)

ours <- solar_zenith(tt, lon, lat)
ref  <- noaa_zenith(tt, lon, lat, refraction = FALSE)
refr <- noaa_zenith(tt, lon, lat, refraction = TRUE)

d_geom <- ours - ref     # algorithm/transcription difference
d_refr <- ours - refr    # total error incl. the omitted physics

band <- function(z) cut(z, c(0, 85, 96, 102, 180),
                        labels = c("day <85", "twilight 85-96",
                                   "late twi 96-102", "night >102"))
b <- band(ref)

cat(sprintf("n = %d, 2021-2023, lat 30-60 N\n\n", n))

cat("=== OURS vs NOAA, both GEOMETRIC (no refraction) ===\n")
cat(sprintf("  overall: bias %+.5f deg | rms %.5f | max|d| %.5f\n",
            mean(d_geom), sqrt(mean(d_geom^2)), max(abs(d_geom))))
s <- tapply(d_geom, b, function(x)
  sprintf("bias %+.5f | rms %.5f | max %.5f", mean(x), sqrt(mean(x^2)), max(abs(x))))
for (k in names(s)) cat(sprintf("  %-16s %s\n", k, s[[k]]))

cat("\n=== OURS (geometric) vs NOAA WITH refraction ===\n")
cat("  i.e. how far our zenith sits from the sun's APPARENT position\n")
cat(sprintf("  overall: bias %+.5f deg | rms %.5f | max|d| %.5f\n",
            mean(d_refr), sqrt(mean(d_refr^2)), max(abs(d_refr))))
s2 <- tapply(d_refr, b, function(x)
  sprintf("bias %+.5f | rms %.5f | max %.5f", mean(x), sqrt(mean(x^2)), max(abs(x))))
for (k in names(s2)) cat(sprintf("  %-16s %s\n", k, s2[[k]]))

cat("\n=== what a zenith offset is worth in latitude ===\n")
cat("  dz/dlat = 1 exactly when the sun is due N or S; the twilight band is\n")
cat("  where geolocation reads latitude, so a systematic offset there maps\n")
cat("  to latitude at roughly 1:1 in the worst case.\n")
tw <- d_refr[b == "twilight 85-96"]
cat(sprintf("  twilight-band mean offset: %+.4f deg -> up to %+.0f km of latitude\n",
            mean(tw), mean(tw) * 111))
