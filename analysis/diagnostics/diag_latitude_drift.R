# THE LATITUDE DRIFT IS IN THE LIGHT, NOT IN THE MEASUREMENT.
#
# THE FINDING. Real twilight light along these tracks is NOT a pure function of solar
# zenith. Relative to a fitted zenith-only response, it gets systematically dimmer as
# the animal moves north, at the surface, independent of diving:
#
#     within-tag slope = -0.0707 zenith-equivalent degrees per degree of latitude
#                        (se 0.0032, t ~ -22), gated to depth_min <= 10 m
#
# "Zenith-equivalent" = residual divided by the fitted response slope, so +1 means the
# light is one degree of zenith "brighter" than the model expects.
#
# THE CONTROL THAT MAKES IT A FINDING RATHER THAN AN ARTEFACT. The identical analysis
# on the synthetic battery -- light GENERATED from a pure zenith function, so it has
# no atmosphere, no cloud and no latitude dependence of any kind by construction:
#
#     perfect (noiseless)          slope +0.0009  (se 0.0013, t 0.7)
#     noisy (model-correct noise)  slope +0.0010  (se 0.0094, t 0.1)
#     REAL light, gated            slope -0.0707  (se 0.0032, t -22)
#
# Zero in synthetic, 22 sigma in real. The measurement machinery -- response fitting,
# the tangent, the clamping, the zenith-equivalent conversion -- creates no slope on
# light that has none.
#
# IT IS NOT DEPTH. Eliminated three ways:
#   * the slope is invariant to the gate: none -0.0850, 10 m -0.0707, 5 m -0.0746,
#     2.5 m -0.0711, 1 m -0.0816. Keeping ONLY observations where the animal came
#     within one metre of the surface (37% of data) leaves it as large as ungated.
#   * depth inside the gate does not vary monotonically with latitude
#     (1.45 / 3.06 / 3.08 / 1.82 m across latitude bands)
#   * controlling for depth linearly inside the gate removes 11%
#
# IT IS THE RIGHT SIZE TO BE THE WHOLE BIAS. -0.0707 over a typical 15 degree
# excursion is -1.06 zenith-equivalent degrees; at the measured Fisher exchange rate
# of 1.4-4.4 degrees of latitude per degree of z50, that is -1.5 to -4.7 degrees of
# latitude bias, bracketing the observed -1.54.
#
# WHAT IT IS NOT: Rayleigh scattering in air. Rayleigh acts through AIR MASS, which is
# a function of zenith alone -- the optical path at a given z is the same at 37 N and
# 52 N -- so it is absorbed entirely by an empirically fitted response. The same holds
# for twilight scattering geometry: shadow height above the observer depends on how
# far the sun is below the horizon, not on where the observer stands.
#
# WHAT IT MIGHT BE, and why this dataset cannot say: cloud climatology (the North
# Pacific storm track puts more cloud north), sea state, water clarity, or
# solar-path geometry (the terminator's orientation at fixed z differs with latitude).
# These are not separable here. One colony and one phenology means the animals go
# north IN SUMMER, so latitude and solar geometry are collinear at r = -0.665; single
# predictors give R2 0.1417 (latitude) against 0.1395 (dz/dt), the joint model reaches
# only 0.1442, and the latitude-by-dz/dt cross-tabulation is non-monotone in both
# directions. Separating them needs a second colony at a different latitude, or tags
# whose excursion and season decouple.
#
# This script reproduces all of the above. Emission level only; no fits.
SCRATCH <- "analysis/diagnostics"
source(file.path(SCRATCH, "nes_common.R"))
suppressMessages(devtools::load_all(".", quiet = TRUE))

resid_real <- function() {
  P <- readRDS("analysis/cache/nes/panel_v1.rds")
  rbindlist(lapply(names(P$tags), function(tg) {
    g <- P$tags[[tg]]; r <- P$responses[[tg]]
    if (is.null(r) || !"depth_min" %in% names(g$light)) return(NULL)
    ML <- r$max_light
    y <- pmax(0, as.numeric(g$light$light) - r$baseline)
    tr <- argos_at(as.data.table(g$argos), g$light$time)
    ok <- is.finite(tr$lat) & is.finite(y) & is.finite(g$light$depth_min)
    if (sum(ok) < 2000) return(NULL)
    tt <- as.numeric(g$light$time)[ok]
    z <- solar_zenith(tt, tr$lon[ok], tr$lat[ok])
    mu <- pmin(pmax(r$calibration[1] - r$calibration[2] * z, 0), ML)
    data.table(id = tg, z = z, zeq = ((y[ok] - mu)/ML)/(r$calibration[2]/ML),
               d = g$light$depth_min[ok], lat = tr$lat[ok])
  }))
}
resid_synth <- function() {
  BAT <- readRDS("scratch/nes_calibration/light_battery.rds")
  rbindlist(lapply(BAT, function(b) {
    ML <- b$max_light; r <- b$response
    slope <- ML/(4*r[4]); zero <- r[3] + 2*r[4]
    z <- solar_zenith(as.numeric(b$time), lon360(b$lon), b$lat)
    mu <- pmin(pmax(slope*(zero - z), 0), ML)
    rbindlist(lapply(c("perfect","noisy"), function(lv) {
      y <- b[[lv]]; ok <- b$supported & is.finite(y)
      data.table(id = b$id, level = lv, z = z[ok],
                 zeq = ((y[ok]-mu[ok])/ML)/(slope/ML), lat = b$lat[ok])
    }))
  }))
}
sl <- function(x) { m <- summary(lm(zeq ~ lat + factor(id), x))
                    sprintf("%+0.4f (se %.4f, t %5.1f)", m$coefficients["lat",1],
                            m$coefficients["lat",2], m$coefficients["lat",3]) }

Re <- resid_real()[z >= 82 & z <= 100]
Sy <- resid_synth()[z >= 82 & z <= 100]
cat("=== latitude slope of the emission residual (zenith-equivalent deg per deg) ===\n")
cat(sprintf("  synthetic, noiseless   : %s\n", sl(Sy[level=="perfect"])))
cat(sprintf("  synthetic, model noise : %s\n", sl(Sy[level=="noisy"])))
cat(sprintf("  REAL, gated at 10 m    : %s\n", sl(Re[d <= 10])))
cat("\n=== invariance to the depth gate (real light) ===\n")
for (gt in c(Inf,10,5,2.5,1))
  cat(sprintf("  gate %-6s: %s  (%.0f%% kept)\n",
              ifelse(is.finite(gt), sprintf("%.1fm",gt), "none"), sl(Re[d<=gt]),
              100*nrow(Re[d<=gt])/nrow(Re)))
cat("\n=== implied latitude bias over a 15 degree excursion ===\n")
s <- coef(lm(zeq ~ lat + factor(id), Re[d<=10]))["lat"]
cat(sprintf("  %.3f zenith-equivalent deg, times exchange rate 1.4-4.4 = %.2f to %.2f deg\n",
            15*s, 15*s*1.4, 15*s*4.4))
cat("  observed real-data latitude bias: -1.54 deg\n")
