# DOES GATING ON MEASURED DEPTH REMOVE THE LATITUDE-STRUCTURED EMISSION DRIFT?
#
# THE FINDING THIS TESTS. `depth_min` (shallowest depth in each 30-minute decimation
# window) is in the cached panel and was never regressed against the emission
# residual at observation level. Verified independently:
#   * it is a GATE, not a gradient -- median residual is flat at +0.028 / +0.022 for
#     depth <2.5 m and 2.5-10 m, then plunges to -0.103 / -0.246 / -0.492 for
#     10-25 / 25-100 / >100 m;
#   * **95.1% of slab-scale events (residual < -0.3 of max_light) occur below 10 m**;
#   * depth explains **51.3%** of within-tag residual variance in the graded band;
#   * 41.9% of graded windows never come shallower than 2.5 m, so the 30-minute
#     MAXIMUM is not recovering surface light.
#
# The last point contradicts the recorded verdict "depth attenuation: dead", whose
# stated reason was that 96% of windows reach within 10 m of the surface. Reaching
# 10 m is not reaching the surface, and the gate sits exactly there.
#
# WHY IT COULD MATTER WHERE GLOBAL KNOBS DID NOT. Only ~4.2% of graded windows are
# below 10 m, but their median residual is -0.1 to -0.5 of max_light. Removing them
# shifts the mean residual by roughly 0.013 of max_light, about 0.35 zenith-equivalent
# degrees, which at the measured exchange rate of 1.4-4.4 deg latitude per degree of
# z50 is **0.5-1.5 degrees of latitude**. Structured in space and season rather than
# global -- which is why every one-parameter sweep returned null.
#
# THE TEST. The bias is latitude-proportional WITHIN track (reported at roughly -0.21
# deg of bias per degree north of the colony). If depth shading is the cause, gating
# on depth should flatten that slope. If the slope survives gating, the depth
# mechanism is dead at the emission level and no fit need be spent on it.
#
# PRE-REGISTERED: gating at 10 m should at least HALVE the latitude slope of the mean
# signed residual. Anything less and this is not the mechanism.
#
# Emission-level only. No fits. Aggregates only.
SCRATCH <- "analysis/diagnostics"
source(file.path(SCRATCH, "nes_common.R"))
suppressMessages(devtools::load_all(".", quiet = TRUE))
P <- readRDS("analysis/cache/nes/panel_v1.rds")

R <- list()
for (tg in names(P$tags)) {
  g <- P$tags[[tg]]; r <- P$responses[[tg]]
  if (is.null(r) || !"depth_min" %in% names(g$light)) next
  ML <- r$max_light
  y <- pmax(0, as.numeric(g$light$light) - r$baseline)
  tr <- argos_at(as.data.table(g$argos), g$light$time)
  ok <- is.finite(tr$lat) & is.finite(y) & is.finite(g$light$depth_min)
  if (sum(ok) < 2000) next
  tt <- as.numeric(g$light$time)[ok]
  z  <- solar_zenith(tt, tr$lon[ok], tr$lat[ok])
  mu <- pmin(pmax(r$calibration[1] - r$calibration[2] * z, 0), ML)
  R[[tg]] <- data.table(id = tg, z = z, res = (y[ok] - mu) / ML,
                        d = g$light$depth_min[ok], lat = tr$lat[ok],
                        decl = solar_declination(as.POSIXct(tt, origin="1970-01-01", tz="UTC")),
                        slope_per_deg = r$calibration[2] / ML)
}
D <- rbindlist(R)
G <- D[z >= 82 & z <= 100]                      # graded band; the rest is clamped
cat(sprintf("%d tags, %d graded-band observations\n", uniqueN(G$id), nrow(G)))
cat(sprintf("gate at 10 m removes %.1f%% of them\n\n", 100 * mean(G$d > 10)))

# residual expressed in ZENITH-EQUIVALENT degrees: how many degrees of zenith the
# light is "off" by, which is the unit the latitude exchange rate applies to
G[, zeq := res / slope_per_deg]

lat_slope <- function(x) {
  f <- lm(zeq ~ lat + factor(id), data = x)      # within-tag slope
  c(slope = unname(coef(f)["lat"]), se = summary(f)$coefficients["lat", 2])
}
a <- lat_slope(G); b <- lat_slope(G[d <= 10])
cat("=== latitude slope of the mean signed residual (zenith-equivalent deg per deg lat) ===\n")
cat(sprintf("  ungated  : %+0.4f (se %.4f)\n", a["slope"], a["se"]))
cat(sprintf("  gated 10m: %+0.4f (se %.4f)\n", b["slope"], b["se"]))
cat(sprintf("  reduction: %.0f%%   -> %s\n", 100*(1 - abs(b["slope"])/abs(a["slope"])),
            if (abs(b["slope"]) <= 0.5*abs(a["slope"])) "MEETS the pre-registered halving"
            else "does NOT meet it"))
cat(sprintf("\n  implied latitude bias slope at exchange rate 1.4-4.4:\n"))
cat(sprintf("    ungated  %+0.3f to %+0.3f deg lat per deg lat\n", a["slope"]*1.4, a["slope"]*4.4))
cat(sprintf("    gated    %+0.3f to %+0.3f\n", b["slope"]*1.4, b["slope"]*4.4))

cat("\n=== mean signed residual (zenith-equivalent deg) by true latitude ===\n")
G[, lb := cut(lat, c(-Inf,38,43,48,Inf), labels=c("<38 (colony)","38-43","43-48",">48"))]
t1 <- G[, .(n=.N, ungated=round(mean(zeq),3)), by=lb]
t2 <- G[d<=10, .(gated=round(mean(zeq),3)), by=lb]
print(as.data.frame(merge(t1,t2,by="lb")[order(lb)]), row.names=FALSE)

cat("\n=== and by season ===\n")
G[, sb := cut(decl, c(-24,-10,10,24), labels=c("NH winter","equinox","NH summer"))]
s1 <- G[!is.na(sb), .(n=.N, ungated=round(mean(zeq),3)), by=sb]
s2 <- G[d<=10 & !is.na(sb), .(gated=round(mean(zeq),3)), by=sb]
print(as.data.frame(merge(s1,s2,by="sb")[order(sb)]), row.names=FALSE)

cat("\n=== does the gate make the spread expressible? ===\n")
emp <- function(x) c(ratio = mean(abs(x[x<0]))/mean(x[x>0]), tail = mean(abs(x) > 0.5))
G[, zb := cut(z, c(82,86,90,94,98,100))]
u <- G[, as.list(emp(res)), by=zb][order(zb)]
gg <- G[d<=10, as.list(emp(res)), by=zb][order(zb)]
m <- merge(u, gg, by="zb", suffixes=c("_ungated","_gated"))
print(as.data.frame(m[, lapply(.SD, function(x) if(is.numeric(x)) round(x,3) else x)]), row.names=FALSE)
cat("  model assumes ratio 2.00 and tail (prob_slab) 0.10 in every band.\n")
