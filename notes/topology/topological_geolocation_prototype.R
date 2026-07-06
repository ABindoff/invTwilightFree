# Topological recasting of light-level geolocation: prototype (R).
#
# Companion to phase_residual_prototype.png / topo_geoloc.py. Reproduces the
# four core fields against the package's own solar model where possible.
#
# Claims demonstrated:
#   A. Each observation at KNOWN time is a circle of position on the sphere
#      (locus of equal solar zenith); several circles intersect at the fix.
#   B. The diurnal-phase residual phi(x) = wrap(H_true - H_x) is circle-valued;
#      its zero-locus is exactly the true meridian. Longitude is a pure phase,
#      latitude- and time-free. Winding of phi = longitude in turns.
#   C. Latitude rides on day length, whose sensitivity d(daylength)/d(lat) -> 0
#      at the equinox (declination 0): the fiber pinches.
#   D. Hemisphere monodromy: a southern track and its northern mirror give
#      near-identical day length, coinciding at the equinoxes (branch points).
#
# Run with the package loaded (uses invTwilightFree::solar_zenith for panel A).

library(invTwilightFree)

J2000 <- 946728000  # unix seconds at 2000-01-01 12:00 UTC

# Inline almanac (matches src/rust/src/lib.rs) for the quantities the package
# does not export (declination, GMST, hour angle, day length).
solar_parts <- function(t) {
  d <- (t - J2000) / 86400
  g <- (357.529 + 0.98560028 * d) * pi / 180
  q <- 280.459 + 0.98564736 * d
  Lsun <- (q + 1.915 * sin(g) + 0.020 * sin(2 * g)) * pi / 180
  e <- (23.439 - 3.6e-7 * d) * pi / 180
  decl <- asin(sin(e) * sin(Lsun))
  ra <- atan2(cos(e) * sin(Lsun), cos(Lsun)) * 180 / pi
  gmst <- (18.697374558 + 24.06570982441908 * d) * 15  # degrees
  list(decl = decl, ra = ra, gmst = gmst)
}

hour_angle_deg <- function(lon, t) {           # diurnal phase, 0 at local noon
  s <- solar_parts(t)
  ((s$gmst + lon - s$ra + 180) %% 360) - 180
}

day_length_h <- function(lat_deg, decl_rad) {
  cosH0 <- -tan(lat_deg * pi / 180) * tan(decl_rad)
  H0 <- acos(pmax(-1, pmin(1, cosH0)))
  24 * H0 / pi
}

# ---- grid and truth --------------------------------------------------------
lons <- seq(-180, 180, by = 1)
lats <- seq(-80, 80, by = 1)
grid <- expand.grid(lon = lons, lat = lats)
true_lon <- 150; true_lat <- -50
t0 <- as.numeric(as.POSIXct("2024-06-21 00:00", tz = "UTC"))  # austral solstice

op <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

# ---- A. circles of position + fix -----------------------------------------
plot(NA, xlim = c(-180, 180), ylim = c(-80, 80), xlab = "longitude",
     ylab = "latitude", main = "A. Circles of position meet at the fix")
cols <- c("#1f77b4", "#2ca02c", "#d62728")
for (i in seq_along(c(2, 9, 16))) {
  t <- t0 + c(2, 9, 16)[i] * 3600
  z <- matrix(solar_zenith(rep(t, nrow(grid)), grid$lon, grid$lat),
              nrow = length(lons))
  z_obs <- solar_zenith(t, true_lon, true_lat)
  contour(lons, lats, z, levels = z_obs, add = TRUE, col = cols[i],
          drawlabels = FALSE, lwd = 2)
}
points(true_lon, true_lat, pch = 21, bg = "yellow", cex = 1.6)

# ---- B. phase residual: zero-locus = true meridian -------------------------
H_true <- hour_angle_deg(true_lon, t0)
phi <- matrix(((H_true - hour_angle_deg(grid$lon, t0)) + 180) %% 360 - 180,
              nrow = length(lons))
image(lons, lats, phi, col = hcl.colors(64, "Blue-Red 3"),
      xlab = "longitude", ylab = "latitude",
      main = "B. Phase residual: zero-locus = true meridian")
contour(lons, lats, phi, levels = 0, add = TRUE, lwd = 2, drawlabels = FALSE)
cat(sprintf("B. sd of phi across latitude at a fixed lon: %.6g (expect 0)\n",
            sd(phi[which.min(abs(lons - 0)), ])))

# ---- C. day length vs latitude: equinox pinch ------------------------------
lat_axis <- seq(-75, 75, length.out = 400)
plot(NA, xlim = c(-75, 75), ylim = c(0, 24), xlab = "latitude",
     ylab = "day length (h)", main = "C. Latitude fiber pinches at equinox")
decls <- c(23.4, 12, 5, 0); dc <- c("#d62728", "#ff7f0e", "#9467bd", "#1f77b4")
for (j in seq_along(decls)) {
  lines(lat_axis, day_length_h(lat_axis, decls[j] * pi / 180), col = dc[j], lwd = 2)
}
abline(v = c(true_lat, -true_lat), lty = 3)
legend("topleft", legend = paste0("decl ", decls), col = dc, lwd = 2, bty = "n", cex = 0.8)

# ---- D. hemisphere monodromy over the year ---------------------------------
days <- 0:364
ty <- as.numeric(as.POSIXct("2024-01-01 12:00", tz = "UTC")) + days * 86400
decl_year <- vapply(ty, function(t) solar_parts(t)$decl, numeric(1))
D_true <- day_length_h(true_lat, decl_year)
D_mirror <- day_length_h(-true_lat, decl_year)
plot(days, D_true, type = "l", col = "#1f77b4", lwd = 2, ylim = c(0, 24),
     xlab = "day of 2024", ylab = "day length (h)",
     main = "D. True vs mirror coincide at equinoxes")
lines(days, D_mirror, col = "#d62728", lwd = 2, lty = 2)
abline(v = which(diff(sign(decl_year)) != 0), lty = 3)
legend("bottom", legend = c("true (S)", "mirror (N)"),
       col = c("#1f77b4", "#d62728"), lwd = 2, lty = c(1, 2), bty = "n", cex = 0.8)

par(op)
