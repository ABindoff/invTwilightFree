# IS THE MOVEMENT KERNEL BIASED POLEWARD?
#
# The grid HMM accumulates, for destination cell i,
#     sum over j of  alpha_j * exp(-dist(i,j)^2 / var2) * exp(log_norm_const)
# and `log_norm_const` is a CONSTANT, not a per-cell normaliser. So the total
# weight flowing into cell i is proportional to how many grid cells lie within
# the kernel's reach of it.
#
# On a grid uniform in DEGREES that count is not uniform. A one-degree cell of
# longitude is 111 km at the equator and 111*cos(lat) km at latitude lat, so
# high-latitude cells are packed closer together in kilometres and MORE of them
# fall inside a given radius. The kernel therefore hands high latitudes more
# incoming mass simply because the mesh is finer there.
#
# Equivalently and more properly: a density on the sphere discretised onto a
# lon-lat grid needs the cell-area Jacobian cos(lat) on the DESTINATION cell.
# Omitting it treats every cell as equal area, which over-weights the small
# high-latitude cells. Prediction: incoming mass proportional to 1/cos(lat).
#
# No data, no fitting. Pure geometry on the grid the analysis actually uses.
LON <- seq(150.5, 249.5, by = 1)
LAT <- seq(20.5, 69.5, by = 1)
SIGMA <- 110                       # diffusion, km per 12-hour step
R_E <- 6371

g <- expand.grid(lon = LON, lat = LAT)
gc <- function(lon1, lat1, lon2, lat2) {
  p1 <- lat1*pi/180; p2 <- lat2*pi/180
  a <- sin((p2-p1)/2)^2 + cos(p1)*cos(p2)*sin((lon2-lon1)*pi/360)^2
  2*R_E*asin(pmin(1, sqrt(a)))
}
# incoming weight for a destination cell, from a UNIFORM distribution over the
# grid: the quantity the forward recursion accumulates when the light says
# nothing.
incoming <- function(lat0, lon0) {
  d <- gc(lon0, lat0, g$lon, g$lat)
  k <- d <= 5*SIGMA
  sum(exp(-(d[k]^2) / (2*SIGMA^2)))
}
probe <- data.frame(lat = seq(22, 68, by = 4))
probe$weight <- vapply(probe$lat, function(la) incoming(la, 200), 0)
probe$rel <- probe$weight / probe$weight[probe$lat == 42]
probe$pred_1_over_cos <- (1/cos(probe$lat*pi/180)) / (1/cos(42*pi/180))
cat("=== incoming kernel mass by destination latitude (uniform source) ===\n")
print(data.frame(lat = probe$lat,
                 relative_mass = round(probe$rel, 3),
                 predicted_1_over_cos = round(probe$pred_1_over_cos, 3)),
      row.names = FALSE)
cat(sprintf("\ncorrelation with 1/cos(lat): %.4f\n", cor(probe$rel, probe$pred_1_over_cos)))
cat(sprintf("mass at 68 N is %.2fx the mass at 22 N\n",
            probe$rel[probe$lat == 68] / probe$rel[probe$lat == 22]))

cat("\n=== what that is worth per step, over the latitudes the seals use ===\n")
lo <- incoming(38, 200); hi <- incoming(50, 200)
cat(sprintf("38 N vs 50 N: %.4f vs %.4f, a factor of %.3f, i.e. %+.3f nats PER STEP\n",
            lo, hi, hi/lo, log(hi/lo)))
cat(sprintf("over 480 knots, unopposed: %+.1f nats of prior preference for 50 N\n",
            480*log(hi/lo)))
cat("\nThe light opposes it, which is why the bias is degrees rather than a\n")
cat("collapse to the pole -- and why WEAKENING the light (tempering) sent the\n")
cat("fitted latitude marching north past the domain centre to 51.\n")

cat("\n=== the fix: add the destination cell's area, ln(cos(lat)) ===\n")
probe$fixed <- probe$weight * cos(probe$lat*pi/180)
probe$fixed_rel <- probe$fixed / probe$fixed[probe$lat == 42]
print(data.frame(lat = probe$lat, before = round(probe$rel, 3),
                 after = round(probe$fixed_rel, 3)), row.names = FALSE)
cat(sprintf("\nspread across 22-68 N: before %.2fx, after %.2fx\n",
            max(probe$rel)/min(probe$rel),
            max(probe$fixed_rel)/min(probe$fixed_rel)))
