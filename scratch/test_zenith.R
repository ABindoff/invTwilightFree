
library(invTwilightFree)
library(SGAT)

# Test locations and times
t <- seq(as.POSIXct("2024-01-01", tz="UTC"), as.POSIXct("2024-01-02", tz="UTC"), by="1 hour")
lon <- rep(158.860, length(t))
lat <- rep(-54.619, length(t))

# 1. invTF zenith
z_inv <- invTwilightFree::solar_zenith(as.numeric(t), lon, lat)

# We'll just compare invTF and R_script, and print SGAT::zenith to see its source.
z_sgat <- rep(NA, length(t))

# 3. generate_seal_light.R zenith
get_z <- function(t, lon, lat) {
  d <- (as.numeric(t) - 946728000) / 86400
  g <- (357.529 + 0.98560028 * d) * pi / 180
  q <- 280.459 + 0.98564736 * d
  l_sun <- (q + 1.915 * sin(g) + 0.020 * sin(2 * g)) * pi / 180
  e <- (23.439 - 0.00000036 * d) * pi / 180
  delta <- asin(sin(e) * sin(l_sun))
  ra <- atan2(cos(e) * sin(l_sun), cos(l_sun))
  gmst <- (18.697374558 + 24.06570982441908 * d) * 15 * pi / 180
  lmst <- gmst + lon * pi / 180
  h <- lmst - ra
  phi <- lat * pi / 180
  cos_z <- sin(phi) * sin(delta) + cos(phi) * cos(delta) * cos(h)
  acos(pmax(-1, pmin(1, cos_z))) * 180 / pi
}

z_R <- mapply(get_z, t, lon, lat)

df <- data.frame(
  time = t,
  invTF = z_inv,
  SGAT = z_sgat,
  R_script = z_R
)

print(head(df))
cat("\nMax diff invTF vs SGAT: ", max(abs(df$invTF - df$SGAT)), "\n")
cat("Max diff invTF vs R_script: ", max(abs(df$invTF - df$R_script)), "\n")

# Let's test a changing longitude to see if there's a difference there
lon2 <- seq(158, 175, length.out=length(t))
lat2 <- seq(-54, -68, length.out=length(t))
z_inv2 <- invTwilightFree::solar_zenith(as.numeric(t), lon2, lat2)
z_sgat2 <- SGAT::zenith(t, cbind(lon2, lat2))
z_R2 <- mapply(get_z, t, lon2, lat2)

cat("\nMoving track test:\n")
cat("Max diff invTF vs SGAT: ", max(abs(z_inv2 - z_sgat2)), "\n")
cat("Max diff invTF vs R_script: ", max(abs(z_inv2 - z_R2)), "\n")
