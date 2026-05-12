
library(SGAT)
cached <- readRDS('c:/Users/bindoffa/antigravity projects/invTwilightFree/vignettes/cache_sgat.rds')
fit_sgat <- cached$fit

cat("Initial path x0 summary:\n")
x0 <- fit_sgat$model$x0
print(head(x0))

# Calculate error relative to truth if possible
# We need the track data for truth
# Let's just look at the longitude trend in x0
lons <- x0[,1]
times <- as.numeric(fit_sgat$model$z0)
cat("First Lon:", lons[1], "at", times[1], "\n")
cat("Last Lon:", lons[length(lons)], "at", times[length(times)], "\n")

# Simple linear regression to check drift
fit <- lm(lons ~ times)
cat("Drift (deg/sec):", coef(fit)[2], "\n")
cat("Drift (deg/day):", coef(fit)[2] * 86400, "\n")
