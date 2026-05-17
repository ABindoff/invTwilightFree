
library(TwilightFree)
library(SGAT)
library(raster)

track_data <- readRDS("scratch/simulated_seal_light_scenarios.rds")

# Check names
print(names(track_data))

df <- data.frame(
  Date = track_data$time, 
  Light = track_data$light_ideal
)
df$Date <- as.POSIXct(df$Date, tz="UTC")

fit <- TwilightFree(df, deployed.at=c(158.86, -54.619), retrieved.at=c(158.86, -54.619))

# Essie needs a grid
grid <- makeGrid(lon=c(130, 180), lat=c(-80, -40), cell.size=1, mask="none")

# Run Essie
cat("Running essie...\n")
ess <- SGAT::essie(fit, grid, epsilon1=0.01)

print(names(ess$lattice[[1]]))
str(ess$lattice[[1]], max.level=1)
saveRDS(ess, "scratch/test_essie.rds")
