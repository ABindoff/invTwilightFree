
library(TwilightFree)
track_data <- readRDS("scratch/simulated_seal_light_scenarios.rds")
df <- data.frame(Date = as.POSIXct(track_data$time, tz="UTC"), Light = track_data$light_ideal)
fit <- TwilightFree(df, deployed.at=c(158.86,-54.619), retrieved.at=c(158.86,-54.619))
str(fit)
saveRDS(fit, "scratch/test_orig_tf_fit.rds")
