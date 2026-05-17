
library(dplyr)
track_data <- readRDS("scratch/simulated_seal_light_scenarios.rds")

# Check solar zenith ranges per day
track_data$date <- as.Date(track_data$time)
daily_stats <- track_data %>%
  group_by(date) %>%
  summarize(
    min_z = min(z),
    max_z = max(z),
    true_lat = mean(true_lat),
    has_night = any(z > 96),
    has_day = any(z < 90)
  )

cat("Solar conditions throughout simulation:\n")
print(head(daily_stats))
cat("\n...\n")
print(tail(daily_stats))

cat("\nDays without night (Midnight Sun):", sum(!daily_stats$has_night), "\n")
cat("Days without day (Polar Night):", sum(!daily_stats$has_day), "\n")

# Save summary to file for inspection
saveRDS(daily_stats, "scratch/daily_solar_stats.rds")
