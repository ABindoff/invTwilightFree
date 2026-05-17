track <- readRDS('scratch/simulated_seal_light_scenarios.rds')
idx <- which(track$time > as.POSIXct('2024-01-18 09:30:00', tz='UTC') & track$time < as.POSIXct('2024-01-18 10:30:00', tz='UTC'))
print(data.frame(Time=track$time[idx], Light=track$light_ideal[idx]))
