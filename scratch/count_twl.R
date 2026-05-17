library(TwGeos)
track <- readRDS('scratch/simulated_seal_light_scenarios.rds')
twl_data <- data.frame(Date = track$time, Light = track$light_ideal)
twl <- findTwilights(twl_data, threshold = 10, include = twl_data$Date)
cat('Num Twilights IDEAL:', nrow(twl), '\n')

twl_data2 <- data.frame(Date = track$time, Light = track$light_cloudy)
twl2 <- findTwilights(twl_data2, threshold = 10, include = twl_data2$Date)
cat('Num Twilights CLOUDY:', nrow(twl2), '\n')
