library(TwGeos)
track <- readRDS('scratch/simulated_seal_light_scenarios.rds')
twl_data <- data.frame(Date = track$time, Light = track$light_ideal)
twl <- findTwilights(twl_data, threshold = 10, include = twl_data$Date)
twl <- twl[c(TRUE, diff(as.numeric(twl$Twilight)) > 7200), ]
cat('Num Twilights:', nrow(twl), '\n')
