
library(invTwilightFree)
library(SGAT)
library(FLightR)
library(TwGeos)

# Load Scenarios
track_all <- readRDS("scratch/simulated_seal_light_scenarios.rds")

run_scenario <- function(light_col, scenario_name) {
  cat("\n--- Processing Scenario:", scenario_name, "---\n")
  prefix <- paste0("_seal_", scenario_name)
  
  light <- track_all[[light_col]]
  date_time <- track_all$time
  
  # 1. Preprocess Twilights (Automated threshold = 10)
  twl_data <- data.frame(Date = date_time, Light = light)
  twl <- try(TwGeos::findTwilights(twl_data, threshold = 10, include = twl_data$Date))
  if (inherits(twl, "try-error") || nrow(twl) < 50) {
    cat("Warning: Twilight detection failed or insufficient twilights for", scenario_name, "\n")
    return(NULL)
  }
  
  # Filter fake high-frequency twilights caused by light noise
  twl <- twl[c(TRUE, diff(as.numeric(twl$Twilight)) > 7200), ]
  cat("Filtered to", nrow(twl), "true twilights.\n")
  
  # 2. FLightR
  cat("Running FLightR...\n")
  time_flightr <- system.time({
    flightr_data <- data.frame(
      datetime = format(date_time, "%Y-%m-%dT%H:%M:%SZ", tz="UTC"),
      light = light,
      twilight = 0,
      interp = FALSE, excluded = FALSE
    )
    for (i in 1:nrow(twl)) {
      idx <- which.min(abs(as.numeric(date_time) - as.numeric(twl$Twilight[i])))
      flightr_data$twilight[idx] <- ifelse(twl$Rise[i], 1, 2)
    }
    flightr_data$interp[1] <- TRUE
    tmp_file <- paste0("scratch/tmp_flightr_", scenario_name, ".csv")
    write.csv(flightr_data, tmp_file, row.names = FALSE, quote = FALSE)
    
    fit_flightr <- tryCatch({
      proc_data <- FLightR::get.tags.data(tmp_file, log.light.borders = c(0, 64), measurement.period = 240)
      Grid <- FLightR::make.grid(left = 130, right = 180, bottom = -80, top = -40, distance.from.land.allowed.to.use = c(-Inf, Inf))
      
      # 14-day calibration
      calib_df <- data.frame(calibration.start = date_time[1], calibration.stop = date_time[1] + 14*86400, lon = track_all$true_lon[1], lat = track_all$true_lat[1])
      calib_flightr <- FLightR::make.calibration(proc_data, Calibration.periods = calib_df)
      prerun <- FLightR::make.prerun.object(proc_data, Grid, start = c(track_all$true_lon[1], track_all$true_lat[1]), Calibration = calib_flightr)
      FLightR::run.particle.filter(prerun, nParticles = 30000, plot = FALSE)
    }, error = function(e) {
      cat("FLightR failed:", e$message, "\n")
      NULL
    })
  })

  saveRDS(list(fit = fit_flightr, time = time_flightr), paste0("vignettes/cache_flightr", prefix, ".rds"))
  
  cat("Running SGAT...\n")
  time_sgat <- system.time({
    fit_sgat <- tryCatch({
      path <- SGAT::thresholdPath(twl$Twilight, twl$Rise, unfold = FALSE)
      x0 <- path$x
      if(any(is.na(x0[,2]))) x0[,2] <- approx(which(!is.na(x0[,2])), x0[!is.na(x0[,2]), 2], xout = 1:nrow(x0), rule = 2)$y
      z0 <- twl$Twilight
      
      alpha_mat <- matrix(c(0, 20), nrow = nrow(twl), ncol = 2, byrow = TRUE)
      model <- SGAT::thresholdModel(twl$Twilight, twl$Rise, twilight.model = "Normal", alpha = alpha_mat, beta = c(2, 2), x0 = x0, z0 = z0, zenith = 94.276)
      proposal <- SGAT::mvnorm(S = diag(c(0.1, 0.1)^2), n = nrow(x0))
      SGAT::stellaMetropolis(model, proposal, x0 = x0, iters = 20000, thin = 20)
    }, error = function(e) {
      cat("SGAT failed:", e$message, "\n")
      NULL
    })
  })
  saveRDS(list(fit = fit_sgat, time = time_sgat), paste0("vignettes/cache_sgat", prefix, ".rds"))
}

# Execute
# run_scenario("light_ideal", "ideal")
# run_scenario("light_cloudy", "cloudy")
run_scenario("light_shaded", "shaded")
run_scenario("light_alan", "alan")
