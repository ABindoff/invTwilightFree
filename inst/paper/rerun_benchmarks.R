# =============================================================================
# Benchmark re-measurement for Tables 1 and 2 (the re-run flagged in MS 2.10).
#
# The original comparator runtimes were not measured: run_benchmark_final.R
# stored `time = c(elapsed = 131*60)` and `134*60` as literal constants, and
# run_benchmark_ideal.R stored `time = 0`. This script measures every method
# back to back in one process on a quiesced machine, repeats the cheap panel,
# pins all thread counts to 1, and records the Monte Carlo effort next to each
# timing so the numbers are comparable in the only sense they can be.
#
# Results are appended to inst/paper/benchmark_rerun.csv after every fit, so a
# partial run is still usable. Fits from the first repetition are cached for
# Figure 1.
#
#   R-4.6.0/bin/Rscript inst/paper/rerun_benchmarks.R
#
# Requires R 4.6.0: SGAT, FLightR, TwGeos and raster are installed there only.
# =============================================================================

for (v in c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
            "RAYON_NUM_THREADS", "NUMEXPR_NUM_THREADS"))
  do.call(Sys.setenv, setNames(list("1"), v))

suppressMessages({
  library(invTwilightFree); library(TwGeos); library(SGAT); library(FLightR)
})

OUT       <- "inst/paper/benchmark_rerun.csv"
FITS      <- "inst/paper/fig1_fits.rds"
REPS_SEAL <- 3L
REPS_GEOM <- 1L
SEAL_TH   <- 10   # twilight detection threshold
SEAL_GAP  <- 4    # h, minimum separation; the only value that lets SGAT sample
                  # on all three seal scenarios (see MS 3.2)

gc_dist <- function(lon1, lat1, lon2, lat2) {
  r <- pi / 180
  acos(pmin(1, sin(lat1 * r) * sin(lat2 * r) +
              cos(lat1 * r) * cos(lat2 * r) * cos((lon2 - lon1) * r))) * 6371
}
rmse_to_truth <- function(time, lon, lat, tt, tlat, tlon) {
  t <- as.numeric(time)
  sqrt(mean(gc_dist(lon, lat,
                    approx(as.numeric(tt), tlon, xout = t)$y,
                    approx(as.numeric(tt), tlat, xout = t)$y)^2, na.rm = TRUE))
}
# A fixed arithmetic probe. Comparing it before and after the run detects
# machine load drifting under us, which is the whole reason for the re-run.
probe <- function() {
  t0 <- proc.time(); x <- 0
  for (i in 1:2e6) x <- x + sqrt(i)
  unname((proc.time() - t0)["elapsed"])
}
emit <- function(...) {
  row <- data.frame(..., stringsAsFactors = FALSE)
  write.table(row, OUT, sep = ",", row.names = FALSE, na = "",
              col.names = !file.exists(OUT), append = file.exists(OUT))
  cat(sprintf("  %-9s %-20s rep%d  %8s km  %8.1f s\n", row$method, row$scenario,
              row$rep, ifelse(is.na(row$rmse_km), "FAIL", sprintf("%.0f", row$rmse_km)),
              row$seconds))
}
if (file.exists(OUT)) file.remove(OUT)
fits_cache <- list()

# ---- twilight preparation shared by SGAT and FLightR -----------------------
prep_twl <- function(time, light, threshold, gap_h) {
  twl <- TwGeos::findTwilights(data.frame(Date = time, Light = light),
                               threshold = threshold, include = time)
  if (is.null(twl) || nrow(twl) < 2) return(NULL)
  if (!is.na(gap_h))
    twl <- twl[c(TRUE, diff(as.numeric(twl$Twilight)) > gap_h * 3600), ]
  twl
}

run_sgat <- function(twl, iters, thin, prop_sd, zenith) {
  fit <- tryCatch({
    path <- SGAT::thresholdPath(twl$Twilight, twl$Rise, unfold = FALSE)
    x0 <- path$x
    if (any(is.na(x0[, 2])))
      x0[, 2] <- approx(which(!is.na(x0[, 2])), x0[!is.na(x0[, 2]), 2],
                        xout = 1:nrow(x0), rule = 2)$y
    am <- matrix(c(0, 45), nrow = nrow(twl), ncol = 2, byrow = TRUE)
    m  <- SGAT::thresholdModel(twl$Twilight, twl$Rise, twilight.model = "Normal",
                               alpha = am, beta = c(1, 0.08), x0 = x0,
                               z0 = twl$Twilight, zenith = zenith)
    p  <- SGAT::mvnorm(S = diag(rep(prop_sd, 2)^2), n = nrow(x0))
    SGAT::stellaMetropolis(m, p, x0 = x0, iters = iters, thin = thin)
  }, error = function(e) {
    cat("    SGAT error:", conditionMessage(e), "\n"); NULL
  })
  if (is.null(fit)) return(NULL)
  st <- if (!is.null(fit$model$time)) fit$model$time else fit$model$z
  list(time = st,
       lon = apply(fit$x[[1]][, 1, ], 1, mean),
       lat = apply(fit$x[[1]][, 2, ], 1, mean))
}

run_flightr <- function(time, light, twl, grid_box, start, meas_period, nPart) {
  fit <- tryCatch({
    fd <- data.frame(datetime = format(time, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
                     light = light, twilight = 0, interp = FALSE, excluded = FALSE)
    for (i in seq_len(nrow(twl))) {
      idx <- which.min(abs(as.numeric(time) - as.numeric(twl$Twilight[i])))
      fd$twilight[idx] <- if (twl$Rise[i]) 1 else 2
    }
    fd$interp[1] <- TRUE
    tmp <- tempfile(fileext = ".csv")
    write.csv(fd, tmp, row.names = FALSE, quote = FALSE)
    pd <- FLightR::get.tags.data(tmp, log.light.borders = c(0, 64),
                                 measurement.period = meas_period)
    G  <- FLightR::make.grid(left = grid_box[1], right = grid_box[2],
                             bottom = grid_box[3], top = grid_box[4],
                             distance.from.land.allowed.to.use = c(-Inf, Inf))
    cal <- FLightR::make.calibration(pd, Calibration.periods = data.frame(
      calibration.start = time[1], calibration.stop = time[1] + 3 * 86400,
      lon = start[1], lat = start[2]))
    pr <- FLightR::make.prerun.object(pd, G, start = start, Calibration = cal)
    FLightR::run.particle.filter(pr, threads = 1, nParticles = nPart, plot = FALSE)
  }, error = function(e) {
    cat("    FLightR error:", conditionMessage(e), "\n"); NULL
  })
  if (is.null(fit)) return(NULL)
  q <- fit$Results$Quantiles
  list(time = q$time, lon = q$Medianlon, lat = q$Medianlat)
}

p_start <- probe()
cat(sprintf("machine probe (start): %.3f s\n\n", p_start))

# =============================================================================
# Table 2: 40-day seal CRW, three noise scenarios, REPS_SEAL repetitions
# =============================================================================
td   <- readRDS("scratch/simulated_seal_light_scenarios.rds")
seal <- list(list(label = "Cloudy",              col = "light_ideal"),
             list(label = "Shaded (ARS diving)", col = "light_shaded"),
             list(label = "ALAN near colony",    col = "light_alan"))
seal_start <- c(td$true_lon[1], td$true_lat[1])

cat("=== Table 2: seal scenarios ===\n")
for (rep in seq_len(REPS_SEAL)) for (sc in seal) {
  lt <- td[[sc$col]]

  tm <- system.time(f <- TwilightFreeSMC(
    date_time = td$time, light = lt,
    start_lat = td$true_lat[1],        start_lon = td$true_lon[1],
    end_lat   = td$true_lat[nrow(td)], end_lon   = td$true_lon[nrow(td)],
    method = "guided", n_particles = 1000, seed = 42))
  kt <- as.POSIXct(f$knot_times, origin = "1970-01-01", tz = "UTC")
  emit(table = "2", scenario = sc$label, method = "invTF", rep = rep,
       rmse_km = rmse_to_truth(kt, f$lon, f$lat, td$time, td$true_lat, td$true_lon),
       seconds = unname(tm["elapsed"]), effort = "1000 particles")
  if (rep == 1L)
    fits_cache[[paste0(sc$label, "|invTF")]] <- data.frame(time = kt, lon = f$lon, lat = f$lat)

  twl   <- prep_twl(td$time, lt, SEAL_TH, SEAL_GAP)
  n_twl <- if (is.null(twl)) 0L else nrow(twl)

  tm <- system.time(r <- if (n_twl >= 50)
    run_sgat(twl, iters = 50000, thin = 20, prop_sd = 0.1, zenith = 94.276) else NULL)
  emit(table = "2", scenario = sc$label, method = "SGAT", rep = rep,
       rmse_km = if (is.null(r)) NA_real_ else
         rmse_to_truth(r$time, r$lon, r$lat, td$time, td$true_lat, td$true_lon),
       seconds = unname(tm["elapsed"]),
       effort = sprintf("50000 iters/thin 20, %d twilights", n_twl))
  if (rep == 1L && !is.null(r))
    fits_cache[[paste0(sc$label, "|SGAT")]] <- data.frame(time = r$time, lon = r$lon, lat = r$lat)

  tm <- system.time(r <- if (n_twl >= 50)
    run_flightr(td$time, lt, twl, c(130, 185, -85, -40), seal_start, 240, 50000) else NULL)
  emit(table = "2", scenario = sc$label, method = "FLightR", rep = rep,
       rmse_km = if (is.null(r)) NA_real_ else
         rmse_to_truth(r$time, r$lon, r$lat, td$time, td$true_lat, td$true_lon),
       seconds = unname(tm["elapsed"]),
       effort = sprintf("50000 particles, %d twilights", n_twl))
  if (rep == 1L && !is.null(r))
    fits_cache[[paste0(sc$label, "|FLightR")]] <- data.frame(time = r$time, lon = r$lon, lat = r$lat)

  if (rep == 1L) {
    fits_cache[[paste0(sc$label, "|truth")]] <-
      data.frame(time = td$time, lon = td$true_lon, lat = td$true_lat)
    saveRDS(fits_cache, FITS)
  }
}
saveRDS(fits_cache, FITS)

# =============================================================================
# Table 1: 180-day geometric track, ideal and adverse
# =============================================================================
simulate_geom <- function(seed, adverse) {
  set.seed(seed)
  n     <- 180 * 24 * 6
  times <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = "10 mins", length.out = n)
  s     <- seq(0, 180, length.out = n)
  lats  <- approx(cumsum(c(0, 30, 60, 60, 30)), c(-45, -60, -65, -45, -45), xout = s)$y
  lons  <- approx(cumsum(c(0, 30, 60, 60, 30)), c(140, 110, 80, 140, 140), xout = s)$y
  z     <- solar_zenith(as.numeric(times), lons, lats)
  l     <- pmin(pmax(558.5 - 5.818 * z, 0), 64)
  if (adverse) {
    hit    <- runif(n) < 0.05
    l[hit] <- pmax(0, l[hit] - runif(sum(hit), 10, 40))
  } else {
    l <- l * runif(n, 0.95, 1.05)
  }
  data.frame(time = times, true_lat = lats, true_lon = lons,
             light = pmin(pmax(l, 0), 64))
}
geom <- list(list(label = "Ideal (clear sky)",    seed = 42,  adverse = FALSE, prop_sd = 0.05),
             list(label = "Adverse (5% shading)", seed = 123, adverse = TRUE,  prop_sd = 0.10))

cat("\n=== Table 1: geometric track ===\n")
for (rep in seq_len(REPS_GEOM)) for (sc in geom) {
  trk <- simulate_geom(sc$seed, sc$adverse)
  st  <- c(trk$true_lon[1], trk$true_lat[1])

  tm <- system.time(f <- TwilightFreeSMC(
    date_time = trk$time, light = trk$light,
    start_lat = trk$true_lat[1],         start_lon = trk$true_lon[1],
    end_lat   = trk$true_lat[nrow(trk)], end_lon   = trk$true_lon[nrow(trk)],
    method = "guided", n_particles = 1000, seed = 42))
  kt <- as.POSIXct(f$knot_times, origin = "1970-01-01", tz = "UTC")
  emit(table = "1", scenario = sc$label, method = "invTF", rep = rep,
       rmse_km = rmse_to_truth(kt, f$lon, f$lat, trk$time, trk$true_lat, trk$true_lon),
       seconds = unname(tm["elapsed"]), effort = "1000 particles")

  twl   <- prep_twl(trk$time, trk$light, 5, NA)
  n_twl <- if (is.null(twl)) 0L else nrow(twl)

  tm <- system.time(r <- run_sgat(twl, iters = 100000, thin = 100,
                                  prop_sd = sc$prop_sd, zenith = 94.276))
  emit(table = "1", scenario = sc$label, method = "SGAT", rep = rep,
       rmse_km = if (is.null(r)) NA_real_ else
         rmse_to_truth(r$time, r$lon, r$lat, trk$time, trk$true_lat, trk$true_lon),
       seconds = unname(tm["elapsed"]),
       effort = sprintf("100000 iters/thin 100, %d twilights", n_twl))

  tm <- system.time(r <- run_flightr(trk$time, trk$light, twl,
                                     c(70, 180, -80, 0), st, 600, 50000))
  emit(table = "1", scenario = sc$label, method = "FLightR", rep = rep,
       rmse_km = if (is.null(r)) NA_real_ else
         rmse_to_truth(r$time, r$lon, r$lat, trk$time, trk$true_lat, trk$true_lon),
       seconds = unname(tm["elapsed"]),
       effort = sprintf("50000 particles, %d twilights", n_twl))
}

p_end <- probe()
cat(sprintf("\nmachine probe: start %.3f s, end %.3f s, drift %+.1f%%\n",
            p_start, p_end, 100 * (p_end - p_start) / p_start))
cat("results:", OUT, "\nfigure-1 fits:", FITS, "\n")
