# =============================================================================
# Completes the §2.10 re-run in two parts, appending to the same results file.
#
# (a) Table 1 SGAT. The geometric drivers use per-scenario settings that differ
#     from the seal panel: threshold 2 or 5 with no thinning, alpha c(0,5) or
#     c(0,20), beta c(2.2,0.08), a weak latitude prior, and a proposal sd of
#     0.05 or 0.1. They also initialise the chain AT THE TRUE TRACK, which is
#     not information an analyst has; that is reproduced here so the numbers
#     are comparable with the cached ones, and reported as a caveat.
#
# (b) Matched Monte Carlo effort. Table 2 compares invTwilightFree at 1000
#     particles with FLightR at 50,000. This sweeps the particle count on the
#     seal scenarios so the accuracy difference can be read against effort
#     rather than confounded with it.
#
#   R-4.6.0/bin/Rscript inst/paper/rerun_finish.R
# =============================================================================

for (v in c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
            "RAYON_NUM_THREADS", "NUMEXPR_NUM_THREADS"))
  do.call(Sys.setenv, setNames(list("1"), v))

suppressMessages({library(invTwilightFree); library(TwGeos); library(SGAT)})

OUT <- "inst/paper/benchmark_rerun.csv"

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
emit <- function(...) {
  row <- data.frame(..., stringsAsFactors = FALSE)
  write.table(row, OUT, sep = ",", row.names = FALSE, na = "",
              col.names = !file.exists(OUT), append = file.exists(OUT))
  cat(sprintf("  %-9s %-22s rep%d  %8s km  %8.1f s  [%s]\n", row$method, row$scenario,
              row$rep, ifelse(is.na(row$rmse_km), "FAIL", sprintf("%.0f", row$rmse_km)),
              row$seconds, row$effort))
}

# -----------------------------------------------------------------------------
# (b) invTwilightFree against particle count, seal scenarios
# -----------------------------------------------------------------------------
td   <- readRDS("scratch/simulated_seal_light_scenarios.rds")
seal <- list(list(label = "Cloudy", col = "light_ideal"),
             list(label = "Shaded (ARS diving)", col = "light_shaded"),
             list(label = "ALAN near colony", col = "light_alan"))

cat("\n=== invTwilightFree against Monte Carlo effort ===\n")
for (np in c(1000L, 5000L, 20000L, 50000L)) for (sc in seal) {
  tm <- system.time(f <- TwilightFreeSMC(
    date_time = td$time, light = td[[sc$col]],
    start_lat = td$true_lat[1],        start_lon = td$true_lon[1],
    end_lat   = td$true_lat[nrow(td)], end_lon   = td$true_lon[nrow(td)],
    method = "guided", n_particles = np, seed = 42))
  kt <- as.POSIXct(f$knot_times, origin = "1970-01-01", tz = "UTC")
  emit(table = "2-effort", scenario = sc$label, method = "invTF", rep = 1L,
       rmse_km = rmse_to_truth(kt, f$lon, f$lat, td$time, td$true_lat, td$true_lon),
       seconds = unname(tm["elapsed"]), effort = sprintf("%d particles", np))
}
cat("\ndone\n")
