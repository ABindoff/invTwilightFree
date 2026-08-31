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
# (a) Table 1 SGAT, original per-scenario settings
# -----------------------------------------------------------------------------
simulate_geom <- function(seed, adverse) {
  set.seed(seed)
  n     <- 180 * 24 * 6
  times <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = "10 mins", length.out = n)
  s     <- seq(0, 180, length.out = n)
  lats  <- approx(cumsum(c(0, 30, 60, 60, 30)), c(-45, -60, -65, -45, -45), xout = s)$y
  lons  <- approx(cumsum(c(0, 30, 60, 60, 30)), c(140, 110, 80, 140, 140), xout = s)$y
  l     <- pmin(pmax(558.5 - 5.818 * solar_zenith(as.numeric(times), lons, lats), 0), 64)
  if (adverse) {
    hit <- runif(n) < 0.05; l[hit] <- pmax(0, l[hit] - runif(sum(hit), 10, 40))
  } else {
    l <- l * runif(n, 0.95, 1.05)
  }
  data.frame(time = times, true_lat = lats, true_lon = lons, light = pmin(pmax(l, 0), 64))
}

geom <- list(
  list(label = "Ideal (clear sky)", seed = 42, adverse = FALSE,
       threshold = 2, alpha = c(0, 5),  prop_sd = 0.05, jitter = FALSE),
  list(label = "Adverse (5% shading)", seed = 123, adverse = TRUE,
       threshold = 5, alpha = c(0, 20), prop_sd = 0.10, jitter = TRUE))

cat("=== Table 1: SGAT, original per-scenario settings ===\n")
for (sc in geom) {
  trk <- simulate_geom(sc$seed, sc$adverse)
  twl <- TwGeos::findTwilights(data.frame(Date = trk$time, Light = trk$light),
                               threshold = sc$threshold, include = trk$time)
  set.seed(1)
  tm <- system.time(fit <- tryCatch({
    # the original drivers seed the chain at the true track
    x0 <- cbind(approx(as.numeric(trk$time), trk$true_lon, xout = as.numeric(twl$Twilight))$y,
                approx(as.numeric(trk$time), trk$true_lat, xout = as.numeric(twl$Twilight))$y)
    if (sc$jitter) x0 <- x0 + matrix(rnorm(length(x0), 0, 0.01), ncol = 2)
    m <- SGAT::thresholdModel(twl$Twilight, twl$Rise, twilight.model = "Normal",
           alpha = sc$alpha, beta = c(2.2, 0.08),
           logp.x = function(x) dnorm(x[, 2], 0, 90, log = TRUE),
           x0 = x0, z0 = twl$Twilight)
    p <- SGAT::mvnorm(S = diag(rep(sc$prop_sd, 2)^2), n = nrow(x0))
    SGAT::stellaMetropolis(m, p, x0 = x0, iters = 100000, thin = 100)
  }, error = function(e) {cat("    SGAT error:", conditionMessage(e), "\n"); NULL}))
  r <- if (is.null(fit)) NULL else {
    st <- if (!is.null(fit$model$time)) fit$model$time else fit$model$z
    list(time = st, lon = apply(fit$x[[1]][, 1, ], 1, mean),
         lat = apply(fit$x[[1]][, 2, ], 1, mean))
  }
  emit(table = "1", scenario = sc$label, method = "SGAT", rep = 2L,
       rmse_km = if (is.null(r)) NA_real_ else
         rmse_to_truth(r$time, r$lon, r$lat, trk$time, trk$true_lat, trk$true_lon),
       seconds = unname(tm["elapsed"]),
       effort = sprintf("100000 iters/thin 100, %d twilights, truth-initialised", nrow(twl)))
}

cat("
done
")
