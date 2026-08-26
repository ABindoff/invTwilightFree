# IS THE 1-DEGREE GRID BUYING ANYTHING? Arm A at CELL = 2.
#
# The grid HMM's transition step is a dense all-pairs double loop, so cost is
# O(n_cells^2) -- measured exponent 1.99, with 325/1250/5000 cells taking
# 0.24/2.70/42.58 s on the same window. Halving the resolution is therefore 16x
# faster for free, with no code change and no numerical-validity risk.
#
# The arithmetic says the resolution is probably wasted. Grid quantisation error is
# h/sqrt(12), so 1 deg contributes 32 km and 2 deg contributes 64 km, against a
# measured latitude error of 250 km. Added in quadrature that is 252 vs 258 km, so
# about 3% for a 16x speedup. This measures it rather than assuming it.
#
# Identical to arm A_shipped of fit_final_29.R in every respect except CELL, and
# reads the same panel_v1, so the difference is attributable to resolution alone.
#
# RESUMABLE: checkpointed on id.
SCRATCH <- "analysis/diagnostics"
source(file.path(SCRATCH, "nes_common.R"))
suppressMessages(devtools::load_all(".", quiet = TRUE))

P <- readRDS("analysis/cache/nes/panel_v1.rds")
tags <- P$tags; responses <- P$responses
OUT_TAG  <- "scratch/nes_calibration/cell2A_tags.csv"
OUT_KNOT <- "scratch/nes_calibration/cell2A_knots.csv"
STEP_HOURS <- 12; DIFFUSION <- 110; CELL <- 2
PSLAB <- 0.10; LAM <- 1.00

grid <- rast(xmin = 150, xmax = 250, ymin = 20, ymax = 70,
             resolution = CELL, crs = "EPSG:4326")
values(grid) <- 1
cat(sprintf("grid: %d cells at %.0f deg (arm A used %d at 1 deg)\n\n",
            terra::ncell(grid), CELL, 5000L))

done <- if (file.exists(OUT_TAG)) {
  d0 <- fread(OUT_TAG, colClasses = list(character = "id")); d0$id
} else character(0)
cat(sprintf("%d tags; %d done\n\n", length(tags), length(done)))

for (tg in names(tags)) {
  if (tg %in% done) next
  g <- tags[[tg]]; r <- responses[[tg]]
  if (is.null(r)) { message(tg, ": no response"); next }
  L <- as.data.table(g$light)
  message("fit: ", tg)
  el <- system.time(
    invisible(capture.output(
      f <- TwilightFreeGrid(L$time, pmax(0, L$light - r$baseline), grid = grid,
        start_lon = unname(g$p0[1]), start_lat = unname(g$p0[2]),
        end_lon = unname(g$p1[1]), end_lat = unname(g$p1[2]),
        step_hours = STEP_HOURS, diffusion = DIFFUSION,
        calibration = r$calibration,
        likelihood_params = c(LAM / (r$max_light * 0.5), r$max_light, PSLAB)))))[["elapsed"]]
  gp <- grid_posterior(f)
  sdp <- numeric(nrow(gp$P))
  for (i in seq_len(nrow(gp$P))) {
    w <- gp$P[i, ]; s <- sum(w)
    if (!is.finite(s) || s <= 0) { sdp[i] <- NA; next }
    w <- w / s; sdp[i] <- sqrt(sum(w * (gp$lat - sum(w * gp$lat))^2))
  }
  tm <- as.POSIXct(f$fit$time, origin = "1970-01-01", tz = "UTC")
  tr <- argos_at(as.data.table(g$argos), tm)
  ep <- f$fit$lat - tr$lat
  e  <- gc_km(lon360(f$fit$lon), f$fit$lat, tr$lon, tr$lat)
  ok <- is.finite(ep) & is.finite(sdp) & sdp > 0
  fwrite(data.table(id = tg, arm = "A_cell2", season = g$season, family = g$family,
                    n = sum(ok), median_km = median(e[ok]),
                    bias_lat = mean(ep[ok]), rmse_lat = sqrt(mean(ep[ok]^2)),
                    lat_sd = mean(sdp[ok]),
                    cover_lat = mean(abs(ep[ok]) <= 1.96 * sdp[ok]),
                    log_z = f$log_z, secs = el), OUT_TAG, append = file.exists(OUT_TAG))
  fwrite(data.table(id = tg, arm = "A_cell2", time = tm[ok], true_lat = tr$lat[ok],
                    err_lat = ep[ok], err_km = e[ok], lat_sd = sdp[ok],
                    decl = solar_declination(tm[ok])),
         OUT_KNOT, append = file.exists(OUT_KNOT))
}
cat("\ndone\n")
