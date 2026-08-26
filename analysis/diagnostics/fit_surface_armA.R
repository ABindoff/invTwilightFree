# ARM A on the surface-conditioned decimation. The new control.
#
# Identical to arm A_shipped of fit_final_29.R in every respect except the panel:
# same grid, same STEP_HOURS/DIFFUSION/CELL, same likelihood parameters
# (lambda 1.0, prob_slab 0.10), same scorer, same outputs. The ONLY difference is
# that light comes from panel_surface_v1.rds instead of panel_v1.rds.
#
# So the comparison against arm A_shipped isolates the decimation, which is the
# whole point: the shipped `max` over 30 minutes against a median of the surface
# samples on a drift-corrected depth zero.
#
# WHAT TO EXPECT, pre-registered so it cannot be rationalised afterwards:
#
#   * The calibration already moved in the predicted direction (pooled z50 91.92
#     -> 91.25 for Mk9_219), so the bias is real. What is NOT established is that
#     removing it improves POSITIONS.
#   * Arm B today is the cautionary case: removing shading made accuracy, bias and
#     coverage all worse, dose-dependently (rho -0.783 between fraction dropped and
#     change in bias), because the shaded observations had been partly cancelling
#     a bias of opposite sign. The max bias and the shading bias have opposite
#     signs, so removing the max bias could unmask the other one just as the gate
#     did.
#   * Therefore headline median km may WORSEN even though the input is more
#     honest. The informative quantities are the per-knot structure -- the
#     seasonal slope and the latitude slope -- not the median.
#   * The one thing that would falsify the whole account: if bias, seasonal
#     structure AND latitude slope are all unchanged, then the decimation was not
#     carrying the error and the search moves elsewhere.
#
# RESUMABLE: checkpointed on id.
SCRATCH <- "analysis/diagnostics"
source(file.path(SCRATCH, "nes_common.R"))
suppressMessages(devtools::load_all(".", quiet = TRUE))

P <- readRDS("analysis/cache/nes/panel_surface_v1.rds")
tags <- P$tags; responses <- P$responses
OUT_TAG  <- "scratch/nes_calibration/surfaceA_tags.csv"
OUT_KNOT <- "scratch/nes_calibration/surfaceA_knots.csv"
STEP_HOURS <- 12; DIFFUSION <- 110; CELL <- 1
PSLAB <- 0.10; LAM <- 1.00

grid <- rast(xmin = 150, xmax = 250, ymin = 20, ymax = 70,
             resolution = CELL, crs = "EPSG:4326")
values(grid) <- 1

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
  invisible(capture.output(
    f <- TwilightFreeGrid(L$time, pmax(0, L$light - r$baseline), grid = grid,
      start_lon = unname(g$p0[1]), start_lat = unname(g$p0[2]),
      end_lon = unname(g$p1[1]), end_lat = unname(g$p1[2]),
      step_hours = STEP_HOURS, diffusion = DIFFUSION,
      calibration = r$calibration,
      likelihood_params = c(LAM / (r$max_light * 0.5), r$max_light, PSLAB))))
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
  fwrite(data.table(id = tg, arm = "A_surface", season = g$season, family = g$family,
                    n = sum(ok), median_km = median(e[ok]),
                    bias_lat = mean(ep[ok]), rmse_lat = sqrt(mean(ep[ok]^2)),
                    lat_sd = mean(sdp[ok]),
                    cover_lat = mean(abs(ep[ok]) <= 1.96 * sdp[ok]),
                    log_z = f$log_z), OUT_TAG, append = file.exists(OUT_TAG))
  fwrite(data.table(id = tg, arm = "A_surface", time = tm[ok], true_lat = tr$lat[ok],
                    err_lat = ep[ok], err_km = e[ok], lat_sd = sdp[ok],
                    decl = solar_declination(tm[ok])),
         OUT_KNOT, append = file.exists(OUT_KNOT))
}
cat("\ndone\n")
