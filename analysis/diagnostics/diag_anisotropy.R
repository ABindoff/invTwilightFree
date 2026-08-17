# Does a per-axis movement prior calibrate both axes at once?
#
# THE PREDICTION, made before running. The prior fraction
# pi = sigma^-2/(sigma^-2 + tau^-2) is the share of what the posterior knows
# that comes from the prior. Measured under the isotropic prior it is 0.87 in
# latitude and 0.30 in longitude, because the light likelihood is 6.2x wider in
# latitude (333 km) than in longitude (54 km) while the prior is the same on
# both. That is why latitude intervals under-cover (0.66) and longitude
# intervals over-cover (~1.00), and why no single scale fixes both: they need
# opposite adjustments.
#
# Equalising the two prior fractions requires sigma_lat/sigma_lon = tau_lat/tau_lon
# = 6.2. Holding the geometric mean at the current 110 gives sigma_lat ~ 275,
# sigma_lon ~ 44. So the sweep brackets that ratio.
#
# What counts as success: latitude coverage RISES from 0.66 toward 0.95 while
# longitude coverage FALLS from ~1.00 toward 0.95, without median error
# degrading. If instead accuracy degrades as coverage improves, the anisotropy
# is just buying intervals with bias, which is what widening an isotropic prior
# already does and would not be worth having.
#
# RESUMABLE: appends per fit, skips what is done.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))
OUT <- "scratch/nes_calibration/anisotropy_results.csv"

STEP_HOURS <- 12; CELL <- 1; CAL_DAYS <- 15
CONFIGS <- list(
  list(nm = "iso 110/110 (current)", ns = 110, ew = 110),
  list(nm = "aniso 220/55  (4:1)",   ns = 220, ew = 55),
  list(nm = "aniso 275/44  (6:1)",   ns = 275, ew = 44),
  list(nm = "aniso 400/40 (10:1)",   ns = 400, ew = 40))

arch <- lapply(readRDS("analysis/cache/nes/archives_v3_30min.rds"), `[[`, "main")
meta <- nes_meta(); ar <- nes_argos()
tags <- list()
for (id in names(arch)) {
  m <- meta[match(id, meta$id), ]
  if (is.na(m$ptt) || !nzchar(m$ptt)) next
  df <- ar$deploy[Ptt == m$ptt]; a <- ar$fixes[Ptt == m$ptt]
  if (!nrow(df) || nrow(a) < 50) next
  x <- as.data.table(arch[[id]])
  t0 <- max(min(x$time), as.POSIXct(max(m$deploy_date, df$user_date), tz = "UTC"))
  t1 <- min(max(x$time), max(a$time))
  d <- x[time >= t0 & time <= t1]
  if (nrow(d) < 1000) next
  tf <- a[time >= max(time) - 72*3600]
  tags[[id]] <- list(id = id, light = d, argos = a,
                     p0 = c(df$deploy_lon, df$deploy_lat),
                     p1 = c(median(tf$lon), median(tf$lat)))
}
grid <- rast(xmin = 150, xmax = 250, ymin = 20, ymax = 70,
             resolution = CELL, crs = "EPSG:4326")
values(grid) <- 1
KM_DEG <- 111.32

score <- function(id, r, cfg) {
  g <- tags[[id]]
  t <- system.time(invisible(capture.output(
    f <- TwilightFreeGrid(g$light$time, pmax(0, g$light$light - r$baseline), grid = grid,
      start_lon = unname(g$p0[1]), start_lat = unname(g$p0[2]),
      end_lon = unname(g$p1[1]), end_lat = unname(g$p1[2]),
      step_hours = STEP_HOURS, diffusion = cfg$ns, diffusion_lon = cfg$ew,
      calibration = r$calibration,
      likelihood_params = c(1/(r$max_light*0.5), r$max_light, 0.10)))))[3]
  gp <- grid_posterior(f)
  sdp <- sdl <- numeric(nrow(gp$P))
  for (i in seq_len(nrow(gp$P))) {
    w <- gp$P[i, ]; s <- sum(w)
    if (!is.finite(s) || s <= 0) { sdp[i] <- sdl[i] <- NA; next }
    w <- w / s
    sdl[i] <- sqrt(sum(w * (gp$lon - sum(w * gp$lon))^2))
    sdp[i] <- sqrt(sum(w * (gp$lat - sum(w * gp$lat))^2))
  }
  tm <- as.POSIXct(f$fit$time, origin = "1970-01-01", tz = "UTC")
  tr <- argos_at(g$argos, tm)
  ep <- f$fit$lat - tr$lat; el <- dlon(lon360(f$fit$lon), tr$lon)
  e <- gc_km(lon360(f$fit$lon), f$fit$lat, tr$lon, tr$lat)
  ok <- is.finite(ep) & is.finite(sdp) & sdp > 0 & is.finite(sdl) & sdl > 0
  # realised prior fractions, from the errors this configuration produced
  tau_lat <- sd(ep[ok]) * KM_DEG
  tau_lon <- sd(el[ok]) * KM_DEG * cos(mean(tr$lat[ok], na.rm = TRUE) * pi/180)
  s_ns <- cfg$ns * sqrt(STEP_HOURS/24); s_ew <- cfg$ew * sqrt(STEP_HOURS/24)
  data.frame(id = id, config = cfg$nm, ns = cfg$ns, ew = cfg$ew,
    median_km = round(median(e[ok])),
    rmse_lat = round(sqrt(mean(ep[ok]^2)), 2),
    rmse_lon = round(sqrt(mean(el[ok]^2)), 2),
    cover_lat = round(mean(abs(ep[ok]) <= 1.96*sdp[ok]), 2),
    cover_lon = round(mean(abs(el[ok]) <= 1.96*sdl[ok]), 2),
    pi_lat = round((1/s_ns^2)/((1/s_ns^2)+(1/tau_lat^2)), 3),
    pi_lon = round((1/s_ew^2)/((1/s_ew^2)+(1/tau_lon^2)), 3),
    log_z = round(f$log_z, 1), secs = round(t), row.names = NULL)
}

scale_fits <- lapply(tags, function(g) {
  k <- g$light$time <= min(g$light$time) + CAL_DAYS * 86400
  fit_light_response(g$light$time[k], g$light$light[k], g$p0[1], g$p0[2])
})
geom_fits <- lapply(tags, function(g) {
  dep <- departure_time(g$light)
  k <- is.finite(dep) & as.numeric(g$light$time) < dep
  if (sum(k) < 200) return(NULL)
  fit_light_response(g$light$time[k], g$light$light[k], g$p0[1], g$p0[2])
})
responses <- pool_light_responses(geom_fits, scale_fits)

done <- if (file.exists(OUT)) {
  d <- data.table::fread(OUT); paste(d$id, d$config) } else character(0)
cat(sprintf("%d fits done, %d to go\n", length(done),
            length(tags)*length(CONFIGS) - length(done)))
for (cfg in CONFIGS) for (id in names(tags)) {
  if (paste(id, cfg$nm) %in% done) next
  message(cfg$nm, " : ", id)
  data.table::fwrite(score(id, responses[[id]], cfg), OUT, append = file.exists(OUT))
}

o <- data.table::fread(OUT)
cat("\n=== pooled over tags ===\n")
print(as.data.frame(o[, .(tags = .N, median_km = round(mean(median_km)),
  rmse_lat = round(mean(rmse_lat), 2), rmse_lon = round(mean(rmse_lon), 2),
  cover_lat = round(mean(cover_lat), 2), cover_lon = round(mean(cover_lon), 2),
  pi_lat = round(mean(pi_lat), 3), pi_lon = round(mean(pi_lon), 3),
  mean_log_z = round(mean(log_z), 1)), by = config][order(-cover_lat)]),
  row.names = FALSE)
cat("\nSuccess = latitude coverage up toward 0.95 AND longitude coverage down\n")
cat("toward 0.95, with median error not degrading. Watch pi_lat and pi_lon\n")
cat("converge; that is the mechanism, the coverages are the consequence.\n")
