# Is behavioural state-switching worth having here, before anyone wires dive
# effort into it?
#
# The grid HMM already supports it: pass a VECTOR of diffusions and a transition
# matrix and it runs a latent state chain, marginalising the state exactly. So
# the first-order question costs nothing to answer.
#
# What to expect, and why the answer is not obvious. Accuracy is FLAT across an
# eightfold change in a single diffusion (282 km at D=55 against 223 at D=440),
# so switching is very unlikely to move the point estimate. But COVERAGE is not
# flat -- it ran 0.21 to 0.93 over that same sweep, and latitude coverage is the
# open problem (0.66 against a nominal 0.95). A model that can be slow where the
# animal is slow and fast where it is fast should be better CALIBRATED even if
# it is no more accurate, because a single scale has to be wrong in both places
# at once. So a result of "no accuracy change, better coverage" is success.
#
# `log_z` is recorded as well, which asks a different question: do the DATA
# prefer the switching model, regardless of what Argos says about accuracy.
#
# RESUMABLE. Every fit is appended to states_results.csv as it finishes, and a
# restart skips what is already there. A five-hour run should not be hostage to
# a closing laptop.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))
OUT <- "scratch/nes_calibration/states_results.csv"

STEP_HOURS <- 12; CELL <- 1; CAL_DAYS <- 15
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
# the truth-free domain the report now uses
grid <- rast(xmin = 150, xmax = 250, ymin = 20, ymax = 70,
             resolution = CELL, crs = "EPSG:4326")
values(grid) <- 1

# Argos says the true one-step scale is about 31 km / 12 h, i.e. D ~ 44
# km/sqrt(day); the report runs a single D = 110. The pairs below bracket that:
# a slow state near the ARS scale and a fast one near transit.
CONFIGS <- list(
  "1-state D=110 (report)"   = list(d = 110,        tp = NULL),
  "1-state D=44 (true step)" = list(d = 44,         tp = NULL),
  "2-state 44/220"           = list(d = c(44, 220), tp = c(0.9, 0.1, 0.1, 0.9)),
  "2-state 30/300"           = list(d = c(30, 300), tp = c(0.9, 0.1, 0.1, 0.9)),
  "2-state 44/220 sticky"    = list(d = c(44, 220), tp = c(0.97, 0.03, 0.03, 0.97)),
  "3-state 20/110/300"       = list(d = c(20, 110, 300),
                                    tp = c(0.9,0.05,0.05, 0.05,0.9,0.05, 0.05,0.05,0.9)))

score <- function(id, r, cfg, label) {
  g <- tags[[id]]
  t <- system.time(invisible(capture.output(
    f <- TwilightFreeGrid(g$light$time, pmax(0, g$light$light - r$baseline), grid = grid,
      start_lon = unname(g$p0[1]), start_lat = unname(g$p0[2]),
      end_lon = unname(g$p1[1]), end_lat = unname(g$p1[2]),
      step_hours = STEP_HOURS, diffusion = cfg$d, trans_prob = cfg$tp,
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
  ok <- is.finite(ep) & is.finite(sdp) & sdp > 0
  data.frame(id = id, config = label, median_km = round(median(e[ok])),
    rmse_lat = round(sqrt(mean(ep[ok]^2)), 2),
    post_sd_lat = round(median(sdp[ok]), 2),
    post_sd_lon = round(median(sdl[ok]), 2),
    cover_lat = round(mean(abs(ep[ok]) <= 1.96 * sdp[ok]), 2),
    cover_lon = round(mean(abs(el[ok]) <= 1.96 * sdl[ok]), 2),
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
cat(sprintf("%d fits already on disk; %d to do\n",
            length(done), length(tags) * length(CONFIGS) - length(done)))

for (nm in names(CONFIGS)) for (id in names(tags)) {
  if (paste(id, nm) %in% done) next
  message(nm, " : ", id)
  row <- score(id, responses[[id]], CONFIGS[[nm]], nm)
  data.table::fwrite(row, OUT, append = file.exists(OUT))
}

o <- data.table::fread(OUT)
cat("\n=== pooled ===\n")
print(as.data.frame(o[, .(tags = .N, median_km = round(mean(median_km)),
  worst_km = max(median_km), rmse_lat = round(mean(rmse_lat), 2),
  post_sd_lat = round(mean(post_sd_lat), 2),
  cover_lat = round(mean(cover_lat), 2), cover_lon = round(mean(cover_lon), 2),
  mean_log_z = round(mean(log_z), 1), mins = round(sum(secs)/60, 1)),
  by = config][order(config)]), row.names = FALSE)
cat("\nlog_z is the marginal likelihood: whether the DATA prefer switching,\n")
cat("independently of whether Argos says it is more accurate.\n")
