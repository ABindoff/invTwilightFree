# Is the reported accuracy partly an artefact of the search domain?
#
# The report builds its grid from the Argos envelope plus 8/4/6/8 degrees of
# padding. That is ground truth choosing where the answer may live, and the
# latitude bounds are the ones that matter: a light-based estimator's
# characteristic failure is a wild latitude excursion near the equinox, and a
# box drawn round the true track clips exactly that.
#
# Compare against a domain that uses NO per-animal truth. The legitimate prior
# knowledge is: a northern elephant seal from Ano Nuevo, in the North Pacific,
# bounded by continents. That is a statement about the species and the ocean,
# not about these ten tracks.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))

STEP_HOURS <- 12; DIFFUSION <- 110; CELL <- 1; CAL_DAYS <- 15

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

all_lon <- unlist(lapply(tags, function(g) g$argos$lon))
all_lat <- unlist(lapply(tags, function(g) g$argos$lat))
cat(sprintf("Argos envelope: lon %.1f..%.1f  lat %.1f..%.1f\n",
            min(all_lon), max(all_lon), min(all_lat), max(all_lat)))

DOMAINS <- list(
  # what the report does: the truth envelope, padded
  "A argos envelope + pad (report)" = c(
    floor(min(all_lon)) - 8, ceiling(max(all_lon)) + 4,
    floor(min(all_lat)) - 6, ceiling(max(all_lat)) + 8),
  # no per-animal truth: the North Pacific, from the dateline region to the
  # American coast, subtropics to the sub-Arctic
  "B north pacific (no truth)"      = c(150, 250, 20, 70),
  # wider still, to see whether the estimator is merely being fenced in
  "C wide north pacific"            = c(130, 255, 10, 75))

score <- function(id, r, grid, label) {
  g <- tags[[id]]
  t <- system.time(invisible(capture.output(
    f <- TwilightFreeGrid(g$light$time, pmax(0, g$light$light - r$baseline), grid = grid,
      start_lon = unname(g$p0[1]), start_lat = unname(g$p0[2]),
      end_lon = unname(g$p1[1]), end_lat = unname(g$p1[2]),
      step_hours = STEP_HOURS, diffusion = DIFFUSION, calibration = r$calibration,
      likelihood_params = c(1/(r$max_light*0.5), r$max_light, 0.10)))))[3]
  tm <- as.POSIXct(f$fit$time, origin = "1970-01-01", tz = "UTC")
  tr <- argos_at(g$argos, tm)
  e <- gc_km(lon360(f$fit$lon), f$fit$lat, tr$lon, tr$lat); ok <- is.finite(e)
  # how often does the estimate sit within one cell of an edge? a track pressed
  # against the boundary is being held there by the box, not by the light
  edge <- mean(f$fit$lat <= ymin_g + CELL | f$fit$lat >= ymax_g - CELL, na.rm = TRUE)
  data.frame(id = id, domain = label, median_km = round(median(e[ok])),
    q90_km = round(unname(quantile(e[ok], .9))),
    bias_lat = round(mean((f$fit$lat - tr$lat)[ok]), 2),
    rmse_lat = round(sqrt(mean(((f$fit$lat - tr$lat)[ok])^2)), 2),
    lat_min = round(min(f$fit$lat, na.rm = TRUE), 1),
    lat_max = round(max(f$fit$lat, na.rm = TRUE), 1),
    at_lat_edge = round(edge, 3), secs = round(t), row.names = NULL)
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

res <- list()
for (nm in names(DOMAINS)) {
  b <- DOMAINS[[nm]]
  ymin_g <<- b[3]; ymax_g <<- b[4]
  grid <- rast(xmin = b[1], xmax = b[2], ymin = b[3], ymax = b[4],
               resolution = CELL, crs = "EPSG:4326")
  values(grid) <- 1
  cat(sprintf("\n%s: lon %.0f..%.0f lat %.0f..%.0f = %d cells\n",
              nm, b[1], b[2], b[3], b[4], ncell(grid)))
  for (id in names(tags)) {
    message(nm, " : ", id)
    res[[length(res)+1]] <- score(id, responses[[id]], grid, nm)
  }
}
o <- as.data.table(do.call(rbind, res))
cat("\n=== per tag ===\n"); print(as.data.frame(o[order(id, domain)]), row.names = FALSE)
cat("\n=== pooled ===\n")
print(as.data.frame(o[, .(tags = .N, median_km = round(mean(median_km)),
  worst_km = max(median_km), q90_km = round(mean(q90_km)),
  bias_lat = round(mean(bias_lat), 2), rmse_lat = round(mean(rmse_lat), 2),
  frac_at_lat_edge = round(mean(at_lat_edge), 3),
  mins = round(sum(secs)/60, 1)), by = domain][order(domain)]), row.names = FALSE)
saveRDS(o, file.path(SCRATCH, "domain.rds"))
