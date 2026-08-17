# What is a "guided Brownian" worth, before anyone writes the sampler?
#
# The idea: pull the track gently toward where the animals are known to spend
# time. As a diffusion that is the Langevin model, whose stationary distribution
# IS the occupancy field. The sampler is the expensive part, but the CEILING can
# be measured first with machinery the package already has: put the occupancy
# field in as a spatial prior via prior_raster() + custom_rule().
#
# The trap this design avoids: an occupancy field estimated from the SAME track
# it guides is a feedback loop, and this project has two measured examples of
# that topology diverging (the light-response refit ran 321 -> 391 -> 755 km).
# So each animal's field is built from the OTHER NINE only. That is leave-one-
# out, it breaks the loop, and it is what a real study has: previous
# deployments. No animal ever sees its own truth.
#
# Read the result as an UPPER BOUND: nine Argos tracks are a far better
# occupancy estimate than nine light-based ones would be.
#
# `lat_span` is reported next to `rmse_lat` on purpose. An attraction term can
# lower latitude error two ways -- by adding information, or by shrinking
# everything toward a blob. The second would look like success and be the same
# mistake as buying coverage by widening intervals. If lat_span collapses while
# rmse_lat falls, it is shrinkage, not information.
#
# RESUMABLE: appends per fit to residence_results.csv, skips what is done.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))
OUT <- "scratch/nes_calibration/residence_results.csv"

STEP_HOURS <- 12; DIFFUSION <- 110; CELL <- 1; CAL_DAYS <- 15
BANDWIDTH <- 3
WEIGHTS   <- c(0, 0.25, 0.5, 1)

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
proto <- rast(xmin = 150, xmax = 250, ymin = 20, ymax = 70,
              resolution = CELL, crs = "EPSG:4326")
values(proto) <- 1

# Occupancy field from a set of tracks: a Gaussian kernel density on the grid,
# returned as a LOG density so it adds to the light log-likelihood. Floored at a
# low quantile so an animal is never forbidden from anywhere outright.
occupancy_log <- function(lon, lat, r, bw = BANDWIDTH, floor_q = 0.02) {
  xy <- xyFromCell(r, seq_len(ncell(r)))
  d <- numeric(nrow(xy))
  for (i in seq_along(lon)) {
    dx <- ((xy[, 1] - lon[i] + 180) %% 360 - 180) * cos(xy[, 2] * pi/180)
    dy <- xy[, 2] - lat[i]
    d <- d + exp(-0.5 * (dx^2 + dy^2) / bw^2)
  }
  d <- d / sum(d)
  fl <- quantile(d[d > 0], floor_q)
  d[d < fl] <- fl
  out <- r; values(out) <- log(d / sum(d))
  out
}

score <- function(id, r_resp, term, label) {
  g <- tags[[id]]
  t <- system.time(invisible(capture.output(
    f <- TwilightFreeGrid(g$light$time, pmax(0, g$light$light - r_resp$baseline),
      grid = proto,
      start_lon = unname(g$p0[1]), start_lat = unname(g$p0[2]),
      end_lon = unname(g$p1[1]), end_lat = unname(g$p1[2]),
      step_hours = STEP_HOURS, diffusion = DIFFUSION,
      terms = if (is.null(term)) list() else list(term),
      calibration = r_resp$calibration,
      likelihood_params = c(1/(r_resp$max_light*0.5), r_resp$max_light, 0.10))))) [3]
  gp <- grid_posterior(f)
  sdp <- numeric(nrow(gp$P))
  for (i in seq_len(nrow(gp$P))) {
    w <- gp$P[i, ]; s <- sum(w)
    if (!is.finite(s) || s <= 0) { sdp[i] <- NA; next }
    w <- w / s
    sdp[i] <- sqrt(sum(w * (gp$lat - sum(w * gp$lat))^2))
  }
  tm <- as.POSIXct(f$fit$time, origin = "1970-01-01", tz = "UTC")
  tr <- argos_at(g$argos, tm)
  ep <- f$fit$lat - tr$lat
  e <- gc_km(lon360(f$fit$lon), f$fit$lat, tr$lon, tr$lat)
  ok <- is.finite(ep) & is.finite(sdp) & sdp > 0
  data.frame(id = id, weight = label, median_km = round(median(e[ok])),
    q90_km = round(unname(quantile(e[ok], .9))),
    bias_lat = round(mean(ep[ok]), 2),
    rmse_lat = round(sqrt(mean(ep[ok]^2)), 2),
    rmse_lon = round(sqrt(mean((dlon(lon360(f$fit$lon), tr$lon)[ok])^2)), 2),
    post_sd_lat = round(median(sdp[ok]), 2),
    cover_lat = round(mean(abs(ep[ok]) <= 1.96 * sdp[ok]), 2),
    lat_span = round(diff(range(f$fit$lat, na.rm = TRUE)), 1),
    true_span = round(diff(range(tr$lat, na.rm = TRUE)), 1),
    secs = round(t), row.names = NULL)
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
  d <- data.table::fread(OUT); paste(d$id, d$weight) } else character(0)
cat(sprintf("%d fits already on disk; %d to do\n",
            length(done), length(tags) * length(WEIGHTS) - length(done)))

for (id in names(tags)) {
  todo <- WEIGHTS[!(paste(id, sprintf("%.2f", WEIGHTS)) %in% done)]
  if (!length(todo)) next
  others <- setdiff(names(tags), id)
  ol <- unlist(lapply(others, function(o) tags[[o]]$argos$lon))
  oa <- unlist(lapply(others, function(o) tags[[o]]$argos$lat))
  keep <- seq(1, length(ol), by = max(1, floor(length(ol) / 4000)))
  fld <- occupancy_log(ol[keep], oa[keep], proto)
  for (w in todo) {
    message(id, " weight ", w)
    term <- if (w == 0) NULL else location_term(
      name = "residence", source = prior_raster(fld, is_log = TRUE),
      rule = custom_rule(function(obs, expected) expected), weight = w)
    row <- score(id, responses[[id]], term, sprintf("%.2f", w))
    data.table::fwrite(row, OUT, append = file.exists(OUT))
  }
}

o <- data.table::fread(OUT)
cat("\n=== what a leave-one-out occupancy prior is worth (upper bound) ===\n")
print(as.data.frame(o[, .(tags = .N, median_km = round(mean(median_km)),
  worst_km = max(median_km), q90_km = round(mean(q90_km)),
  bias_lat = round(mean(bias_lat), 2), rmse_lat = round(mean(rmse_lat), 2),
  rmse_lon = round(mean(rmse_lon), 2), cover_lat = round(mean(cover_lat), 2),
  lat_span = round(mean(lat_span), 1), true_span = round(mean(true_span), 1),
  mins = round(sum(secs)/60, 1)), by = weight][order(weight)]), row.names = FALSE)
cat("\nIf rmse_lat falls while lat_span collapses below true_span, that is\n")
cat("shrinkage to a blob rather than information.\n")
