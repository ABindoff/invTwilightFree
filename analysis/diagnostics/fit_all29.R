# Fit all 29 double-tagged deployments (2021 + the 2022/2023 delivery) with the
# validated recipe, and write per-knot errors. This is the gateway run: every
# question that stalled at n = 10 and p ~ 0.1 becomes answerable from its output.
#
# THE RECIPE, unchanged from the 2021 analysis and validated there:
#   geometry  from the pre-departure haul-out (departure dated from the dive
#             record, no ground truth), POOLED across deployments
#   scale     from each deployment's own first 15 days
#   tangent   amplitude tied to max_light
#   domain    150-250 E, 20-70 N, defined from the species range not the tracks
#
# ONE CHANGE, forced by the new data: pooling is now WITHIN TAG FAMILY. The 2021
# panel was all Mk9 219xxxx, so "pool across tags" and "pool within model" were
# the same thing. The delivery adds four 18A deployments, and response geometry
# is a property of the light channel, so pooling across families would impose one
# family's curve on the other. With 25 Mk9 and 4 18A there is enough to pool each
# separately.
#
# RESUMABLE: appends per deployment.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))
OUT_KNOT <- "scratch/nes_calibration/all29_knots.csv"
OUT_TAG  <- "scratch/nes_calibration/all29_tags.csv"

STEP_HOURS <- 12; DIFFUSION <- 110; CELL <- 1; CAL_DAYS <- 15

# ---- assemble the panel -----------------------------------------------------
old_arch <- lapply(readRDS("analysis/cache/nes/archives_v3_30min.rds"), `[[`, "main")
new_arch <- lapply(readRDS("analysis/cache/nes/fvilches_archives_v1_30min.rds"), `[[`, "main")
new_argos <- readRDS("analysis/cache/nes/fvilches_argos_v1.rds")
new_man <- fread("scratch/nes_calibration/fvilches_manifest.csv", colClasses = list(character = "topp"))
meta <- nes_meta(); ar <- nes_argos()

# tag serial per deployment, for family-wise pooling
ser <- rbind(
  data.table(topp = names(old_arch),
             serial = sub("^[0-9]+_(.+)\\.csv$", "\\1",
                          basename(list.files("data/nes_untracked", pattern = "^[0-9]+_.*\\.csv$")))[
                            match(names(old_arch),
                                  sub("^([0-9]+)_.*$", "\\1",
                                      basename(list.files("data/nes_untracked",
                                                          pattern = "^[0-9]+_.*\\.csv$"))))]),
  data.table(topp = sub("^([0-9]+)_.*$", "\\1",
                        list.files("fvilches/extracted/TDR raw", pattern = "Archive[.]csv$")),
             serial = sub("^[0-9]+_(.+)-out-Archive[.]csv$", "\\1",
                          list.files("fvilches/extracted/TDR raw", pattern = "Archive[.]csv$"))),
  fill = TRUE)
ser <- unique(ser, by = "topp")
family <- function(s) fifelse(is.na(s), "unknown",
                       fifelse(grepl("^219", s), "Mk9_219",
                        fifelse(grepl("^18A", s), "F18A", "other")))

tags <- list()
for (tg in names(old_arch)) {                       # 2021, from the first delivery
  m <- meta[match(tg, meta$id), ]
  if (is.na(m$ptt) || !nzchar(m$ptt)) next
  df <- ar$deploy[Ptt == m$ptt]; a <- ar$fixes[Ptt == m$ptt]
  if (!nrow(df) || nrow(a) < 50) next
  x <- as.data.table(old_arch[[tg]])
  t0 <- max(min(x$time), as.POSIXct(max(m$deploy_date, df$user_date), tz = "UTC"))
  t1 <- min(max(x$time), max(a$time))
  d <- x[time >= t0 & time <= t1]
  if (nrow(d) < 1000) next
  tf <- a[time >= max(time) - 72*3600]
  tags[[tg]] <- list(id = tg, season = "2021", light = d, argos = a,
                     p0 = c(df$deploy_lon, df$deploy_lat),
                     p1 = c(median(tf$lon), median(tf$lat)))
}
for (tg in names(new_arch)) {                       # 2022/2023, already clipped
  mm <- new_man[topp == tg]
  if (!nrow(mm)) next
  a <- as.data.table(new_argos[[tg]])
  d <- as.data.table(new_arch[[tg]])
  d <- d[time >= min(a$time) & time <= max(a$time)]
  if (nrow(d) < 1000 || nrow(a) < 50) next
  # fread returns season as integer while the 2021 branch sets it as character;
  # vapply over the mixed panel then fails on type. Coerce at the source.
  tags[[tg]] <- list(id = tg, season = as.character(mm$season[1]), light = d, argos = a,
                     p0 = c(mm$deploy_lon[1], mm$deploy_lat[1]),
                     p1 = c(mm$recover_lon[1], mm$recover_lat[1]))
}
for (tg in names(tags)) tags[[tg]]$family <- family(ser[topp == tg]$serial[1])
cat(sprintf("panel: %d deployments\n", length(tags)))
print(table(vapply(tags, function(g) g$season, ""),
            vapply(tags, function(g) g$family, "")))

# ---- calibration ------------------------------------------------------------
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
cat(sprintf("\nhaul-out geometry available for %d of %d deployments\n",
            sum(!vapply(geom_fits, is.null, logical(1))), length(tags)))

# pool within family
responses <- list()
for (fam in unique(vapply(tags, function(g) g$family, ""))) {
  ids <- names(tags)[vapply(tags, function(g) g$family, "") == fam]
  gf <- geom_fits[ids]; sf <- scale_fits[ids]
  if (!sum(!vapply(gf, is.null, logical(1)))) {
    message("family ", fam, ": no haul-out geometry, skipped"); next
  }
  pooled <- pool_light_responses(gf, sf)
  p <- attr(pooled, "pooled")
  cat(sprintf("family %-8s: %d deployments, geometry pooled from %d -> z50 %.2f, width %.1f\n",
              fam, length(ids), p[["n"]], p[["z50"]], 4.394 * p[["scale"]]))
  for (i in ids) responses[[i]] <- pooled[[i]]
}

grid <- rast(xmin = 150, xmax = 250, ymin = 20, ymax = 70,
             resolution = CELL, crs = "EPSG:4326")
values(grid) <- 1

done <- if (file.exists(OUT_TAG)) as.character(fread(OUT_TAG)$id) else character(0)
cat(sprintf("\n%d fitted, %d to go\n\n", length(done), length(tags) - length(done)))

for (tg in names(tags)) {
  if (tg %in% done) next
  g <- tags[[tg]]; r <- responses[[tg]]
  if (is.null(r)) { message(tg, ": no response, skipped"); next }
  message("fit: ", tg)
  t <- system.time(invisible(capture.output(
    f <- TwilightFreeGrid(g$light$time, pmax(0, g$light$light - r$baseline), grid = grid,
      start_lon = unname(g$p0[1]), start_lat = unname(g$p0[2]),
      end_lon = unname(g$p1[1]), end_lat = unname(g$p1[2]),
      step_hours = STEP_HOURS, diffusion = DIFFUSION, calibration = r$calibration,
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
  fwrite(data.table(id = tg, season = g$season, family = g$family, time = tm,
                    est_lon = lon360(f$fit$lon), est_lat = f$fit$lat,
                    true_lon = tr$lon, true_lat = tr$lat,
                    err_km = e, err_lat = ep, err_lon = el,
                    lat_sd = sdp, lon_sd = sdl),
         OUT_KNOT, append = file.exists(OUT_KNOT))
  ok <- is.finite(ep) & is.finite(sdp) & sdp > 0
  fwrite(data.table(id = tg, season = g$season, family = g$family,
                    n_knots = length(tm), n_scored = sum(ok),
                    own_geometry = !is.null(geom_fits[[tg]]),
                    median_km = round(median(e[ok])),
                    rmse_lat = round(sqrt(mean(ep[ok]^2)), 3),
                    rmse_lon = round(sqrt(mean(el[ok]^2)), 3),
                    bias_lat = round(mean(ep[ok]), 3),
                    cover_lat = round(mean(abs(ep[ok]) <= 1.96*sdp[ok]), 3),
                    cover_lon = round(mean(abs(el[ok]) <= 1.96*sdl[ok]), 3),
                    log_z = round(f$log_z, 1), secs = round(t)),
         OUT_TAG, append = file.exists(OUT_TAG))
}

o <- fread(OUT_TAG)
cat("\n=== by season ===\n")
print(as.data.frame(o[, .(n = .N, median_km = round(mean(median_km)),
  rmse_lat = round(mean(rmse_lat), 2), cover_lat = round(mean(cover_lat), 2),
  cover_lon = round(mean(cover_lon), 2)), by = season][order(season)]), row.names = FALSE)
cat("\n=== by tag family ===\n")
print(as.data.frame(o[, .(n = .N, median_km = round(mean(median_km)),
  rmse_lat = round(mean(rmse_lat), 2)), by = family]), row.names = FALSE)
cat(sprintf("\nPOOLED OVER ALL %d: median %.0f km, rmse_lat %.2f, cover_lat %.2f\n",
            nrow(o), mean(o$median_km), mean(o$rmse_lat), mean(o$cover_lat)))
