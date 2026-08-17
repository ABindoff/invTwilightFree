# Now that `calibration_logistic` is on the engine's scale, is the four-parameter
# logistic actually better than the tangent?
#
# The historical answer was no, 833 km against 354, and that is why the tangent
# is what the engine gets. But that comparison was run with the floor on the RAW
# scale while the engine was handed baselined light, so the logistic sat 11% of
# the response range too high at every zenith. It was never a fair test.
#
# The measured response says the tangent's single largest error is that it forces
# expected light to zero at night while the tag actually reads 39-44 (25% of the
# range). A correctly-scaled logistic can represent that and the tangent cannot,
# so if the floor matters at all, this is where it shows.
#
# Same 29 deployments, same everything else. Two arms per tag.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))
OUT  <- "scratch/nes_calibration/logistic_v2_results.csv"
OUTK <- "scratch/nes_calibration/logistic_v2_knots.csv"

STEP_HOURS <- 12; DIFFUSION <- 110; CELL <- 1; CAL_DAYS <- 15
ARMS <- c("tangent", "logistic")

old <- lapply(readRDS("analysis/cache/nes/archives_v3_30min.rds"), `[[`, "main")
new <- lapply(readRDS("analysis/cache/nes/fvilches_archives_v1_30min.rds"), `[[`, "main")
new_argos <- readRDS("analysis/cache/nes/fvilches_argos_v1.rds")
new_man <- fread("scratch/nes_calibration/fvilches_manifest.csv",
                 colClasses = list(character = "topp"))
meta <- nes_meta(); ar <- nes_argos()

# Tag serial per deployment. Pooling must be WITHIN FAMILY, exactly as
# fit_all29.R does it, or the tangent arm here is not the same tangent that
# produced the 244 km headline and the comparison loses its anchor.
ser <- rbind(
  data.table(topp = names(old),
             serial = sub("^[0-9]+_(.+)\\.csv$", "\\1",
                          basename(list.files("data/nes_untracked", pattern = "^[0-9]+_.*\\.csv$")))[
                            match(names(old),
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
for (tg in names(old)) {
  m <- meta[match(tg, meta$id), ]
  if (is.na(m$ptt) || !nzchar(m$ptt)) next
  df <- ar$deploy[Ptt == m$ptt]; a <- ar$fixes[Ptt == m$ptt]
  if (!nrow(df) || nrow(a) < 50) next
  x <- as.data.table(old[[tg]])
  t0 <- max(min(x$time), as.POSIXct(max(m$deploy_date, df$user_date), tz = "UTC"))
  t1 <- min(max(x$time), max(a$time))
  d <- x[time >= t0 & time <= t1]; if (nrow(d) < 1000) next
  tf <- a[time >= max(time) - 72*3600]
  tags[[tg]] <- list(id = tg, light = d, argos = a,
                     p0 = c(df$deploy_lon, df$deploy_lat),
                     p1 = c(median(tf$lon), median(tf$lat)))
}
for (tg in names(new)) {
  mm <- new_man[topp == tg]; if (!nrow(mm)) next
  a <- as.data.table(new_argos[[tg]]); d <- as.data.table(new[[tg]])
  d <- d[time >= min(a$time) & time <= max(a$time)]; if (nrow(d) < 1000) next
  tags[[tg]] <- list(id = tg, light = d, argos = a,
                     p0 = c(mm$deploy_lon[1], mm$deploy_lat[1]),
                     p1 = c(mm$recover_lon[1], mm$recover_lat[1]))
}
for (tg in names(tags)) tags[[tg]]$family <- family(ser[topp == tg]$serial[1])

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

cat("\n=== the correction, per tag ===\n")
chk <- rbindlist(lapply(names(responses), function(tg) {
  r <- responses[[tg]]; if (is.null(r)) return(NULL)
  data.table(id = tg, dark = round(r$dark_level, 1), base = round(r$baseline, 1),
             floor_used = round(r$calibration_logistic[1], 1),
             max_light = round(r$max_light, 1),
             pct = round(100*r$calibration_logistic[1]/r$max_light, 1))
}))
print(as.data.frame(chk), row.names = FALSE)
cat(sprintf("\nmean floor now supplied: %.1f units = %.1f%% of range (was the raw %.1f)\n\n",
            mean(chk$floor_used), mean(chk$pct), mean(chk$dark)))

grid <- rast(xmin = 150, xmax = 250, ymin = 20, ymax = 70,
             resolution = CELL, crs = "EPSG:4326")
values(grid) <- 1
truth_at <- function(a, tm) {
  setorder(a, time)
  list(lon = approx(as.numeric(a$time), a$lon, as.numeric(tm), rule = 2)$y,
       lat = approx(as.numeric(a$time), a$lat, as.numeric(tm), rule = 2)$y)
}
done <- if (file.exists(OUT)) paste(fread(OUT)$id, fread(OUT)$arm) else character(0)
cat(sprintf("%d of %d fits already done\n\n", length(done), length(tags)*length(ARMS)))

for (arm in ARMS) for (tg in names(tags)) {
  if (paste(tg, arm) %in% done) next
  g <- tags[[tg]]; r <- responses[[tg]]; if (is.null(r)) next
  message(arm, " : ", tg)
  cal <- if (arm == "tangent") r$calibration else r$calibration_logistic
  t <- system.time(invisible(capture.output(
    f <- TwilightFreeGrid(g$light$time, pmax(0, g$light$light - r$baseline),
      grid = grid,
      start_lon = unname(g$p0[1]), start_lat = unname(g$p0[2]),
      end_lon = unname(g$p1[1]), end_lat = unname(g$p1[2]),
      step_hours = STEP_HOURS, diffusion = DIFFUSION, calibration = cal,
      likelihood_params = c(1/(r$max_light*0.5), r$max_light, 0.10)))))[3]
  gp <- grid_posterior(f)
  sdp <- sdl <- numeric(nrow(gp$P))
  for (i in seq_len(nrow(gp$P))) {
    w <- gp$P[i, ]; s <- sum(w)
    if (!is.finite(s) || s <= 0) { sdp[i] <- sdl[i] <- NA; next }
    w <- w/s
    sdl[i] <- sqrt(sum(w*(gp$lon - sum(w*gp$lon))^2))
    sdp[i] <- sqrt(sum(w*(gp$lat - sum(w*gp$lat))^2))
  }
  tm <- as.POSIXct(f$fit$time, origin = "1970-01-01", tz = "UTC")
  tt <- truth_at(as.data.table(g$argos), tm)
  ep <- f$fit$lat - tt$lat; el <- dlon(lon360(f$fit$lon), tt$lon)
  e <- gc_km(lon360(f$fit$lon), f$fit$lat, tt$lon, tt$lat)
  ok <- is.finite(ep) & is.finite(sdp) & sdp > 0
  fwrite(data.table(id = tg, arm = arm, time = tm, err_lat = ep, err_lon = el,
                    err_km = e, lat_sd = sdp),
         OUTK, append = file.exists(OUTK))
  fwrite(data.table(id = tg, arm = arm, n = sum(ok),
                    median_km = round(median(e[ok])),
                    bias_lat = round(mean(ep[ok]), 3),
                    rmse_lat = round(sqrt(mean(ep[ok]^2)), 3),
                    rmse_lon = round(sqrt(mean(el[ok]^2)), 3),
                    cover_lat = round(mean(abs(ep[ok]) <= 1.96*sdp[ok]), 3),
                    cover_lon = round(mean(abs(el[ok]) <= 1.96*sdl[ok]), 3),
                    secs = round(t)),
         OUT, append = file.exists(OUT))
}
o <- fread(OUT)
cat("\n=== tangent vs a correctly-scaled logistic ===\n")
print(as.data.frame(o[, .(tags = .N, median_km = round(mean(median_km)),
  bias_lat = round(mean(bias_lat), 3), rmse_lat = round(mean(rmse_lat), 3),
  rmse_lon = round(mean(rmse_lon), 3), cover_lat = round(mean(cover_lat), 3)),
  by = arm]), row.names = FALSE)
