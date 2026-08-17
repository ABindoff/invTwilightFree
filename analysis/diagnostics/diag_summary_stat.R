# Are we summarising the posterior badly?
#
# The grid HMM reports `best_i`, the argmax of the smoothed marginal at each
# knot -- the posterior MODE. That choice has never been tested. The posterior is
# skewed even where it is not folded (the latitude likelihood flattens near the
# equinox rather than going bimodal, because the hemisphere mirror lies outside
# the search domain), so mode, mean and median differ and one may be
# systematically better. If it is, that is an improvement with no modelling
# change at all: the same posterior, summarised differently.
#
# Four summaries per knot per axis:
#   mode    the current behaviour, argmax of the marginal
#   mean    the probability-weighted centroid
#   median  the 50th percentile of the marginal
#   trimmed the centroid over the cells holding the central 80% of the mass,
#           which is the mean's robustness to a long tail without the mode's
#           sensitivity to a single cell
#
# Scored per axis separately, because the two axes are in different regimes:
# longitude is likelihood-dominated (pi ~ 0.30) and latitude prior-dominated
# (pi ~ 0.87), so the right summary need not be the same for both.
#
# RESUMABLE: appends per tag.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))
OUT <- "scratch/nes_calibration/summary_stat_results.csv"

STEP_HOURS <- 12; DIFFUSION <- 110; CELL <- 1; CAL_DAYS <- 15

arch <- lapply(readRDS("analysis/cache/nes/archives_v3_30min.rds"), `[[`, "main")
meta <- nes_meta(); ar <- nes_argos()
tags <- list()
for (tg in names(arch)) {
  m <- meta[match(tg, meta$id), ]
  if (is.na(m$ptt) || !nzchar(m$ptt)) next
  df <- ar$deploy[Ptt == m$ptt]; a <- ar$fixes[Ptt == m$ptt]
  if (!nrow(df) || nrow(a) < 50) next
  x <- as.data.table(arch[[tg]])
  t0 <- max(min(x$time), as.POSIXct(max(m$deploy_date, df$user_date), tz = "UTC"))
  t1 <- min(max(x$time), max(a$time))
  d <- x[time >= t0 & time <= t1]
  if (nrow(d) < 1000) next
  tf <- a[time >= max(time) - 72*3600]
  tags[[tg]] <- list(id = tg, light = d, argos = a,
                     p0 = c(df$deploy_lon, df$deploy_lat),
                     p1 = c(median(tf$lon), median(tf$lat)))
}
grid <- rast(xmin = 150, xmax = 250, ymin = 20, ymax = 70,
             resolution = CELL, crs = "EPSG:4326")
values(grid) <- 1

# weighted quantile over a discrete marginal
wquant <- function(x, w, p = 0.5) {
  o <- order(x); x <- x[o]; w <- w[o]
  cw <- cumsum(w) / sum(w)
  x[which(cw >= p)[1]]
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

done <- if (file.exists(OUT)) unique(data.table::fread(OUT)$id) else character(0)
cat(sprintf("%d tags done, %d to go\n", length(done), length(tags) - length(done)))

for (tg in names(tags)) {
  if (tg %in% done) next
  message("summaries: ", tg)
  g <- tags[[tg]]; r <- responses[[tg]]
  invisible(capture.output(
    f <- TwilightFreeGrid(g$light$time, pmax(0, g$light$light - r$baseline), grid = grid,
      start_lon = unname(g$p0[1]), start_lat = unname(g$p0[2]),
      end_lon = unname(g$p1[1]), end_lat = unname(g$p1[2]),
      step_hours = STEP_HOURS, diffusion = DIFFUSION,
      calibration = r$calibration,
      likelihood_params = c(1/(r$max_light*0.5), r$max_light, 0.10))))
  gp <- grid_posterior(f)
  K <- nrow(gp$P)
  est <- data.table(mode_lon = f$fit$lon, mode_lat = f$fit$lat,
                    mean_lon = NA_real_, mean_lat = NA_real_,
                    med_lon = NA_real_, med_lat = NA_real_,
                    trim_lon = NA_real_, trim_lat = NA_real_)
  for (i in seq_len(K)) {
    w <- gp$P[i, ]; s <- sum(w)
    if (!is.finite(s) || s <= 0) next
    w <- w / s
    est$mean_lon[i] <- sum(w * gp$lon); est$mean_lat[i] <- sum(w * gp$lat)
    est$med_lon[i] <- wquant(gp$lon, w); est$med_lat[i] <- wquant(gp$lat, w)
    keep <- w >= quantile(w[w > 0], 0.2)          # central mass, drop the tail
    est$trim_lon[i] <- sum(w[keep] * gp$lon[keep]) / sum(w[keep])
    est$trim_lat[i] <- sum(w[keep] * gp$lat[keep]) / sum(w[keep])
  }
  tm <- as.POSIXct(f$fit$time, origin = "1970-01-01", tz = "UTC")
  tr <- argos_at(g$argos, tm)
  rows <- rbindlist(lapply(c("mode", "mean", "med", "trim"), function(s) {
    lo <- lon360(est[[paste0(s, "_lon")]]); la <- est[[paste0(s, "_lat")]]
    e <- gc_km(lo, la, tr$lon, tr$lat)
    ep <- la - tr$lat; el <- dlon(lo, tr$lon)
    ok <- is.finite(e)
    data.table(id = tg, summary = s, n = sum(ok),
               median_km = round(median(e[ok])),
               q90_km = round(unname(quantile(e[ok], .9))),
               bias_lat = round(mean(ep[ok]), 3),
               rmse_lat = round(sqrt(mean(ep[ok]^2)), 3),
               rmse_lon = round(sqrt(mean(el[ok]^2)), 3))
  }))
  data.table::fwrite(rows, OUT, append = file.exists(OUT))
}

o <- data.table::fread(OUT)
cat("\n=== pooled over tags ===\n")
print(as.data.frame(o[, .(tags = .N, median_km = round(mean(median_km)),
  q90_km = round(mean(q90_km)), bias_lat = round(mean(bias_lat), 3),
  rmse_lat = round(mean(rmse_lat), 3), rmse_lon = round(mean(rmse_lon), 3)),
  by = summary][order(median_km)]), row.names = FALSE)

cat("\n=== paired against the mode, per tag ===\n")
base <- o[summary == "mode", .(id, km0 = median_km, la0 = rmse_lat, lo0 = rmse_lon)]
for (s in c("mean", "med", "trim")) {
  m <- merge(base, o[summary == s, .(id, km = median_km, la = rmse_lat, lo = rmse_lon)],
             by = "id")
  if (!nrow(m)) next
  cat(sprintf("  %-5s  d_median %+5.1f km  better on %d/%d   d_rmse_lat %+.3f  d_rmse_lon %+.3f  p=%.3f\n",
              s, mean(m$km - m$km0), sum(m$km < m$km0), nrow(m),
              mean(m$la - m$la0), mean(m$lo - m$lo0),
              tryCatch(wilcox.test(m$km0, m$km, paired = TRUE)$p.value, error = function(e) NA)))
}
cat("\nNegative deltas favour the alternative over the mode. The two axes are in\n")
cat("different regimes (longitude likelihood-dominated, latitude prior-dominated)\n")
cat("so the best summary need not be the same for both.\n")
