# Recipe L on all ten animals, over each animal's whole deployment.
#
# On three tags and a fixed 2.5-month window it took the pooled median error
# from 321 km to 104 km with no ground truth, which is a large enough jump to
# be suspicious of. This is the honest test: every animal, the full record, the
# same endpoints the analysis uses.
#
# One trap guarded here. The Mk9 is switched on before it goes on the animal,
# so the "not diving" period at the head of the archive can include time on a
# bench under a roof. The haul-out window is therefore clipped to start at the
# deployment fix, not at the first sample.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))

arch <- lapply(readRDS("analysis/cache/nes/archives_v3_30min.rds"), `[[`, "main")
meta <- nes_meta(); ar <- nes_argos()
grid <- rast(xmin = 150, xmax = 250, ymin = 20, ymax = 66, resolution = 2, crs = "EPSG:4326")
values(grid) <- 1

fit_and_score <- function(id, r, label, d, a, p0, p1) {
  if (is.null(r)) return(data.frame(id = id, recipe = label, z50 = NA_real_,
    width = NA_real_, median_km = NA_real_, bias_lat = NA_real_, rmse_lat = NA_real_,
    post_sd_lat = NA_real_, cover_lat = NA_real_, secs = NA_real_))
  t <- system.time(invisible(capture.output(
    f <- TwilightFreeGrid(d$time, pmax(0, d$light - r$baseline), grid = grid,
      start_lon = p0[1], start_lat = p0[2], end_lon = p1[1], end_lat = p1[2],
      step_hours = 12, diffusion = 250, calibration = r$calibration,
      likelihood_params = c(1/(r$max_light*0.5), r$max_light, 0.10)))))[3]
  gp <- grid_posterior(f)
  sdp <- vapply(seq_len(nrow(gp$P)), function(i) { w <- gp$P[i, ]; s <- sum(w)
    if (!is.finite(s) || s <= 0) return(NA_real_); w <- w/s
    sqrt(sum(w*(gp$lat - sum(w*gp$lat))^2)) }, 0)
  tm <- as.POSIXct(f$fit$time, origin = "1970-01-01", tz = "UTC")
  tr <- argos_at(a, tm)
  ep <- f$fit$lat - tr$lat
  e <- gc_km(lon360(f$fit$lon), f$fit$lat, tr$lon, tr$lat)
  g <- is.finite(ep) & is.finite(sdp) & sdp > 0
  data.frame(id = id, recipe = label, z50 = round(r$z50, 1),
    width = round(4.394*r$scale, 1), median_km = round(median(e[g])),
    bias_lat = round(mean(ep[g]), 2), rmse_lat = round(sqrt(mean(ep[g]^2)), 2),
    post_sd_lat = round(median(sdp[g]), 2),
    cover_lat = round(mean(abs(ep[g]) <= 1.96*sdp[g]), 2),
    secs = round(t), row.names = NULL)
}

res <- list(); win <- list()
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
  p0 <- c(df$deploy_lon, df$deploy_lat)
  tail_fix <- a[time >= max(time) - 72*3600]
  p1 <- c(median(tail_fix$lon), median(tail_fix$lat))

  # windows: the standard 15 days, and the pre-departure haul-out clipped to
  # start at the deployment fix
  k15 <- d$time <= t0 + 15*86400
  dep <- departure_time(x)
  # Preferred: haul-out after the deployment fix, so no bench time can creep in.
  # Where the animal left within a day of being tagged that window is empty, and
  # the only haul-out on record is before the deployment DATE (which is recorded
  # to the day, so it clips real on-animal time). Fall back to it and say so.
  kpre <- as.numeric(d$time) < dep
  src <- "post-fix"
  pre_t <- d$time[kpre]; pre_l <- d$light[kpre]
  if (length(pre_t) < 200) {
    j <- as.numeric(x$time) < dep
    if (sum(j) >= 200) { pre_t <- x$time[j]; pre_l <- x$light[j]; src <- "full-archive" }
    else src <- "none"
  }
  win[[length(win)+1]] <- data.frame(id = id, ptt = m$ptt,
    days = round(as.numeric(difftime(t1, t0, units = "days")), 1),
    n_15d = sum(k15), n_haulout_postfix = sum(kpre), n_haulout_used = length(pre_t),
    haulout_days = round(length(pre_t)*0.5/24, 1), window = src, row.names = NULL)

  r15  <- fit_light_response(d$time[k15], d$light[k15], p0[1], p0[2])
  rpre <- if (length(pre_t) >= 200)
    fit_light_response(pre_t, pre_l, p0[1], p0[2]) else NULL

  res[[length(res)+1]] <- fit_and_score(id, r15, "A 15d deploy-fix (current)", d, a, p0, p1)
  res[[length(res)+1]] <- fit_and_score(id, graft_response(rpre, r15),
                                        "L haulout geom + 15d scale", d, a, p0, p1)
}
W <- do.call(rbind, win)
cat("=== calibration windows ===\n"); print(W, row.names = FALSE)
o <- as.data.table(do.call(rbind, res))
cat("\n=== per tag ===\n"); print(as.data.frame(o[order(id, recipe)]), row.names = FALSE)
cat("\n=== pooled over tags with a result for both recipes ===\n")
both <- o[!is.na(median_km), .N, by = id][N == 2, id]
print(as.data.frame(o[id %in% both, .(tags = .N,
  median_km = round(mean(median_km)), worst_km = max(median_km),
  bias_lat = round(mean(bias_lat), 2), rmse_lat = round(mean(rmse_lat), 2),
  post_sd_lat = round(mean(post_sd_lat), 2), cover_lat = round(mean(cover_lat), 2),
  mins = round(sum(secs)/60, 1)), by = recipe][order(recipe)]), row.names = FALSE)
saveRDS(o, file.path(SCRATCH, "all10.rds"))
