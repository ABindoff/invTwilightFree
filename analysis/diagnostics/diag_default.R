# Should the DEFAULT path change too?
#
# The tangent's amplitude decides where the expected curve tops out. Tying it to
# `max_light` (the range the observations are expressed on) beat tying it to the
# envelope's `amp` by 259 km against 411 km once geometry and scale came from
# different windows. But the two are not equal even when they come from the SAME
# window -- on these tags `amp` runs 126 to 152 against a `max_light` of 149 to
# 163, a gap of up to 23 per cent -- so the plain call is likely wrong the same
# way. If it is, the conditional in fit_light_response() disappears and the
# function behaves identically whether or not `scale_light` is passed the light
# it already has, which is the least surprising thing it can do.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))

arch <- lapply(readRDS("analysis/cache/nes/archives_v3_30min.rds"), `[[`, "main")
meta <- nes_meta(); ar <- nes_argos()
grid <- rast(xmin = 150, xmax = 250, ymin = 20, ymax = 66, resolution = 2, crs = "EPSG:4326")
values(grid) <- 1

use_rng_slope <- function(r) {
  if (is.null(r)) return(NULL)
  slope <- r$max_light / (4 * r$scale)
  r$calibration <- c(slope * (r$z50 + 2*r$scale), slope)
  r$slope <- slope
  r
}

score <- function(id, r, label, d, a, p0, p1) {
  if (is.null(r)) return(NULL)
  invisible(capture.output(
    f <- TwilightFreeGrid(d$time, pmax(0, d$light - r$baseline), grid = grid,
      start_lon = p0[1], start_lat = p0[2], end_lon = p1[1], end_lat = p1[2],
      step_hours = 12, diffusion = 250, calibration = r$calibration,
      likelihood_params = c(1/(r$max_light*0.5), r$max_light, 0.10))))
  gp <- grid_posterior(f)
  sdp <- vapply(seq_len(nrow(gp$P)), function(i) { w <- gp$P[i, ]; s <- sum(w)
    if (!is.finite(s) || s <= 0) return(NA_real_); w <- w/s
    sqrt(sum(w*(gp$lat - sum(w*gp$lat))^2)) }, 0)
  tm <- as.POSIXct(f$fit$time, origin = "1970-01-01", tz = "UTC")
  tr <- argos_at(a, tm)
  ep <- f$fit$lat - tr$lat
  e <- gc_km(lon360(f$fit$lon), f$fit$lat, tr$lon, tr$lat)
  g <- is.finite(ep) & is.finite(sdp) & sdp > 0
  data.frame(id = id, recipe = label, amp = round(r$amp, 1),
    max_light = round(r$max_light, 1), median_km = round(median(e[g])),
    bias_lat = round(mean(ep[g]), 2), rmse_lat = round(sqrt(mean(ep[g]^2)), 2),
    cover_lat = round(mean(abs(ep[g]) <= 1.96*sdp[g]), 2), row.names = NULL)
}

res <- list()
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
  tail_fix <- a[time >= max(time) - 72*3600]
  p0 <- c(df$deploy_lon, df$deploy_lat); p1 <- c(median(tail_fix$lon), median(tail_fix$lat))
  k15 <- d$time <= t0 + 15*86400
  r15 <- fit_light_response(d$time[k15], d$light[k15], p0[1], p0[2])
  res[[length(res)+1]] <- score(id, r15, "current (amp slope)", d, a, p0, p1)
  res[[length(res)+1]] <- score(id, use_rng_slope(r15), "rng slope", d, a, p0, p1)
}
o <- as.data.table(do.call(rbind, res))
cat("=== default path: same window for geometry and scale ===\n")
print(as.data.frame(o[order(id, recipe)]), row.names = FALSE)
cat("\n=== pooled ===\n")
print(as.data.frame(o[, .(tags = .N, median_km = round(mean(median_km)),
  worst_km = max(median_km), bias_lat = round(mean(bias_lat), 2),
  rmse_lat = round(mean(rmse_lat), 2), cover_lat = round(mean(cover_lat), 2)),
  by = recipe]), row.names = FALSE)
