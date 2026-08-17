# Settling an API wart before it ships.
#
# The tangent handed to the engine has slope amp/(4*scale). Two candidates for
# `amp` once geometry and scale come from different windows:
#
#   amp  -- the clear-sky amplitude measured off the geometry window's envelope.
#           A property of the light channel, so it should transfer between
#           windows unchanged. This is what the function uses today.
#   rng  -- the scale window's q95 - q05, which is also `max_light`. Ties the
#           top of the expected curve to the top of the range the observations
#           are expressed on. This is what recipe L happened to use.
#
# If they score the same, use `amp` for both: the function then behaves
# identically whether or not `scale_light` is passed the same light it already
# has, which is the least surprising thing it can do. If `rng` really is better,
# the conditional stays and the docs have to earn it.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))

arch <- lapply(readRDS("analysis/cache/nes/archives_v3_30min.rds"), `[[`, "main")
meta <- nes_meta(); ar <- nes_argos()
grid <- rast(xmin = 150, xmax = 250, ymin = 20, ymax = 66, resolution = 2, crs = "EPSG:4326")
values(grid) <- 1

# geometry from rg, scale from rs, with the tangent's amplitude chosen explicitly
graft2 <- function(rg, rs, amp_from = c("rng", "amp")) {
  if (is.null(rg) || is.null(rs)) return(NULL)
  amp_from <- match.arg(amp_from)
  amp_lin <- if (amp_from == "rng") rs$max_light else rg$amp
  slope <- amp_lin / (4 * rg$scale)
  rg$calibration <- c(slope * (rg$z50 + 2*rg$scale), slope)
  rg$baseline <- rs$baseline; rg$max_light <- rs$max_light
  rg$amp_lin <- amp_lin
  rg
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
  data.frame(id = id, recipe = label, amp_lin = round(r$amp_lin, 1),
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
  kpre <- as.numeric(d$time) < departure_time(x)
  if (sum(kpre) < 200) next
  tail_fix <- a[time >= max(time) - 72*3600]
  p0 <- c(df$deploy_lon, df$deploy_lat); p1 <- c(median(tail_fix$lon), median(tail_fix$lat))
  k15 <- d$time <= t0 + 15*86400
  r15  <- fit_light_response(d$time[k15], d$light[k15], p0[1], p0[2])
  rpre <- fit_light_response(d$time[kpre], d$light[kpre], p0[1], p0[2])
  res[[length(res)+1]] <- score(id, graft2(rpre, r15, "rng"), "rng slope", d, a, p0, p1)
  res[[length(res)+1]] <- score(id, graft2(rpre, r15, "amp"), "amp slope", d, a, p0, p1)
}
o <- as.data.table(do.call(rbind, res))
cat("=== tangent amplitude: scale window's range against the envelope's ===\n")
print(as.data.frame(o[order(id, recipe)]), row.names = FALSE)
cat("\n=== pooled ===\n")
print(as.data.frame(o[, .(tags = .N, median_km = round(mean(median_km)),
  worst_km = max(median_km), bias_lat = round(mean(bias_lat), 2),
  rmse_lat = round(mean(rmse_lat), 2), cover_lat = round(mean(cover_lat), 2)),
  by = recipe]), row.names = FALSE)
