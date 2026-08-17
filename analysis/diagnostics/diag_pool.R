# Completing the recipe.
#
# Haul-out geometry works where there is a haul-out: on the seven animals that
# have one, pooled median error falls from 508 km to 259 km over the full
# deployments, and the worst tag from 1255 km to 150 km. Three animals left
# within a day of being tagged and have almost no pre-departure record at all.
#
# But the response geometry is a property of the LIGHT CHANNEL, and all ten
# animals carry the same tag model. So the seven good haul-outs are seven
# measurements of nearly the same thing, and their median should serve the three
# that have none. Two questions, and the second one matters more:
#
#   1. Does pooled geometry rescue the three tags with no haul-out?
#   2. Is pooled geometry just as good for the seven that DO have one? If it is,
#      the recipe collapses to something far simpler and far more portable.
#
# The intensity scale stays per tag and per record either way -- that is the one
# thing the earlier tests showed does not transfer.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))

arch <- lapply(readRDS("analysis/cache/nes/archives_v3_30min.rds"), `[[`, "main")
meta <- nes_meta(); ar <- nes_argos()
grid <- rast(xmin = 150, xmax = 250, ymin = 20, ymax = 66, resolution = 2, crs = "EPSG:4326")
values(grid) <- 1

NA_ROW <- function(id, label) data.frame(id = id, recipe = label, z50 = NA_real_,
  width = NA_real_, median_km = NA_real_, bias_lat = NA_real_, rmse_lat = NA_real_,
  post_sd_lat = NA_real_, cover_lat = NA_real_, secs = NA_real_)

fit_and_score <- function(id, r, label, d, a, p0, p1) {
  if (is.null(r)) return(NA_ROW(id, label))
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

# ---- pass 1: assemble every tag's windows and its haul-out response ----------
P <- list()
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
  k15 <- d$time <= t0 + 15*86400
  kpre <- as.numeric(d$time) < departure_time(x)
  P[[id]] <- list(id = id, d = d, a = a,
    p0 = c(df$deploy_lon, df$deploy_lat),
    p1 = c(median(tail_fix$lon), median(tail_fix$lat)),
    r15 = fit_light_response(d$time[k15], d$light[k15], df$deploy_lon, df$deploy_lat),
    rpre = if (sum(kpre) >= 200)
      fit_light_response(d$time[kpre], d$light[kpre], df$deploy_lon, df$deploy_lat)
      else NULL)
}

have <- names(P)[vapply(P, function(p) !is.null(p$rpre), logical(1))]
z50_p <- median(vapply(P[have], function(p) p$rpre$z50, 0))
sc_p  <- median(vapply(P[have], function(p) p$rpre$scale, 0))
cat(sprintf("pooled geometry from %d haul-outs: z50 = %.2f, scale = %.2f (width %.1f deg)\n",
            length(have), z50_p, sc_p, 4.394*sc_p))
cat("per-tag haul-out geometry:\n")
print(data.frame(id = have,
  z50 = round(vapply(P[have], function(p) p$rpre$z50, 0), 2),
  width = round(vapply(P[have], function(p) 4.394*p$rpre$scale, 0), 1),
  row.names = NULL))

pooled_response <- function(rs) {
  if (is.null(rs)) return(NULL)
  r <- rs; r$z50 <- z50_p; r$scale <- sc_p
  slope <- rs$max_light / (4 * sc_p)
  r$calibration <- c(slope * (z50_p + 2*sc_p), slope)
  r$baseline <- rs$baseline; r$max_light <- rs$max_light
  r
}

# ---- pass 2: score ----------------------------------------------------------
res <- list()
for (id in names(P)) {
  p <- P[[id]]
  res[[length(res)+1]] <- fit_and_score(id, p$r15, "A 15d deploy-fix (current)",
                                        p$d, p$a, p$p0, p$p1)
  res[[length(res)+1]] <- fit_and_score(id, graft_response(p$rpre, p$r15),
                                        "L haulout geom", p$d, p$a, p$p0, p$p1)
  res[[length(res)+1]] <- fit_and_score(id, pooled_response(p$r15),
                                        "P pooled geom", p$d, p$a, p$p0, p$p1)
}
o <- as.data.table(do.call(rbind, res))
o[, has_haulout := id %in% have]
cat("\n=== per tag ===\n")
print(as.data.frame(o[order(id, recipe)]), row.names = FALSE)
cat("\n=== pooled, all ten tags ===\n")
print(as.data.frame(o[!is.na(median_km), .(tags = .N,
  median_km = round(mean(median_km)), worst_km = max(median_km),
  bias_lat = round(mean(bias_lat), 2), rmse_lat = round(mean(rmse_lat), 2),
  cover_lat = round(mean(cover_lat), 2)), by = recipe][order(recipe)]), row.names = FALSE)
cat("\n=== split by whether the tag has a usable haul-out ===\n")
print(as.data.frame(o[!is.na(median_km), .(tags = .N,
  median_km = round(mean(median_km)), worst_km = max(median_km),
  rmse_lat = round(mean(rmse_lat), 2), cover_lat = round(mean(cover_lat), 2)),
  by = .(has_haulout, recipe)][order(-has_haulout, recipe)]), row.names = FALSE)
cat("\n=== fit times (watching for the 5913 s outlier in the last run) ===\n")
print(as.data.frame(o[!is.na(secs), .(id, recipe, secs)][order(-secs)][1:6]), row.names = FALSE)
saveRDS(o, file.path(SCRATCH, "pool.rds"))
