# Run the rewritten light-response chunk's logic exactly as the Rmd will, on the
# real records, before committing to a re-run of every batch. Cheap: no fits.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))

arch <- lapply(readRDS("analysis/cache/nes/archives_v3_30min.rds"), `[[`, "main")
meta <- nes_meta(); ar <- nes_argos()
CAL_DAYS <- 15; DECIMATE <- 30

tags <- list()
for (id in names(arch)) {
  m <- meta[match(id, meta$id), ]
  if (is.na(m$ptt) || !nzchar(m$ptt)) next
  df <- ar$deploy[Ptt == m$ptt]; a <- ar$fixes[Ptt == m$ptt]
  if (!nrow(df) || nrow(a) < 50) next
  x <- as.data.frame(arch[[id]])
  t0 <- max(min(x$time), as.POSIXct(max(m$deploy_date, df$user_date), tz = "UTC"))
  t1 <- min(max(x$time), max(a$time))
  x <- x[x$time >= t0 & x$time <= t1, , drop = FALSE]
  if (nrow(x) < 1000) next
  tags[[id]] <- list(id = id, light = x,
                     deploy = c(lon = df$deploy_lon, lat = df$deploy_lat))
}

# ---- verbatim from the chunk ------------------------------------------------
departure_time <- function(d, dive_m = 50) {
  b <- as.integer(floor(as.numeric(d$time) / 86400))
  deep <- tapply(d$depth_max > dive_m, b, mean, na.rm = TRUE)
  day <- as.integer(names(deep))
  at_sea <- deep > 0.5
  if (!any(at_sea)) return(NA_real_)
  last_ashore <- if (any(!at_sea)) max(which(!at_sea)) else 0L
  i <- if (last_ashore >= length(at_sea)) which(at_sea)[1] else last_ashore + 1L
  day[i] * 86400
}

scale_fits <- lapply(tags, function(g) {
  k <- g$light$time <= min(g$light$time) + CAL_DAYS * 86400
  fit_light_response(g$light$time[k], g$light$light[k],
                     lon = g$deploy[["lon"]], lat = g$deploy[["lat"]])
})

geom_fits <- lapply(tags, function(g) {
  dep <- departure_time(g$light)
  k <- is.finite(dep) & as.numeric(g$light$time) < dep
  if (sum(k) < 200) return(NULL)
  fit_light_response(g$light$time[k], g$light$light[k],
                     lon = g$deploy[["lon"]], lat = g$deploy[["lat"]])
})

haulout_table <- do.call(rbind, lapply(names(tags), function(id) {
  dep <- departure_time(tags[[id]]$light)
  k <- is.finite(dep) & as.numeric(tags[[id]]$light$time) < dep
  data.frame(id = id,
             departs = if (is.finite(dep))
               format(as.POSIXct(dep, origin = "1970-01-01", tz = "UTC"), "%Y-%m-%d")
               else NA_character_,
             haulout_days = round(sum(k) * DECIMATE / 1440, 1),
             n = sum(k), used = !is.null(geom_fits[[id]]), row.names = NULL)
}))
print(haulout_table, row.names = FALSE)

responses <- pool_light_responses(geom_fits, scale_fits)
pooled_geom <- attr(responses, "pooled")
cat(sprintf("\npooled geometry from %d of %d tags: z50 = %.2f, transition %.1f degrees\n",
            pooled_geom[["n"]], length(tags), pooled_geom[["z50"]],
            4.394 * pooled_geom[["scale"]]))

stopifnot(pooled_geom[["z50"]] > 80, pooled_geom[["z50"]] < 105,
          4.394 * pooled_geom[["scale"]] > 5, 4.394 * pooled_geom[["scale"]] < 60)

resp_table <- do.call(rbind, lapply(names(responses), function(id) {
  r <- responses[[id]]
  if (is.null(r)) return(data.frame(id = id, slope = NA_real_, width_deg = NA_real_,
    saturate_at = NA_real_, zero_at = NA_real_, dark_level = NA_real_,
    baseline = NA_real_, own_geometry = NA))
  data.frame(id = id, slope = round(r$slope, 2), width_deg = round(r$width_deg, 1),
             saturate_at = round(r$saturate_at, 1), zero_at = round(r$zero_at, 1),
             dark_level = if (is.finite(r$dark_level)) round(r$dark_level) else NA_real_,
             baseline = round(r$baseline),
             own_geometry = !is.null(geom_fits[[id]]), row.names = NULL)
}))
cat("\n"); print(resp_table, row.names = FALSE)

# the panel calibration the hierarchical engine will be handed
ok <- Filter(Negate(is.null), responses)
cat(sprintf("\ncal_panel = c(%.2f, %.3f), max_light = %.1f  (over %d tags)\n",
            stats::median(vapply(ok, function(r) r$calibration[1], 0)),
            stats::median(vapply(ok, function(r) r$calibration[2], 0)),
            stats::median(vapply(ok, function(r) r$max_light, 0)), length(ok)))
stopifnot(all(vapply(responses, function(r) is.null(r) ||
  (is.finite(r$calibration[1]) && r$calibration[2] > 0 && r$max_light > 0), logical(1))))
cat("all responses well formed\n")
