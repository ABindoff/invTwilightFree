# The full render disagrees with my harness, so one of them is measuring the
# wrong thing. My harness said naive 473 km, pooled 237. The report says naive
# 704, pooled 730 -- both far worse, and the ranking reversed.
#
# Two candidate explanations, and they are separable:
#   (a) the CLOCK. The report applies estimated offsets; my harness used raw
#       timestamps. The recomputed offsets ramp from 9 min at deployment to 41.7
#       at recovery, a 33-minute drift that is about 8 degrees of longitude by
#       the end of the trip. That is not credible for ten tags of one model, and
#       it is new: before the response changed, every estimate ran to the search
#       bound and was rejected unused.
#   (b) the CONFIGURATION. The report grids at 1 degree over the study bbox with
#       diffusion 110; my harness used 2 degrees over a smaller box at 250.
#
# So: report configuration, RAW timestamps, and the two halves of my change
# crossed against each other. If (a) is the story, everything here lands near my
# harness's numbers and the ranking is restored.
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

# the report's study domain and grid, to the letter
all_lon <- unlist(lapply(tags, function(g) g$argos$lon))
all_lat <- unlist(lapply(tags, function(g) g$argos$lat))
bbox <- c(xmin = floor(min(all_lon)) - 8, xmax = ceiling(max(all_lon)) + 4,
          ymin = floor(min(all_lat)) - 6, ymax = ceiling(max(all_lat)) + 8)
grid <- rast(xmin = floor(bbox[["xmin"]]), xmax = ceiling(bbox[["xmax"]]),
             ymin = floor(bbox[["ymin"]]), ymax = ceiling(bbox[["ymax"]]),
             resolution = CELL, crs = "EPSG:4326")
values(grid) <- 1
cat(sprintf("grid: %d x %d = %d cells at %.1f deg | lon %.0f..%.0f lat %.0f..%.0f\n",
            ncol(grid), nrow(grid), ncell(grid), CELL,
            bbox[["xmin"]], bbox[["xmax"]], bbox[["ymin"]], bbox[["ymax"]]))

# The installed package now ties the tangent to max_light. Put it back on the
# envelope's amp to recover the previous release's behaviour.
amp_slope <- function(r) {
  if (is.null(r)) return(NULL)
  slope <- r$amp / (4 * r$scale)
  r$calibration <- c(slope * (r$z50 + 2*r$scale), slope); r$slope <- slope
  r
}

score <- function(id, r, label) {
  if (is.null(r)) return(NULL)
  g <- tags[[id]]
  invisible(capture.output(
    f <- TwilightFreeGrid(g$light$time, pmax(0, g$light$light - r$baseline), grid = grid,
      start_lon = unname(g$p0[1]), start_lat = unname(g$p0[2]),
      end_lon = unname(g$p1[1]), end_lat = unname(g$p1[2]),
      step_hours = STEP_HOURS, diffusion = DIFFUSION,
      calibration = r$calibration,
      likelihood_params = c(1/(r$max_light*0.5), r$max_light, 0.10))))
  tm <- as.POSIXct(f$fit$time, origin = "1970-01-01", tz = "UTC")
  tr <- argos_at(g$argos, tm)
  e <- gc_km(lon360(f$fit$lon), f$fit$lat, tr$lon, tr$lat)
  ok <- is.finite(e)
  data.frame(id = id, recipe = label, median_km = round(median(e[ok])),
    q90_km = round(unname(quantile(e[ok], .9))),
    bias_lon = round(mean(dlon(lon360(f$fit$lon)[ok], tr$lon[ok])), 2),
    bias_lat = round(mean((f$fit$lat - tr$lat)[ok]), 2),
    rmse_lat = round(sqrt(mean(((f$fit$lat - tr$lat)[ok])^2)), 2), row.names = NULL)
}

# ---- responses --------------------------------------------------------------
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
pooled <- pool_light_responses(geom_fits, scale_fits)
cat(sprintf("pooled: z50 %.2f, width %.1f\n", attr(pooled, "pooled")[["z50"]],
            4.394 * attr(pooled, "pooled")[["scale"]]))

RECIPES <- list(
  "a 15d geom, amp slope (previous release)" = lapply(scale_fits, amp_slope),
  "b 15d geom, max_light slope"              = scale_fits,
  "c pooled geom, amp slope"                 = lapply(pooled, amp_slope),
  "d pooled geom, max_light slope (current)" = pooled)

res <- list()
for (nm in names(RECIPES)) for (id in names(tags)) {
  message(nm, " : ", id)
  res[[length(res)+1]] <- score(id, RECIPES[[nm]][[id]], nm)
}
o <- as.data.table(do.call(rbind, res))
cat("\n=== per tag ===\n"); print(as.data.frame(o[order(id, recipe)]), row.names = FALSE)
cat("\n=== pooled over tags, RAW timestamps, report configuration ===\n")
print(as.data.frame(o[, .(tags = .N, median_km = round(mean(median_km)),
  worst_km = max(median_km), q90_km = round(mean(q90_km)),
  bias_lon = round(mean(bias_lon), 2), bias_lat = round(mean(bias_lat), 2),
  rmse_lat = round(mean(rmse_lat), 2)), by = recipe][order(recipe)]), row.names = FALSE)
saveRDS(o, file.path(SCRATCH, "report_config.rds"))
