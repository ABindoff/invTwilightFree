# Assemble the panel and the pooled calibration ONCE, and cache them.
#
# The panel-plus-calibration block has now been copied verbatim three times
# (fit_all29.R, fit_rescore29.R, fit_family_control.R). Copying it again to answer
# the next question would be asking for drift -- and drift here is expensive,
# because a harness that silently diverges from the pipeline produces numbers that
# look fine and mean nothing. This caches the result so downstream diagnostics load
# it instead of rebuilding it.
#
# It writes `analysis/cache/nes/panel_v1.rds`, which is gitignored (analysis/cache/)
# and contains real positions, so it stays local.
#
# The block below is still a verbatim copy of fit_all29.R -- this script does not
# improve it, it just stops the copying. Verified by fit_family_control.R
# reproducing rescore29's tangent arm to 0.999/1.002 on two tags, with the residual
# fully attributed to the argos_at/truth_at scorer difference.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))
CAL_DAYS <- 15; CELL <- 1

# ================= copied verbatim from fit_all29.R: panel ====================
old_arch <- lapply(readRDS("analysis/cache/nes/archives_v3_30min.rds"), `[[`, "main")
new_arch <- lapply(readRDS("analysis/cache/nes/fvilches_archives_v1_30min.rds"), `[[`, "main")
new_argos <- readRDS("analysis/cache/nes/fvilches_argos_v1.rds")
new_man <- fread("scratch/nes_calibration/fvilches_manifest.csv",
                 colClasses = list(character = "topp"))
meta <- nes_meta(); ar <- nes_argos()

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
for (tg in names(old_arch)) {
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
for (tg in names(new_arch)) {
  mm <- new_man[topp == tg]
  if (!nrow(mm)) next
  a <- as.data.table(new_argos[[tg]])
  d <- as.data.table(new_arch[[tg]])
  d <- d[time >= min(a$time) & time <= max(a$time)]
  if (nrow(d) < 1000 || nrow(a) < 50) next
  tags[[tg]] <- list(id = tg, season = as.character(mm$season[1]), light = d, argos = a,
                     p0 = c(mm$deploy_lon[1], mm$deploy_lat[1]),
                     p1 = c(mm$recover_lon[1], mm$recover_lat[1]))
}
for (tg in names(tags)) tags[[tg]]$family <- family(ser[topp == tg]$serial[1])
cat(sprintf("panel: %d deployments\n", length(tags)))

# ================= copied verbatim from fit_all29.R: calibration ==============
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
cat(sprintf("haul-out geometry available for %d of %d deployments\n",
            sum(!vapply(geom_fits, is.null, logical(1))), length(tags)))

responses <- list()
for (fam in unique(vapply(tags, function(g) g$family, ""))) {
  ids <- names(tags)[vapply(tags, function(g) g$family, "") == fam]
  gf <- geom_fits[ids]; sf <- scale_fits[ids]
  if (!sum(!vapply(gf, is.null, logical(1)))) next
  pooled <- pool_light_responses(gf, sf)
  p <- attr(pooled, "pooled")
  cat(sprintf("family %-8s: %d deployments, pooled from %d -> z50 %.2f, width %.1f\n",
              fam, length(ids), p[["n"]], p[["z50"]], 4.394 * p[["scale"]]))
  for (i in ids) responses[[i]] <- pooled[[i]]
}

grid <- rast(xmin = 150, xmax = 250, ymin = 20, ymax = 70,
             resolution = CELL, crs = "EPSG:4326")
values(grid) <- 1

saveRDS(list(tags = tags, responses = responses, geom_fits = geom_fits,
             scale_fits = scale_fits, grid_extent = c(150, 250, 20, 70)),
        "analysis/cache/nes/panel_v1.rds")
cat(sprintf("\ncached: %d tags, %d with responses -> analysis/cache/nes/panel_v1.rds\n",
            length(tags), length(responses)))
