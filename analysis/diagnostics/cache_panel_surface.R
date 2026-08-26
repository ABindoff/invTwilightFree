# Panel and pooled calibration built on the SURFACE-CONDITIONED decimation.
#
# This is cache_panel.R with exactly one thing changed: where the light comes
# from. Argos provenance, the deployment/recovery metadata, the inclusion
# criteria, the calibration window and the pooling are all left as they are, so a
# difference between this panel and panel_v1 is attributable to the decimation
# and nothing else.
#
# Note that light in panel_v1 came from TWO sources (archives_v3_30min.rds for
# the 2021 tags, fvilches_archives_v1_30min.rds for the rest) which were built by
# different runs of the same decimator. Here all 29 come from one file produced
# by one pass over the raw 4 s archives, which removes a provenance split that
# was never itself checked.
#
# THE FALSIFIABLE CHECK, before any position is scored: the shipped decimator
# takes a MAXIMUM, which biases light upward, and the fitted response absorbs
# that bias into z50. Replacing it with a symmetric statistic should therefore
# LOWER z50, by roughly the ~1.1 deg of zenith measured within surfacing bouts.
# panel_v1 has z50 = 91.9. If the new pooled z50 does not fall by about a degree,
# the bias account is wrong and re-running arm A would be premature.
#
# Writes analysis/cache/nes/panel_surface_v1.rds (gitignored, contains real
# positions, stays local).
SCRATCH <- "analysis/diagnostics"
source(file.path(SCRATCH, "nes_common.R"))
CAL_DAYS <- 15; CELL <- 1
SURF <- "analysis/cache/nes/archives_surface_v1_30min.rds"

stopifnot(file.exists(SURF))
surf_arch <- lapply(readRDS(SURF), `[[`, "main")
atten <- vapply(readRDS(SURF), function(z) z$atten, 0)

# Which tags belonged to which delivery: preserved so the metadata joins below
# stay exactly as they were, even though the light now has a single source.
old_ids <- names(readRDS("analysis/cache/nes/archives_v3_30min.rds"))
new_ids <- names(readRDS("analysis/cache/nes/fvilches_archives_v1_30min.rds"))
old_arch <- surf_arch[intersect(old_ids, names(surf_arch))]
new_arch <- surf_arch[intersect(new_ids, names(surf_arch))]
cat(sprintf("light source: %s (%d tags; %d old-delivery, %d new-delivery)\n",
            basename(SURF), length(surf_arch), length(old_arch), length(new_arch)))

new_argos <- readRDS("analysis/cache/nes/fvilches_argos_v1.rds")
new_man <- fread("scratch/nes_calibration/fvilches_manifest.csv",
                 colClasses = list(character = "topp"))
meta <- nes_meta(); ar <- nes_argos()

ser <- rbind(
  data.table(topp = old_ids,
             serial = sub("^[0-9]+_(.+)\\.csv$", "\\1",
                          basename(list.files("data/nes_untracked", pattern = "^[0-9]+_.*\\.csv$")))[
                            match(old_ids,
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

old_panel <- readRDS("analysis/cache/nes/panel_v1.rds")
miss <- setdiff(names(old_panel$tags), names(tags))
extra <- setdiff(names(tags), names(old_panel$tags))
cat(sprintf("composition vs panel_v1: %d shared, %d dropped%s, %d added%s\n",
            length(intersect(names(tags), names(old_panel$tags))),
            length(miss), if (length(miss)) paste0(" (", paste(miss, collapse = ","), ")") else "",
            length(extra), if (length(extra)) paste0(" (", paste(extra, collapse = ","), ")") else ""))

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

responses <- list(); pooled_z50 <- c()
for (fam in unique(vapply(tags, function(g) g$family, ""))) {
  ids <- names(tags)[vapply(tags, function(g) g$family, "") == fam]
  gf <- geom_fits[ids]; sf <- scale_fits[ids]
  if (!sum(!vapply(gf, is.null, logical(1)))) next
  pooled <- pool_light_responses(gf, sf)
  p <- attr(pooled, "pooled")
  cat(sprintf("family %-8s: %d deployments, pooled from %d -> z50 %.2f, width %.1f\n",
              fam, length(ids), p[["n"]], p[["z50"]], 4.394 * p[["scale"]]))
  pooled_z50[fam] <- p[["z50"]]
  for (i in ids) responses[[i]] <- pooled[[i]]
}

if (!length(responses))
  stop("no responses were pooled: every geom_fit was NULL. Check that depth_max in\n",
       "  the decimated light still exceeds 50 m, since departure_time() classifies\n",
       "  at-sea days by that and a decimator that only keeps surface samples will\n",
       "  cap it below the threshold.", call. = FALSE)

cat("\n=== THE CHECK: pooled z50, shipped max-decimation vs surface decimation ===\n")
for (fam in names(pooled_z50)) {
  oldz <- unique(vapply(old_panel$responses[
    names(old_panel$tags)[vapply(old_panel$tags, function(g) g$family, "") == fam]],
    function(r) r$z50, 0))
  cat(sprintf("  %-8s  old z50 %.2f -> new z50 %.2f   change %+.2f deg\n",
              fam, oldz[1], pooled_z50[fam], pooled_z50[fam] - oldz[1]))
}
cat("  predicted: a fall of roughly 1 deg if the max-selection bias account is right\n")

sh <- intersect(names(tags), names(old_panel$tags))
bl <- vapply(sh, function(i) responses[[i]]$baseline, 0)
bo <- vapply(sh, function(i) old_panel$responses[[i]]$baseline, 0)
ml <- vapply(sh, function(i) responses[[i]]$max_light, 0)
mo <- vapply(sh, function(i) old_panel$responses[[i]]$max_light, 0)
cat(sprintf("\n  baseline  mean %.1f -> %.1f (%+.1f light units)\n", mean(bo), mean(bl), mean(bl - bo)))
cat(sprintf("  max_light mean %.1f -> %.1f (%+.1f light units)\n", mean(mo), mean(ml), mean(ml - mo)))
nn <- vapply(sh, function(i) nrow(tags[[i]]$light), 0)
no <- vapply(sh, function(i) nrow(old_panel$tags[[i]]$light), 0)
cat(sprintf("  observations per tag mean %.0f -> %.0f (%+.1f%%)\n",
            mean(no), mean(nn), 100 * (mean(nn) / mean(no) - 1)))

grid <- rast(xmin = 150, xmax = 250, ymin = 20, ymax = 70,
             resolution = CELL, crs = "EPSG:4326")
values(grid) <- 1
saveRDS(list(tags = tags, responses = responses, geom_fits = geom_fits,
             scale_fits = scale_fits, grid_extent = c(150, 250, 20, 70),
             atten = atten, decimation = "surface_median_v1"),
        "analysis/cache/nes/panel_surface_v1.rds")
cat(sprintf("\ncached: %d tags, %d with responses -> analysis/cache/nes/panel_surface_v1.rds\n",
            length(tags), length(responses)))
