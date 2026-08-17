# Does light_response_table() recover the response the diagnostic measured?
#
# empirical_response.csv holds the clear-sky envelope pooled across deployments
# at their TRUE Argos positions, which is the best available statement of what
# the channel actually does. The table is built the other way round: from each
# tag's own fitted envelope, without truth. If the two agree, the builder is
# recovering the real shape from data a user would actually have.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
.libPaths(c(file.path(SCRATCH, "testlib"), .libPaths()))
source(file.path(SCRATCH, "nes_common.R"))
CAL_DAYS <- 15

old <- lapply(readRDS("analysis/cache/nes/archives_v3_30min.rds"), `[[`, "main")
new <- lapply(readRDS("analysis/cache/nes/fvilches_archives_v1_30min.rds"), `[[`, "main")
new_argos <- readRDS("analysis/cache/nes/fvilches_argos_v1.rds")
new_man <- fread("scratch/nes_calibration/fvilches_manifest.csv",
                 colClasses = list(character = "topp"))
meta <- nes_meta(); ar <- nes_argos()

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
  tags[[tg]] <- list(id = tg, light = d, p0 = c(df$deploy_lon, df$deploy_lat),
                     fam = "Mk9_219")
}
for (tg in names(new)) {
  mm <- new_man[topp == tg]; if (!nrow(mm)) next
  a <- as.data.table(new_argos[[tg]]); d <- as.data.table(new[[tg]])
  d <- d[time >= min(a$time) & time <= max(a$time)]; if (nrow(d) < 1000) next
  ser <- sub("^[0-9]+_(.+)-out-Archive[.]csv$", "\1",
             list.files("fvilches/extracted/TDR raw",
                        pattern = paste0("^", tg, "_.*Archive[.]csv$")))[1]
  tags[[tg]] <- list(id = tg, light = d, p0 = c(mm$deploy_lon[1], mm$deploy_lat[1]),
                     fam = if (grepl("^18A", ser)) "F18A" else "Mk9_219")
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

e <- fread("scratch/nes_calibration/empirical_response.csv")
tabs <- list()
for (fm in unique(vapply(tags, function(g) g$fam, ""))) {
  ids <- names(tags)[vapply(tags, function(g) g$fam, "") == fm]
  gf <- Filter(Negate(is.null), geom_fits[ids])
  if (!length(gf)) { message(fm, ": no haul-out geometry"); next }
  tab <- light_response_table(gf, from = 30, to = 140, by = 1)
  tabs[[fm]] <- tab
  zz <- seq(30, 140, 1)
  ml <- median(vapply(scale_fits[ids], function(r) if (is.null(r)) NA else r$max_light, 0),
               na.rm = TRUE)
  cmp <- merge(data.table(zb = zz, table = tab[-(1:3)]),
               e[fam == fm, .(zb, measured = env95, tangent, logistic)], by = "zb")
  cat(sprintf("\n=== %s: %d fits pooled, max_light %.0f ===\n", fm, length(gf), ml))
  print(as.data.frame(cmp[zb %in% seq(40, 130, 10),
        .(zenith = zb, measured = round(measured), table = round(table),
          gap_table = round(table - measured, 1),
          tangent = round(tangent), gap_tangent = round(tangent - measured, 1))]),
        row.names = FALSE)
  k <- cmp[zb >= 40 & zb <= 130]
  cat(sprintf("rmse to the measured envelope:  table %.1f   tangent %.1f   logistic(floor 0) %.1f\n",
              sqrt(mean((k$table - k$measured)^2)),
              sqrt(mean((k$tangent - k$measured)^2)),
              sqrt(mean((k$logistic - k$measured)^2))))
  tw <- cmp[zb >= 84 & zb <= 100]
  cat(sprintf("through twilight (84-100):      table %.1f   tangent %.1f\n",
              sqrt(mean((tw$table - tw$measured)^2)),
              sqrt(mean((tw$tangent - tw$measured)^2))))
  nt <- cmp[zb >= 110]
  cat(sprintf("at night (>110):                table %.1f   tangent %.1f\n",
              sqrt(mean((nt$table - nt$measured)^2)),
              sqrt(mean((nt$tangent - nt$measured)^2))))
}
saveRDS(tabs, "scratch/nes_calibration/response_tables.rds")
cat("\ntables written to scratch/nes_calibration/response_tables.rds\n")
