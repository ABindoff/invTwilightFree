# Why does the corrected-floor response still run +1.4 degrees NORTH?
#
# North, in a northern-summer record, means the model infers LONGER days than
# the truth. The darkness regime reads day length from whether observations
# exceed the dark envelope, so if that envelope is set too LOW for a given tag,
# its ordinary night readings look like "not dark", the night looks shorter, and
# the position runs poleward. Too high and the reverse.
#
# The table pools its shape across the family and rescales by each tag's
# max_light, so the envelope a tag receives is not the envelope that tag
# measured. Test whether that discrepancy predicts the tag's bias. If it does,
# the fix is a per-tag floor rather than a pooled one, and it is a one-line
# change to how the table is assembled.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
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
  tags[[tg]] <- list(id = tg, light = d, p0 = c(df$deploy_lon, df$deploy_lat))
}
for (tg in names(new)) {
  mm <- new_man[topp == tg]; if (!nrow(mm)) next
  a <- as.data.table(new_argos[[tg]]); d <- as.data.table(new[[tg]])
  d <- d[time >= min(a$time) & time <= max(a$time)]; if (nrow(d) < 1000) next
  tags[[tg]] <- list(id = tg, light = d, p0 = c(mm$deploy_lon[1], mm$deploy_lat[1]))
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
src <- Filter(Negate(is.null), geom_fits)

R <- rbindlist(lapply(names(tags), function(tg) {
  r <- responses[[tg]]; if (is.null(r)) return(NULL)
  tab <- light_response_table(src, from = 30, to = 140, by = 1, max_light = r$max_light)
  own <- max(r$dark_level - r$baseline, 0)      # this tag's own measured floor
  data.table(id = tg, table_floor = tail(tab, 1), own_floor = own,
             gap = tail(tab, 1) - own, max_light = r$max_light)
}))
b <- fread("scratch/nes_calibration/table_results.csv")[arm == "table_flat"]
b[, id := as.character(id)]
M <- merge(R, b[, .(id, bias_lat, median_km, cover_lat)], by = "id")
cat(sprintf("%d tags\n", nrow(M)))
cat(sprintf("table floor mean %.1f, own floor mean %.1f, gap mean %+.1f (sd %.1f)\n\n",
            mean(M$table_floor), mean(M$own_floor), mean(M$gap), sd(M$gap)))

cat("=== does the floor a tag is GIVEN minus the floor it MEASURED predict bias? ===\n")
for (v in c("bias_lat", "median_km", "cover_lat")) {
  ct <- cor.test(M$gap, M[[v]], method = "spearman", exact = FALSE)
  cat(sprintf("  gap vs %-10s rho = %+.3f (p = %.3f)\n", v, ct$estimate, ct$p.value))
}
cat("\nThe mechanism predicts a NEGATIVE rho against bias: a table floor set BELOW\n")
cat("what the tag measures makes its nights look short and drives it poleward.\n")
cat("A null here means the pooled floor is not the source of the bias, and the\n")
cat("per-tag floor is not the fix.\n")
cat(sprintf("\nfor scale: bias runs %+.2f to %+.2f, gap runs %+.1f to %+.1f units\n",
            min(M$bias_lat), max(M$bias_lat), min(M$gap), max(M$gap)))
fwrite(M, "scratch/nes_calibration/floor_bias.csv")
