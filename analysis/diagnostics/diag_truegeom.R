# The ALAN hypothesis is dead: 2021032 is the DARKEST tag at night, not the
# brightest. But the zenith plot shows something else -- its twilight transition
# sits deeper and looser than a clean tag's, while the pooled calibration hands
# every animal z50 = 92.09.
#
# So: fit each tag's response against its TRUE Argos positions over the whole
# at-sea record, and compare with the pooled value the model actually used. This
# is the response the animal really has, which is not obtainable in practice --
# it is a diagnostic, not a recipe.
#
# If the tags that were fitted worst are the tags whose true response sits
# furthest from the pool, then pooling is the cost, and the fix is partial
# pooling: shrink each tag toward the panel rather than replacing it. If they are
# not, pooling is exonerated and 2021032's problem lies elsewhere again.
#
# Note what makes this different from the earlier haul-out comparison: that used
# a few days ashore, this uses 230 days at sea. If the two disagree, the
# haul-out curve is not the curve that applies during the trip -- which would be
# a real limitation of the shipped recipe and worth knowing.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))
CAL_DAYS <- 15

arch <- lapply(readRDS("analysis/cache/nes/archives_v3_30min.rds"), `[[`, "main")
meta <- nes_meta(); ar <- nes_argos()

rows <- list()
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

  dep <- departure_time(d)
  kh <- is.finite(dep) & as.numeric(d$time) < dep
  r_haul <- if (sum(kh) >= 200)
    fit_light_response(d$time[kh], d$light[kh], df$deploy_lon, df$deploy_lat) else NULL

  # the response along the TRUE track, at sea only
  at_sea <- if (is.finite(dep)) d[as.numeric(time) >= dep] else d
  tr <- argos_at(a, at_sea$time)
  ok <- is.finite(tr$lon) & is.finite(tr$lat)
  r_true <- if (sum(ok) >= 200)
    fit_light_response(at_sea$time[ok], at_sea$light[ok], tr$lon[ok], tr$lat[ok]) else NULL

  rows[[length(rows)+1]] <- data.frame(id = id,
    haul_z50  = if (is.null(r_haul)) NA_real_ else round(r_haul$z50, 2),
    haul_wid  = if (is.null(r_haul)) NA_real_ else round(r_haul$width_deg, 1),
    true_z50  = if (is.null(r_true)) NA_real_ else round(r_true$z50, 2),
    true_wid  = if (is.null(r_true)) NA_real_ else round(r_true$width_deg, 1),
    row.names = NULL)
}
G <- rbindlist(rows)
z_pool <- median(G$haul_z50, na.rm = TRUE)
w_pool <- median(G$haul_wid, na.rm = TRUE)
G[, `:=`(err_z50 = round(true_z50 - z_pool, 2), err_wid = round(true_wid - w_pool, 1))]

k <- as.data.table(read.csv("analysis/output/nes/per_knot_errors.csv"))
k <- k[engine == "hier" & batch == "light" & is.finite(err_lat)]
k[, id := as.character(id)]
sc <- k[, .(offset = round(mean(err_lat), 2), scatter = round(sd(err_lat), 2),
            cover = round(mean(abs(err_lat) <= 1.96*lat_sd, na.rm = TRUE), 2)), by = id]
M <- merge(G, sc, by = "id")

cat(sprintf("=== pooled geometry actually used: z50 %.2f, width %.1f ===\n\n",
            z_pool, w_pool))
cat("haul_* is the calibration the recipe uses; true_* is the response the animal\n")
cat("really had at sea; err_* is how far the pooled value sits from the truth.\n\n")
print(as.data.frame(M[order(-abs(err_z50))]), row.names = FALSE)

hz <- M[is.finite(err_z50)]
cat("\n=== does the pooled calibration's error explain which tags fit badly? ===\n")
for (v in c("err_z50", "err_wid")) {
  for (y in c("scatter", "cover", "offset")) {
    ct <- suppressWarnings(cor.test(abs(hz[[v]]), hz[[y]]))
    cat(sprintf("  |%s| vs %-8s r = %+.3f  (p = %.3f)\n", v, y, ct$estimate, ct$p.value))
  }
}
cat("\n=== does the haul-out curve predict the at-sea curve? ===\n")
hh <- M[is.finite(haul_z50) & is.finite(true_z50)]
cat(sprintf("  z50:   haul-out vs at-sea  r = %+.3f (n = %d)\n",
            cor(hh$haul_z50, hh$true_z50), nrow(hh)))
cat(sprintf("  width: haul-out vs at-sea  r = %+.3f\n",
            cor(hh$haul_wid, hh$true_wid)))
cat(sprintf("  mean |z50 haul-out minus at-sea| = %.2f deg\n",
            mean(abs(hh$haul_z50 - hh$true_z50))))
cat("\nIf that correlation is weak, the haul-out curve is NOT the curve that\n")
cat("applies at sea, and the shipped recipe is calibrating on the wrong regime.\n")
