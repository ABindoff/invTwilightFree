# "Per-tag latitude offset", made concrete.
#
# Split each animal's latitude error into a CONSTANT part (the offset: it sits
# north or south of truth all deployment) and the SCATTER around it. If the
# offset dominates, one number per tag would fix most of the error -- and the
# question becomes whether that number can be got from the tag itself.
#
# Then the hypothesis that matters, because it needs no positional tag at all:
# section 5.2 POOLS the response geometry across animals. If the true geometry
# differs slightly between tags, pooling imposes the panel's curve on each one,
# and a tag whose own curve sits furthest from the pool should be pushed
# furthest off in latitude. That is testable: correlate each tag's own haul-out
# z50 (minus the pooled value) against its latitude offset.
#
# If it holds, the fix is partial pooling -- shrink each tag toward the panel
# instead of replacing it -- which uses only the tag's own light.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))
CAL_DAYS <- 15

k <- as.data.table(read.csv("analysis/output/nes/per_knot_errors.csv"))
k <- k[engine == "hier" & batch == "light" & is.finite(err_lat)]
k[, id := as.character(id)]

cat("=== how much of each tag's latitude error is a constant offset? ===\n")
s <- k[, .(knots = .N,
           offset = round(mean(err_lat), 2),
           scatter = round(sd(err_lat), 2),
           rmse = round(sqrt(mean(err_lat^2)), 2)), by = id]
s[, rmse_if_offset_removed := round(scatter, 2)]
s[, pct_of_error_that_is_offset := round(100 * offset^2 / (offset^2 + scatter^2))]
print(as.data.frame(s[order(-abs(offset))]), row.names = FALSE)
cat(sprintf("\npooled RMSE now %.2f deg; removing a per-tag constant would give %.2f\n",
            sqrt(mean(k$err_lat^2)),
            sqrt(mean(k[, .(v = mean((err_lat - mean(err_lat))^2)), by = id]$v))))

# ---- does the offset come from POOLING the response geometry? ---------------
arch <- lapply(readRDS("analysis/cache/nes/archives_v3_30min.rds"), `[[`, "main")
meta <- nes_meta(); ar <- nes_argos()
geo <- list()
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
  kk <- is.finite(dep) & as.numeric(d$time) < dep
  r <- if (sum(kk) >= 200)
    fit_light_response(d$time[kk], d$light[kk], df$deploy_lon, df$deploy_lat) else NULL
  geo[[id]] <- data.frame(id = id,
    own_z50 = if (is.null(r)) NA_real_ else round(r$z50, 2),
    own_width = if (is.null(r)) NA_real_ else round(r$width_deg, 1),
    has_own = !is.null(r), row.names = NULL)
}
G <- rbindlist(geo)
z50_p <- median(G$own_z50, na.rm = TRUE); w_p <- median(G$own_width, na.rm = TRUE)
G[, `:=`(d_z50 = round(own_z50 - z50_p, 2), d_width = round(own_width - w_p, 1))]
M <- merge(s[, .(id, offset, scatter)], G, by = "id")
cat(sprintf("\n=== pooled geometry: z50 %.2f, width %.1f ===\n", z50_p, w_p))
print(as.data.frame(M[order(d_z50)]), row.names = FALSE)

hz <- M[!is.na(d_z50)]
if (nrow(hz) >= 4) {
  cat(sprintf("\ncorrelation of latitude offset with (own z50 - pooled z50): %+.3f (n = %d)\n",
              cor(hz$d_z50, hz$offset), nrow(hz)))
  cat(sprintf("correlation with (own width - pooled width):               %+.3f\n",
              cor(hz$d_width, hz$offset)))
  f <- lm(offset ~ d_z50, data = hz)
  cat(sprintf("slope %.2f deg latitude per deg of z50 (p = %.3g, R2 = %.2f)\n",
              coef(f)[2], summary(f)$coefficients[2, 4], summary(f)$r.squared))
}
cat("\n=== do the three tags with NO haul-out fare worse? ===\n")
print(as.data.frame(M[, .(tags = .N, mean_abs_offset = round(mean(abs(offset)), 2),
                          mean_scatter = round(mean(scatter), 2)), by = has_own]),
      row.names = FALSE)
