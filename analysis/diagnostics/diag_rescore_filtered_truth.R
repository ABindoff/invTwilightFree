# Does speed-filtering the Argos reference change the campaign's headline numbers?
#
# The FITS do not change -- only the reference they are scored against does. So this
# needs no refitting: take the stored per-knot fitted positions, recompute the true
# position from the FILTERED Argos at the same knot times, and compare.
#
# What is at stake: 10 of 29 deployments (the whole 2021 delivery) were scored against
# raw Argos with 15-25% of steps implying speeds above 5 m/s and a maximum of
# 159,722 m/s. The other 19 arrived pre-filtered at 3 m/s. If the headline moves, then
# every accuracy and bias figure in the campaign is partly a measurement of the
# reference. If it does not, the reference was noisy but unbiased and the figures
# stand -- which is itself worth knowing rather than assuming.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))
source("scratch/nes_calibration/argos_filter.R")

K <- fread("scratch/nes_calibration/all29_knots.csv")
K[, id := as.character(id)][, time := as.POSIXct(time, tz = "UTC")]
cat(sprintf("%d knots across %d deployments\n", nrow(K), uniqueN(K$id)))

new_argos <- readRDS("analysis/cache/nes/fvilches_argos_v1.rds")
meta <- nes_meta(); ar <- nes_argos()
argos_for <- function(tg) {
  if (!is.null(new_argos[[tg]])) return(list(src = "fvilches", a = as.data.table(new_argos[[tg]])))
  m <- meta[match(tg, meta$id), ]
  if (is.na(m$ptt) || !nzchar(m$ptt)) return(NULL)
  a <- ar$fixes[Ptt == m$ptt]
  if (nrow(a) < 50) NULL else list(src = "2021", a = a)
}

R <- rbindlist(lapply(unique(K$id), function(tg) {
  g <- argos_for(tg); if (is.null(g)) return(NULL)
  k <- K[id == tg][order(time)]
  a <- g$a[order(time)]
  f <- argos_speed_filter(a, vmax = 3)
  tt <- as.numeric(k$time)
  tr <- function(x) list(lon = approx(as.numeric(x$time), lon360(x$lon), tt, rule = 2)$y,
                         lat = approx(as.numeric(x$time), x$lat, tt, rule = 2)$y)
  t0 <- tr(a); t1 <- tr(f)
  e0 <- gc_km(lon360(k$est_lon), k$est_lat, t0$lon, t0$lat)
  e1 <- gc_km(lon360(k$est_lon), k$est_lat, t1$lon, t1$lat)
  data.table(id = tg, src = g$src, n = nrow(k),
             removed = attr(f, "frac_removed"),
             km_raw = median(e0, na.rm = TRUE), km_filt = median(e1, na.rm = TRUE),
             bias_raw = mean(k$est_lat - t0$lat, na.rm = TRUE),
             bias_filt = mean(k$est_lat - t1$lat, na.rm = TRUE),
             rmse_raw = sqrt(mean((k$est_lat - t0$lat)^2, na.rm = TRUE)),
             rmse_filt = sqrt(mean((k$est_lat - t1$lat)^2, na.rm = TRUE)))
}))

cat("\n=== per deployment ===\n")
cat(sprintf("%-9s %-9s %7s %9s %9s %9s %9s\n", "tag", "delivery", "%drop",
            "km raw", "km filt", "bias raw", "bias filt"))
for (i in seq_len(nrow(R)))
  cat(sprintf("%-9s %-9s %6.1f%% %9.0f %9.0f %+9.3f %+9.3f\n", R$id[i], R$src[i],
              100*R$removed[i], R$km_raw[i], R$km_filt[i], R$bias_raw[i], R$bias_filt[i]))

cat("\n=== headline, all 29 ===\n")
cat(sprintf("  median km   : raw %.0f -> filtered %.0f  (change %+.0f)\n",
            median(R$km_raw), median(R$km_filt), median(R$km_filt) - median(R$km_raw)))
cat(sprintf("  mean  km    : raw %.0f -> filtered %.0f  (change %+.0f)\n",
            mean(R$km_raw), mean(R$km_filt), mean(R$km_filt) - mean(R$km_raw)))
cat(sprintf("  mean bias   : raw %+.3f -> filtered %+.3f  (change %+.3f deg)\n",
            mean(R$bias_raw), mean(R$bias_filt), mean(R$bias_filt) - mean(R$bias_raw)))
cat(sprintf("  mean rmse   : raw %.3f -> filtered %.3f\n",
            mean(R$rmse_raw), mean(R$rmse_filt)))

o <- R[src == "2021"]; f <- R[src == "fvilches"]
cat("\n=== split by delivery (only 2021 was unfiltered) ===\n")
cat(sprintf("  2021     (n=%2d): km %.0f -> %.0f | bias %+.3f -> %+.3f\n", nrow(o),
            median(o$km_raw), median(o$km_filt), mean(o$bias_raw), mean(o$bias_filt)))
cat(sprintf("  fvilches (n=%2d): km %.0f -> %.0f | bias %+.3f -> %+.3f  (must be unchanged)\n",
            nrow(f), median(f$km_raw), median(f$km_filt), mean(f$bias_raw), mean(f$bias_filt)))
cat(sprintf("\n  CONTROL: fvilches unchanged? %s\n",
            if (max(abs(f$km_raw - f$km_filt)) < 0.5 &&
                max(abs(f$bias_raw - f$bias_filt)) < 0.005) "YES" else "*** NO"))

p <- suppressWarnings(wilcox.test(o$km_raw, o$km_filt, paired = TRUE)$p.value)
cat(sprintf("  2021 paired km change: median %+.1f km, p = %.4f\n",
            median(o$km_filt - o$km_raw), p))
pb <- suppressWarnings(wilcox.test(o$bias_raw, o$bias_filt, paired = TRUE)$p.value)
cat(sprintf("  2021 paired bias change: median %+.4f deg, p = %.4f\n",
            median(o$bias_filt - o$bias_raw), pb))
fwrite(R, "scratch/nes_calibration/rescore_filtered_truth.csv")
