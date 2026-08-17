# What does the speed filter remove, and how far does it move the "truth"?
#
# Two things have to be true for the filter to be worth applying:
#  1. it must bring the 2021 delivery into line with the already-filtered fvilches
#     one, otherwise the reference stays inconsistent between deliveries;
#  2. it must be INERT on the fvilches tags, which are already filtered at 3 m/s --
#     if it removes fixes there too, it is deleting real movement, not errors.
# The second is the control.
#
# Then: how far does the interpolated ground truth actually move? That is the number
# that says whether the campaign's error and bias figures need recomputing.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))
source("scratch/nes_calibration/argos_filter.R")

old <- lapply(readRDS("analysis/cache/nes/archives_v3_30min.rds"), `[[`, "main")
new <- lapply(readRDS("analysis/cache/nes/fvilches_archives_v1_30min.rds"), `[[`, "main")
new_argos <- readRDS("analysis/cache/nes/fvilches_argos_v1.rds")
meta <- nes_meta(); ar <- nes_argos()

REC <- list()
for (tg in names(old)) {
  m <- meta[match(tg, meta$id), ]
  if (is.na(m$ptt) || !nzchar(m$ptt)) next
  a <- ar$fixes[Ptt == m$ptt]; if (nrow(a) < 50) next
  REC[[tg]] <- list(src = "2021", a = a[order(time)], l = as.data.table(old[[tg]]))
}
for (tg in names(new)) {
  a <- as.data.table(new_argos[[tg]]); if (is.null(a) || !nrow(a)) next
  REC[[tg]] <- list(src = "fvilches", a = a[order(time)], l = as.data.table(new[[tg]]))
}

spd <- function(a) {
  a <- a[order(a$time), ]; n <- nrow(a); if (n < 3) return(numeric(0))
  d <- gc_km(lon360(a$lon[-n]), a$lat[-n], lon360(a$lon[-1]), a$lat[-1])
  dt <- as.numeric(difftime(a$time[-1], a$time[-n], units = "hours"))
  v <- d[dt > 0] / dt[dt > 0] * 1000/3600; v[is.finite(v)]
}

R <- rbindlist(lapply(names(REC), function(tg) {
  r <- REC[[tg]]
  f <- argos_speed_filter(r$a, vmax = 3)
  v0 <- spd(as.data.frame(r$a)); v1 <- spd(f)
  # how far does the interpolated truth move, at the light record's own times?
  tt <- as.numeric(r$l$time)
  la0 <- approx(as.numeric(r$a$time), r$a$lat, tt, rule = 2)$y
  lo0 <- approx(as.numeric(r$a$time), lon360(r$a$lon), tt, rule = 2)$y
  la1 <- approx(as.numeric(f$time), f$lat, tt, rule = 2)$y
  lo1 <- approx(as.numeric(f$time), lon360(f$lon), tt, rule = 2)$y
  shift <- gc_km(lo0, la0, lo1, la1)
  data.table(id = tg, src = r$src, n0 = nrow(r$a), removed = attr(f, "n_removed"),
             frac = attr(f, "frac_removed"), iters = attr(f, "iterations"),
             vmax0 = max(v0), vmax1 = if (length(v1)) max(v1) else NA_real_,
             shift_med = median(shift), shift_q95 = quantile(shift, .95),
             shift_max = max(shift), lat_shift_med = median(abs(la1 - la0)))
}))

cat("=== what the filter removes ===\n")
cat(sprintf("%-9s %-9s %6s %8s %7s %10s %10s\n", "tag", "delivery", "n_fix",
            "removed", "%", "v_max before", "v_max after"))
for (i in seq_len(nrow(R)))
  cat(sprintf("%-9s %-9s %6d %8d %6.1f%% %10.0f %10.2f\n", R$id[i], R$src[i], R$n0[i],
              R$removed[i], 100*R$frac[i], R$vmax0[i], R$vmax1[i]))

cat("\n=== CONTROL: the filter must be inert on the already-filtered delivery ===\n")
fv <- R[src == "fvilches"]
cat(sprintf("  fvilches: removed %d fixes across %d tags (%.3f%%)  -> %s\n",
            sum(fv$removed), nrow(fv), 100*mean(fv$frac),
            if (sum(fv$removed) == 0) "INERT, as required"
            else "*** removes real fixes: threshold too aggressive"))
o <- R[src == "2021"]
cat(sprintf("  2021:     removed %d fixes across %d tags (%.1f%% median per tag)\n",
            sum(o$removed), nrow(o), 100*median(o$frac)))
cat(sprintf("  v_max after filtering: 2021 %.2f | fvilches %.2f  -> %s\n",
            max(o$vmax1), max(fv$vmax1),
            if (abs(max(o$vmax1) - max(fv$vmax1)) < 0.5) "deliveries now consistent"
            else "still inconsistent"))

cat("\n=== how far does the ground truth move? ===\n")
cat(sprintf("%-9s %-9s %10s %10s %10s %12s\n", "tag", "delivery", "med km", "q95 km",
            "max km", "med |dlat| deg"))
for (i in seq_len(nrow(R)))
  cat(sprintf("%-9s %-9s %10.1f %10.1f %10.0f %12.3f\n", R$id[i], R$src[i],
              R$shift_med[i], R$shift_q95[i], R$shift_max[i], R$lat_shift_med[i]))
cat(sprintf("\n  2021 tags: median shift %.1f km, q95 %.1f km, median |dlat| %.3f deg\n",
            median(o$shift_med), median(o$shift_q95), median(o$lat_shift_med)))
cat("  For scale, the campaign reports 276 km accuracy and -1.54 deg latitude bias.\n")
cat("  A shift comparable to either means those figures were measured against a\n")
cat("  reference that itself moves by that much, and must be recomputed.\n")
fwrite(R, "scratch/nes_calibration/argos_filter_effect.csv")
