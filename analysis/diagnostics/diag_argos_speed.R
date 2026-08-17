# How bad is the unfiltered Argos reference?
#
# Every accuracy and bias figure in this campaign is measured against Argos fixes
# interpolated to knot times, and those fixes have never been speed-filtered. If a
# meaningful fraction imply impossible movement then the reference is degenerate:
# the "truth" wanders somewhere the animal never was, the clear-sky light computed at
# it is wrong, and the -1.54 degree bias is partly a measurement of the reference
# rather than of the method.
#
# Measure it first, filter second. A northern elephant seal transits at roughly
# 0.7-1.0 m/s and can sustain about 2 m/s; anything far above that is a location
# error, not an animal.
#
# Aggregates only; no per-fix rows printed.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))

new_argos <- readRDS("analysis/cache/nes/fvilches_argos_v1.rds")
ar <- nes_argos()
meta <- nes_meta()

# assemble every deployment's Argos series, both deliveries
A <- list()
for (tg in names(new_argos)) A[[tg]] <- as.data.table(new_argos[[tg]])
for (i in seq_len(nrow(meta))) {
  p <- meta$ptt[i]; if (is.na(p) || !nzchar(p)) next
  f <- ar$fixes[Ptt == p]
  if (nrow(f) >= 50) A[[meta$id[i]]] <- f
}
cat(sprintf("%d deployments with Argos\n\n", length(A)))

step_speed <- function(a) {
  a <- unique(a[order(time)], by = "time")          # duplicate timestamps -> Inf speed
  n <- nrow(a); if (n < 3) return(NULL)
  d <- gc_km(lon360(a$lon[-n]), a$lat[-n], lon360(a$lon[-1]), a$lat[-1])
  dt <- as.numeric(difftime(a$time[-1], a$time[-n], units = "hours"))
  ok <- dt > 0
  list(v = d[ok] / dt[ok] * 1000 / 3600, dt = dt[ok], n = n)   # m/s
}

S <- rbindlist(lapply(names(A), function(tg) {
  s <- step_speed(A[[tg]]); if (is.null(s)) return(NULL)
  data.table(id = tg, n_fix = s$n, med_gap_h = median(s$dt),
             v_med = median(s$v), v_q95 = quantile(s$v, .95), v_max = max(s$v),
             f_gt2 = mean(s$v > 2), f_gt5 = mean(s$v > 5), f_gt10 = mean(s$v > 10))
}))
cat("implied step speeds (m/s) between consecutive Argos fixes\n")
cat(sprintf("%-9s %6s %8s %7s %7s %9s %8s %8s %8s\n", "tag", "n_fix", "gap_h",
            "v_med", "v_q95", "v_max", ">2m/s", ">5m/s", ">10m/s"))
for (i in seq_len(nrow(S)))
  cat(sprintf("%-9s %6d %8.2f %7.2f %7.2f %9.0f %8.3f %8.3f %8.3f\n",
              S$id[i], S$n_fix[i], S$med_gap_h[i], S$v_med[i], S$v_q95[i],
              S$v_max[i], S$f_gt2[i], S$f_gt5[i], S$f_gt10[i]))

allv <- unlist(lapply(A, function(a) { s <- step_speed(a); if (is.null(s)) NULL else s$v }))
cat(sprintf("\npooled over %d steps:\n", length(allv)))
cat(sprintf("  median %.2f | q90 %.2f | q99 %.2f | max %.0f m/s\n",
            median(allv), quantile(allv, .9), quantile(allv, .99), max(allv)))
for (v in c(2, 3, 5, 10, 20))
  cat(sprintf("  above %2d m/s: %6.3f%% of steps\n", v, 100 * mean(allv > v)))
cat(sprintf("\nfor scale: a northern elephant seal transits at 0.7-1.0 m/s and sustains ~2 m/s.\n"))
cat(sprintf("%.0f m/s is %.0f km/h.\n", max(allv), max(allv) * 3.6))

# what does an implausible fix cost in position terms?
cat("\nhow far do the worst fixes displace the reference?\n")
D <- rbindlist(lapply(names(A), function(tg) {
  a <- unique(A[[tg]][order(time)], by = "time"); n <- nrow(a); if (n < 3) return(NULL)
  d <- gc_km(lon360(a$lon[-n]), a$lat[-n], lon360(a$lon[-1]), a$lat[-1])
  dt <- as.numeric(difftime(a$time[-1], a$time[-n], units = "hours"))
  v <- d / pmax(dt, 1e-6) * 1000/3600
  data.table(id = tg, km = d[is.finite(v) & v > 5])
}))
if (nrow(D)) cat(sprintf("  steps faster than 5 m/s: n = %d, median jump %.0f km, max %.0f km\n",
                         nrow(D), median(D$km), max(D$km)))
