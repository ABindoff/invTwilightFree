# Are there out-and-back excursions a speed filter cannot catch?
#
# A speed filter permits vmax TIMES THE GAP. On a sparse track that is a lot of rope:
# at 3 m/s a 10-hour gap allows 108 km, so a single fix displaced 100 km passes the
# speed test and produces a triangular detour -- out and straight back. The signature
# is a large DETOUR RATIO,
#      (d(i-1,i) + d(i,i+1)) / d(i-1,i+1)
# which is ~1 for travel along a path and large for a spike. It is scale-free, so it
# works where a speed threshold does not.
#
# Reported per deployment, with the gap distribution, because the two interact: dense
# tracks are protected by the speed filter, sparse ones are not.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))
source("scratch/nes_calibration/argos_filter.R")
new_argos <- readRDS("analysis/cache/nes/fvilches_argos_v1.rds")
meta <- nes_meta(); ar <- nes_argos()
old <- lapply(readRDS("analysis/cache/nes/archives_v3_30min.rds"), `[[`, "main")

argos_for <- function(tg) {
  if (!is.null(new_argos[[tg]])) return(list(src = "fvilches", a = as.data.table(new_argos[[tg]])))
  m <- meta[match(tg, meta$id), ]
  if (is.na(m$ptt) || !nzchar(m$ptt)) return(NULL)
  a <- ar$fixes[Ptt == m$ptt]; if (nrow(a) < 50) NULL else list(src = "2021", a = a)
}
ids <- unique(c(names(old), names(new_argos)))

detour <- function(a) {
  a <- a[order(a$time), ]; n <- nrow(a); if (n < 3) return(NULL)
  lo <- lon360(a$lon); la <- a$lat
  d_in  <- gc_km(lo[1:(n-2)], la[1:(n-2)], lo[2:(n-1)], la[2:(n-1)])
  d_out <- gc_km(lo[2:(n-1)], la[2:(n-1)], lo[3:n],     la[3:n])
  d_dir <- gc_km(lo[1:(n-2)], la[1:(n-2)], lo[3:n],     la[3:n])
  list(r = (d_in + d_out) / pmax(d_dir, 0.5),      # 0.5 km floor: Argos noise
       excursion = pmin(d_in, d_out),               # how far out the spike goes
       gap = as.numeric(difftime(a$time[3:n], a$time[1:(n-2)], units = "hours")))
}

R <- rbindlist(lapply(ids, function(tg) {
  g <- argos_for(tg); if (is.null(g)) return(NULL)
  f <- argos_speed_filter(g$a, vmax = 3)          # after the speed filter
  d <- detour(f); if (is.null(d)) return(NULL)
  data.table(id = tg, src = g$src, n = nrow(f),
             gap_med = median(d$gap, na.rm = TRUE),
             r_med = median(d$r), r_q99 = quantile(d$r, .99), r_max = max(d$r),
             f_r5 = mean(d$r > 5), f_r10 = mean(d$r > 10),
             worst_km = max(d$excursion[d$r > 10], na.rm = TRUE))
}))
R[!is.finite(worst_km), worst_km := NA]
setorder(R, -f_r10)

cat("AFTER the 3 m/s speed filter: out-and-back excursions that survive it\n")
cat("detour ratio r = (in + out) / direct;  r ~ 1 is travel, r >> 1 is a spike\n\n")
cat(sprintf("%-9s %-9s %6s %8s %7s %8s %8s %8s %10s\n", "tag", "delivery", "n_fix",
            "gap_h", "r_med", "r_q99", "%r>5", "%r>10", "worst_km"))
for (i in seq_len(nrow(R)))
  cat(sprintf("%-9s %-9s %6d %8.2f %7.2f %8.1f %7.1f%% %7.1f%% %10.0f\n",
              R$id[i], R$src[i], R$n[i], R$gap_med[i], R$r_med[i], R$r_q99[i],
              100*R$f_r5[i], 100*R$f_r10[i], R$worst_km[i]))

cat(sprintf("\nsparse tracks are the exposed ones: cor(median gap, %% with r>10) = %+.3f\n",
            cor(R$gap_med, R$f_r10)))
b <- R[id == "2023041"]
if (nrow(b)) cat(sprintf("\n2023041 specifically: %d fixes, median gap %.1f h, %.1f%% of triples have r>10, worst excursion %.0f km\n",
                         b$n, b$gap_med, 100*b$f_r10, b$worst_km))
fwrite(R, "scratch/nes_calibration/detour_stats.csv")
