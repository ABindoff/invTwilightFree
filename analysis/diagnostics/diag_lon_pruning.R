# How much computation would longitude pruning actually save, and how wide does the
# retained band have to be?
#
# The premise: longitude is the carrier PHASE and is well identified (rmse ~1.15 deg,
# coverage 0.96-1.00), while latitude is the hard, weakly identified axis. So most of
# the grid is longitude columns that never hold posterior mass, and evaluating the
# emission there is wasted. Once longitude is pinned the problem is nearly 1-D.
#
# WHAT HAS TO BE MEASURED BEFORE BUILDING IT:
#  1. how many longitude columns actually hold the mass, per knot -- this sets the
#     saving AND the margin. Sizing the band from the marginal posterior is the only
#     honest way; guessing a window risks clipping the truth on tags where longitude
#     is worse than average.
#  2. the WORST knot, not the median. A band that works 99% of the time and clips a
#     handful of knots will corrupt those knots silently.
#  3. whether the required band is stable across tags, or whether some tags need a
#     much wider one (a clock offset would do that: longitude is set by the clock, and
#     the campaign's clock correction was tried and EXCLUDED at +863 km, so residual
#     clock error is real and unmodelled).
#
# Reported as the fraction of columns needed for 1 - eps of the marginal mass, for
# several eps.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))
CAL_DAYS <- 15; STEP_HOURS <- 12; DIFFUSION <- 110; CELL <- 1

old <- lapply(readRDS("analysis/cache/nes/archives_v3_30min.rds"), `[[`, "main")
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
  tf <- a[time >= max(time) - 72*3600]
  tags[[tg]] <- list(id = tg, light = d, argos = a,
                     p0 = c(df$deploy_lon, df$deploy_lat),
                     p1 = c(median(tf$lon), median(tf$lat)))
}
tags <- tags[seq_len(min(5, length(tags)))]
resp <- lapply(tags, function(g) {
  k <- g$light$time <= min(g$light$time) + CAL_DAYS * 86400
  fit_light_response(g$light$time[k], g$light$light[k], g$p0[1], g$p0[2])
})

grid <- rast(xmin = 150, xmax = 250, ymin = 20, ymax = 70,
             resolution = CELL, crs = "EPSG:4326")
values(grid) <- 1
EPS <- c(1e-2, 1e-4, 1e-6, 1e-9)

cat(sprintf("grid: %d lon columns x %d lat rows = %d cells\n\n",
            length(seq(150.5, 249.5, CELL)), length(seq(20.5, 69.5, CELL)),
            length(seq(150.5, 249.5, CELL)) * length(seq(20.5, 69.5, CELL))))
cat(sprintf("%-10s %7s %s\n", "tag", "knots",
            paste(sprintf("%18s", sprintf("eps=%.0e med/max", EPS)), collapse = "")))

ALL <- list()
for (tg in names(tags)) {
  g <- tags[[tg]]; r <- resp[[tg]]
  invisible(capture.output(
    f <- TwilightFreeGrid(g$light$time, pmax(0, g$light$light - r$baseline), grid = grid,
      start_lon = unname(g$p0[1]), start_lat = unname(g$p0[2]),
      end_lon = unname(g$p1[1]), end_lat = unname(g$p1[2]),
      step_hours = STEP_HOURS, diffusion = DIFFUSION,
      calibration = r$calibration,
      likelihood_params = c(1/(r$max_light*0.5), r$max_light, 0.10),
      area_correction = TRUE)))
  gp <- grid_posterior(f)
  ulon <- sort(unique(gp$lon)); ncol_tot <- length(ulon)
  # marginal longitude posterior per knot, then the smallest CONTIGUOUS band holding
  # 1-eps of it (contiguous, because a pruning rule has to keep an interval)
  need <- matrix(NA_real_, nrow(gp$P), length(EPS))
  for (k in seq_len(nrow(gp$P))) {
    w <- gp$P[k, ]; s <- sum(w); if (!is.finite(s) || s <= 0) next
    m <- vapply(ulon, function(L) sum(w[gp$lon == L]), 0) / s
    ord <- order(m, decreasing = TRUE)
    cs <- cumsum(m[ord])
    for (e in seq_along(EPS)) {
      kk <- which(cs >= 1 - EPS[e])[1]
      sel <- ulon[ord[seq_len(kk)]]
      need[k, e] <- (max(sel) - min(sel)) / CELL + 1   # contiguous span
    }
  }
  ALL[[tg]] <- need
  cat(sprintf("%-10s %7d %s\n", tg, nrow(gp$P),
              paste(sprintf("%18s", sprintf("%.0f / %.0f",
                    apply(need, 2, median, na.rm = TRUE),
                    apply(need, 2, max, na.rm = TRUE))), collapse = "")))
}

N <- do.call(rbind, ALL)
ncol_tot <- length(seq(150.5, 249.5, CELL))
cat(sprintf("\npooled over %d knots, columns needed out of %d:\n", nrow(N), ncol_tot))
for (e in seq_along(EPS))
  cat(sprintf("  eps = %.0e : median %3.0f | q99 %3.0f | MAX %3.0f  -> speedup %.1fx (median) %.1fx (worst-case band)\n",
              EPS[e], median(N[, e], na.rm = TRUE), quantile(N[, e], .99, na.rm = TRUE),
              max(N[, e], na.rm = TRUE),
              ncol_tot / median(N[, e], na.rm = TRUE),
              ncol_tot / max(N[, e], na.rm = TRUE)))
cat("\nThe MAX column is the one that matters for correctness: a band sized on the\n")
cat("median clips the worst knots silently. The speedup that can actually be claimed\n")
cat("is the worst-case one, unless the band is sized per knot.\n")
