# NOON versus MIDNIGHT as the longitude estimator, scored against Argos.
#
# Two estimates per day from the same light record:
#   noon_d     = midpoint of the LIT run  (bounded by dawn_d and dusk_d)
#   midnight_d = midpoint of the DARK run (bounded by dusk_d and dawn_{d+1})
#
# They are not independent -- each shares one edge with the other -- but they are
# different combinations of edges, so a corrupted edge hits them differently. Both
# are built from threshold crossings rather than from the peak, which matters
# because midday light is CLIPPED in these tags (daytime median 181-186 against a
# fitted max_light of 159.6), so the daytime maximum carries almost no timing
# information. The night minimum sits inside the dynamic range.
#
# Reported:
#   1. longitude error of each estimator against Argos
#   2. availability -- near solstice at high latitude the dark run can vanish
#   3. error after a 3-day rolling median, which is the question that actually
#      matters for a smoothing step
#   4. the identity residual on real data:
#        (noon_d - midnight_d) - (dawn_{d+1} - dawn_d)/2  should equal  dL/4
#      which is a self-consistency test of the estimated edges
#
# Runs on panel_v1, i.e. the SHIPPED max decimation.
suppressMessages({ library(data.table) })
D2R <- pi/180
lon180 <- function(x) ((x + 180) %% 360) - 180
eot_h <- function(t) { n <- as.numeric(format(t, "%j")); B <- 2*pi*(n - 81)/364
  (9.87*sin(2*B) - 7.53*cos(B) - 1.5*sin(B)) / 60 }
gc_km <- function(l1, p1, l2, p2, R = 6371.0088) { d <- D2R
  2*R*asin(pmin(1, sqrt(sin((p2-p1)*d/2)^2 + cos(p1*d)*cos(p2*d)*sin((l2-l1)*d/2)^2))) }

P <- readRDS("analysis/cache/nes/panel_v1.rds")
MATCH_H <- 6      # Argos fix must be within this of the estimate, so noon and
                  # midnight are not scored against the same fix

rows <- list()
for (tg in names(P$tags)) {
  g <- P$tags[[tg]]; r <- P$responses[[tg]]
  L <- as.data.table(g$light)[, .(time, light)]; setorder(L, time)
  A <- as.data.table(g$argos)[is.finite(lon) & is.finite(lat)]; setorder(A, time)
  if (!nrow(A)) next
  thr <- r$baseline + (r$max_light - r$baseline) / 2
  L[, lit := light >= thr]
  rr <- rle(L$lit)
  ends <- cumsum(rr$lengths); starts <- ends - rr$lengths + 1L
  seg <- data.table(lit = rr$values, i0 = starts, i1 = ends)
  seg[, `:=`(t0 = as.numeric(L$time[i0]), t1 = as.numeric(L$time[i1]))]
  seg[, dur_h := (t1 - t0) / 3600]
  seg <- seg[dur_h >= 4 & dur_h <= 20]
  if (nrow(seg) < 4) next
  seg[, mid := (t0 + t1) / 2]
  atn <- as.numeric(A$time)                  # numeric throughout: abs() on POSIXt errors
  at <- function(tt) {                       # Argos truth at a time
    j <- pmax(1L, pmin(length(atn), findInterval(tt, atn)))
    j2 <- pmin(length(atn), j + 1L)
    pick <- ifelse(abs(atn[j] - tt) <= abs(atn[j2] - tt), j, j2)
    ok <- abs(atn[pick] - tt) <= MATCH_H * 3600
    list(lon = ifelse(ok, A$lon[pick], NA_real_),
         lat = ifelse(ok, A$lat[pick], NA_real_))
  }
  for (k in seq_len(nrow(seg))) {
    tt <- seg$mid[k]
    tc <- as.POSIXct(tt, origin = "1970-01-01", tz = "UTC")
    e <- eot_h(tc)
    uh <- (tt %% 86400) / 3600
    est <- if (seg$lit[k]) lon180(15 * (12 - uh - e)) else lon180(15 * (0 - uh - e))
    tr <- at(tt)
    if (!is.finite(tr$lon)) next
    rows[[length(rows) + 1L]] <- data.table(
      id = tg, kind = if (seg$lit[k]) "noon" else "midnight",
      t = tt, est_lon = est, true_lon = tr$lon, true_lat = tr$lat,
      dur_h = seg$dur_h[k])
  }
}
R <- rbindlist(rows)
R[, err := lon180(est_lon - true_lon)]
R[, km := abs(err) * 111 * cos(true_lat * D2R)]
R <- R[abs(err) < 30]        # drop gross failures symmetric between estimators

cat("=== longitude error by estimator, all 29 tags ===\n")
s <- R[, .(n = .N, tags = uniqueN(id), med_err = median(err),
           robust_sd = IQR(err) / 1.349,
           med_abs_km = median(km)), by = kind]
print(as.data.frame(s[, .(kind, n, tags, med_err = round(med_err, 2),
                          robust_sd_deg = round(robust_sd, 2),
                          med_abs_km = round(med_abs_km))]), row.names = FALSE)

cat("\n=== availability: estimates per tag ===\n")
av <- dcast(R[, .N, by = .(id, kind)], id ~ kind, value.var = "N", fill = 0)
cat("noon estimates per tag    (min|med|max):",
    paste(round(quantile(av$noon, c(0, .5, 1))), collapse = " | "), "\n")
cat("midnight estimates per tag(min|med|max):",
    paste(round(quantile(av$midnight, c(0, .5, 1))), collapse = " | "), "\n")
cat("tags where midnight is scarcer than noon:", sum(av$midnight < av$noon), "of", nrow(av), "\n")

cat("\n=== after smoothing (the question that matters for a smoothing step) ===\n")
cat("Window is TIME-based, +/- HALF_D days, not a fixed point count: a 7-point\n")
cat("window spans 7 days for a once-daily series but only 3.5 for the interleaved\n")
cat("one, which would have made BOTH look worse purely by being less smoothed.\n\n")
HALF_D <- 1.5
sm <- function(d) { setorder(d, t)
  tt <- d$t; v <- d$est_lon
  d[, sm := vapply(seq_along(tt), function(i) {
      k <- abs(tt - tt[i]) <= HALF_D * 86400
      if (sum(k) < 3) NA_real_ else median(v[k]) }, 0)]
  d[is.finite(sm), .(n = .N, robust_sd = IQR(lon180(sm - true_lon)) / 1.349,
                     med_abs_km = median(abs(lon180(sm - true_lon)) * 111 * cos(true_lat * D2R)))] }
out <- list()
for (k in c("noon", "midnight")) {
  z <- rbindlist(lapply(split(R[kind == k], R[kind == k]$id), sm))
  out[[k]] <- c(sd = weighted.mean(z$robust_sd, z$n),
                km = weighted.mean(z$med_abs_km, z$n), n = sum(z$n))
}
zb <- rbindlist(lapply(split(R, R$id), sm))
out[["BOTH"]] <- c(sd = weighted.mean(zb$robust_sd, zb$n),
                   km = weighted.mean(zb$med_abs_km, zb$n), n = sum(zb$n))
for (k in names(out))
  cat(sprintf("  %-9s robust SD %5.2f deg | median |err| %4.0f km  (n=%d)\n",
              k, out[[k]]["sd"], out[[k]]["km"], out[[k]]["n"]))

cat("\n=== identity check on real data ===\n")
cat("(noon_d - midnight_d) - (dawn_d+1 - dawn_d)/2 should equal dL/4\n")
idc <- list()
for (tg in unique(R$id)) {
  g <- P$tags[[tg]]; r <- P$responses[[tg]]
  L <- as.data.table(g$light)[, .(time, light)]; setorder(L, time)
  thr <- r$baseline + (r$max_light - r$baseline) / 2
  L[, lit := light >= thr]; rr <- rle(L$lit)
  ends <- cumsum(rr$lengths); starts <- ends - rr$lengths + 1L
  s <- data.table(lit = rr$values, t0 = as.numeric(L$time[starts]),
                  t1 = as.numeric(L$time[ends]))
  s[, dur := (t1 - t0) / 3600]; s <- s[dur >= 4 & dur <= 20]
  # `lit` is also a COLUMN of s, and data.table evaluates i-expressions in the
  # frame of the table, so a variable of that name gets shadowed by the column.
  LIT <- s[lit == TRUE]; DRK <- s[lit == FALSE]
  if (nrow(LIT) < 3 || nrow(DRK) < 3) next
  for (i in seq_len(nrow(LIT) - 1L)) {
    dawn0 <- LIT$t0[i]; dawn1 <- LIT$t0[i + 1L]; dusk0 <- LIT$t1[i]
    if (dawn1 - dawn0 > 30*3600 || dawn1 - dawn0 < 18*3600) next
    dk <- DRK[t0 >= dusk0 & t1 <= dawn1]
    if (nrow(dk) != 1L) next
    noon <- (LIT$t0[i] + LIT$t1[i]) / 2; midn <- (dk$t0 + dk$t1) / 2
    dL <- (LIT$dur[i + 1L] - LIT$dur[i]) * 3600
    idc[[length(idc) + 1L]] <- data.table(
      id = tg, lhs_min = (noon - midn) / 60 + 720,   # midnight follows noon here
      rhs_min = (dawn1 - dawn0) / 2 / 60, dL4_min = dL / 4 / 60)
  }
}
I <- rbindlist(idc)
I[, resid := lhs_min - rhs_min]
cat(sprintf("  n = %d day-pairs across %d tags\n", nrow(I), uniqueN(I$id)))
cat(sprintf("  residual (lhs - rhs), min : median %+6.2f | IQR %5.2f\n",
            median(I$resid), IQR(I$resid)))
cat(sprintf("  predicted dL/4,       min : median %+6.2f | IQR %5.2f\n",
            median(I$dL4_min), IQR(I$dL4_min)))
cat(sprintf("  correlation(residual, dL/4) = %+.3f\n", cor(I$resid, I$dL4_min)))
cat(sprintf("  SD of (residual - dL/4) = %.2f min  <- edge-estimation noise\n",
            IQR(I$resid - I$dL4_min) / 1.349))
