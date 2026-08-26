# SURFACE-CONDITIONED DECIMATION. Replaces `decimate_dive()` in ingest_fvilches.R.
#
# WHY -------------------------------------------------------------------------
# The shipped decimator takes the MAXIMUM light over each 30-minute window. That
# was a deliberate de-shading choice and it works, because the max recovers the
# surface value whenever the animal surfaced. But a maximum is an extreme order
# statistic, so it is biased upward, and the bias is measurable:
#
#   * within surfacing bouts, max minus median is 4.0-4.5 light units at twilight,
#     which at 3.74 units per degree of zenith is ~1.1 deg
#   * on a day-length estimator the run-extension inflates day length by +26.6 min,
#     which at ~0.20 deg latitude per minute predicts +5.3 deg; the measured median
#     latitude error of that estimator is +5.23 deg
#
# The shipped decimator did get one thing exactly right, which the first version of
# THIS script broke: it stamped the retained value with the time of the sample it
# came from, so the (time, light) pair was always a real observation.
#
# THE BUG IN v1, AND THE FIX --------------------------------------------------
# v1 pooled every surface sample in a 30-minute window, took the MEDIAN light but
# the MEAN time, and emitted one observation per window. Those two estimators
# weight samples differently. Surface samples inside one window span a median of
# 16 minutes (51% over 10 min, 32% over 20 min) because a window usually holds
# more than one surfacing bout, so the light came from whichever bout had more
# samples while the timestamp landed between the bouts, at a moment the animal was
# at depth. That is a 2.8 deg zenith smear at the median and 5.1 deg at the 90th
# percentile, against the 1.1 deg of max bias it was meant to remove. Arm A on v1
# was 471 km against 177 km shipped, with the seasonal swing 6.0 -> 15.7 deg.
#
# v2 fixes it two ways:
#
#   A. ONE OBSERVATION PER SURFACING, not per window. Contiguous surface runs are
#      split into chunks of at most MAX_SPAN_MIN, so no observation ever blends
#      light across more than that. Bounds the smear at ~0.9 deg of zenith and
#      removes the two-bouts-in-one-window failure entirely.
#
#   B. CONSISTENT ESTIMATORS. Light and time are both MEDIANS over the same
#      samples. If light is monotone in time across the chunk, which it is over a
#      few minutes, then median(light) = f(median(time)) exactly, because a
#      monotone map preserves order statistics. The timestamp is also a real
#      sample time rather than a constructed midpoint. Mean/mean would only be
#      unbiased for a linear ramp; median/median holds for any monotone one.
#
# ALSO FIXED, and it is not optional ------------------------------------------
# The depth baseline drifts through the deployment: median +1.5 m, up to +5.5 m,
# almost always upward, and 10 of 29 tags move more than 2.5 m. A FIXED threshold
# therefore changes meaning over time, and deployment phase correlates with
# latitude because the animals migrate out and back. Without this correction a
# 2.5 m threshold retains 99% of windows on one tag and 5.8% on another. The zero
# is re-estimated as a slowly varying low quantile.
#
# NOT APPLIED -----------------------------------------------------------------
# Residual attenuation. Mk9 `Light Level` is already log-compressed irradiance, so
# attenuation is additive and correctable with one coefficient, but retained
# samples sit ~1 m below the corrected zero, worth ~0.6 light units or 0.16 deg of
# zenith. The coefficient is also poorly identified (direct fit 0.941 units/m,
# this estimator 0.591, the gap being saturation flattening shallow bins). k is
# COMPUTED and STORED per tag but not applied; ATTEN_ON runs it as its own arm.
#
# CONSEQUENCE FOR THE GATE ARMS -----------------------------------------------
# Every observation is now a surface observation, so `depth_min` is ~0 throughout
# and the shipped `depth_min > 10 m` gate is a no-op on this decimation. Surface
# conditioning subsumes it. The columns are kept so downstream scripts still run.
#
# Usage:  Rscript analysis/diagnostics/decimate_surface.R [tag ...]
#         no arguments -> all archives -> analysis/cache/nes/archives_surface_v1_30min.rds
#         with tags    -> those tags only, diagnostics printed, nothing written
suppressMessages({ library(data.table) })
Sys.setlocale("LC_TIME", "C")

TDR          <- "fvilches/extracted/TDR raw"
OUT          <- "analysis/cache/nes/archives_surface_v1_30min.rds"
SURF_M       <- 2.5     # metres above the drift-corrected zero
MAX_SPAN_MIN <- 5       # no observation blends light across more than this
MIN_N        <- 5       # samples required per observation (20 s at 4 s sampling)
THIN_WINDOW  <- TRUE    # keep one observation per 30-min window (see header)
WINDOW_SEC   <- 1800
ATTEN_ON     <- FALSE   # apply the per-tag attenuation correction (see header)
ZERO_Q       <- 0.005   # quantile of daily depth taken as the tag's zero
ZERO_SMOOTH  <- 7       # days, running median applied to the daily zero
ATTEN_MAX_M  <- 60      # depth range over which light is linear in depth
LIGHT_NEED   <- c("Time", "Depth", "Light Level")
TEMP_NAMES   <- c("External Temp", "External Temperature")

# ---- the drifting zero -------------------------------------------------------
depth_zero <- function(time, depth) {
  d <- data.table(day = as.integer(as.numeric(time) %/% 86400), depth = depth)
  z <- d[is.finite(depth), .(z = as.numeric(quantile(depth, ZERO_Q))), by = day][order(day)]
  if (nrow(z) < 3) return(rep(median(z$z), length(time)))
  k <- max(3L, ZERO_SMOOTH) %/% 2L
  zs <- vapply(seq_len(nrow(z)), function(i)
    median(z$z[max(1L, i - k):min(nrow(z), i + k)]), 0)
  stats::approx(z$day + 0.5, zs, xout = as.numeric(time) / 86400, rule = 2)$y
}

# ---- the per-tag attenuation coefficient -------------------------------------
# Daylight must be flagged by something independent of the sample's own depth, or
# the estimate collapses: selecting bright samples keeps only the brightest of
# each deep bin and flattens the slope. The flag comes from the enclosing
# window's SURFACE light, which is depth-independent by construction.
atten_coef <- function(light, depth_c, bright) {
  ok <- bright & is.finite(light) & is.finite(depth_c) &
        depth_c >= 0 & depth_c <= ATTEN_MAX_M
  if (sum(ok) < 2000) return(NA_real_)
  b <- data.table(bin = round(depth_c[ok] / 5) * 5, li = light[ok])
  s <- b[, .(n = .N, med = median(li)), by = bin][n >= 200][order(bin)]
  if (nrow(s) < 4) return(NA_real_)
  as.numeric(-coef(lm(med ~ bin, data = s))[2])
}

decimate_surface <- function(d) {
  d <- d[is.finite(depth) & is.finite(light) & !is.na(time)]
  if (nrow(d) < 1000) return(NULL)
  setorder(d, time)
  d[, z0 := depth_zero(time, depth)]
  d[, depth_c := depth - z0]

  # attenuation coefficient, from windows whose surface light is high
  d[, .w := as.numeric(time) %/% 1800]
  s1 <- d[depth_c < SURF_M, .(sl = as.numeric(median(light))), by = .w]
  d[, sl := s1$sl[match(.w, s1$.w)]]
  cut_bright <- as.numeric(quantile(s1$sl, 0.6, na.rm = TRUE))
  k <- atten_coef(d$light, d$depth_c, !is.na(d$sl) & d$sl >= cut_bright)
  if (!is.finite(k) || k < 0 || k > 10) k <- 0
  d[, light_c := if (ATTEN_ON) light + k * pmax(0, depth_c) else as.numeric(light)]

  # ---- one observation per surfacing, chunked so no observation spans long ----
  s <- d[depth_c < SURF_M]
  if (nrow(s) < MIN_N) return(NULL)
  ts <- as.numeric(s$time)
  # a new run starts wherever the surface samples are not consecutive in time
  gap <- c(TRUE, diff(ts) > 2 * median(diff(as.numeric(d$time))))
  run <- cumsum(gap)
  # split any run that lasts longer than MAX_SPAN_MIN into equal-time chunks
  s[, run := run]
  s[, chunk := paste(run, (ts - ave(ts, run, FUN = min)) %/% (MAX_SPAN_MIN * 60))]
  out <- s[, .(time = as.POSIXct(median(as.numeric(time)), origin = "1970-01-01", tz = "UTC"),
               light = as.numeric(median(light_c)),
               temp_surf = as.numeric(temp[which.min(depth_c)][1]),
               n_surf = .N,
               span_min = (max(as.numeric(time)) - min(as.numeric(time))) / 60,
               light_sd = if (.N > 2) as.numeric(sd(light_c)) else NA_real_),
           by = chunk]
  out <- out[n_surf >= MIN_N]
  out[, chunk := NULL]
  setorder(out, time)

  # depth_max/depth_min must describe the ENCLOSING WINDOW over ALL samples, not
  # the surface chunk. Taking them from the chunk caps depth_max below SURF_M,
  # which silently destroys the dive record -- and `departure_time()` in
  # nes_common.R classifies at-sea days by `depth_max > 50 m`, so every day looked
  # like a haul-out day, every geom_fit came back NULL, and the pooled calibration
  # vanished (0 of 29 instead of 20 of 29). The dive record is load-bearing
  # elsewhere even though this decimation no longer uses it to select light.
  wd <- d[, .(wmax = max(depth_c, na.rm = TRUE),
              wmin = min(depth_c, na.rm = TRUE)), by = .w]
  out[, .w := as.numeric(time) %/% WINDOW_SEC]
  out[, depth_max := wd$wmax[match(.w, wd$.w)]]
  out[, depth_min := wd$wmin[match(.w, wd$.w)]]
  # Per-surfacing emission gives ~93 obs/day against the shipped 48, because
  # haul-out periods chunk into many pieces. That would double the likelihood
  # weight and sharpen the posterior on its own, confounding the statistic change
  # with a count change, and lambda was tuned at 48/day. Keep the best-determined
  # observation per window so the ONLY change is the statistic and its timestamp.
  if (THIN_WINDOW) {
    setorder(out, .w, -n_surf)
    out <- out[, .SD[1L], by = .w]
    setorder(out, time)
  }
  out[, .w := NULL]
  # departure_time() needs deep dives to be visible; assert they survived.
  if (max(out$depth_max, na.rm = TRUE) < 50)
    warning("depth_max never exceeds 50 m: departure_time() will fail", call. = FALSE)
  out <- out[is.finite(light) & !is.na(time)]
  if (nrow(out) < 1000) return(NULL)
  list(main = as.data.frame(out), atten = k, atten_applied = ATTEN_ON,
       zero_drift = as.numeric(diff(quantile(d$z0, c(0.02, 0.98)))),
       surf_m = SURF_M, max_span_min = MAX_SPAN_MIN, n_raw = nrow(d))
}

read_archive <- function(path) {
  hdr <- names(fread(path, nrows = 0L))
  if (!all(LIGHT_NEED %in% hdr)) return(NULL)
  tcol <- intersect(TEMP_NAMES, hdr)
  cols <- c("Time", "Depth", if (length(tcol)) tcol[1] else NULL, "Light Level")
  d <- fread(path, select = cols)
  setnames(d, cols, c("time", "depth", if (length(tcol)) "temp" else NULL, "light"))
  if (!length(tcol)) d[, temp := NA_real_]
  d[, time := as.POSIXct(time, format = "%H:%M:%S %d-%b-%Y", tz = "UTC")]
  if (all(is.na(d$time))) d[, time := as.POSIXct(time, tz = "UTC")]
  d
}

# Sourcing with options(decimate.lib = TRUE) gives the functions and constants
# without running the driver, so the validator can exercise decimate_surface()
# directly instead of duplicating it.
if (isTRUE(getOption("decimate.lib", FALSE))) {
  message("decimate_surface.R loaded as a library")
} else {

args <- commandArgs(trailingOnly = TRUE)
files <- list.files(TDR, pattern = "Archive", full.names = TRUE)
files <- files[endsWith(files, ".csv")]
ids <- substr(basename(files), 1, 7)
if (length(args)) { keep <- ids %in% args; files <- files[keep]; ids <- ids[keep] }
cat(sprintf("%d archive(s)\n\n", length(files)))

res <- list()
for (i in seq_along(files)) {
  raw <- read_archive(files[i])
  if (is.null(raw)) { cat(ids[i], ": no light channel\n"); next }
  r <- decimate_surface(raw)
  if (is.null(r)) { cat(ids[i], ": too few samples\n"); next }
  res[[ids[i]]] <- r
  m <- r$main
  days <- as.numeric(difftime(max(m$time), min(m$time), units = "days"))
  cat(sprintf("%s  raw %7d -> %6d obs (%4.1f/day) | span med %4.1f max %4.1f min | n med %3.0f | atten %.3f | drift %+.1f m\n",
              ids[i], r$n_raw, nrow(m), nrow(m) / days,
              median(m$span_min), max(m$span_min), median(m$n_surf),
              r$atten, r$zero_drift))
}
if (!length(args)) {
  saveRDS(res, OUT)
  cat(sprintf("\nwritten: %s (%d tags)\n", OUT, length(res)))
} else {
  cat("\ntest mode: nothing written\n")
}

}  # end driver guard
