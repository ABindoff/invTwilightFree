# DOES EACH EMITTED (time, light) PAIR TELL THE TRUTH?
#
# The v1 decimator failed because light and timestamp described different moments.
# That is checkable against the raw record without any simulation and without any
# fit: for every emitted observation, go back to the 4 s archive, take the SURFACE
# samples lying within +/- TOL seconds of the emitted TIMESTAMP, and ask whether
# their light matches the emitted LIGHT.
#
# An honest pair agrees. A pair whose light came from a different moment does not.
# Discrepancy is reported in light units and in zenith-equivalent degrees at
# 3.74 units per degree, so it is directly comparable to the 1.1 deg of max bias
# this work exists to remove.
#
# The v1 estimator is reconstructed alongside as a POSITIVE CONTROL. If the check
# cannot detect the bug that is known to be there, the check is worthless.
#
# Usage: Rscript analysis/diagnostics/validate_decimator.R [tag]
suppressMessages({ library(data.table) })
Sys.setlocale("LC_TIME", "C")
options(decimate.lib = TRUE)
source("analysis/diagnostics/decimate_surface.R")

UNITS_PER_DEG <- 3.74
TOL <- 45          # seconds either side of the emitted timestamp

args <- commandArgs(trailingOnly = TRUE)
tag <- if (length(args)) args[1] else "2021023"
f <- list.files(TDR, pattern = paste0("^", tag), full.names = TRUE)
f <- f[endsWith(f, ".csv")][1]
stopifnot(!is.na(f))
cat("tag", tag, "\n\n")

raw <- read_archive(f)
raw <- raw[is.finite(depth) & is.finite(light) & !is.na(time)]
setorder(raw, time)
raw[, z0 := depth_zero(time, depth)][, depth_c := depth - z0]
surf <- raw[depth_c < SURF_M, .(ts = as.numeric(time), light = as.numeric(light))]
setorder(surf, ts)

# FIRST VERSION OF THIS CHECK WAS BROKEN, and in the same way as everything else
# today: it compared the emitted light to surface samples within +/- TOL of the
# emitted TIME, and returned NA when there were none. But "no surface samples at
# the emitted time" is not missing data, it is THE FAILURE -- the decimator has
# stamped an observation at a moment the animal was underwater. Dropping those
# rows made the known-broken v1 estimator look clean (63% checkable, tiny
# discrepancy) because the 37% it dropped were precisely its bad pairs.
#
# So now: truth is surface light INTERPOLATED across the gap, which is defined
# everywhere, and the distance to the nearest real surface sample is reported
# separately as the orphan rate. An honest decimator has an orphan rate near zero
# because its timestamps are real surface times.
# TWO metrics, because they answer different questions and one alone is blind.
#
#  ORPHAN RATE, from `near_s`: does surface truth exist at the emitted time? This
#  is the only thing that catches v1, whose timestamps fall between bouts.
#
#  DISCREPANCY vs the NEIGHBOURHOOD MEDIAN, from `truth_med`: is the emitted value
#  a representative surface reading or an extreme one? This is the only thing that
#  catches the shipped max. Comparing against light INTERPOLATED AT THE EMITTED
#  INSTANT cannot catch it, because the shipped timestamp IS the argmax sample's
#  own time, so interpolation returns that sample and the comparison is circular.
#  That circularity is why an earlier version of this script scored the shipped
#  decimator at +0.049 deg when the neighbourhood comparison gives +1.475.
#
# Orphans are never dropped: the neighbourhood is widened to WIDE_TOL for them, so
# a broken pairing is measured rather than excluded.
WIDE_TOL <- 600
truth_med <- function(t) {
  band <- function(t, tol) {
    lo <- findInterval(t - tol, surf$ts) + 1L
    hi <- findInterval(t + tol, surf$ts)
    vapply(seq_along(t), function(i)
      if (hi[i] >= lo[i]) median(surf$light[lo[i]:hi[i]]) else NA_real_, 0)
  }
  m <- band(t, TOL)
  bad <- !is.finite(m)
  if (any(bad)) m[bad] <- band(t[bad], WIDE_TOL)
  m
}
near_s <- function(t) {
  i <- findInterval(t, surf$ts)
  before <- ifelse(i >= 1L, t - surf$ts[pmax(1L, i)], Inf)
  after <- ifelse(i < nrow(surf), surf$ts[pmin(nrow(surf), i + 1L)] - t, Inf)
  pmin(before, after)
}

report <- function(lab, tm, li) {
  t <- as.numeric(tm)
  keep <- is.finite(t) & is.finite(li)
  t <- t[keep]; li <- li[keep]
  gap <- near_s(t)
  orph <- gap > TOL
  d <- li - truth_med(t)
  ok <- is.finite(d)
  d <- d[ok]; orph2 <- orph[ok]
  f <- function(x) if (!length(x)) "        n/a          " else
    sprintf("bias %+6.2f   90th |d| %5.2f", mean(x), quantile(abs(x), 0.90))
  cat(sprintf("%-28s n=%5d\n", lab, length(t)))
  cat(sprintf("   [pairing]   orphan timestamps %5.1f%%  | gap to nearest surface sample: median %4.0f s, 90th pct %5.0f s\n",
              100 * mean(orph), median(gap), quantile(gap, 0.90)))
  cat(sprintf("   [statistic] vs neighbourhood median, light units: %s\n", f(d)))
  cat(sprintf("               orphans only                       : %s\n", f(d[orph2])))
  cat(sprintf("   [statistic] ZENITH-EQUIVALENT: bias %+6.3f deg | 90th pct |d| %5.3f deg\n\n",
              mean(d) / UNITS_PER_DEG, quantile(abs(d), 0.90) / UNITS_PER_DEG))
  invisible(d)
}

# ---- v2, as shipped in decimate_surface.R -----------------------------------
r <- decimate_surface(copy(raw[, .(time, depth, light, temp)]))
m <- as.data.table(r$main)
cat("=== v2: one observation per surfacing, median light at median time ===\n")
report("v2 (current)", m$time, m$light)

# ---- v1, reconstructed: median light, MEAN time, one per 30-min window ------
v1 <- raw[, .b := as.numeric(time) %/% 1800][
  , { s <- depth_c < SURF_M
      if (!any(s)) .(time = as.POSIXct(NA), light = NA_real_)
      else .(time = as.POSIXct(mean(as.numeric(time[s])), origin = "1970-01-01", tz = "UTC"),
             light = as.numeric(median(light[s]))) }, by = .b][is.finite(light)]
cat("=== v1: POSITIVE CONTROL, the known-broken pairing ===\n")
report("v1 (median light, mean time)", v1$time, v1$light)

# ---- the shipped decimator: max over the whole window, argmax timestamp -----
sh <- raw[, { j <- which.max(light)
              .(time = time[j], light = as.numeric(light[j])) }, by = .b][is.finite(light)]
cat("=== shipped: max over window at the argmax sample's own time ===\n")
cat("   (its pairing is honest by construction; the discrepancy here IS the\n")
cat("    upward order-statistic bias, which is the thing being removed)\n")
report("shipped (max)", sh$time, sh$light)

cat("HOW TO READ IT\n")
cat("  ORPHAN RATE is the primary discriminator, not the discrepancy. An honest\n")
cat("  decimator stamps observations at real surface times, so its orphan rate is\n")
cat("  ~0. v1's rate is the bug made visible; the first version of this script\n")
cat("  hid it by treating orphans as missing data.\n")
cat("  The shipped max should show a positive BIAS with a low orphan rate: its\n")
cat("  pairing is honest, its statistic is not.\n")
