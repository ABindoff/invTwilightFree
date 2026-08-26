# DOES THE DELIVERED RESPONSE TABLE BEAT THE PANEL'S CLAMPED LINEAR RAMP?
#
# Tested WITHOUT running the HMM. `eval_logpk_grid()` is exported and takes the
# calibration directly, returning the summed per-cell log-likelihood, so the
# emission model can be compared on its own: over a window where the animal is
# roughly stationary, does the likelihood peak nearer the truth under the table
# than under the ramp?
#
# THE LEVEL PROBLEM, which is not optional. The table is in the units of a
# per-minute median of surface-conditioned samples: plateau 172, floor 46 for
# Mk9. The engine compares against pmax(0, light - baseline) on the 30-minute
# max-decimated panel series, and clamps expected light to [0, max_light]. Fed raw,
# the table would predict 46 units of light at night where the panel reads near
# zero. Measured on tag 2021023 the max sits +5.32 units above a surface median
# overall, but band-dependent (about 6.0 in daylight, 4.0-4.5 at twilight), so a
# fitted affine leaves a structured residual worst at twilight.
#
# So the map is ANCHORED, not fitted. The anchors matter and I got them wrong
# first time, which cost 2641 km of peak error against the ramp's 245, consistently
# ~22 deg too far south on every tag:
#
#   WRONG: their floor -> r$floor, their plateau -> r$floor + r$amp.
#     r$floor (74.2) and r$amp (144.9) are RAW-ENVELOPE parameters -- night level
#     74.2 and daytime peak 219.1 on the raw scale. Mapping onto them predicted 81
#     units of light at night, where the baseline-subtracted observations read 0-3,
#     and the fit fled south to accommodate the excess.
#
#   RIGHT: their base -> 0, their plateau (base + amp) -> max_light.
#     The engine's ramp runs 0 to max_light and is applied to
#     pmax(0, light - baseline); subtracting `baseline` (a q05) has ALREADY removed
#     the floor, so the target series bottoms out at zero by construction. Their
#     table has a floor at 0.27 of plateau precisely because they did not subtract
#     one. Verified against the observed series: baseline-subtracted median light
#     runs 144 at zenith 20 down to 0-3 beyond 110 deg.
suppressMessages({ library(data.table) })
source("analysis/diagnostics/nes_common.R")
suppressMessages(devtools::load_all(".", quiet = TRUE))

GEO <- "../wavscat/scratch/geo/out"
source(file.path(GEO, "calibration_tables.R"))
PAR <- fread(file.path(GEO, "response_params.csv"), colClasses = list(character = "topp"))
P <- readRDS("analysis/cache/nes/panel_v1.rds")

TBL <- list(Mk9_219 = calib_Mk9_2190, F18A = calib_18A_19A)
cat("=== table sanity ===\n")
for (nm in names(TBL)) {
  v <- TBL[[nm]]
  ny <- length(v) - 3L
  cat(sprintf("  %-8s length %d (3 + %d bins) | z %.2f to %.2f by %.2f | y range %.0f to %.0f | monotone-decreasing fraction %.2f\n",
              nm, length(v), ny, v[2], v[2] + (ny - 1) * v[3], v[3],
              min(v[-(1:3)]), max(v[-(1:3)]), mean(diff(v[-(1:3)]) <= 0)))
}

# anchored affine: their (base, base+amp) -> (0, max_light)
map_table <- function(v, their_base, their_amp, max_light) {
  b <- max_light / their_amp
  a <- -b * their_base
  c(v[1:3], a + b * v[-(1:3)])
}

score_window <- function(tg, days = 7) {
  g <- P$tags[[tg]]; r <- P$responses[[tg]]
  fam <- g$family
  if (is.null(TBL[[fam]]) || is.null(r)) return(NULL)
  pr <- PAR[topp == tg]
  if (!nrow(pr)) return(NULL)
  L <- as.data.table(g$light); A <- as.data.table(g$argos)[is.finite(lon) & is.finite(lat)]
  t0 <- min(as.numeric(L$time)) + 90 * 86400          # mid-deployment, away from the colony
  W <- L[as.numeric(time) >= t0 & as.numeric(time) < t0 + days * 86400]
  Aw <- A[as.numeric(time) >= t0 & as.numeric(time) < t0 + days * 86400]
  if (nrow(W) < 100 || nrow(Aw) < 5) return(NULL)
  true_lon <- lon360(median(Aw$lon)); true_lat <- median(Aw$lat)
  gr <- expand.grid(lon = seq(150.5, 249.5, by = 1), lat = seq(20.5, 69.5, by = 1))
  obs <- pmax(0, W$light - r$baseline); tt <- as.numeric(W$time)
  lp <- c(1.0 / (r$max_light * 0.5), r$max_light, 0.10)
  tab <- map_table(TBL[[fam]], pr$base[1], pr$amp[1], r$max_light)
  out <- list()
  for (nm in c("ramp", "table")) {
    cal <- if (nm == "ramp") r$calibration else tab
    ll <- eval_logpk_grid(gr$lon, gr$lat, tt, obs, cal, lp, 2)
    i <- which.max(ll)
    out[[nm]] <- data.table(
      id = tg, family = fam, cal = nm,
      peak_km = gc_km(gr$lon[i], gr$lat[i], true_lon, true_lat),
      d_lat = gr$lat[i] - true_lat,
      finite = mean(is.finite(ll)),
      ll_at_truth = ll[which.min(gc_km(gr$lon, gr$lat, true_lon, true_lat))],
      ll_at_peak = ll[i])
  }
  rbindlist(out)
}

sel <- names(P$tags)[vapply(P$tags, function(g) g$family, "") %in% names(TBL)]
sel <- c(head(sel[vapply(P$tags[sel], function(g) g$family, "") == "Mk9_219"], 6),
         head(sel[vapply(P$tags[sel], function(g) g$family, "") == "F18A"], 2))
R <- rbindlist(lapply(sel, score_window))
cat("\n=== emission-only peak location, 7-day windows at day 90 ===\n")
print(as.data.frame(R[, .(id, family, cal, peak_km = round(peak_km),
                          d_lat = round(d_lat, 1), finite = round(finite, 3),
                          gap_ll = round(ll_at_peak - ll_at_truth, 1))]), row.names = FALSE)
cat("\n=== summary ===\n")
s <- R[, .(n = .N, med_peak_km = round(median(peak_km)),
           med_gap = round(median(ll_at_peak - ll_at_truth), 1)), by = cal]
print(as.data.frame(s), row.names = FALSE)
cat("\n  peak_km: distance from the likelihood peak to the Argos truth (lower better)\n")
cat("  gap_ll : log-likelihood at the peak minus at the truth (lower better; 0 = truth is the MLE)\n")
