# PER-AXIS MOVEMENT SCALE, AND HOW IT GROWS WITH LAG.
#
# The grid HMM's movement prior is isotropic Brownian with one scale, `diffusion`,
# in km/sqrt(day). Two facts already measured make that suspect:
#   - the true 12 h step is 17.2 km north-south and 22.4 east-west, against a prior
#     step of 77.8 km on BOTH axes: the latitude prior is 4.5x too wide;
#   - fitted bias scales with prior width (perfect arm: -0.567 at D=110, -1.670 at
#     D=440), because bias ~ tilt x sum-of-covariances and the cell-area tilt is the
#     thing being amplified.
#
# But tightening is not free, and the reason is that these animals are DIRECTED. For
# Brownian motion sd(displacement over lag t) = D sqrt(t), so D(lag) is FLAT. For a
# directed animal it GROWS with lag. A memoryless prior has to pick one number, and
# the one-step value gives honest steps but far too narrow marginal intervals
# (documented: D=44 gives latitude coverage 0.21), while the displacement-equivalent
# value gives honest intervals at a movement scale no seal ever swam.
#
# So measure D(lag) SEPARATELY PER AXIS. The interesting possibility is that the two
# axes are directed to different degrees -- the step ratio favours east-west (1.34)
# but the total span favours north-south (~2130 km vs ~1660) -- in which case the
# right anisotropy for honest INTERVALS may point the opposite way from the step
# ratio, and the current isotropic prior is wrong on both axes in opposite
# directions.
#
# Argos truth from the battery's offset-0 records: already interpolated and clipped
# per deployment, so no PTT-reuse contamination. Aggregates only.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))
suppressMessages(devtools::load_all(".", quiet = TRUE))
BAT <- readRDS("scratch/nes_calibration/light_battery.rds")
REC <- Filter(function(z) z$offset == 0, BAT)

LAGS <- c(0.5, 1, 2, 5, 10, 20, 40)   # days

scales <- rbindlist(lapply(REC, function(b) {
  tt <- as.numeric(b$time)
  g  <- seq(min(tt), max(tt), by = 6 * 3600)          # 6 h lattice
  la <- approx(tt, b$lat, g)$y
  lo <- approx(tt, lon360(b$lon), g)$y
  rbindlist(lapply(LAGS, function(L) {
    step <- round(L * 24 / 6)
    if (step >= length(g) - 5) return(NULL)
    i <- seq_len(length(g) - step)
    dlat <- (la[i + step] - la[i]) * 111.32
    dl   <- lo[i + step] - lo[i]
    dl[dl >  180] <- dl[dl >  180] - 360
    dl[dl < -180] <- dl[dl < -180] + 360
    dlon <- dl * 111.32 * cos(la[i] * pi / 180)
    data.table(id = b$id, lag = L,
               D_ns = sd(dlat) / sqrt(L), D_ew = sd(dlon) / sqrt(L))
  }))
}))

cat("=== per-axis diffusion scale D(lag) = sd(displacement) / sqrt(lag), km/sqrt(day) ===\n")
cat("    flat in lag => Brownian.  growing => directed.\n\n")
S <- scales[, .(D_ns = round(median(D_ns), 1), D_ew = round(median(D_ew), 1),
                ratio_ew_ns = round(median(D_ew) / median(D_ns), 2)), by = lag][order(lag)]
print(as.data.frame(S), row.names = FALSE)

cat(sprintf("\n  growth factor 0.5 d -> 20 d :  NS %.2fx   EW %.2fx\n",
            S[lag == 20]$D_ns / S[lag == 0.5]$D_ns,
            S[lag == 20]$D_ew / S[lag == 0.5]$D_ew))
cat("  a larger growth factor means MORE directed on that axis, so a memoryless\n")
cat("  prior must be inflated further above the step scale to cover it honestly.\n")

cat(sprintf("\n  shipped prior: D = 110 on BOTH axes.\n"))
for (L in c(0.5, 5, 20)) {
  r <- S[lag == L]
  cat(sprintf("    vs lag %4.1f d : NS %.1f (%.2fx too wide)   EW %.1f (%.2fx)\n",
              L, r$D_ns, 110 / r$D_ns, r$D_ew, 110 / r$D_ew))
}

cat("\n=== candidate settings for the fit test ===\n")
d1 <- S[lag == 0.5]; d5 <- S[lag == 5]
cat(sprintf("  one-step (0.5 d)      : NS %.0f  EW %.0f\n", d1$D_ns, d1$D_ew))
cat(sprintf("  mid-lag  (5 d)        : NS %.0f  EW %.0f\n", d5$D_ns, d5$D_ew))
cat(sprintf("  shipped               : NS 110 EW 110\n"))
cat("\n  The mid-lag value is the compromise a memoryless model is implicitly making:\n")
cat("  big enough to cover multi-day displacement, small enough to be a real animal.\n")
