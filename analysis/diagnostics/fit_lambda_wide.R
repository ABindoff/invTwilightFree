# STAGE 1: a 2-D grid search for the source of the latitude bias.
#
# THE HYPOTHESIS BEING TESTED. The bias must come from one of two places:
#   (a) CALIBRATION -- what light level we take to correspond to civil twilight,
#       i.e. where the modelled clear-sky curve sits along the zenith axis; or
#   (b) SHADING     -- given movement, priors and light, how much shading the model
#       must accept to believe it is at a given location.
# Both are tunable against Argos truth, so a grid over the two settles it.
#
# THE TWO KNOBS, and why they are parametrised this way.
#
#   delta (degrees of zenith).  The tangent response is
#         expected(z) = clamp(a - b*z, 0, max_light)
#   and "what light do we call civil twilight" is expected(96). A RIGID SHIFT along
#   the zenith axis,
#         expected(z) = clamp(a - b*(z - delta), 0, max_light) = c(a + b*delta, b)
#   changes exactly that and nothing else: delta maps linearly to the light level at
#   96 deg, and the daytime slope is untouched, so it cannot be confounded with a
#   gain or scale change. Positive delta = light persists to higher zenith =
#   a brighter, later twilight.
#
#   lam_mult (dimensionless).  lambda is the spike's lower-arm rate, i.e. the price
#   of observing light BELOW the expected curve, which is what shading is. Swept as a
#   multiplier on the current 1/(max_light*0.5). Swept in BOTH directions on purpose:
#   the measured residual scale says the default is 3.1x too loose, while a
#   simulation in this campaign said tightening it makes the shading-induced latitude
#   bias worse. Those cannot both be right and the grid will say which.
#   shade_ratio is held at 2 (TwilightFreeGrid's default) so this is a clean 2-D
#   experiment rather than a 3-D one.
#
# WHAT IS MEASURED. Emission only: the light log-likelihood profiled over latitude at
# the TRUE Argos longitude, per knot window. No HMM, no movement prior. This isolates
# the two knobs from the estimator entirely -- which is the point, since the movement
# kernel has already been corrected and should no longer be a suspect.
#
# CONTROLS, all of which must hold or the run is void:
#   1. delta = 0, lam_mult = 1 IS the shipped tangent. Its band biases must reproduce
#      the known values (day approx +0.36, twilight approx -1.25, night approx -1.68).
#      If they do not, the harness differs from fit_rescore29.R and nothing else here
#      can be read.
#   2. The DAY band is nearly uninformative (profile sd approx 12.7 vs 6.4 at
#      twilight). Any config that appears to fix the bias by flattening the profile is
#      buying it with information, not accuracy, so prof_sd is carried everywhere and
#      a widening config is disqualified however good its bias looks.
#   3. eval_logpk_grid's 7th argument is SHADE_RATIO, not a cell count (cells =
#      length(lon)). Passing a count silently sets an ~80x asymmetric spike and
#      produces a spurious poleward bias. It is 2 here.
#
# VALIDATION. Tuning against Argos on the same 29 deployments is in-sample. 29
# deployments come from only 15 distinct loggers (86% reuse; 42% of the bias variance
# is between loggers), so leaving out a DEPLOYMENT would leak through the shared
# sensor. Leave-one-LOGGER-out is the honest version and is reported alongside the
# in-sample optimum. The gap between them is the overfitting.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))
OUT <- "scratch/nes_calibration/bias_lambda_wide.csv"

STEP_HOURS <- 12; CAL_DAYS <- 15
LATS   <- seq(20, 70, by = 0.5)
MAXWIN <- 25
SHADE_RATIO <- 2

DELTAS <- 0                              # settled as noise in the 2-D grid: weak, monotone in the unhelpful direction, and its LOLO choice was unstable across all 15 folds
LAMS   <- 2^seq(-4, 5, by = 0.5)         # 0.0625x to 32x; 10 of 15 loggers optimised OUTSIDE the old 0.25-4x window

# ---- tags: all 29, assembled exactly as fit_rescore29.R does ----------------
old <- lapply(readRDS("analysis/cache/nes/archives_v3_30min.rds"), `[[`, "main")
new <- lapply(readRDS("analysis/cache/nes/fvilches_archives_v1_30min.rds"), `[[`, "main")
new_argos <- readRDS("analysis/cache/nes/fvilches_argos_v1.rds")
new_man <- fread("scratch/nes_calibration/fvilches_manifest.csv",
                 colClasses = list(character = "topp"))
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
  tags[[tg]] <- list(id = tg, light = d, argos = a, p0 = c(df$deploy_lon, df$deploy_lat))
}
for (tg in names(new)) {
  mm <- new_man[topp == tg]; if (!nrow(mm)) next
  a <- as.data.table(new_argos[[tg]]); d <- as.data.table(new[[tg]])
  d <- d[time >= min(a$time) & time <= max(a$time)]; if (nrow(d) < 1000) next
  tags[[tg]] <- list(id = tg, light = d, argos = a,
                     p0 = c(mm$deploy_lon[1], mm$deploy_lat[1]))
}

# ---- serial (logger) per deployment, for clustering and for LOLO -----------
a1 <- data.table(file = list.files("data/nes_untracked", pattern = "^[0-9]+_.*\\.csv$"))
a1[, `:=`(topp = sub("^([0-9]+)_.*$", "\\1", file),
          serial = sub("^[0-9]+_(.+)\\.csv$", "\\1", file))]
b1 <- data.table(file = list.files("fvilches/extracted/TDR raw", pattern = "Archive[.]csv$"))
b1[, `:=`(topp = sub("^([0-9]+)_.*$", "\\1", file),
          serial = sub("^[0-9]+_(.+)-out-Archive[.]csv$", "\\1", file))]
SER <- unique(rbind(a1[, .(topp, serial)], b1[, .(topp, serial)]), by = c("topp", "serial"))
family <- function(s) fifelse(is.na(s), "unknown",
                       fifelse(grepl("^219", s), "Mk9_219",
                        fifelse(grepl("^18A", s), "F18A", "other")))
for (tg in names(tags)) tags[[tg]]$family <- family(SER[topp == tg]$serial[1])

# ---- responses: pooled WITHIN FAMILY, exactly as fit_rescore29.R ------------
# Pooling within family is what took the campaign 473 -> 208 km. Pooling across all
# 29 gives a different tangent, and then control 1 cannot pass.
scale_fits <- lapply(tags, function(g) {
  k <- g$light$time <= min(g$light$time) + CAL_DAYS * 86400
  fit_light_response(g$light$time[k], g$light$light[k], g$p0[1], g$p0[2])
})
geom_fits <- lapply(tags, function(g) {
  dep <- departure_time(g$light)
  k <- is.finite(dep) & as.numeric(g$light$time) < dep
  if (sum(k) < 200) return(NULL)
  fit_light_response(g$light$time[k], g$light$light[k], g$p0[1], g$p0[2])
})
responses <- list()
for (fam in unique(vapply(tags, function(g) g$family, ""))) {
  ids <- names(tags)[vapply(tags, function(g) g$family, "") == fam]
  gf <- geom_fits[ids]; sf <- scale_fits[ids]
  if (!sum(!vapply(gf, is.null, logical(1)))) next
  pooled <- pool_light_responses(gf, sf)
  for (i in ids) responses[[i]] <- pooled[[i]]
}
cat(sprintf("%d deployments, %d loggers, %d responses\n",
            length(tags), uniqueN(SER[topp %in% names(tags)]$serial), length(responses)))
cat(sprintf("grid: %d deltas x %d lambda multipliers = %d configs\n\n",
            length(DELTAS), length(LAMS), length(DELTAS) * length(LAMS)))

profile_stats <- function(ll, truth) {
  ll <- ll - max(ll); w <- exp(ll); s <- sum(w)
  if (!is.finite(s) || s <= 0) return(c(NA_real_, NA_real_))
  w <- w / s; mu <- sum(w * LATS)
  c(mu - truth, sqrt(sum(w * (LATS - mu)^2)))
}

R <- list()
for (tg in names(tags)) {
  g <- tags[[tg]]; r <- responses[[tg]]; if (is.null(r)) next
  a <- as.data.table(g$argos); setorder(a, time)
  tt  <- as.numeric(g$light$time)
  obs <- pmax(0, g$light$light - r$baseline)
  brk <- seq(min(tt), max(tt), by = STEP_HOURS * 3600)
  kk   <- findInterval(tt, brk)
  tlon <- approx(as.numeric(a$time), a$lon, brk, rule = 2)$y
  tlat <- approx(as.numeric(a$time), a$lat, brk, rule = 2)$y
  lam0 <- 1 / (r$max_light * 0.5)
  use <- which(tabulate(kk, length(brk)) >= 16)
  if (!length(use)) next
  use <- use[seq(1, length(use), length.out = min(MAXWIN, length(use)))]
  message(tg, " (", length(use), " windows)")
  for (i in use) {
    j <- which(kk == i); if (length(j) < 16) next
    z <- solar_zenith(tt[j], rep(tlon[i], length(j)), rep(tlat[i], length(j)))
    bands <- list(all = rep(TRUE, length(j)), day = z < 85,
                  twilight = z >= 85 & z <= 105, night = z > 105)
    for (dl in DELTAS) {
      cal <- c(r$calibration[1] + r$calibration[2] * dl, r$calibration[2])
      for (lm in LAMS) {
        lp <- c(lam0 * lm, r$max_light, 0.10)
        for (bn in names(bands)) {
          sel <- bands[[bn]]; if (sum(sel) < 6) next
          ll <- eval_logpk_grid(rep(tlon[i], length(LATS)), LATS, tt[j][sel],
                                obs[j][sel], cal, lp, SHADE_RATIO)
          st <- profile_stats(ll, tlat[i])
          R[[length(R) + 1]] <- data.table(id = tg, delta = dl, lam_mult = lm,
                                           band = bn, bias = st[1], sd = st[2])
        }
      }
    }
  }
}
D <- rbindlist(R)
D <- merge(D, SER, by.x = "id", by.y = "topp", all.x = TRUE)
fwrite(D, OUT)
cat(sprintf("\nwrote %d rows to %s\n", nrow(D), OUT))
