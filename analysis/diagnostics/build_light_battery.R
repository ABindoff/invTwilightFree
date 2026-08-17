# A battery of light records with KNOWN truth, built from REAL movement and REAL
# shading.
#
# WHY, and what went wrong before. The simulation run earlier today was a wash: the
# tracks were free random walks with an invented shading process, and the resulting
# fits were two to three times worse than anything real (|bias| 3.8 vs 1.5, coverage
# 0.24 vs 0.66). A schedule that redistributes a factor of 3.7 in lambda has nothing
# to work with when the fit is that badly determined, so the comparison could not
# have detected a real effect. The regime, not the method, defeated it.
#
# THE DESIGN HERE keeps everything real except the one thing being controlled:
#   * MOVEMENT is real -- true positions come from Argos, interpolated to knots. No
#     random walk, no invented step-length distribution.
#   * SHADING is real -- the per-observation attenuation is measured from that tag's
#     own light record against its own clear-sky expectation, and re-applied. Dive
#     bouts, cloud, haul-outs and whatever else is in there comes along with it.
#   * TIME OF YEAR is the only thing varied, by shifting the DATE in whole days.
#
# Whole days matters. Shifting by a whole number of days preserves local time of day
# exactly, so the diel structure of the shading stays aligned with the sun -- a seal
# that dives at dawn still dives at dawn. Shifting by a fraction of a day would slide
# the animal's behaviour against the solar cycle and manufacture an artefact.
#
# WHAT IS KNOWN EXACTLY: the true position at every knot, the clear-sky light that
# position implies, and the attenuation applied to it. That is the point of the
# exercise -- the bias can be chased against a known answer.
#
# HONEST LIMITS, stated up front:
#  * The date shift moves the animal's real track to a season it never experienced.
#    The MOVEMENT is real but the combination of movement and season is synthetic;
#    a real seal in December is not where this track puts it. That is deliberate --
#    it is what breaks the phenology confound -- but it is not a real animal.
#  * The shading is transplanted onto a different solar geometry, so an attenuation
#    that occurred at (say) local noon still occurs at local noon, but the sun is at
#    a different elevation. Diel alignment is preserved; solar-elevation alignment
#    is not, and cannot be for a date shift.
#  * Argos itself has error (tens of km) and is sampled irregularly, so "truth" is
#    an interpolation of a noisy reference, not a perfect one.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))
source("scratch/nes_calibration/argos_filter.R")
suppressMessages(devtools::load_all(".", quiet = TRUE))
OUT <- "scratch/nes_calibration/light_battery.rds"

CAL_DAYS  <- 15
OFFSETS   <- c(0, 90, 180, 270)      # whole days: preserves local time of day
N_TRACK   <- 6
MAX_TRUTH_GAP_H <- 12   # do not claim a true position further than this from a fix
SKY_MIN   <- 0.05                    # fraction of max_light below which "sky" is nil
S_CAP     <- 3
GEN_PSLAB <- 0.10        # slab weight used to GENERATE the noisy level                       # cap on the attenuation factor (see below)

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
  d <- x[time >= t0 & time <= t1]; if (nrow(d) < 2000) next
  tags[[tg]] <- list(id = tg, light = d, argos = a, p0 = c(df$deploy_lon, df$deploy_lat))
}
for (tg in names(new)) {
  mm <- new_man[topp == tg]; if (!nrow(mm)) next
  a <- as.data.table(new_argos[[tg]]); d <- as.data.table(new[[tg]])
  d <- d[time >= min(a$time) & time <= max(a$time)]; if (nrow(d) < 2000) next
  tags[[tg]] <- list(id = tg, light = d, argos = a,
                     p0 = c(mm$deploy_lon[1], mm$deploy_lat[1]))
}

# Pick tracks by TRUTH SUPPORT first, then latitude range.
#
# The first version scored latitude range x fix count, which rewards how MANY fixes a
# tag has and says nothing about WHEN they are. It selected 2023041, whose 406 fixes
# leave 34-day gaps: after honest filtering only 58% of its samples have a fix within
# 12 hours, and interpolating the rest produced a 1800 km triangular excursion that
# was visible by eye. A battery track is only useful where the truth is known, so
# supported fraction is the first criterion and latitude range breaks ties.
support_of <- function(g) {
  a <- as.data.table(argos_clean2(as.data.table(g$argos), vmax = 3, max_detour = 10,
                                  min_km = 15, max_isolation_h = 48))
  if (nrow(a) < 50) return(c(supp = 0, latr = 0))
  at <- as.numeric(a$time); tt <- as.numeric(g$light$time)
  c(supp = mean(vapply(tt, function(x) min(abs(at - x)), 0) / 3600 <= MAX_TRUTH_GAP_H),
    latr = diff(range(a$lat)))
}
sup <- vapply(tags, support_of, c(supp = 0, latr = 0))
cat("truth support by deployment (fraction of samples within 12 h of a fix):
")
o <- order(sup["supp", ], decreasing = TRUE)
for (i in o) cat(sprintf("  %-9s supported %.0f%%  lat range %.1f deg
",
                         names(tags)[i], 100 * sup["supp", i], sup["latr", i]))
elig <- names(tags)[sup["supp", ] >= 0.90]
picked <- elig[order(sup["latr", elig], decreasing = TRUE)][seq_len(min(N_TRACK, length(elig)))]
cat(sprintf("
%d deployments available, %d with >=90%% support; using %s

",
            length(tags), length(elig), paste(picked, collapse = ", ")))

resp <- lapply(tags[picked], function(g) {
  k <- g$light$time <= min(g$light$time) + CAL_DAYS * 86400
  fit_light_response(g$light$time[k], g$light$light[k], g$p0[1], g$p0[2])
})

BAT <- list()
for (tg in picked) {
  g <- tags[[tg]]; r <- resp[[tg]]
  a <- as.data.table(g$argos); setorder(a, time)
  # SPEED-FILTER the reference first. The 2021 delivery is raw: 15-25% of its steps
  # imply speeds above 5 m/s and one reaches 159,722 m/s, while the fvilches delivery
  # arrived already filtered at 3 m/s. Interpolating an unfiltered series to the
  # 30-minute observation times puts the "true" position somewhere the animal never
  # was, and the clear-sky light computed there is then not clear-sky light for this
  # animal at all.
  #
  # Note this matters MORE here than for the campaign's own error metrics. Those are
  # scored at 12-hourly knots, where interpolation steps over the brief spikes: the
  # headline moved by 0 km and 0.000 deg when rescored on filtered truth. At 30-minute
  # resolution the spikes are not stepped over -- truth shifts by 8-26 km at q95 and
  # up to 3634 km.
  # Two stages. The speed filter alone is not enough on a sparse track: it permits
  # vmax TIMES THE GAP, so at 3 m/s a 13-hour gap allows 140 km and a fix displaced
  # 65 km passes cleanly, leaving a triangular out-and-back detour. The detour filter
  # catches those with a scale-free ratio and removes only 70 fixes across all 29
  # deployments, against 8479 for the speed stage.
  af <- argos_clean2(a, vmax = 3, max_detour = 10, min_km = 15, max_isolation_h = 48)
  cat(sprintf("%-9s Argos %d -> %d (speed %d, detour %d, isolated %d)",
              tg, nrow(a), nrow(af), attr(af, "n_speed"), attr(af, "n_detour"),
              attr(af, "n_isolated")))
  a <- as.data.table(af)
  tt <- as.numeric(g$light$time)
  # TRUTH: filtered Argos positions interpolated to the observation times
  tlon <- approx(as.numeric(a$time), lon360(a$lon), tt, rule = 2)$y
  tlat <- approx(as.numeric(a$time), a$lat, tt, rule = 2)$y

  # ...but only where a fix is actually near. Removing an unverifiable fix does not
  # make the truth known across the gap it sat in -- it only stops pretending. Argos
  # gaps here reach 34 days even after filtering, and a straight line drawn across
  # one is a guess, not a reference. A battery whose whole purpose is exact truth
  # must say where it does not have any.
  at <- as.numeric(a$time)
  near <- vapply(tt, function(x) min(abs(at - x)), 0) / 3600
  supported <- near <= MAX_TRUTH_GAP_H
  cat(sprintf(" | truth supported on %.0f%% of samples
", 100 * mean(supported)))

  # clear-sky expectation at the TRUE position on the REAL dates
  ml   <- r$max_light
  flo  <- r$calibration_logistic[1]
  expc <- function(times, lon, lat) {
    z <- solar_zenith(times, lon, lat)
    pmin(pmax(r$calibration_logistic[1] +
                r$calibration_logistic[2] /
                (1 + exp((z - r$calibration_logistic[3]) / r$calibration_logistic[4])),
              0), ml)
  }
  perfect0 <- expc(tt, tlon, tlat)
  sky0 <- pmax(perfect0 - flo, 0)
  obs_sky <- pmax(as.numeric(g$light$light) - r$baseline - flo, 0)

  # MEASURED attenuation. Defined only where the sky actually contributes: where the
  # clear-sky curve is essentially at the floor there is no sky to attenuate, and a
  # ratio there would be dividing noise by noise. Those samples keep the tag's REAL
  # reading instead, which carries the true night behaviour (dark current, and any
  # artificial light or moonlight that is in there).
  lit <- sky0 > SKY_MIN * ml
  s <- rep(NA_real_, length(tt))
  s[lit] <- pmin(obs_sky[lit] / sky0[lit], S_CAP)

  # Attenuation is a property of the ANIMAL's behaviour, not of the sun, so it has
  # to be defined at every sample -- including ones that were dark on the original
  # dates but are lit after a date shift. Fill the gaps with a local median over a
  # +/- 1 day window of lit samples, which follows the dive schedule rather than
  # imposing a global constant.
  fill_local <- function(v, halfwidth = 48) {
    out <- v; na <- which(is.na(v))
    for (i in na) {
      w <- max(1, i - halfwidth):min(length(v), i + halfwidth)
      m <- stats::median(v[w], na.rm = TRUE)
      out[i] <- if (is.finite(m)) m else stats::median(v, na.rm = TRUE)
    }
    out
  }
  s_fill <- fill_local(s)

  # Night contamination -- dark current above the fitted floor, plus whatever real
  # artificial light or moonlight is in the record. Kept as an ADDITIVE term so it
  # can never replace a daylight value, which is what the first version did: it
  # classified night from the ORIGINAL geometry, so a date shift transplanted real
  # night readings into what had become daylight and collapsed the medians.
  alan <- ifelse(lit, 0, pmax(obs_sky, 0))
  cat(sprintf("%-9s  lit %.0f%% of samples | attenuation: median %.2f, q10 %.2f, q90 %.2f | capped %.2f%%\n",
              tg, 100 * mean(lit), median(s[lit]), quantile(s[lit], .1),
              quantile(s[lit], .9), 100 * mean(obs_sky[lit] / sky0[lit] > S_CAP)))

  for (off in OFFSETS) {
    tt2 <- tt + off * 86400            # whole days: local time of day preserved
    perfect <- expc(tt2, tlon, tlat)
    sky <- pmax(perfect - flo, 0)
    # clear sky for the NEW geometry, attenuated by the REAL measured factor, plus
    # the real night contamination as an additive term
    degraded <- flo + sky * s_fill + alan
    degraded <- pmin(pmax(degraded, 0), ml)

    # THIRD LEVEL: clear sky plus the observation noise the model itself assumes,
    # but NO attenuation. Without this rung, perfect-vs-shaded confounds two things
    # at once -- the observation model and the shading -- and a bias appearing
    # between them could not be attributed to either.
    #
    # Generated from the engine's own spike-and-slab at lambda = 1/(max_light/2),
    # which is the package default, so fitting this level at lambda multiplier 1 is a
    # second null: the model is then exactly correct and any bias is the estimator's.
    lam_gen <- 1 / (ml * 0.5)
    m_lo <- 1 - exp(-lam_gen * perfect)
    m_hi <- 0.5 * (1 - exp(-2 * lam_gen * (ml - perfect)))
    u1 <- runif(length(perfect)); u2 <- runif(length(perfect)); u3 <- runif(length(perfect))
    lower <- u1 < m_lo / (m_lo + m_hi)
    noisy <- ifelse(lower,
                    perfect + log(1 - u2 * m_lo) / lam_gen,
                    perfect - log(1 - u2 * (2 * m_hi)) / (2 * lam_gen))
    noisy <- ifelse(u3 < GEN_PSLAB, runif(length(perfect), 0, ml), noisy)
    noisy <- pmin(pmax(noisy, 0), ml)
    BAT[[length(BAT) + 1]] <- list(
      id = tg, offset = off,
      time = as.POSIXct(tt2, origin = "1970-01-01", tz = "UTC"),
      lon = tlon, lat = tlat,
      perfect = pmin(pmax(perfect, 0), ml),
      noisy = noisy,
      degraded = degraded,
      atten = s, supported = supported, max_light = ml, floor = flo,
      response = r$calibration_logistic, baseline = r$baseline)
  }
}
# ---- ACCEPTANCE CONTROL -----------------------------------------------------
# At offset 0 the construction must reproduce the tag's REAL record: the same
# geometry, the same measured attenuation, the same night term. If it does not, the
# degradation model is not a description of these data and nothing built on it can
# be trusted. This is the check whose absence let the first version ship a bug.
cat("
=== ACCEPTANCE CONTROL: offset 0 must reproduce the real record ===
")
worst <- 0
for (tg in picked) {
  b <- BAT[[which(vapply(BAT, function(z) z$id == tg && z$offset == 0, TRUE))[1]]]
  real <- pmin(pmax(as.numeric(tags[[tg]]$light$light) - resp[[tg]]$baseline, 0),
               resp[[tg]]$max_light)
  d <- b$degraded - real
  rel <- stats::median(abs(d)) / resp[[tg]]$max_light
  worst <- max(worst, rel)
  cat(sprintf("  %-9s median |degraded - real| = %.2f units (%.2f%% of max_light) | cor %.4f
",
              tg, stats::median(abs(d)), 100 * rel, stats::cor(b$degraded, real)))
}
cat(sprintf("  verdict: %s
", if (worst < 0.02) "REPRODUCES -- battery is valid"
            else "*** DOES NOT REPRODUCE: the degradation model does not describe these data"))

saveRDS(BAT, OUT)
cat(sprintf("\n%d records written to %s (%d tracks x %d date offsets)\n",
            length(BAT), OUT, length(picked), length(OFFSETS)))
d <- rbindlist(lapply(BAT, function(b) data.table(id = b$id, offset = b$offset,
  n = length(b$time), start = min(b$time), end = max(b$time),
  lat_min = min(b$lat), lat_max = max(b$lat),
  dec_min = min(abs(solar_declination(b$time))),
  dec_max = max(abs(solar_declination(b$time))),
  perfect_max = max(b$perfect), degraded_med = median(b$degraded))))
print(as.data.frame(d), row.names = FALSE, digits = 4)
