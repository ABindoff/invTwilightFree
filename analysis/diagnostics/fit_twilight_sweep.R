# Re-tune the darkness regime against TWILIGHT-BAND bias, not median km.
#
# WHY. The band decomposition says daylight is unbiased in every response
# (+0.15 to +0.75) and all the negative bias sits in twilight and night. Twilight is
# also 2-3x sharper than either other band (profile sd 4.4-8.0 vs 11-14), so it is
# the only band that is both informative and wrong -- it is where latitude is
# decided. The darkness regime already halves the twilight bias (-0.695 vs the
# tangent's -1.313) AND tightens the profile (4.41 vs 6.60), which is why it wins on
# bias while tying on distance.
#
# But its three constants were tuned on MEDIAN KM, and the earlier sweep showed
# bias_lat RISING with `ratio` -- so that tuning actively selected against the thing
# we care about. `dark_frac` was never swept at all. Retune on the right objective.
#
# COST. Emission only: profile the light likelihood over latitude at the true Argos
# position, per knot window. No HMM, no movement prior, no fits. Seconds per config.
#
# OBJECTIVE. Drive |twilight bias| toward 0 WITHOUT widening the twilight profile.
# Both matter: a config that flattens the likelihood will show less bias simply by
# saying less, so any candidate whose prof_sd grows is disqualified.
#
# CONTROLS.
#  - DAY band must be unchanged across configs. The darkness regime only acts where
#    the expected curve is below dark_at, so if daylight moves, the regime is
#    leaking outside its support and the sweep is measuring something else.
#  - `tangent` (no darkness regime) is carried as the reference the -1.54 came from.
#  - eval_logpk_grid's 7th argument is SHADE_RATIO, not a cell count. Passing a cell
#    count silently sets an ~80x asymmetric spike; that error produced a spurious
#    poleward bias earlier this week. It is 2 here, matching TwilightFreeGrid.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))
OUT <- "scratch/nes_calibration/twilight_sweep.csv"

STEP_HOURS <- 12; CAL_DAYS <- 15
LATS <- seq(20, 70, by = 0.5)
MAXWIN <- 25            # knot windows per tag
SHADE_RATIO <- 2        # 7th arg of eval_logpk_grid; NOT a cell count

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
cat(sprintf("%d deployments\n", length(tags)))

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
responses <- pool_light_responses(geom_fits, scale_fits)
src <- Filter(Negate(is.null), geom_fits)

profile_stats <- function(ll, truth) {
  ll <- ll - max(ll); w <- exp(ll); s <- sum(w)
  if (!is.finite(s) || s <= 0) return(c(NA, NA))
  w <- w / s; mu <- sum(w * LATS)
  c(mu - truth, sqrt(sum(w * (LATS - mu)^2)))
}

# configs: the current setting and the tangent are carried as references
# Round 1 (dark_frac .20-.50, ratio 20-200, pslab .15-.55) put the optimum on the
# EDGE at f0.50_r200_p0.15 -- the same failure the original km-based sweep had, where
# the answer was outside the grid rather than inside it. Extend all three axes in the
# direction the surface points, and keep the round-1 corner so the two are comparable.
# There is a natural ceiling on `ratio`: the cost of a sharp arm saturates at the
# slab, -log(pslab/max_light), so a very large ratio cannot run away, it just stops
# improving. If the optimum is STILL on the edge, that saturation is the reason.
# Round 2 put the apparent optimum at dark_frac = 0.80 -- and the DAY control caught
# it: day bias moved +0.741 -> +0.285. dark_frac is the threshold as a fraction of
# max_light, so at 0.80 the "darkness" regime has swallowed most of the daylight
# curve and is no longer a darkness regime. Those configs are disqualified however
# good their twilight number looks.
#
# Round 3 holds dark_frac at or below its current 0.35, where day does not move, and
# pushes the two axes that are still improving. ratio 400 -> 800 bought only
# -0.245 -> -0.230, which is the slab saturation predicted above, so this should
# converge rather than run to the edge again.
G <- expand.grid(dark_frac = c(0.25, 0.35),
                 ratio = c(200, 400, 800, 1600),
                 pslab_dark = c(0.01, 0.03, 0.05, 0.10))
# The shipped setting must always be in the grid: it is the reference every
# admissibility test is measured against, and round 3 dropped it (ratio 100), which
# left cur_sd empty and killed the report.
G <- unique(rbind(G, data.frame(dark_frac = 0.35, ratio = 100, pslab_dark = 0.35)))
cat(sprintf("%d darkness configs + tangent reference\n", nrow(G)))

R <- list()
for (tg in names(tags)) {
  g <- tags[[tg]]; r <- responses[[tg]]; if (is.null(r)) next
  a <- as.data.table(g$argos); setorder(a, time)
  tt <- as.numeric(g$light$time)
  obs <- pmax(0, g$light$light - r$baseline)
  tabf <- light_response_table(src, from = 30, to = 140, by = 1,
                               max_light = r$max_light,
                               flat_outside = c(r$saturate_at, r$zero_at))
  brk <- seq(min(tt), max(tt), by = STEP_HOURS * 3600)
  kk <- findInterval(tt, brk)
  tlon <- approx(as.numeric(a$time), a$lon, brk, rule = 2)$y
  tlat <- approx(as.numeric(a$time), a$lat, brk, rule = 2)$y
  lam <- 1 / (r$max_light * 0.5)
  use <- which(tabulate(kk, length(brk)) >= 16)
  if (!length(use)) next
  use <- use[seq(1, length(use), length.out = min(MAXWIN, length(use)))]
  message(tg, " (", length(use), " windows)")
  for (i in use) {
    j <- which(kk == i); if (length(j) < 16) next
    z <- solar_zenith(tt[j], rep(tlon[i], length(j)), rep(tlat[i], length(j)))
    bands <- list(day = z < 85, twilight = z >= 85 & z <= 105, night = z > 105)
    run1 <- function(cal, lp, cfg) {
      for (bn in names(bands)) {
        sel <- bands[[bn]]; if (sum(sel) < 6) next
        ll <- eval_logpk_grid(rep(tlon[i], length(LATS)), LATS, tt[j][sel],
                              obs[j][sel], cal, lp, SHADE_RATIO)
        st <- profile_stats(ll, tlat[i])
        R[[length(R) + 1]] <<- data.table(id = tg, cfg = cfg, band = bn,
                                          bias = st[1], sd = st[2])
      }
    }
    run1(r$calibration, c(lam, r$max_light, 0.10), "tangent")
    for (q in seq_len(nrow(G)))
      run1(tabf, c(lam, r$max_light, 0.10, G$dark_frac[q], G$ratio[q], G$pslab_dark[q]),
           sprintf("f%.2f_r%d_p%.2f", G$dark_frac[q], G$ratio[q], G$pslab_dark[q]))
  }
}
D <- rbindlist(R)
fwrite(D, OUT)

# ---- logger clustering: 15 loggers behind 29 deployments -------------------
a1 <- data.table(file = list.files("data/nes_untracked", pattern = "^[0-9]+_.*\\.csv$"))
a1[, `:=`(topp = sub("^([0-9]+)_.*$", "\\1", file), serial = sub("^[0-9]+_(.+)\\.csv$", "\\1", file))]
b1 <- data.table(file = list.files("fvilches/extracted/TDR raw", pattern = "Archive[.]csv$"))
b1[, `:=`(topp = sub("^([0-9]+)_.*$", "\\1", file),
          serial = sub("^[0-9]+_(.+)-out-Archive[.]csv$", "\\1", file))]
ser <- unique(rbind(a1[, .(topp, serial)], b1[, .(topp, serial)]), by = c("topp", "serial"))
D <- merge(D, ser, by.x = "id", by.y = "topp", all.x = TRUE)

per_log <- D[, .(bias = mean(bias, na.rm = TRUE), sd = mean(sd, na.rm = TRUE)),
             by = .(cfg, band, serial)]
S <- per_log[, .(bias = round(mean(bias, na.rm = TRUE), 3),
                 sd = round(mean(sd, na.rm = TRUE), 2),
                 loggers = .N), by = .(cfg, band)]

ref <- S[cfg == "tangent"]
cat(sprintf("\nREFERENCE tangent: day %+0.3f (sd %.2f) | twilight %+0.3f (sd %.2f) | night %+0.3f\n",
            ref[band=="day"]$bias, ref[band=="day"]$sd,
            ref[band=="twilight"]$bias, ref[band=="twilight"]$sd, ref[band=="night"]$bias))

tw <- S[band == "twilight" & cfg != "tangent"][order(abs(bias))]
dy <- S[band == "day"]; ni <- S[band == "night"]
cur <- "f0.35_r100_p0.35"
cat(sprintf("CURRENT  %s: twilight %+0.3f (sd %.2f)\n\n", cur,
            tw[cfg == cur]$bias, tw[cfg == cur]$sd))
cat("top 12 by |twilight bias|. A config is only admissible if BOTH the twilight\n")
cat("profile does not widen AND the day band does not move: dark_frac too large\n")
cat("makes the darkness regime swallow daylight, which flatters twilight for the\n")
cat("wrong reason. day_ref is the day bias at the smallest dark_frac tested.\n\n")
cat(sprintf("%-20s %9s %7s %9s %9s %9s %10s\n",
            "config", "twi_bias", "twi_sd", "night", "day", "day_shift", "verdict"))
cur_sd <- tw[cfg == cur]$sd
day_ref <- dy[cfg != "tangent"][which.max(dy[cfg != "tangent"]$bias)]$bias
for (i in seq_len(min(12, nrow(tw)))) {
  c0 <- tw$cfg[i]
  dsh <- dy[cfg == c0]$bias - day_ref
  bad <- character(0)
  if (tw$sd[i] > cur_sd * 1.02) bad <- c(bad, "WIDER")
  if (abs(dsh) > 0.05) bad <- c(bad, "DAY-LEAK")
  cat(sprintf("%-20s %+9.3f %7.2f %+9.3f %+9.3f %+9.3f %10s\n", c0, tw$bias[i], tw$sd[i],
              ni[cfg == c0]$bias, dy[cfg == c0]$bias, dsh,
              if (!length(bad)) "admissible" else paste(bad, collapse = "+")))
}
adm <- tw[vapply(tw$cfg, function(c0)
            abs(dy[cfg == c0]$bias - day_ref) <= 0.05, logical(1)) & tw$sd <= cur_sd * 1.02]
if (nrow(adm)) cat(sprintf("\nBEST ADMISSIBLE: %s  twilight %+0.3f (sd %.2f) vs current %+0.3f (sd %.2f)\n",
                           adm$cfg[1], adm$bias[1], adm$sd[1],
                           tw[cfg == cur]$bias, tw[cfg == cur]$sd))
# CONTROL, computed over the darkness configs ONLY. Including `tangent` here was
# wrong: it uses a different RESPONSE (clamped linear vs table), so its day bias
# differs for a reason that has nothing to do with the regime, and the spread looked
# like a leak when it was a response difference.
dyd <- dy[cfg != "tangent"]
cat(sprintf("\nCONTROL day-band spread across DARKNESS configs: %.4f\n",
            diff(range(dyd$bias, na.rm = TRUE))))
cat(sprintf("  (tangent's day bias %+0.3f differs because of the RESPONSE, not the regime)\n",
            dy[cfg == "tangent"]$bias))
cat("  the regime acts only below dark_at, so daylight must not move across configs\n")
cat(sprintf("  verdict: %s\n", if (diff(range(dyd$bias, na.rm = TRUE)) < 0.1) "clean"
            else "*** the regime is leaking into daylight"))
edge <- tw$cfg[1]
cat(sprintf("\nbest config %s -- check whether it is on the grid edge again\n", edge))
cat(sprintf("  grid: dark_frac %s | ratio %s | pslab_dark %s\n",
            paste(unique(G$dark_frac), collapse = "/"),
            paste(unique(G$ratio), collapse = "/"),
            paste(unique(G$pslab_dark), collapse = "/")))
