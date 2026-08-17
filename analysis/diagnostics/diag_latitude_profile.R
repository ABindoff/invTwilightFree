# THE EMISSION'S OWN OPINION ABOUT LATITUDE.
#
# Every bias test so far has fitted whole tracks and read the outcome, which
# confounds the emission with the movement prior, the endpoint bridge and the
# smoother. Ten mechanisms have died in that loop. This instrument removes the
# confound: for each 12-hour window, evaluate the log-likelihood across a
# latitude grid at the TRUE longitude, and read the profile directly. That is
# the emission's entire opinion about latitude, with nothing else in it.
#
# Two things come out.
#
# 1. EMISSION-ONLY BIAS. Posterior mean of the profile minus the truth. If this
#    is small while track bias is large, the response is not the problem and
#    every response we have tried has been beside the point.
#
# 2. WHERE IT IS MADE. Recompute the profile from daytime observations only,
#    then twilight only, then night only. Whichever band's profile is displaced
#    is where the bias is manufactured. That localises it regardless of which
#    mechanism is responsible, which no hypothesis-driven test could do.
#
# Also reported: the profile's own width against the spread of its peak. That is
# the honest information content of a single window, and it speaks to the
# session's central finding that error does not track available information.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))
CAL_DAYS <- 15; STEP_HOURS <- 12
LATS <- seq(20, 70, by = 0.25)

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

# The four responses that have been scored against Argos, so the profile can be
# checked against a known outcome for each.
ARMS <- c("tangent", "logistic", "table_flat", "dark")

profile_stats <- function(ll, truth) {
  ll <- ll - max(ll)
  w <- exp(ll); s <- sum(w)
  if (!is.finite(s) || s <= 0) return(c(NA, NA, NA))
  w <- w / s
  mu <- sum(w * LATS)
  sd <- sqrt(sum(w * (LATS - mu)^2))
  c(mu - truth, LATS[which.max(ll)] - truth, sd)
}

R <- list()
for (tg in names(tags)) {
  g <- tags[[tg]]; r <- responses[[tg]]; if (is.null(r)) next
  a <- as.data.table(g$argos); setorder(a, time)
  tt <- as.numeric(g$light$time)
  obs <- pmax(0, g$light$light - r$baseline)
  tab <- light_response_table(src, from = 30, to = 140, by = 1, max_light = r$max_light)
  win <- c(r$saturate_at, r$zero_at)
  tabf <- light_response_table(src, from = 30, to = 140, by = 1,
                               max_light = r$max_light, flat_outside = win)
  brk <- seq(min(tt), max(tt), by = STEP_HOURS * 3600)
  kk <- findInterval(tt, brk)
  tlon <- approx(as.numeric(a$time), a$lon, brk, rule = 2)$y
  tlat <- approx(as.numeric(a$time), a$lat, brk, rule = 2)$y
  lam <- 1 / (r$max_light * 0.5)
  cfg <- list(
    tangent    = list(cal = r$calibration,          lp = c(lam, r$max_light, 0.10)),
    logistic   = list(cal = r$calibration_logistic, lp = c(lam, r$max_light, 0.10)),
    table_flat = list(cal = tabf,                   lp = c(lam, r$max_light, 0.10)),
    dark       = list(cal = tabf, lp = c(lam, r$max_light, 0.10, 0.35, 100, 0.35)))

  use <- which(tabulate(kk, length(brk)) >= 16)
  use <- use[seq(1, length(use), length.out = min(40, length(use)))]
  for (i in use) {
    j <- which(kk == i); if (length(j) < 16) next
    z <- solar_zenith(tt[j], rep(tlon[i], length(j)), rep(tlat[i], length(j)))
    bands <- list(all = rep(TRUE, length(j)),
                  day = z < 85, twilight = z >= 85 & z <= 105, night = z > 105)
    for (arm in ARMS) for (bn in names(bands)) {
      sel <- bands[[bn]]
      if (sum(sel) < 6) next
      ll <- eval_logpk_grid(rep(tlon[i], length(LATS)), LATS, tt[j][sel],
                            obs[j][sel], cfg[[arm]]$cal, cfg[[arm]]$lp, 2)
      st <- profile_stats(ll, tlat[i])
      R[[length(R) + 1]] <- data.table(id = tg, arm = arm, band = bn,
                                       n_obs = sum(sel), true_lat = tlat[i],
                                       bias_mean = st[1], bias_mode = st[2],
                                       prof_sd = st[3])
    }
  }
}
R <- rbindlist(R)
cat(sprintf("%d profiles: %d tags, %d windows\n\n", nrow(R), uniqueN(R$id),
            uniqueN(paste(R$id, R$true_lat))))

trk <- data.table(arm = c("tangent","logistic","table_flat","dark"),
                  track_bias = c(-0.14, -2.24, +1.96, +1.40),
                  track_km = c(242, 615, 451, 289))
cat("=== EMISSION-ONLY BIAS (all bands), against the fitted-track bias ===\n")
s <- R[band == "all" & is.finite(bias_mean),
       .(windows = .N, emis_bias = round(mean(bias_mean), 3),
         emis_mode = round(mean(bias_mode), 3),
         prof_sd = round(mean(prof_sd), 2)), by = arm]
print(as.data.frame(merge(s, trk, by = "arm")), row.names = FALSE)
cat("\nIf emis_bias tracks track_bias, the response IS the problem and step 2 is\n")
cat("well posed. If it is near zero for every arm, the bias is made downstream.\n")

cat("\n=== WHERE IS IT MADE? emission bias by zenith band ===\n")
print(as.data.frame(dcast(R[is.finite(bias_mean),
      .(b = round(mean(bias_mean), 2)), by = .(arm, band)],
      arm ~ band, value.var = "b")), row.names = FALSE)
cat("\n=== and how much does each band think it knows? (profile sd) ===\n")
print(as.data.frame(dcast(R[is.finite(prof_sd), .(s = round(mean(prof_sd), 1)),
      by = .(arm, band)], arm ~ band, value.var = "s")), row.names = FALSE)

cat("\n=== information content: profile sd vs the spread of the profile mean ===\n")
for (aa in ARMS) {
  x <- R[arm == aa & band == "all" & is.finite(bias_mean)]
  cat(sprintf("  %-11s says sd %.2f, actually scatters %.2f  (ratio %.2f)\n",
              aa, mean(x$prof_sd), sd(x$bias_mean), mean(x$prof_sd)/sd(x$bias_mean)))
}
fwrite(R, "scratch/nes_calibration/latitude_profiles.csv")
