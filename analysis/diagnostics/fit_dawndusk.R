# A GROUND-TRUTH-FREE criterion for lambda: dawn/dusk self-consistency.
#
# WHY THIS AND NOT LIKELIHOOD. Per-tag lambda by marginal likelihood was tried and
# failed with 0/10 interior optima -- log_z is monotone in lambda, because the light
# residuals are autocorrelated (ACF 0.26-0.33 over 1-4 h) so the likelihood
# over-counts independent observations and always rewards a sharper emission. Any
# criterion treating observations as independent inherits that. This one does not:
# it compares two halves of the SAME day, so shared slow structure is differenced out
# rather than counted twice.
#
# THE IDEA. Each day contains two transitions. Estimate latitude from the rising half
# alone and from the setting half alone; both estimate the same quantity, so their
# difference needs no truth. Form
#     z = (lat_dawn - lat_dusk) / sqrt(sd_dawn^2 + sd_dusk^2)
# If the emission model is calibrated, sd(z) = 1 across days. Too tight a lambda makes
# the model over-confident: claimed sds shrink faster than the real disagreement, so
# sd(z) > 1. Too loose and sd(z) < 1. Choosing lambda so that sd(z) = 1 therefore has
# an interior optimum BY CONSTRUCTION, which is exactly what log_z lacked.
#
# mean(z) is reported too and is a different diagnostic: a systematic dawn-vs-dusk
# disagreement is not a dispersion problem at all, it is an ASYMMETRY -- a response
# that is wrong on one limb, or a clock offset. If mean(z) is large, lambda is not
# the whole story and no amount of tuning it will help.
#
# TRUTH-FREE THROUGHOUT. The only external quantity used is the tag's own light.
#  - longitude is estimated per window from the data (2-D coarse pass, both halves
#    together), not taken from Argos. It is estimated ONCE at lambda = 1 and reused
#    across lambda, because longitude is the carrier PHASE and barely depends on the
#    shading rate -- checked below by reporting its spread.
#  - the dawn/dusk split is the sign of d(zenith)/dt at the estimated position, i.e.
#    the physical rising/setting split, not a clock convention.
#
# HONEST CAVEAT, stated before the numbers: the two halves of one day share that
# day's weather and the animal's behavioural state, so they are NOT fully
# independent. Shared conditions make them agree MORE than chance, which deflates
# sd(z) and biases the chosen lambda TIGHT. The size of that is unknown; it is
# bounded below by using halves of the same day rather than different days, and it is
# why the criterion must be validated against the Argos-optimal lambda rather than
# trusted on its own.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))
OUT <- "scratch/nes_calibration/dawndusk.csv"

CAL_DAYS <- 15
LATS <- seq(20, 70, by = 0.5)
LONS <- seq(150, 250, by = 2)        # coarse: only needed to place the window
LAMS <- 2^seq(-4, 5, by = 0.5)       # same ladder as the wide Argos sweep
MAXWIN <- 20                         # 24 h windows per tag
SHADE_RATIO <- 2
MIN_OBS <- 12                        # per half

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
  tags[[tg]] <- list(id = tg, light = d, p0 = c(df$deploy_lon, df$deploy_lat))
}
for (tg in names(new)) {
  mm <- new_man[topp == tg]; if (!nrow(mm)) next
  a <- as.data.table(new_argos[[tg]]); d <- as.data.table(new[[tg]])
  d <- d[time >= min(a$time) & time <= max(a$time)]; if (nrow(d) < 1000) next
  tags[[tg]] <- list(id = tg, light = d, p0 = c(mm$deploy_lon[1], mm$deploy_lat[1]))
}

a1 <- data.table(file = list.files("data/nes_untracked", pattern = "^[0-9]+_.*\\.csv$"))
a1[, `:=`(topp = sub("^([0-9]+)_.*$", "\\1", file), serial = sub("^[0-9]+_(.+)\\.csv$", "\\1", file))]
b1 <- data.table(file = list.files("fvilches/extracted/TDR raw", pattern = "Archive[.]csv$"))
b1[, `:=`(topp = sub("^([0-9]+)_.*$", "\\1", file),
          serial = sub("^[0-9]+_(.+)-out-Archive[.]csv$", "\\1", file))]
SER <- unique(rbind(a1[, .(topp, serial)], b1[, .(topp, serial)]), by = c("topp", "serial"))
family <- function(s) fifelse(is.na(s), "unknown",
                       fifelse(grepl("^219", s), "Mk9_219",
                        fifelse(grepl("^18A", s), "F18A", "other")))
for (tg in names(tags)) tags[[tg]]$family <- family(SER[topp == tg]$serial[1])

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
cat(sprintf("%d deployments; %d lambdas\n\n", length(tags), length(LAMS)))

pstats <- function(ll, grid) {
  ll <- ll - max(ll); w <- exp(ll); s <- sum(w)
  if (!is.finite(s) || s <= 0) return(c(NA_real_, NA_real_))
  w <- w / s; mu <- sum(w * grid)
  c(mu, sqrt(sum(w * (grid - mu)^2)))
}

R <- list(); LONCHK <- list()
for (tg in names(tags)) {
  g <- tags[[tg]]; r <- responses[[tg]]; if (is.null(r)) next
  tt  <- as.numeric(g$light$time)
  obs <- pmax(0, g$light$light - r$baseline)
  cal <- r$calibration; lam0 <- 1 / (r$max_light * 0.5)
  brk <- seq(min(tt), max(tt), by = 86400)          # 24 h windows: one dawn, one dusk
  kk  <- findInterval(tt, brk)
  use <- which(tabulate(kk, length(brk)) >= 2 * MIN_OBS)
  if (!length(use)) next
  use <- use[seq(1, length(use), length.out = min(MAXWIN, length(use)))]
  message(tg, " (", length(use), " days)")
  for (i in use) {
    j <- which(kk == i); if (length(j) < 2 * MIN_OBS) next
    lp1 <- c(lam0, r$max_light, 0.10)
    # --- longitude and a reference latitude, from the data alone ---------------
    gr <- expand.grid(lon = LONS, lat = LATS)
    ll2 <- eval_logpk_grid(gr$lon, gr$lat, tt[j], obs[j], cal, lp1, SHADE_RATIO)
    bi <- which.max(ll2); lon_hat <- gr$lon[bi]; lat_ref <- gr$lat[bi]
    LONCHK[[length(LONCHK) + 1]] <- data.table(id = tg, win = i, lon_hat = lon_hat)
    # --- physical rising/setting split at the estimated position ---------------
    z <- solar_zenith(tt[j], rep(lon_hat, length(j)), rep(lat_ref, length(j)))
    if (diff(range(z)) < 20) next            # no usable transition (polar day/night)
    dz <- c(diff(z), 0)
    dawn <- dz < 0; dusk <- dz > 0           # zenith falling = rising sun
    if (sum(dawn) < MIN_OBS || sum(dusk) < MIN_OBS) next
    for (lm in LAMS) {
      lp <- c(lam0 * lm, r$max_light, 0.10)
      sa <- pstats(eval_logpk_grid(rep(lon_hat, length(LATS)), LATS, tt[j][dawn],
                                   obs[j][dawn], cal, lp, SHADE_RATIO), LATS)
      sb <- pstats(eval_logpk_grid(rep(lon_hat, length(LATS)), LATS, tt[j][dusk],
                                   obs[j][dusk], cal, lp, SHADE_RATIO), LATS)
      if (any(!is.finite(c(sa, sb))) || sa[2] <= 0 || sb[2] <= 0) next
      R[[length(R) + 1]] <- data.table(
        id = tg, win = i, lam_mult = lm,
        lat_dawn = sa[1], sd_dawn = sa[2], lat_dusk = sb[1], sd_dusk = sb[2],
        z = (sa[1] - sb[1]) / sqrt(sa[2]^2 + sb[2]^2))
    }
  }
}
D <- rbindlist(R)
D <- merge(D, SER, by.x = "id", by.y = "topp", all.x = TRUE)
fwrite(D, OUT)
L <- rbindlist(LONCHK)
cat(sprintf("\nwrote %d rows (%d tags, %d day-windows)\n",
            nrow(D), uniqueN(D$id), uniqueN(paste(D$id, D$win))))
cat(sprintf("longitude estimated once per window at lambda = 1; per-tag sd of lon_hat: median %.2f deg\n",
            median(L[, sd(lon_hat), by = id]$V1, na.rm = TRUE)))
cat("  (a small spread supports reusing it across lambda; the carrier phase should\n")
cat("   not care about the shading rate)\n")
