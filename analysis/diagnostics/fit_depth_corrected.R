# Does correcting light for measured depth move the tracks?
#
# THE TEST (user's framing, and it is sharper than treating this as a bias fix):
# the spike's lower arm exists to absorb shading. If it is doing its job, removing
# the shading from the DATA should leave the fitted tracks essentially where they
# were. Divergence is the finding -- it says the shading model is not absorbing
# what it was built to absorb, and that the likelihood is paying for attenuation
# with position instead.
#
# Light is decimated by MAXIMUM within each 30-min bin, so the retained sample was
# taken at that bin's shallowest point, `depth_min`, and is attenuated by
# Beer-Lambert from there:
#
#     reading = dark + sky * exp(-k * depth_min)
#  => sky_at_surface = (reading - dark) * exp(+k * depth_min)
#
# Only the SKY component is attenuated; the dark current is the sensor and does not
# care how deep it is. Correcting the raw reading instead would inflate the night
# floor into fictitious light, which is exactly the artefact the darkness regime was
# invented to absorb, so the two would become entangled.
#
# k = 0.07 /m, mid-range for clear open ocean in the blue-green band (Jerlov I-II).
# No sweep: the question is whether the tracks move at all, and a single physically
# defensible k answers that. If they move, the size of k becomes worth arguing about.
#
# WHAT MAKES THE COMPARISON CLEAN: both arms use the same calibration, the same
# likelihood parameters, the same knots, the same endpoints and the same seed. The
# ONLY difference is the light column. TwilightFreeHier takes one calibration for
# all tags while max_light runs 147-175, so each tag's light is normalised by its
# own max_light onto a common 0-100 scale first -- applied identically to both arms,
# so it cannot manufacture a difference.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))
OUT <- "scratch/nes_calibration/depth_corrected_results.csv"

K_ATTEN   <- 0.07      # /m
# Beer-Lambert inversion amplifies NOISE exponentially in depth, and at depth the
# reading is mostly dark current, so a generous cap manufactures fake daylight out
# of sensor noise: at a 30 m cap (8.2x) the amplified deep samples became the
# brightest in the record and the fits collapsed to -40 deg of latitude bias.
# 10 m keeps the multiplier under 2x and still covers 96% of samples untouched.
DEPTH_CAP <- 10        # m; multiplier caps at exp(0.07*10) = 2.01x
STEP_HOURS <- 12; CAL_DAYS <- 15; DIFFUSION <- 110; CELL <- 1
HEADROOM  <- 1.4       # widen the likelihood support so the corrected arm is not clipped
SWEEPS <- 2500L; BURN <- 800L; THIN <- 5L; SEED <- 4L

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
  tf <- a[time >= max(time) - 72*3600]
  tags[[tg]] <- list(id = tg, light = d, argos = a,
                     p0 = c(df$deploy_lon, df$deploy_lat),
                     p1 = c(median(tf$lon), median(tf$lat)))
}
for (tg in names(new)) {
  mm <- new_man[topp == tg]; if (!nrow(mm)) next
  a <- as.data.table(new_argos[[tg]]); d <- as.data.table(new[[tg]])
  d <- d[time >= min(a$time) & time <= max(a$time)]; if (nrow(d) < 1000) next
  tags[[tg]] <- list(id = tg, light = d, argos = a,
                     p0 = c(mm$deploy_lon[1], mm$deploy_lat[1]),
                     p1 = c(mm$recover_lon[1], mm$recover_lat[1]))
}
MAXTAGS <- as.integer(Sys.getenv("NES_MAX_TAGS", "0"))
if (MAXTAGS > 0 && length(tags) > MAXTAGS) {
  tags <- tags[seq_len(MAXTAGS)]
  cat(sprintf("SMOKE TEST: limited to %d deployments\n", MAXTAGS))
}
cat(sprintf("%d deployments\n", length(tags)))

# per-tag response, only for its baseline and max_light (the scale, not the shape)
resp <- lapply(tags, function(g) {
  k <- g$light$time <= min(g$light$time) + CAL_DAYS * 86400
  fit_light_response(g$light$time[k], g$light$light[k], g$p0[1], g$p0[2])
})

# ---- build the two light columns -------------------------------------------
bind <- 0; tot <- 0; clip <- 0
mk <- function(tg, corrected) {
  g <- tags[[tg]]; r <- resp[[tg]]
  sky <- pmax(as.numeric(g$light$light) - r$baseline, 0)
  if (corrected) {
    dm <- suppressWarnings(as.numeric(g$light$depth_min))
    dm[!is.finite(dm)] <- 0
    dm <- pmin(pmax(dm, 0), DEPTH_CAP)
    bind <<- bind + sum(dm >= DEPTH_CAP); tot <<- tot + length(dm)
    sky <- sky * exp(K_ATTEN * dm)
  }
  # Raw baseline-subtracted units, per tag, exactly as fit_rescore29.R feeds the
  # grid engine -- so the uncorrected arm must reproduce its 276 km / -1.54 result.
  # That reproduction is the control: without it the corrected arm means nothing.
  #
  # Two earlier attempts failed here and both were the harness, not the physics:
  # normalising by each arm's own upper quantile let the correction shrink
  # everything (amplified deep samples set the quantile), and clamping at the old
  # max_light flattened the corrected daytime curve against the ceiling. The
  # likelihood's support is therefore widened per tag to give the corrected arm
  # headroom, identically in both arms.
  clip <<- clip + sum(sky > r$max_light * HEADROOM)
  data.frame(time = g$light$time, light = pmin(sky, r$max_light * HEADROOM))
}
d0 <- lapply(names(tags), function(t) mk(t, FALSE)); names(d0) <- names(tags)
d1 <- lapply(names(tags), function(t) mk(t, TRUE));  names(d1) <- names(tags)
cat(sprintf("depth cap binds on %.3f%% of samples; support clips %.3f%% (both arms)\n",
            100*bind/tot, 100*clip/(2*tot)))
cat(sprintf("mean per-tag ratio of corrected to uncorrected daytime light: %.3f\n\n",
            mean(mapply(function(a, b) {
              k <- a$light > 0.2 * max(a$light)
              if (!any(k)) NA else mean(b$light[k]) / mean(a$light[k])
            }, d0, d1), na.rm = TRUE)))

grid <- rast(xmin = 150, xmax = 250, ymin = 20, ymax = 70,
             resolution = CELL, crs = "EPSG:4326")
values(grid) <- 1

fit_one <- function(tg, dat) {
  g <- tags[[tg]]; r <- resp[[tg]]
  ml <- r$max_light * HEADROOM
  invisible(capture.output(
    f <- TwilightFreeGrid(dat$time, dat$light, grid = grid,
      start_lon = unname(g$p0[1]), start_lat = unname(g$p0[2]),
      end_lon = unname(g$p1[1]), end_lat = unname(g$p1[2]),
      step_hours = STEP_HOURS, diffusion = DIFFUSION,
      calibration = r$calibration, likelihood_params = c(1/(ml*0.5), ml, 0.10),
      area_correction = TRUE)))
  f
}

res <- rbindlist(lapply(names(tags), function(tg) {
  message("fitting ", tg)
  f0 <- fit_one(tg, d0[[tg]]); f1 <- fit_one(tg, d1[[tg]])
  tm <- as.POSIXct(f0$fit$time, origin = "1970-01-01", tz = "UTC")
  tr <- argos_at(as.data.table(tags[[tg]]$argos), tm)
  sep <- gc_km(lon360(f0$fit$lon), f0$fit$lat, lon360(f1$fit$lon), f1$fit$lat)
  e0 <- gc_km(lon360(f0$fit$lon), f0$fit$lat, tr$lon, tr$lat)
  e1 <- gc_km(lon360(f1$fit$lon), f1$fit$lat, tr$lon, tr$lat)
  data.table(id = tg, n = length(e0),
             sep_km = round(median(sep, na.rm = TRUE), 1),
             km0 = round(median(e0, na.rm = TRUE)), km1 = round(median(e1, na.rm = TRUE)),
             bias0 = round(mean(f0$fit$lat - tr$lat, na.rm = TRUE), 3),
             bias1 = round(mean(f1$fit$lat - tr$lat, na.rm = TRUE), 3))
}))
fwrite(res, OUT)
cat("\n=== per-tag ===\n"); print(as.data.frame(res), row.names = FALSE)
cat(sprintf("\nmedian separation between the two tracks: %.1f km\n", median(res$sep_km)))
cat(sprintf("accuracy   uncorrected %.0f km -> corrected %.0f km  (better on %d/%d)\n",
            median(res$km0), median(res$km1), sum(res$km1 < res$km0), nrow(res)))
cat(sprintf("lat bias   uncorrected %+.3f -> corrected %+.3f  (|bias| smaller on %d/%d)\n",
            mean(res$bias0), mean(res$bias1),
            sum(abs(res$bias1) < abs(res$bias0)), nrow(res)))
cat(sprintf("paired Wilcoxon on |bias|: p = %.4f\n",
            suppressWarnings(wilcox.test(abs(res$bias0), abs(res$bias1), paired = TRUE)$p.value)))
cat("\nSMALL separation => the shading arm is absorbing the attenuation, as designed.\n")
cat("LARGE separation => it is not, and position is paying for depth.\n")
