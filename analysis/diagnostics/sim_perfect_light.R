# THE DECISIVE EXPERIMENT: fit PERFECT light generated at the known positions.
#
# Every diagnostic so far has asked what is wrong with the light. This asks
# whether anything need be wrong with it at all. Generate light exactly as the
# model believes it is generated -- the fitted response evaluated at the true
# Argos positions, at the real observation times -- and fit it with the same
# engine, the same calibration, the same prior, the same endpoints.
#
# The model is then EXACTLY correct for the data. Any residual bias or
# under-coverage cannot be observational; it must come from the estimator: the
# grid discretisation, the movement prior, the endpoint bridge, or the
# non-linearity of the map from timing to latitude.
#
# WHY I EXPECT SOME. Near the equinox dH/dphi -> 0, so dphi/dH -> infinity: a
# symmetric error in twilight timing maps to a strongly asymmetric error in
# latitude. That is a transformation-of-variables effect and would bias the
# estimate with perfect data and a perfectly specified model. No test on the
# real light could have found it, because it is not in the light.
#
# THREE ARMS, each isolating one thing:
#   A  perfect, noiseless          -> the estimator alone
#   B  perfect + spike/slab noise  -> estimator + the model's own noise
#   C  perfect + OBSERVED shading  -> estimator + real attenuation structure,
#      by transplanting each tag's own observed-minus-theoretical residuals
#
# If A is unbiased and C is not, the residuals in C are the thing to target and
# the diagnostic script says what they look like.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))
OUT <- "scratch/nes_calibration/sim_perfect_results.csv"
OUTK <- "scratch/nes_calibration/sim_perfect_knots.csv"
set.seed(42)

STEP_HOURS <- 12; DIFFUSION <- 110; CELL <- 1; CAL_DAYS <- 15
ARMS <- c("A_perfect", "B_noise", "C_observed_shading")

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
  d <- x[time >= t0 & time <= t1]
  if (nrow(d) < 1000) next
  tf <- a[time >= max(time) - 72*3600]
  tags[[tg]] <- list(id = tg, light = d, argos = a,
                     p0 = c(df$deploy_lon, df$deploy_lat),
                     p1 = c(median(tf$lon), median(tf$lat)))
}
for (tg in names(new)) {
  mm <- new_man[topp == tg]; if (!nrow(mm)) next
  a <- as.data.table(new_argos[[tg]]); d <- as.data.table(new[[tg]])
  d <- d[time >= min(a$time) & time <= max(a$time)]
  if (nrow(d) < 1000) next
  tags[[tg]] <- list(id = tg, light = d, argos = a,
                     p0 = c(mm$deploy_lon[1], mm$deploy_lat[1]),
                     p1 = c(mm$recover_lon[1], mm$recover_lat[1]))
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

grid <- rast(xmin = 150, xmax = 250, ymin = 20, ymax = 70,
             resolution = CELL, crs = "EPSG:4326")
values(grid) <- 1

# the model's own expected-light curve: clamped linear in zenith
expected <- function(z, cal, maxl) pmin(pmax(cal[1] - cal[2]*z, 0), maxl)

# truth at arbitrary times, interpolated and extended to the ends so the
# simulated series has no gaps (a gap would change the sampling, not the physics)
truth_at <- function(a, tm) {
  setorder(a, time)
  list(lon = approx(as.numeric(a$time), a$lon, as.numeric(tm), rule = 2)$y,
       lat = approx(as.numeric(a$time), a$lat, as.numeric(tm), rule = 2)$y)
}

done <- if (file.exists(OUT)) paste(fread(OUT)$id, fread(OUT)$arm) else character(0)
cat(sprintf("%d of %d fits already done\n\n", length(done), length(tags)*length(ARMS)))

for (arm in ARMS) for (tg in names(tags)) {
  if (paste(tg, arm) %in% done) next
  g <- tags[[tg]]; r <- responses[[tg]]
  if (is.null(r)) next
  message(arm, " : ", tg)
  tr <- truth_at(as.data.table(g$argos), g$light$time)
  z <- solar_zenith(as.numeric(g$light$time), tr$lon, tr$lat)
  mu <- expected(z, r$calibration, r$max_light)

  if (arm == "A_perfect") {
    sim <- mu
  } else if (arm == "B_noise") {
    # the model's own error model: a slab fraction uniform over the range, the
    # rest an asymmetric exponential spike about mu
    lam <- 1/(r$max_light*0.5); n <- length(mu)
    slab <- runif(n) < 0.10
    e <- rexp(n, rate = lam) * ifelse(runif(n) < 0.5, -1, 0.5)
    sim <- ifelse(slab, runif(n, 0, r$max_light), pmin(pmax(mu + e, 0), r$max_light))
  } else {
    # transplant this tag's OWN observed-minus-theoretical residuals, keeping
    # their dependence on zenith by resampling within zenith bins
    obs <- pmax(0, g$light$light - r$baseline)
    res <- obs - mu
    zb <- cut(z, breaks = seq(0, 180, by = 5))
    sim <- mu
    for (b in unique(zb[!is.na(zb)])) {
      i <- which(zb == b)
      sim[i] <- pmin(pmax(mu[i] + sample(res[i], length(i), replace = TRUE), 0),
                     r$max_light)
    }
  }

  t <- system.time(invisible(capture.output(
    f <- TwilightFreeGrid(g$light$time, sim, grid = grid,
      start_lon = unname(g$p0[1]), start_lat = unname(g$p0[2]),
      end_lon = unname(g$p1[1]), end_lat = unname(g$p1[2]),
      step_hours = STEP_HOURS, diffusion = DIFFUSION, calibration = r$calibration,
      likelihood_params = c(1/(r$max_light*0.5), r$max_light, 0.10)))))[3]
  gp <- grid_posterior(f)
  sdp <- sdl <- numeric(nrow(gp$P))
  for (i in seq_len(nrow(gp$P))) {
    w <- gp$P[i, ]; s <- sum(w)
    if (!is.finite(s) || s <= 0) { sdp[i] <- sdl[i] <- NA; next }
    w <- w/s
    sdl[i] <- sqrt(sum(w*(gp$lon - sum(w*gp$lon))^2))
    sdp[i] <- sqrt(sum(w*(gp$lat - sum(w*gp$lat))^2))
  }
  tm <- as.POSIXct(f$fit$time, origin = "1970-01-01", tz = "UTC")
  tt <- truth_at(as.data.table(g$argos), tm)
  ep <- f$fit$lat - tt$lat; el <- dlon(lon360(f$fit$lon), tt$lon)
  e <- gc_km(lon360(f$fit$lon), f$fit$lat, tt$lon, tt$lat)
  ok <- is.finite(ep) & is.finite(sdp) & sdp > 0
  fwrite(data.table(id = tg, arm = arm, time = tm, err_lat = ep, err_lon = el,
                    err_km = e, lat_sd = sdp, true_lat = tt$lat),
         OUTK, append = file.exists(OUTK))
  fwrite(data.table(id = tg, arm = arm, n = sum(ok),
                    median_km = round(median(e[ok])),
                    bias_lat = round(mean(ep[ok]), 3),
                    rmse_lat = round(sqrt(mean(ep[ok]^2)), 3),
                    rmse_lon = round(sqrt(mean(el[ok]^2)), 3),
                    cover_lat = round(mean(abs(ep[ok]) <= 1.96*sdp[ok]), 3),
                    cover_lon = round(mean(abs(el[ok]) <= 1.96*sdl[ok]), 3),
                    secs = round(t)),
         OUT, append = file.exists(OUT))
}

o <- fread(OUT)
cat("\n=== simulated light, fitted with the model that generated it ===\n")
print(as.data.frame(o[, .(tags = .N, median_km = round(mean(median_km)),
  bias_lat = round(mean(bias_lat), 3), rmse_lat = round(mean(rmse_lat), 3),
  rmse_lon = round(mean(rmse_lon), 3),
  cover_lat = round(mean(cover_lat), 3), cover_lon = round(mean(cover_lon), 3)),
  by = arm][order(arm)]), row.names = FALSE)
cat("\nREAL light, for comparison: 244 km, rmse_lat 2.53, cover_lat 0.72\n")
cat("\nArm A is the estimator with nothing wrong with the data. If it is unbiased\n")
cat("and well covered, the bias is observational and arm C localises it. If A is\n")
cat("biased, the estimator is and no amount of work on the light will help.\n")
