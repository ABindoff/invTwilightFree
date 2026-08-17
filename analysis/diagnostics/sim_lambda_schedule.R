# Does a declination-scheduled lambda beat a constant one? Tested on SIMULATION,
# where the truth is known and the confounds of the real data are absent.
#
# WHY SIMULATE. All 29 real deployments share one phenology -- depart May/June,
# return Jan/Feb -- so time of year, the animal's latitude and its behaviour are
# completely confounded, and the measured slope of lambda_opt on |declination|
# (within-deployment r = -0.31) cannot be attributed to declination from those data.
# Simulation breaks the confound outright: the shading process here is CONSTANT in
# time by construction, so if the best lambda still varies with declination it can
# only be the information geometry of day length, which is the theoretical claim.
#
# PART A -- does the simulation even reproduce the phenomenon? If lambda_opt does not
# vary with |declination| on simulated data with constant shading, then either the
# theory is wrong or the simulator is too clean, and PART B would be untestable. This
# has to pass first.
#
# PART B -- the fit-level comparison. A schedule is only worth having if it beats the
# BEST CONSTANT lambda, not merely the shipped one. `declination_lambda_scale()` is
# normalised to geometric mean 1, so it redistributes the rate without changing its
# level; the level is swept identically for both arms and best-of-each is compared.
#
# SCENARIOS, chosen to break the real data's confound:
#   north_span  northern, solstice into equinox (what the seals do)
#   south_span  the mirror image: southern hemisphere, opposite phase
#   equator     crossing the equator, where day length carries least information
#   solstice_only  a short track sitting near a solstice (control: little declination
#                  variation, so the schedule should do NOTHING here)
# The last is the important control. A schedule that "helps" even where declination
# barely moves is not doing what it claims.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
suppressMessages({ library(data.table); library(terra) })
# run from the package root (setwd removed for portability)
suppressMessages(devtools::load_all(".", quiet = TRUE))
set.seed(20260815)
OUT_A <- "scratch/nes_calibration/sim_lambda_partA.csv"
OUT_B <- "scratch/nes_calibration/sim_lambda_partB.csv"

MAXL <- 100; FLOOR <- 0; AMP <- 100; Z50 <- 96; SCALE <- 5
CAL  <- c(FLOOR, AMP, Z50, SCALE)          # logistic response, used to SIMULATE and to FIT
LAM0 <- 1 / (MAXL * 0.5); PSLAB <- 0.10
STEP_H <- 12; DIFF <- 90; CELL <- 2
P_SHADE <- 0.45                            # CONSTANT in time -- the whole point
N_TRACK <- 8

etrue <- function(z) pmin(pmax(FLOOR + AMP / (1 + exp((z - Z50) / SCALE)), 0), MAXL)
sample_spike <- function(mu, lam = LAM0, maxl = MAXL) {
  m_lo <- 1 - exp(-lam * mu); m_hi <- 0.5 * (1 - exp(-2 * lam * (maxl - mu)))
  if (runif(1) < m_lo / (m_lo + m_hi)) mu + log(1 - runif(1) * m_lo) / lam
  else mu - log(1 - runif(1) * (2 * m_hi)) / (2 * lam)
}
sphere_step <- function(lon, lat, km) {
  b <- runif(1, 0, 2 * pi); d <- abs(rnorm(1, 0, km)) / 6371
  p1 <- lat * pi/180; l1 <- lon * pi/180
  p2 <- asin(sin(p1) * cos(d) + cos(p1) * sin(d) * cos(b))
  l2 <- l1 + atan2(sin(b) * sin(d) * cos(p1), cos(d) - sin(p1) * sin(p2))
  c(((l2 * 180/pi + 180) %% 360) - 180, p2 * 180/pi)
}

# one simulated deployment: truth, then light from it
sim_one <- function(start, lon0, lat0, days) {
  t0 <- as.numeric(as.POSIXct(start, tz = "UTC"))
  knots <- t0 + seq(0, days * 86400, by = STEP_H * 3600)
  K <- length(knots)
  lon <- numeric(K); lat <- numeric(K); lon[1] <- lon0; lat[1] <- lat0
  for (k in 2:K) {
    s <- sphere_step(lon[k-1], lat[k-1], DIFF * sqrt(STEP_H / 24))
    lon[k] <- s[1]; lat[k] <- s[2]
  }
  tt <- seq(min(knots), max(knots), by = 1800)
  ki <- pmin(pmax(findInterval(tt, knots, left.open = TRUE) + 1, 1), K)
  z  <- solar_zenith(tt, lon[ki], lat[ki])
  mu <- etrue(z)
  obs <- vapply(mu, function(m) if (runif(1) < PSLAB) runif(1, 0, MAXL)
                                else sample_spike(m), 0)
  # dive-bout shading, stationary in time: attenuate the sky component in runs
  dive <- logical(length(obs)); st <- FALSE
  for (i in seq_along(dive)) {
    st <- if (st) runif(1) > (1 - P_SHADE)/6 else runif(1) < P_SHADE/6
    dive[i] <- st
  }
  obs[dive] <- obs[dive] * runif(sum(dive))
  list(knots = knots, lon = lon, lat = lat,
       time = as.POSIXct(tt, origin = "1970-01-01", tz = "UTC"),
       light = pmin(pmax(obs, 0), MAXL), ki = ki)
}

SCEN <- list(
  north_span    = list(start = "2021-06-10", lon = 200, lat =  40, days = 210),
  south_span    = list(start = "2021-12-10", lon = 200, lat = -40, days = 210),
  equator       = list(start = "2021-08-01", lon = 200, lat =   5, days = 180),
  solstice_only = list(start = "2021-05-20", lon = 200, lat =  40, days =  60))

cat("simulating", N_TRACK, "tracks per scenario\n")
TR <- list()
for (nm in names(SCEN)) for (r in seq_len(N_TRACK)) {
  s <- SCEN[[nm]]
  TR[[length(TR)+1]] <- c(list(scen = nm, rep = r),
                          sim_one(s$start, s$lon, s$lat, s$days))
}
cat(sprintf("%d tracks; declination range per scenario:\n", length(TR)))
for (nm in names(SCEN)) {
  d <- abs(solar_declination(TR[[which(vapply(TR, function(x) x$scen, "") == nm)[1]]]$knots))
  cat(sprintf("  %-14s |dec| %.1f to %.1f deg\n", nm, min(d), max(d)))
}

# ---------------- PART A: does the simulation reproduce lambda_opt(|dec|)? -------
cat("\n=== PART A: is lambda_opt related to |declination| in SIMULATED data? ===\n")
cat("shading is constant in time by construction, so any relationship is geometry\n")
LATS <- seq(-70, 70, by = 0.5); LAMS <- 2^seq(-3, 5, by = 0.5)
A <- list()
for (x in TR[vapply(TR, function(z) z$rep <= 3, TRUE)]) {
  K <- length(x$knots)
  use <- seq(2, K - 1, length.out = min(40, K - 2))
  for (i in round(use)) {
    j <- which(x$ki == i); if (length(j) < 12) next
    for (lm in LAMS) {
      ll <- eval_logpk_grid(rep(x$lon[i], length(LATS)), LATS, as.numeric(x$time[j]),
                            x$light[j], CAL, c(LAM0 * lm, MAXL, PSLAB), 2)
      w <- exp(ll - max(ll)); s <- sum(w); if (!is.finite(s) || s <= 0) next
      w <- w / s
      A[[length(A)+1]] <- data.table(scen = x$scen, rep = x$rep, win = i, lam_mult = lm,
                                     absdec = abs(solar_declination(x$knots[i])),
                                     bias = sum(w * LATS) - x$lat[i])
    }
  }
}
DA <- rbindlist(A); fwrite(DA, OUT_A)
xing <- function(lam, y) { xx <- log2(lam); o <- order(xx); xx <- xx[o]; y <- y[o]
  k <- which(y[-1] * y[-length(y)] < 0); if (!length(k)) return(NA_real_)
  i <- k[1]; 2^(xx[i] - y[i] * (xx[i+1] - xx[i]) / (y[i+1] - y[i])) }
LO <- DA[, .(lam_opt = xing(lam_mult, bias), absdec = absdec[1]), by = .(scen, rep, win)]
LO <- LO[is.finite(lam_opt)]
cat(sprintf("  %d windows with an identifiable lambda_opt\n", nrow(LO)))
if (nrow(LO) > 20) {
  ct <- suppressWarnings(cor.test(LO$absdec, log2(LO$lam_opt)))
  cat(sprintf("  cor(|declination|, log2 lambda_opt) = %+.3f, p = %.3g, slope %+.3f per degree\n",
              ct$estimate, ct$p.value, coef(lm(log2(lam_opt) ~ absdec, data = LO))[2]))
  cat(sprintf("  (real data gave r = -0.31, slope -0.081)\n"))
  cat("  VERDICT: ")
  cat(if (ct$estimate < 0 && ct$p.value < 0.05)
        "reproduced -- the simulation is a valid test bed for PART B\n"
      else "NOT reproduced -- PART B cannot test the schedule on this simulator\n")
}

# ---------------- PART B: scheduled vs constant, at fit level -------------------
cat("\n=== PART B: scheduled vs constant lambda, at matched levels ===\n")
LEVELS <- c(0.5, 1, 2, 4)
B <- list()
for (x in TR) {
  gg <- rast(xmin = x$lon[1] - 45, xmax = x$lon[1] + 45,
             ymin = max(-85, min(x$lat) - 25), ymax = min(85, max(x$lat) + 25),
             resolution = CELL, crs = "EPSG:4326")
  values(gg) <- 1
  for (lev in LEVELS) for (sched in c("constant", "declination")) {
    ls <- if (sched == "constant") NULL else
      declination_lambda_scale(as.POSIXct(x$knots, origin = "1970-01-01", tz = "UTC"))
    f <- try(suppressWarnings(invisible(capture.output(
      fit <- TwilightFreeGrid(x$time, x$light, grid = gg,
        start_lon = x$lon[1], start_lat = x$lat[1],
        end_lon = x$lon[length(x$lon)], end_lat = x$lat[length(x$lat)],
        step_hours = STEP_H, diffusion = DIFF, calibration = CAL,
        likelihood_params = c(LAM0 * lev, MAXL, PSLAB),
        lambda_scale = ls)))), silent = TRUE)
    if (inherits(f, "try-error")) next
    gp <- grid_posterior(fit); sdp <- numeric(nrow(gp$P))
    for (k in seq_len(nrow(gp$P))) {
      w <- gp$P[k, ]; s <- sum(w)
      if (!is.finite(s) || s <= 0) { sdp[k] <- NA; next }
      w <- w / s; m <- sum(w * gp$lat); sdp[k] <- sqrt(sum(w * (gp$lat - m)^2))
    }
    n <- min(length(fit$fit$lat), length(x$lat))
    e <- fit$fit$lat[seq_len(n)] - x$lat[seq_len(n)]
    ok <- is.finite(e) & is.finite(sdp[seq_len(n)]) & sdp[seq_len(n)] > 0
    B[[length(B)+1]] <- data.table(scen = x$scen, rep = x$rep, level = lev, sched = sched,
      bias = mean(e[ok]), absbias = abs(mean(e[ok])), rmse = sqrt(mean(e[ok]^2)),
      cover = mean(abs(e[ok]) <= 1.96 * sdp[seq_len(n)][ok]),
      km = median(6371 * acos(pmin(1, sin(fit$fit$lat[seq_len(n)]*pi/180)*sin(x$lat[seq_len(n)]*pi/180) +
             cos(fit$fit$lat[seq_len(n)]*pi/180)*cos(x$lat[seq_len(n)]*pi/180)*
             cos((fit$fit$lon[seq_len(n)]-x$lon[seq_len(n)])*pi/180)))[ok]))
  }
  message("track ", x$scen, " rep ", x$rep, " done")
}
DB <- rbindlist(B); fwrite(DB, OUT_B)
cat(sprintf("\n%d fits\n", nrow(DB)))
S <- DB[, .(bias = round(mean(bias), 3), absbias = round(mean(absbias), 3),
            rmse = round(mean(rmse), 3), cover = round(mean(cover), 3),
            km = round(mean(km))), by = .(scen, sched, level)]
for (sc in unique(S$scen)) {
  cat(sprintf("\n--- %s ---\n%8s %13s %9s %9s %9s %8s\n", sc, "level", "schedule",
              "bias", "|bias|", "rmse_lat", "cover"))
  for (i in which(S$scen == sc))
    cat(sprintf("%8.2f %13s %+9.3f %9.3f %9.3f %8.3f\n", S$level[i], S$sched[i],
                S$bias[i], S$absbias[i], S$rmse[i], S$cover[i]))
  best <- S[scen == sc][, .SD[which.min(absbias)], by = sched]
  cat(sprintf("  BEST-OF-EACH |bias|: constant %.3f (level %.2f) | scheduled %.3f (level %.2f)\n",
              best[sched == "constant"]$absbias, best[sched == "constant"]$level,
              best[sched == "declination"]$absbias, best[sched == "declination"]$level))
}
cat("\nA schedule earns its place only by beating the BEST CONSTANT level, and only if\n")
cat("it does nothing in solstice_only, where declination barely moves.\n")
