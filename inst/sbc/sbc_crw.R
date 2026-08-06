# =====================================================================
# SBC for the correlated random walk: the per-tag persistence rho and the
# movement variance sig2 that is sampled alongside it.
#
# Self-consistency SBC. Data are simulated from the inference model's own
# prior and likelihood, refit, and the rank of each true value within its
# posterior accumulated. Correct posteriors give uniform ranks, tested
# against the simultaneous ECDF bands of Sailynoja, Burkner & Vehtari
# (2022) so that the WHOLE curve is covered at the stated confidence, not
# each point separately.
#
#   rho    ~ Uniform(-0.9, 0.9)            [matches the sampler's support]
#   sig2   ~ InvGamma(a0, b0)
#   x_1    = deploy (fixed)
#   x_2    = x_1 + N(0, sig2 * Cinv)
#   x_t    = x_{t-1} + rho*(x_{t-1} - x_{t-2}) + N(0, sig2 * Cinv)
#   light  | x : normalised spike-and-slab at the knot's zenith
#
# What this can and cannot show: SBC validates the sampler against its own
# model. It cannot detect a likelihood that is wrong about the world -- for
# that see the real-data interval coverage in analysis/. Both are needed.
#
# Usage:  Rscript -e 'REPS <- 100L; source("inst/sbc/sbc_crw.R")'
# =====================================================================
suppressPackageStartupMessages(library(invTwilightFree))
bands <- system.file("sbc", "ecdf_bands.R", package = "invTwilightFree")
if (!nzchar(bands)) bands <- file.path("inst", "sbc", "ecdf_bands.R")
source(bands)

if (!exists("REPS"))   REPS   <- 100L
if (!exists("SEED"))   SEED   <- 1L
if (!exists("NTAGS"))  NTAGS  <- 3L      # small panel per replicate
if (!exists("NDAYS"))  NDAYS  <- 14
# rho mixes slowly because it is coupled to the latent track: changing it
# changes the prior's shape, which changes the track, which changes rho.
# Measured on a single dataset, 1500 sweeps thinned by 3 gave 334 draws with an
# effective sample size of only 49. SBC ranks a truth against those draws, so
# residual autocorrelation shows up as a rank distribution with heavy tails and
# reads as a calibration failure. Sweep long and thin hard enough that the kept
# draws are close to independent.
if (!exists("SWEEPS")) SWEEPS <- 6000L
if (!exists("BURN"))   BURN   <- 2000L
if (!exists("THIN"))   THIN   <- 12L
if (!exists("OUT"))    OUT    <- "sbc_crw.png"

DEPLOY_LON <- 150; DEPLOY_LAT <- -50
DATE   <- "2024-06-15"          # austral winter: latitude is well identified,
STEP_H <- 12                    # so the movement parameters are identifiable
CAL    <- c(558.5, 5.818)
LP     <- c(0.03, 64, 0.10)
A0     <- 5; B0 <- 3200         # sig2 ~ InvGamma(5, 3200): mean 800 km^2
if (!exists("RHO_MAX")) RHO_MAX <- 0.8
RHO_LO <- -RHO_MAX; RHO_HI <- RHO_MAX
# Persistence approaching 1 is near-non-stationary: displacement grows like
# n^1.5 rather than sqrt(n), so simulated tracks wander far outside the proposal
# surrogate's mesh and cannot be fitted. That shows up as an SBC failure of the
# sampler when it is really a failure of the test regime. At 0.8 the 90th
# percentile excursion is about 14 degrees, against a 30 degree mesh.

km_lon <- 111.0 * cos(DEPLOY_LAT * pi / 180); km_lat <- 111.0
CINV   <- c(km_lon^2, 0, 0, km_lat^2)       # movement metric, km^2 per degree^2

# ---- simulate one tag from the CRW prior ----------------------------
sim_tag <- function(sig2, rho, ndays, t0) {
  k  <- as.integer(ndays * 24 / STEP_H) + 1L
  kt <- as.numeric(t0) + (0:(k - 1)) * STEP_H * 3600
  sdl <- sqrt(sig2) / km_lon      # per-step sd in degrees
  sdt <- sqrt(sig2) / km_lat
  lon <- numeric(k); lat <- numeric(k)
  lon[1] <- DEPLOY_LON; lat[1] <- DEPLOY_LAT
  lon[2] <- lon[1] + rnorm(1, 0, sdl); lat[2] <- lat[1] + rnorm(1, 0, sdt)
  for (t in 3:k) {
    lon[t] <- lon[t-1] + rho * (lon[t-1] - lon[t-2]) + rnorm(1, 0, sdl)
    lat[t] <- lat[t-1] + rho * (lat[t-1] - lat[t-2]) + rnorm(1, 0, sdt)
  }
  # observations: piecewise-constant location within a knot, as the model assumes
  per <- as.integer(STEP_H * 60 / 10)          # 10-minute light
  ot <- numeric(0); ol <- numeric(0)
  for (t in 2:k) {
    tt <- kt[t-1] + (1:per) * 600
    tt <- tt[tt <= kt[t]]
    if (!length(tt)) next
    z <- solar_zenith(tt, rep(lon[t], length(tt)), rep(lat[t], length(tt)))
    mu <- pmin(pmax(CAL[1] - CAL[2] * z, 0), LP[2])
    # draw from the spike-and-slab the engine assumes
    u <- runif(length(mu))
    slab <- u < LP[3]
    y <- numeric(length(mu))
    y[slab] <- runif(sum(slab), 0, LP[2])
    ns <- !slab
    if (any(ns)) {
      # Asymmetric exponential about mu, TRUNCATED to [0, max_light] by inverse
      # CDF. Drawing an untruncated exponential and clamping would pile point
      # masses at 0 and max_light that the model's density does not have, and
      # SBC is sensitive to exactly that kind of generator/model mismatch.
      lam <- LP[1]; mx <- LP[2]
      m <- mu[ns]
      lo_mass <- 1 - exp(-lam * m)                       # arm below mu
      hi_mass <- 0.5 * (1 - exp(-2 * lam * (mx - m)))    # arm above mu
      pick_hi <- runif(length(m)) < hi_mass / (lo_mass + hi_mass)
      u <- runif(length(m))
      y_lo <- m + log(1 - u * (1 - exp(-lam * m))) / lam
      y_hi <- m - log(1 - u * (1 - exp(-2 * lam * (mx - m)))) / (2 * lam)
      y[ns] <- ifelse(pick_hi, y_hi, y_lo)
    }
    ot <- c(ot, tt); ol <- c(ol, y)
  }
  list(df = data.frame(Date = as.POSIXct(ot, origin = "1970-01-01", tz = "UTC"),
                       Light = ol),
       lon = lon, lat = lat, k = k)
}

set.seed(SEED)
t0 <- as.POSIXct(DATE, tz = "UTC")
rank_rho <- integer(0); rank_sig <- integer(0); L <- NA_integer_

for (rep in seq_len(REPS)) {
  sig2_true <- 1 / rgamma(NTAGS, shape = A0, rate = B0)
  rho_true  <- runif(NTAGS, RHO_LO, RHO_HI)
  tags <- lapply(seq_len(NTAGS), function(i) sim_tag(sig2_true[i], rho_true[i], NDAYS, t0))
  data <- stats::setNames(lapply(tags, `[[`, "df"), paste0("T", seq_len(NTAGS)))
  locs <- data.frame(id = names(data),
                     deploy_lon = DEPLOY_LON, deploy_lat = DEPLOY_LAT,
                     retrieve_lon = vapply(tags, function(g) g$lon[g$k], 0),
                     retrieve_lat = vapply(tags, function(g) g$lat[g$k], 0))
  fit <- try(TwilightFreeHier(data, locs, step_hours = STEP_H,
                              calibration = CAL, likelihood_params = LP,
                              a_pop = A0,
                              # pin beta at B0: with these the hyperprior has sd
                              # ~0.3 about B0, so each tag's prior is the
                              # InvGamma(A0, B0) the truth was drawn from
                              hyperprior = c(1e8, 1e8 / B0),
                              surrogate_diffusion = 60, mesh_pad = 30, coarse_res = 1.5,
                              sweeps = SWEEPS, burn = BURN, thin = THIN,
                              metric = "flat", movement = "crw",
                              rho_prior_sd = 1e6, rho_max = RHO_MAX,   # flat on (-RHO_MAX, RHO_MAX)
                              seed = rep), silent = TRUE)
  if (inherits(fit, "try-error")) { message("rep ", rep, " failed"); next }

  nk <- fit$draws$rho
  nkeep <- length(nk) / NTAGS
  R <- matrix(nk, nrow = nkeep, ncol = NTAGS, byrow = TRUE)
  S <- matrix(fit$draws$sig2_kmday, nrow = nkeep, ncol = NTAGS)^2 * (STEP_H / 24)
  if (is.na(L)) L <- nkeep
  for (i in seq_len(NTAGS)) {
    rank_rho <- c(rank_rho, sum(R[, i] < rho_true[i]))
    rank_sig <- c(rank_sig, sum(S[, i] < sig2_true[i]))
  }
  if (rep %% 10 == 0) message("rep ", rep, "/", REPS)
}

cat(sprintf("%sSBC complete: %d ranks, %d posterior draws per rank%s",
            "
", length(rank_rho), L, "
"))
df <- data.frame(rho = rank_rho, sig2 = rank_sig)
out <- sbc_ecdf_report(df,
                       list(rho = "rho (persistence)", sig2 = "sigma^2 (movement)"),
                       L, conf = 0.95, fig = OUT, tag = "SBC CRW")
saveRDS(list(ranks = df, L = L, report = out, reps = REPS, ntags = NTAGS),
        sub("[.]png$", ".rds", OUT))
invisible(out)
