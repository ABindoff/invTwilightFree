# =====================================================================
# SBC Stage B: calibration of the FULL hierarchy (population level beta +
# per-individual sig2_i with pooling). Covers the conjugate beta update and
# the shared-beta coupling that Stage A (fixed prior) does not.
#
# Generator matches the inference model exactly:
#   beta   ~ Gamma(g0, h0)
#   sig2_i ~ InvGamma(a_pop, beta)     i = 1..N
#   track_i, light_i  as in Stage A (knot RW + normalised spike-slab)
# Fit: pooled hierarchical block+polish sampler; rank true beta and true sig2_i.
# =====================================================================
suppressPackageStartupMessages(library(invTwilightFree))

DEPLOY_LON <- 150; DEPLOY_LAT <- -50
DATE   <- "2024-06-15"; NDAYS <- 8; STEP_H <- 12
CAL    <- c(558.5, 5.818); LP <- c(0.03, 64, 0.10)
A_POP  <- 3; G0 <- 4; H0 <- 0.0025      # beta ~ Gamma(4, 0.0025): mean 1600 -> sig2 mean ~800
N_IND  <- 4

CFG <- list(lib_only = TRUE, step_hours = STEP_H,
            calibration = CAL, likpar = LP, lat_ref = DEPLOY_LAT,
            a_pop = A_POP, g0 = G0, h0 = H0, coarse_res = 1.0, mesh_pad = 15.0,
            inflate = 1.5, surrogate_diffusion = 45.0, block_len = 5L, polish = TRUE)
source("da_hier.R")
source("../../inst/sbc/ecdf_bands.R")   # reuse the package's SBC ECDF-band machinery

times <- seq(as.POSIXct(DATE, tz = "UTC"), by = "10 min", length.out = NDAYS * 24 * 6)
ut <- as.numeric(times); t0 <- ut[1]; t1 <- ut[length(ut)]
K  <- ceiling((t1 - t0) / (STEP_H * 3600)) + 1; tstep <- (t1 - t0) / (K - 1)
KNOTS <- t0 + (0:(K - 1)) * tstep
kidx <- pmin(pmax(findInterval(ut, KNOTS, left.open = TRUE) + 1L, 1L), K)
km_lat <- 111.0; km_lon <- 111.0 * cos(DEPLOY_LAT * pi / 180)

sample_spike <- function(mu) {
  lam <- LP[1]; maxl <- LP[2]
  m_lo <- 1 - exp(-lam * mu); m_hi <- 0.5 * (1 - exp(-2 * lam * (maxl - mu)))
  if (runif(1) < m_lo / (m_lo + m_hi)) mu - (-log(1 - runif(1) * m_lo) / lam)
  else                                 mu + (-log(1 - runif(1) * (2 * m_hi)) / (2 * lam))
}
gen_track_light <- function(sig2) {
  s_lon <- sqrt(sig2) / km_lon; s_lat <- sqrt(sig2) / km_lat
  lon_k <- numeric(K); lat_k <- numeric(K); lon_k[1] <- DEPLOY_LON; lat_k[1] <- DEPLOY_LAT
  for (k in 2:K) { lon_k[k] <- lon_k[k-1] + rnorm(1, 0, s_lon); lat_k[k] <- lat_k[k-1] + rnorm(1, 0, s_lat) }
  lon_o <- lon_k[kidx]; lat_o <- lat_k[kidx]
  z <- solar_zenith(ut, lon_o, lat_o); mu <- pmin(pmax(CAL[1] - CAL[2] * z, 0), LP[2])
  light <- vapply(seq_along(ut), function(j)
    if (runif(1) < LP[3]) runif(1, 0, LP[2]) else sample_spike(mu[j]), numeric(1))
  data.frame(time = times, light = light, true_lat = lat_o, true_lon = lon_o)
}

CTR$calls <- 0L; CTR$emit <- 0L
fit_hier <- function(IND, sweeps, burn, thin) {
  N <- length(IND)
  X <- lapply(IND, function(E) E$init())
  beta <- G0 / H0                                  # start at the prior mean (neutral)
  sig2 <- rep(beta / (A_POP - 1), N)               # prior-mean sig2, avoids tiny-init transient
  kb <- numeric(0); ks <- matrix(0, 0, N)
  for (s in seq_len(sweeps)) {
    for (i in seq_len(N)) {
      X[[i]] <- update_block(IND[[i]], X[[i]], sig2[i], cfg$block_len)
      SS <- ss_track(IND[[i]], X[[i]])
      sig2[i] <- 1 / rgamma(1, shape = A_POP + (IND[[i]]$K - 1), rate = beta + 0.5 * SS)
    }
    beta <- rgamma(1, shape = G0 + N * A_POP, rate = H0 + sum(1 / sig2))
    if (s > burn && (s - burn) %% thin == 0) { kb <- c(kb, beta); ks <- rbind(ks, sig2) }
  }
  list(beta = kb, sig = ks)
}

REPS  <- if (exists("REPS_B")) REPS_B else 60L
SWEEPS <- 4000L; BURN <- 1500L; THIN <- 25L    # L = (4000-1500)/25 = 100
rank_beta <- integer(0); rank_sig <- integer(0); t_start <- Sys.time()
for (r in seq_len(REPS)) {
  beta_t <- rgamma(1, G0, H0)
  sig2_t <- 1 / rgamma(N_IND, A_POP, beta_t)
  IND <- lapply(seq_len(N_IND), function(i)
    make_individual(gen_track_light(sig2_t[i]), DEPLOY_LON, DEPLOY_LAT))
  fit <- fit_hier(IND, SWEEPS, BURN, THIN)
  rank_beta <- c(rank_beta, sum(fit$beta < beta_t))
  for (i in seq_len(N_IND)) rank_sig <- c(rank_sig, sum(fit$sig[, i] < sig2_t[i]))
  if (r %% 5L == 0L) cat(sprintf("  SBC-B %d/%d (%.0fs)\n", r, REPS,
                                 as.numeric(Sys.time() - t_start, units = "secs")))
}
L <- (SWEEPS - BURN) %/% THIN
cat(sprintf("\nStage B SBC: %d replicates, L=%d ; %d beta ranks, %d sig2 ranks\n",
            REPS, L, length(rank_beta), length(rank_sig)))
sbc_ecdf_report(data.frame(rank_beta = rank_beta), c(rank_beta = "beta (population scale)"),
                L, conf = 0.95, fig = "sbc_stageB_beta.png", tag = "hierarchy | beta")
sbc_ecdf_report(data.frame(rank_sig2 = rank_sig), c(rank_sig2 = "sigma^2_i (pooled)"),
                L, conf = 0.95, fig = "sbc_stageB_sig.png", tag = "hierarchy | sig2_i")
cat("RANK_BETA:", paste(rank_beta, collapse = ","), "\n")
