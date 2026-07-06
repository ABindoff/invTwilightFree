# A-B: hierarchical partial-pooling sampler, R prototype vs native Rust.
# Same panel, same surrogates, same hyperprior, same sweeps -> the beta and
# sig2_i posteriors must coincide (RNG differs); compare wall time.
suppressPackageStartupMessages(devtools::load_all(
  "C:/Users/bindoffa/antigravity_projects/invtwilightfree", quiet = TRUE))

CFG <- list(lib_only = TRUE, step_hours = 12, coarse_res = 1.0, mesh_pad = 12.0,
            inflate = 1.5, surrogate_diffusion = 45.0, block_len = 5L, polish = TRUE,
            a_pop = 3.0, g0 = 1e-3, h0 = 1e-3)
source("da_hier.R")

N <- 6L; SWEEPS <- if (exists("SW")) SW else 4000L; BURN <- 1500L; THIN <- 5L
PAN <- simulate_panel(N = N, days = 14, step_min = 10, pop_km_per_day = 45)
IND <- lapply(seq_len(N), function(i) make_individual(PAN$ind[[i]], PAN$start$lon[i], PAN$start$lat[i]))
dt_day <- IND[[1]]$S$tstep / 86400
cat(sprintf("panel: N=%d, knots=%s\n", N, paste(sapply(IND, function(E) E$K), collapse=",")))

# ---- flatten all individuals into global arrays for Rust ----
knots_per_ind <- as.integer(sapply(IND, function(E) E$K))
obs_times <- numeric(0); obs_light <- numeric(0)
kob_start <- integer(0); kob_len <- integer(0)
smu <- numeric(0); sp <- numeric(0); cal <- numeric(0); lp <- numeric(0)
cinv <- numeric(0); slon <- numeric(0); slat <- numeric(0)
obs_off <- 0L
for (i in seq_len(N)) {
  E <- IND[[i]]; K <- E$K
  obs_times <- c(obs_times, E$S$ot); obs_light <- c(obs_light, E$S$ol)
  for (k in seq_len(K)) {
    j <- E$S$obin[[k]]
    if (length(j)) { kob_start <- c(kob_start, obs_off + min(j) - 1L); kob_len <- c(kob_len, length(j)) }
    else           { kob_start <- c(kob_start, obs_off); kob_len <- c(kob_len, 0L) }
    if (E$sur$info[k]) { smu <- c(smu, E$sur$mu[k, ]); sp <- c(sp, as.numeric(t(E$sur$P[[k]]))) }
    else               { smu <- c(smu, c(0, 0)); sp <- c(sp, c(0, 0, 0, 0)) }
  }
  cal <- c(cal, E$S$calibration); lp <- c(lp, E$S$likpar)
  cinv <- c(cinv, as.numeric(t(E$Cinv)))
  slon <- c(slon, E$start[1]); slat <- c(slat, E$start[2])
  obs_off <- obs_off + length(E$S$ot)
}

# ---- B: Rust native hierarchy ----
t0 <- Sys.time()
rust <- run_block_hier(N, knots_per_ind, kob_start, kob_len, obs_times, obs_light,
                       cal, lp, smu, sp, cinv, slon, slat,
                       cfg$a_pop, cfg$g0, cfg$h0, cfg$block_len, SWEEPS, BURN, THIN, TRUE, 42)
t_rust <- as.numeric(Sys.time() - t0, units = "secs")
nk <- rust$n_kept
rust_sig <- matrix(rust$sig2, nrow = nk, ncol = N, byrow = TRUE)

# ---- A: R prototype hierarchy (same logic as run_hier, block mode, pooled) ----
set.seed(42)
t0 <- Sys.time()
X <- lapply(IND, function(E) E$init())
sig2R <- sapply(IND, function(E) ss_track(E, E$init()) / (E$K - 1)); betaR <- 1
kb <- numeric(0); ksg <- matrix(0, 0, N)
for (s in seq_len(SWEEPS)) {
  for (i in seq_len(N)) {
    X[[i]] <- update_block(IND[[i]], X[[i]], sig2R[i], cfg$block_len)
    SS <- ss_track(IND[[i]], X[[i]])
    sig2R[i] <- 1 / rgamma(1, shape = cfg$a_pop + (IND[[i]]$K - 1), rate = betaR + 0.5*SS)
  }
  betaR <- rgamma(1, shape = cfg$g0 + N*cfg$a_pop, rate = cfg$h0 + sum(1/sig2R))
  if (s > BURN && (s - BURN) %% THIN == 0) { kb <- c(kb, betaR); ksg <- rbind(ksg, sig2R) }
}
t_r <- as.numeric(Sys.time() - t0, units = "secs")

# ---- compare posteriors (km/day) ----
kmday <- function(v) sqrt(v) / sqrt(dt_day)
popR <- kmday(kb) / sqrt(cfg$a_pop - 1) * sqrt(cfg$a_pop - 1)  # keep simple: pop scale
popR <- sqrt(kb/(cfg$a_pop-1)) / sqrt(dt_day)
popU <- sqrt(rust$beta/(cfg$a_pop-1)) / sqrt(dt_day)
cat("\n================= VALIDATION: Rust vs R hierarchy =================\n")
cat(sprintf("pop scale km/day:  R %.1f [%.1f,%.1f] | Rust %.1f [%.1f,%.1f]\n",
    mean(popR), quantile(popR,.025), quantile(popR,.975),
    mean(popU), quantile(popU,.025), quantile(popU,.975)))
for (i in seq_len(N)) {
  sR <- kmday(ksg[,i]); sU <- kmday(rust_sig[,i])
  dz <- abs(mean(sR)-mean(sU))/((sd(sR)+sd(sU))/2+1e-9)
  cat(sprintf("sig ind %d km/day: R %.1f [%.1f,%.1f] | Rust %.1f [%.1f,%.1f] | dmu/sd %.2f | truth %.1f\n",
      i, mean(sR), quantile(sR,.025), quantile(sR,.975),
      mean(sU), quantile(sU,.025), quantile(sU,.975), dz, PAN$sd_day_true[i]))
}
cat("\n================= A-B: wall time =================\n")
cat(sprintf("N=%d, sweeps=%d ; R = %.1fs | Rust = %.3fs | speedup = %.0fx\n",
    N, SWEEPS, t_r, t_rust, t_r / t_rust))
