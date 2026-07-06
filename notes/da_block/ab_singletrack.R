# A-B: single-track block+polish sampler, R prototype vs native Rust kernel.
# Same track, same surrogate, same movement precision, same sweeps -> the two
# posteriors must coincide (RNG streams differ) and we compare wall time.
suppressPackageStartupMessages(devtools::load_all(
  "C:/Users/bindoffa/antigravity_projects/invtwilightfree", quiet = TRUE))

CFG <- list(lib_only = TRUE, step_hours = 12, coarse_res = 1.0, mesh_pad = 15.0,
            inflate = 1.5, surrogate_diffusion = 45.0, block_len = 5L, polish = TRUE)
source("da_hier.R")   # make_individual, update_block, ss_track, ... (lib_only: no run)

SWEEPS <- if (exists("SW")) SW else 6000L; BURN <- 2000L; THIN <- 3L
sig2 <- 1000                                    # fixed per-knot km-variance

PAN <- simulate_panel(N = 1, days = 14, step_min = 10, pop_km_per_day = 45)
df <- PAN$ind[[1]]; start <- c(PAN$start$lon[1], PAN$start$lat[1])
E <- make_individual(df, start[1], start[2]); K <- E$K
Pmove <- E$Cinv / sig2
cat(sprintf("track: %d knots, %d obs ; sig2=%.0f\n", K, length(E$S$ot), sig2))

# ---- flatten surrogate + obs ranges for Rust ----
mu_flat <- as.numeric(t(E$sur$mu)); mu_flat[is.na(mu_flat)] <- 0
P_flat <- numeric(4 * K)
for (k in seq_len(K)) if (E$sur$info[k]) P_flat[(4*k-3):(4*k)] <- as.numeric(t(E$sur$P[[k]]))
obs_start <- integer(K); obs_len <- integer(K)
for (k in seq_len(K)) { j <- E$S$obin[[k]]
  if (length(j)) { obs_start[k] <- min(j) - 1L; obs_len[k] <- length(j) } }

# ---- B: Rust kernel ----
t0 <- Sys.time()
rust <- run_block_track(obs_start, obs_len, E$S$ot, E$S$ol, E$S$calibration, E$S$likpar,
                        mu_flat, P_flat, as.numeric(t(Pmove)), start[1], start[2],
                        cfg$block_len, SWEEPS, BURN, THIN, TRUE, 42)
t_rust <- as.numeric(Sys.time() - t0, units = "secs")

# ---- A: R prototype ----
set.seed(42)
t0 <- Sys.time()
x <- E$init(); nkeep <- (SWEEPS - BURN) %/% THIN
keep <- array(NA_real_, c(nkeep, K, 2)); s <- 0L
for (sw in seq_len(SWEEPS)) {
  x <- update_block(E, x, sig2, cfg$block_len)
  if (sw > BURN && (sw - BURN) %% THIN == 0) { s <- s + 1L; keep[s, , ] <- x }
}
t_r <- as.numeric(Sys.time() - t0, units = "secs")
r_mlon <- apply(keep[,,1], 2, mean); r_slon <- apply(keep[,,1], 2, sd)
r_mlat <- apply(keep[,,2], 2, mean); r_slat <- apply(keep[,,2], 2, sd)

# ---- compare posteriors (validation) ----
zmean <- function(dm, ds1, ds2) abs(dm) / ((ds1 + ds2)/2 + 1e-9)
zlon <- zmean(r_mlon - rust$mean_lon, r_slon, rust$sd_lon)
zlat <- zmean(r_mlat - rust$mean_lat, r_slat, rust$sd_lat)
cat("\n================= VALIDATION: Rust vs R posterior =================\n")
cat(sprintf("mean gap (deg)   lon: max=%.3f mean=%.3f | lat: max=%.3f mean=%.3f\n",
    max(abs(r_mlon-rust$mean_lon)), mean(abs(r_mlon-rust$mean_lon)),
    max(abs(r_mlat-rust$mean_lat)), mean(abs(r_mlat-rust$mean_lat))))
cat(sprintf("mean gap |dmu|/sd lon: max=%.3f mean=%.3f | lat: max=%.3f mean=%.3f\n",
    max(zlon[-1]), mean(zlon[-1]), max(zlat[-1]), mean(zlat[-1])))
cat(sprintf("sd ratio Rust/R  lon: [%.2f, %.2f] | lat: [%.2f, %.2f]\n",
    min((rust$sd_lon/r_slon)[-1]), max((rust$sd_lon/r_slon)[-1]),
    min((rust$sd_lat/r_slat)[-1]), max((rust$sd_lat/r_slat)[-1])))
cat(sprintf("acceptance: Rust block=%.2f\n", rust$accept))

cat("\n================= A-B: wall time =================\n")
cat(sprintf("sweeps=%d ; R prototype = %.2fs | Rust kernel = %.3fs | speedup = %.0fx\n",
    SWEEPS, t_r, t_rust, t_r / t_rust))
