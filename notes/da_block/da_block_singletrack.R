# =====================================================================
# Single-track surrogate-posterior BLOCK / delayed-acceptance sampler
# for invTwilightFree, validated against a gold always-exact single-site
# sampler. Milestone 1 (exactness) on the REAL twilight-free likelihood.
#
# ll_exact  := eval_logpk_grid()  (the actual spike-and-slab light model)
# ll_cheap  := a per-knot Gaussian surrogate, moment-matched off a mesh
#              of ll_exact. Broad-in-latitude BY CONSTRUCTION: matching a
#              bimodal (equinox) surface yields a Gaussian spanning both
#              hemispheres, so no mode is silently dropped.
#
# The block is drawn from (RW movement prior x surrogate Gaussians) with
# the block endpoints held fixed, so that Gaussian's normaliser cancels
# and the MH correction collapses to sum(ll_exact - ll_cheap) over the
# block. Change the surrogate (or use a MAP draft) and this breaks.
# =====================================================================

suppressPackageStartupMessages(library(invTwilightFree))

# ---------------------------------------------------------------------
# 0. Config (overridable before source-ing)
# ---------------------------------------------------------------------
if (!exists("CFG")) CFG <- list()
cfg <- utils::modifyList(list(
  dataset    = "sim_short",
  step_hours = 12.0,
  diffusion  = 50.0,      # km / sqrt(day)
  block_len  = 5L,
  inflate    = 1.5,       # surrogate variance inflation (>=1; exactness-neutral)
  sweeps     = 4000L,
  burn       = 1500L,
  thin       = 3L,
  seed       = 42L,
  surrogate  = "moment_match",  # "moment_match" | "coarse_hmm"
  coarse_res = 1.5,       # coarse grid-HMM resolution (deg) for the surrogate
  fallback   = FALSE,     # single-site update on coarse-HMM-multimodal knots
  mesh_res   = 0.75,      # moment-match mesh resolution (deg)
  mesh_lon_pad = 15.0,
  mesh_lat_pad = 22.0     # generous in lat to admit a mirror hemisphere
), CFG)

set.seed(cfg$seed)

# ---------------------------------------------------------------------
# 1. Load a simulated track and replicate TwilightFreeGrid's setup
#    (calibration, light shift, knot construction, per-knot obs bins)
# ---------------------------------------------------------------------
e <- new.env(); utils::data(list = cfg$dataset, package = "invTwilightFree", envir = e)
df <- e[[cfg$dataset]]
df <- df[order(df$time), ]

start_lon <- df$true_lon[1]; start_lat <- df$true_lat[1]

setup_track <- function(df, start_lon, start_lat, step_hours) {
  light <- df$light; date_time <- df$time
  # ---- auto-calibration, copied from TwilightFreeGrid() ----
  min_l <- stats::quantile(light, 0.05, na.rm = TRUE)
  max_l <- stats::quantile(light, 0.95, na.rm = TRUE)
  light_shifted <- pmax(0, light - min_l)
  max_shifted   <- as.numeric(max_l - min_l)

  cal_idx   <- date_time < (date_time[1] + 3 * 24 * 3600)
  cal_times <- date_time[cal_idx]
  cal_light <- light_shifted[cal_idx]
  cal_zenith <- solar_zenith(as.numeric(cal_times),
                             rep(start_lon, length(cal_times)),
                             rep(start_lat, length(cal_times)))
  trans_idx <- which(cal_light > 0 & cal_light < max_shifted * 0.95 &
                     cal_zenith > 85 & cal_zenith < 100)
  if (length(trans_idx) > 10) {
    fit_cal <- stats::lm(cal_light[trans_idx] ~ cal_zenith[trans_idx])
    intercept_est <- stats::coef(fit_cal)[1]
    slope_est <- -stats::coef(fit_cal)[2]
    if (is.na(slope_est) || slope_est <= 0) slope_est <- max_shifted / (96 - 85)
  } else {
    slope_est <- max_shifted / (96 - 85); intercept_est <- slope_est * 96
  }
  calibration <- as.numeric(c(intercept_est, slope_est))
  likelihood_params <- as.numeric(c(1.0 / (max_shifted * 0.5), max_shifted, 0.10))

  unix_times <- as.numeric(date_time)
  t_start <- unix_times[1]; t_end <- unix_times[length(unix_times)]
  k_steps <- ceiling((t_end - t_start) / (step_hours * 3600)) + 1
  t_step  <- if (k_steps > 1) (t_end - t_start) / (k_steps - 1) else 0
  knot_times <- t_start + (0:(k_steps - 1)) * t_step

  valid <- !is.na(light_shifted)
  ot <- unix_times[valid]; ol <- light_shifted[valid]

  # per-knot obs bins matching run_grid_hmm: (t_prev, t_curr]
  obs_bin <- vector("list", k_steps)
  for (k in seq_len(k_steps)) {
    t_curr <- knot_times[k]
    t_prev <- if (k == 1) t_curr - (knot_times[2] - knot_times[1]) else knot_times[k - 1]
    obs_bin[[k]] <- which(ot > t_prev & ot <= t_curr)
  }
  list(calibration = calibration, likelihood_params = likelihood_params,
       knot_times = knot_times, k_steps = k_steps, t_step = t_step,
       obs_times = ot, obs_light = ol, obs_bin = obs_bin)
}

S <- setup_track(df, start_lon, start_lat, cfg$step_hours)
K <- S$k_steps
cat(sprintf("dataset=%s  knots=%d  step=%.0fh  obs=%d  cal=(%.2f,%.3f)  lik=(%.4f,%.1f,%.2f)\n",
            cfg$dataset, K, cfg$step_hours, length(S$obs_times),
            S$calibration[1], S$calibration[2],
            S$likelihood_params[1], S$likelihood_params[2], S$likelihood_params[3]))

# instrumented exact-likelihood counters. `calls` = eval_logpk_grid FFI calls
# (gold + block single-site/fallback); `emit` = solar_zenith FFI calls (batched
# block correction). `pts` = exact point-evals.
CTR <- new.env(); CTR$pts <- 0L; CTR$calls <- 0L; CTR$emit <- 0L

# ll_exact at knot k for vectors of (lon,lat): the REAL likelihood.
ll_exact_knot <- function(k, lon, lat) {
  j <- S$obs_bin[[k]]
  if (length(j) == 0L) return(rep(0, length(lon)))   # no info this knot
  CTR$pts   <- CTR$pts + length(lon)
  CTR$calls <- CTR$calls + 1L
  eval_logpk_grid(as.numeric(lon), as.numeric(lat),
                  S$obs_times[j], S$obs_light[j],
                  S$calibration, S$likelihood_params)
}

# scalar likelihood params for the batched pure-R emitter (M3)
intercept <- S$calibration[1]; slope <- S$calibration[2]
lambda <- S$likelihood_params[1]; max_light <- S$likelihood_params[2]
prob_slab <- if (length(S$likelihood_params) > 3)
  S$likelihood_params[3] / (S$likelihood_params[3] + S$likelihood_params[4]) else S$likelihood_params[3]

# Batched exact emitter: per-(knot,point) light log-likelihood in ONE solar_zenith
# call over all (point, its-own-obs) pairs. Replicates the spike-and-slab model of
# eval_logpk_grid to ~1e-14 (validated in validate_emit.R), so the block stays
# exact while collapsing B per-knot FFI calls into one. kv/lonv/latv are parallel
# length-m vectors: point i sits at knot kv[i].
emit_block_ll <- function(kv, lonv, latv) {
  m <- length(kv)
  counts <- lengths(S$obs_bin[kv])
  ll <- numeric(m)
  if (sum(counts) == 0L) return(ll)
  pt <- rep(seq_len(m), counts)
  oj <- unlist(S$obs_bin[kv], use.names = FALSE)
  z <- solar_zenith(S$obs_times[oj], lonv[pt], latv[pt]); CTR$emit <- CTR$emit + 1L
  CTR$pts <- CTR$pts + m
  expc <- pmin(pmax(intercept - slope * z, 0), max_light)
  ol <- S$obs_light[oj]
  raw <- ifelse(ol <= expc, lambda * exp(-lambda * (expc - ol)),
                            lambda * exp(-lambda * 2 * (ol - expc)))
  norm <- pmax(1 - exp(-lambda * expc) + 0.5 * (1 - exp(-2 * lambda * (max_light - expc))), 1e-12)
  den <- (1 - prob_slab) * (raw / norm) + prob_slab * (1 / max_light)
  rs <- rowsum(log(den), pt)
  ll[as.integer(rownames(rs))] <- rs[, 1]
  ll
}

# ---------------------------------------------------------------------
# 2. Movement model: isotropic RW (km) mapped to lon/lat degrees.
#    Same prior in gold and block, so exactness is metric-independent.
# ---------------------------------------------------------------------
lat_ref <- mean(df$true_lat)
km_per_deg_lat <- 111.0
km_per_deg_lon <- 111.0 * cos(lat_ref * pi / 180)
dt_days <- S$t_step / 86400
var_km  <- cfg$diffusion^2 * dt_days                 # per-step variance (km^2)
Sigma_move <- diag(c(var_km / km_per_deg_lon^2, var_km / km_per_deg_lat^2))
P_move     <- solve(Sigma_move)                       # per-step movement precision

# ---------------------------------------------------------------------
# 3. Surrogate per-knot Gaussian (mu_k, Sigma_k). Two constructions:
#    "moment_match" : moment-match the raw per-knot ll_exact (misleading in
#                     latitude; kept for comparison / M2 reproduction).
#    "coarse_hmm"   : per-knot marginal of a cheap COARSE grid-HMM forward-
#                     backward, which fuses accumulated light + movement +
#                     the fixed start, so the surrogate is correctly located.
#    Exactness is identical either way: the SAME Gaussian is used to build
#    the block and inside the sum(ll_exact - ll_cheap) correction.
# ---------------------------------------------------------------------
gauss_from_weights <- function(w, lon, lat, inflate) {
  w <- w / sum(w)
  mu <- c(sum(w * lon), sum(w * lat))
  dx <- lon - mu[1]; dy <- lat - mu[2]
  Sig <- matrix(c(sum(w*dx*dx), sum(w*dx*dy), sum(w*dx*dy), sum(w*dy*dy)), 2, 2) * inflate + diag(1e-6, 2)
  list(mu = mu, S = Sig, P = solve(Sig))
}

# simple bimodality flag on a latitude marginal: a second lat cluster with
# > frac of the peak mass, separated from the peak by a trough < 0.5*secondary.
lat_multimodal <- function(w, lat, frac = 0.25, gap = 8) {
  m <- tapply(w, round(lat), sum); m <- m / max(m)
  lv <- as.numeric(names(m))
  pk <- lv[which.max(m)]
  cand <- lv[m > frac & abs(lv - pk) > gap]
  if (!length(cand)) return(FALSE)
  # check a genuine trough between pk and the farthest candidate
  far <- cand[which.max(abs(cand - pk))]
  between <- m[lv >= min(pk, far) & lv <= max(pk, far)]
  min(between) < 0.5 * min(m[as.character(pk)], m[as.character(far)])
}

surr_mu   <- matrix(NA_real_, K, 2)
surr_S    <- vector("list", K)
surr_P    <- vector("list", K)
surr_info <- logical(K)
surr_multi <- logical(K)   # coarse-HMM multimodality flag (for fallback)

build_moment_match <- function() {
  mesh_lon <- seq(min(df$true_lon) - cfg$mesh_lon_pad, max(df$true_lon) + cfg$mesh_lon_pad, by = cfg$mesh_res)
  mesh_lat <- seq(min(df$true_lat) - cfg$mesh_lat_pad, max(df$true_lat) + cfg$mesh_lat_pad, by = cfg$mesh_res)
  mg <- expand.grid(lon = mesh_lon, lat = mesh_lat)
  cat(sprintf("moment-match mesh: %d points\n", nrow(mg)))
  for (k in seq_len(K)) {
    if (length(S$obs_bin[[k]]) == 0L) { surr_P[[k]] <<- matrix(0,2,2); next }
    ll <- ll_exact_knot(k, mg$lon, mg$lat)
    g <- gauss_from_weights(exp(ll - max(ll)), mg$lon, mg$lat, cfg$inflate)
    surr_mu[k,] <<- g$mu; surr_S[[k]] <<- g$S; surr_P[[k]] <<- g$P; surr_info[k] <<- TRUE
  }
}

# Self-contained coarse grid-HMM forward-backward in pure R. Emissions come
# from eval_logpk_grid (the real light likelihood); the transition is the SAME
# RW metric (P_move) the block sampler uses. Knot 1 is pinned to the start cell.
# Avoids the installed run_grid_hmm (which segfaults in this build) and keeps the
# surrogate's movement model consistent with the sampler's prior.
build_coarse_hmm <- function() {
  clon <- seq(min(df$true_lon) - cfg$mesh_lon_pad, max(df$true_lon) + cfg$mesh_lon_pad, by = cfg$coarse_res)
  clat <- seq(min(df$true_lat) - cfg$mesh_lat_pad, max(df$true_lat) + cfg$mesh_lat_pad, by = cfg$coarse_res)
  cg <- expand.grid(lon = clon, lat = clat); n <- nrow(cg)
  cat(sprintf("coarse HMM grid: %d x %d = %d cells (res %.1f deg)\n", length(clon), length(clat), n, cfg$coarse_res))

  # emissions e[k, ]: normalised per-knot likelihood over cells (uniform if no obs)
  E <- matrix(1 / n, K, n)
  for (k in seq_len(K)) {
    if (length(S$obs_bin[[k]]) == 0L) next
    ll <- ll_exact_knot(k, cg$lon, cg$lat)
    e <- exp(ll - max(ll)); E[k, ] <- e / sum(e)
  }
  # homogeneous transition T[i, j] = p(i -> j) under the RW metric P_move
  A <- outer(cg$lon, cg$lon, function(li, lj) lj - li)   # A[i,j] = lon_j - lon_i
  B <- outer(cg$lat, cg$lat, function(li, lj) lj - li)   # B[i,j] = lat_j - lat_i
  quad <- P_move[1,1]*A*A + 2*P_move[1,2]*A*B + P_move[2,2]*B*B
  Tm <- exp(-0.5 * (quad - apply(quad, 1, min)))
  Tm <- Tm / rowSums(Tm)
  rm(A, B, quad)

  # knot 1 pinned to the nearest cell to the start
  start_cell <- which.min((cg$lon - start_lon)^2 + (cg$lat - start_lat)^2)
  # forward (with scaling)
  al <- matrix(0, K, n)
  a <- numeric(n); a[start_cell] <- 1; al[1, ] <- a
  for (k in 2:K) { a <- as.numeric(crossprod(Tm, a)) * E[k, ]; a <- a / sum(a); al[k, ] <- a }
  # backward
  be <- matrix(0, K, n); be[K, ] <- 1 / n
  for (k in (K-1):1) { b <- as.numeric(Tm %*% (be[k+1, ] * E[k+1, ])); be[k, ] <- b / sum(b) }
  # smoothed marginals -> per-knot Gaussian surrogate
  for (k in seq_len(K)) {
    w <- al[k, ] * be[k, ]; sw <- sum(w)
    if (!is.finite(sw) || sw <= 0) { surr_P[[k]] <<- matrix(0,2,2); next }
    w <- w / sw
    g <- gauss_from_weights(w, cg$lon, cg$lat, cfg$inflate)
    surr_mu[k,] <<- g$mu; surr_S[[k]] <<- g$S; surr_P[[k]] <<- g$P; surr_info[k] <<- TRUE
    surr_multi[k] <<- lat_multimodal(w, cg$lat)
  }
}

if (identical(cfg$surrogate, "coarse_hmm")) build_coarse_hmm() else build_moment_match()
# knots eligible for block updates; multimodal ones fall back to single-site
blockable <- surr_info & !(isTRUE(cfg$fallback) & surr_multi)
cat(sprintf("surrogate=%s ; informative %d/%d ; multimodal %d ; median surrogate sd (lon,lat)=(%.2f, %.2f)\n",
            if (is.null(cfg$surrogate)) "moment_match" else cfg$surrogate,
            sum(surr_info), K, sum(surr_multi),
            stats::median(sapply(which(surr_info), function(k) sqrt(surr_S[[k]][1,1]))),
            stats::median(sapply(which(surr_info), function(k) sqrt(surr_S[[k]][2,2])))))

# log Gaussian surrogate density (the SAME object used to build the block)
ll_cheap_knot <- function(k, lon, lat) {
  if (!surr_info[k]) return(rep(0, length(lon)))
  P <- surr_P[[k]]; mu <- surr_mu[k, ]
  dx <- lon - mu[1]; dy <- lat - mu[2]
  quad <- P[1,1]*dx*dx + 2*P[1,2]*dx*dy + P[2,2]*dy*dy
  -log(2*pi) - 0.5*log(det(surr_S[[k]])) - 0.5*quad
}

# ---------------------------------------------------------------------
# 4. Block draw from the surrogate's Gaussian block posterior
#    (RW bridge with fixed endpoints + per-knot Gaussian obs)
# ---------------------------------------------------------------------
block_draw <- function(idx, xL, xR) {
  B <- length(idx); d <- 2 * B
  Q <- matrix(0, d, d); b <- numeric(d)
  for (a in seq_len(B)) {
    ia <- (2*a - 1):(2*a)
    Q[ia, ia] <- Q[ia, ia] + 2 * P_move
    if (a > 1) { ip <- (2*a - 3):(2*a - 2); Q[ia, ip] <- Q[ia, ip] - P_move; Q[ip, ia] <- Q[ip, ia] - P_move }
    Pk <- surr_P[[idx[a]]]
    Q[ia, ia] <- Q[ia, ia] + Pk
    if (surr_info[idx[a]]) b[ia] <- b[ia] + Pk %*% surr_mu[idx[a], ]
  }
  b[1:2]       <- b[1:2]       + P_move %*% xL
  b[(d-1):d]   <- b[(d-1):d]   + P_move %*% xR
  R  <- chol(Q); mu <- solve(Q, b)
  v  <- as.numeric(mu + backsolve(R, rnorm(d)))
  matrix(v, B, 2, byrow = TRUE)
}

# ---------------------------------------------------------------------
# 5. Samplers. Knot 1 fixed at start. Knot K free (single-site).
# ---------------------------------------------------------------------
draw2 <- function(mean, prec) as.numeric(mean + t(chol(solve(prec))) %*% rnorm(2))

# gold single-site: propose from movement full-conditional, always exact
gold_site <- function(x, t) {
  if (t == K) { m <- x[t-1, ]; prec <- P_move }
  else        { m <- (x[t-1, ] + x[t+1, ]) / 2; prec <- 2 * P_move }
  prop <- draw2(m, prec); cur <- x[t, ]
  le <- ll_exact_knot(t, c(prop[1], cur[1]), c(prop[2], cur[2]))
  if (log(runif(1)) < le[1] - le[2]) x[t, ] <- prop
  x
}
# initial track. cfg$init in {"start","truth","surrogate"}:
#   start     - constant at the deployment fix (cold start; stress test)
#   truth     - the true track (isolates kernel correctness from reachability)
#   surrogate - the coarse-HMM smoother mean (the natural cheap initialiser)
init_track <- function() {
  mode <- if (!is.null(cfg$init)) cfg$init else if (isTRUE(cfg$init_truth)) "truth" else "start"
  x <- cbind(rep(start_lon, K), rep(start_lat, K))
  if (mode == "truth") {
    for (k in seq_len(K)) {
      idx <- which.min(abs(as.numeric(df$time) - S$knot_times[k]))
      x[k, ] <- c(df$true_lon[idx], df$true_lat[idx])
    }
  } else if (mode == "surrogate") {
    for (k in seq_len(K)) if (surr_info[k]) x[k, ] <- surr_mu[k, ]
  }
  x[1, ] <- c(start_lon, start_lat)
  x
}

run_gold <- function(sweeps = cfg$sweeps, burn = cfg$burn, thin = cfg$thin) {
  x <- init_track()
  keep <- array(NA_real_, c((sweeps - burn) %/% thin, K, 2)); s_keep <- 0L
  for (s in seq_len(sweeps)) {
    for (t in 2:K) x <- gold_site(x, t)
    if (s > burn && (s - burn) %% thin == 0) { s_keep <- s_keep + 1L; keep[s_keep, , ] <- x }
  }
  keep[seq_len(s_keep), , , drop = FALSE]
}

# block: interior 2..K-1 in surrogate-posterior blocks, knot K single-site
run_block <- function(sweeps = cfg$sweeps, burn = cfg$burn, thin = cfg$thin) {
  x <- init_track()
  keep <- array(NA_real_, c((sweeps - burn) %/% thin, K, 2)); s_keep <- 0L
  acc <- 0L; att <- 0L
  for (s in seq_len(sweeps)) {
    # random block phase each sweep so boundaries vary (decorrelates the partition)
    phase <- sample(0:(cfg$block_len - 1L), 1L)
    blk_max <- max(1L, cfg$block_len - phase)
    l <- 2L
    while (l <= K - 1L) {
      if (!blockable[l]) { x <- gold_site(x, l); l <- l + 1L; next }   # single-site fallback
      r <- l
      while (r < min(l + blk_max - 1L, K - 1L) && blockable[r + 1L]) r <- r + 1L
      blk_max <- cfg$block_len
      idx <- l:r
      prop <- block_draw(idx, x[l - 1L, ], x[r + 1L, ]); cur <- x[idx, , drop = FALSE]
      B2 <- length(idx)
      # exact correction vs the surrogate used to build the draw. All 2B exact
      # evaluations (prop + cur over the block) in ONE batched solar_zenith call.
      le   <- emit_block_ll(c(idx, idx), c(prop[, 1], cur[, 1]), c(prop[, 2], cur[, 2]))
      le_p <- le[seq_len(B2)]; le_c <- le[B2 + seq_len(B2)]
      lc_p <- vapply(seq_len(B2), function(a) ll_cheap_knot(idx[a], prop[a, 1], prop[a, 2]), numeric(1))
      lc_c <- vapply(seq_len(B2), function(a) ll_cheap_knot(idx[a], cur[a, 1],  cur[a, 2]),  numeric(1))
      corr <- sum(le_p - lc_p) - sum(le_c - lc_c)
      att <- att + 1L
      if (log(runif(1)) < corr) { x[idx, ] <- prop; acc <- acc + 1L }
      l <- r + 1L
    }
    x <- gold_site(x, K)   # free endpoint, single-site exact
    if (s > burn && (s - burn) %% thin == 0) { s_keep <- s_keep + 1L; keep[s_keep, , ] <- x }
  }
  BLOCK_ACC <<- acc / att
  keep[seq_len(s_keep), , , drop = FALSE]
}

# ---------------------------------------------------------------------
# 6. Run both and compare
# ---------------------------------------------------------------------
gold_sweeps <- if (!is.null(cfg$gold_sweeps)) cfg$gold_sweeps else cfg$sweeps
gold_burn   <- if (!is.null(cfg$gold_burn))   cfg$gold_burn   else cfg$burn
message(sprintf("running gold (always-exact single-site), %d sweeps ...", gold_sweeps))
CTR$pts <- 0L; CTR$calls <- 0L; CTR$emit <- 0L; t0 <- Sys.time()
G <- run_gold(gold_sweeps, gold_burn, cfg$thin); gold_pts <- CTR$pts; gold_calls <- CTR$calls; gold_time <- as.numeric(Sys.time() - t0, units = "secs")

message("running block (surrogate-posterior + DA correction) ...")
CTR$pts <- 0L; CTR$calls <- 0L; CTR$emit <- 0L; t0 <- Sys.time()
B <- run_block(); blk_pts <- CTR$pts; blk_calls <- CTR$calls; blk_emit <- CTR$emit; blk_time <- as.numeric(Sys.time() - t0, units = "secs")
blk_acc <- BLOCK_ACC

# per-knot marginal comparison (lon and lat)
summ <- function(A) list(mlon = apply(A[,,1], 2, mean), slon = apply(A[,,1], 2, sd),
                         mlat = apply(A[,,2], 2, mean), slat = apply(A[,,2], 2, sd))
sg <- summ(G); sb <- summ(B)
ks_lon <- ks_lat <- numeric(K)
zmean_lon <- zmean_lat <- numeric(K)   # |mean gap| / pooled sd
dmu_lon <- dmu_lat <- numeric(K)       # |mean gap| in degrees
rsd_lon <- rsd_lat <- numeric(K)       # sd ratio block/gold
for (k in 2:K) {
  ks_lon[k] <- suppressWarnings(ks.test(G[,k,1], B[,k,1])$p.value)
  ks_lat[k] <- suppressWarnings(ks.test(G[,k,2], B[,k,2])$p.value)
  psd_lon <- (sg$slon[k] + sb$slon[k]) / 2 + 1e-9
  psd_lat <- (sg$slat[k] + sb$slat[k]) / 2 + 1e-9
  dmu_lon[k] <- abs(sg$mlon[k] - sb$mlon[k]); dmu_lat[k] <- abs(sg$mlat[k] - sb$mlat[k])
  zmean_lon[k] <- dmu_lon[k] / psd_lon;       zmean_lat[k] <- dmu_lat[k] / psd_lat
  rsd_lon[k] <- sb$slon[k] / (sg$slon[k] + 1e-9); rsd_lat[k] <- sb$slat[k] / (sg$slat[k] + 1e-9)
}

cat("\n================= EXACTNESS: block vs gold =================\n")
cat(sprintf("kept draws: gold=%d block=%d ; block acceptance=%.2f\n", dim(G)[1], dim(B)[1], blk_acc))
cat(sprintf("KS p-value (interior)   lon: min=%.3f median=%.3f | lat: min=%.3f median=%.3f\n",
            min(ks_lon[2:K]), median(ks_lon[2:K]), min(ks_lat[2:K]), median(ks_lat[2:K])))
cat(sprintf("mean gap |dmu|/sd       lon: max=%.3f mean=%.3f | lat: max=%.3f mean=%.3f\n",
            max(zmean_lon[2:K]), mean(zmean_lon[2:K]), max(zmean_lat[2:K]), mean(zmean_lat[2:K])))
cat(sprintf("mean gap (degrees)      lon: max=%.3f mean=%.3f | lat: max=%.3f mean=%.3f\n",
            max(dmu_lon[2:K]), mean(dmu_lon[2:K]), max(dmu_lat[2:K]), mean(dmu_lat[2:K])))
cat(sprintf("sd ratio block/gold     lon: [%.2f, %.2f] | lat: [%.2f, %.2f]\n",
            min(rsd_lon[2:K]), max(rsd_lon[2:K]), min(rsd_lat[2:K]), max(rsd_lat[2:K])))

cat("\n================= COST =================\n")
gold_ffi_sw <- gold_calls / gold_sweeps
blk_ffi_sw  <- (blk_calls + blk_emit) / cfg$sweeps
cat(sprintf("likelihood FFI calls/sweep: gold=%.1f (eval_logpk_grid)  block=%.1f (%.1f emit + %.1f single-site)\n",
            gold_ffi_sw, blk_ffi_sw, blk_emit / cfg$sweeps, blk_calls / cfg$sweeps))
cat(sprintf("  -> per-sweep FFI batching: %.1fx fewer calls into the likelihood\n", gold_ffi_sw / blk_ffi_sw))
cat(sprintf("exact point-evals/sweep:    gold=%.0f  block=%.0f  (same exact work; batched into fewer calls)\n",
            gold_pts / gold_sweeps, blk_pts / cfg$sweeps))
cat(sprintf("total likelihood FFI calls: gold=%d  block=%d (%d emit + %d single-site)\n",
            gold_calls, blk_calls + blk_emit, blk_emit, blk_calls))
cat(sprintf("wall time (s):      gold=%.1f (%d sw)  block=%.1f (%d sw)\n",
            gold_time, gold_sweeps, blk_time, cfg$sweeps))

# worst-offending knots and endpoint behaviour (endpoint K is single-site in BOTH)
wk_sd <- which.min(rsd_lat[2:K]) + 1L; wk_mu <- which.max(zmean_lat[2:K]) + 1L
cat("\n----- worst knots -----\n")
cat(sprintf("min lat sd-ratio at knot %d/%d (endpoint=%d): ratio=%.2f  gold_sd=%.2f blk_sd=%.2f\n",
            wk_sd, K, K, rsd_lat[wk_sd], sg$slat[wk_sd], sb$slat[wk_sd]))
cat(sprintf("max lat mean-gap at knot %d/%d: z=%.2f  gold_mean=%.2f blk_mean=%.2f\n",
            wk_mu, K, zmean_lat[wk_mu], sg$mlat[wk_mu], sb$mlat[wk_mu]))
cat("last-4 lat sd-ratio (near free endpoint):", sprintf("%.2f", rsd_lat[(K-3):K]), "\n")

# most ambiguous knot = largest surrogate latitude sd (informative knots only)
surr_slat <- sapply(seq_len(K), function(k) if (surr_info[k]) sqrt(surr_S[[k]][2,2]) else NA_real_)
amb_k <- which.max(surr_slat)
cat(sprintf("\nmost ambiguous knot %d: surrogate lat sd=%.1f deg ; gold lat draws range [%.1f, %.1f], block [%.1f, %.1f]\n",
            amb_k, surr_slat[amb_k], min(G[,amb_k,2]), max(G[,amb_k,2]), min(B[,amb_k,2]), max(B[,amb_k,2])))

# comparison plot: mean tracks + per-knot latitude band + ambiguous-knot lat marginal
out_png <- if (exists("OUT_PNG")) OUT_PNG else "da_block_compare.png"
grDevices::png(out_png, width = 1500, height = 460)
graphics::par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))
# (a) mean tracks over true track
plot(df$true_lon, df$true_lat, type = "l", col = "grey60", lwd = 3,
     xlab = "lon", ylab = "lat", main = sprintf("%s: mean tracks", cfg$dataset))
graphics::lines(sg$mlon, sg$mlat, col = "black", lwd = 1.5)
graphics::lines(sb$mlon, sb$mlat, col = "red", lwd = 1.5, lty = 2)
graphics::legend("topright", c("truth", "gold", "block"), col = c("grey60","black","red"),
                 lwd = c(3,1.5,1.5), lty = c(1,1,2), bty = "n")
# (b) latitude vs knot with 95% bands
kk <- 1:K
plot(kk, sg$mlat, type = "n", ylim = range(c(sg$mlat-2*sg$slat, sg$mlat+2*sg$slat), na.rm=TRUE),
     xlab = "knot", ylab = "lat", main = "per-knot latitude +/- 2sd")
graphics::polygon(c(kk, rev(kk)), c(sg$mlat-2*sg$slat, rev(sg$mlat+2*sg$slat)),
                  col = grDevices::rgb(0,0,0,0.12), border = NA)
graphics::polygon(c(kk, rev(kk)), c(sb$mlat-2*sb$slat, rev(sb$mlat+2*sb$slat)),
                  col = grDevices::rgb(1,0,0,0.12), border = NA)
graphics::lines(kk, sg$mlat, col = "black", lwd = 1.5)
graphics::lines(kk, sb$mlat, col = "red", lwd = 1.5, lty = 2)
graphics::abline(v = amb_k, col = "blue", lty = 3)
graphics::legend("topright", c("gold","block"), col = c("black","red"), lwd = 1.5, lty = c(1,2), bty = "n")
# (c) latitude marginal at the most ambiguous knot: does either sampler drop a mode?
rng <- range(c(G[,amb_k,2], B[,amb_k,2]))
brk <- seq(rng[1]-0.5, rng[2]+0.5, length.out = 40)
hg <- graphics::hist(G[,amb_k,2], breaks = brk, plot = FALSE)
hb <- graphics::hist(B[,amb_k,2], breaks = brk, plot = FALSE)
plot(hg$mids, hg$density, type = "s", col = "black", lwd = 1.5,
     ylim = c(0, max(hg$density, hb$density)),
     xlab = "lat", ylab = "density", main = sprintf("lat marginal @ knot %d", amb_k))
graphics::lines(hb$mids, hb$density, type = "s", col = "red", lwd = 1.5, lty = 2)
graphics::abline(v = df$true_lat[which.min(abs(as.numeric(df$time) - S$knot_times[amb_k]))], col = "grey40", lwd = 2)
graphics::legend("topright", c("gold","block","truth"), col = c("black","red","grey40"),
                 lwd = c(1.5,1.5,2), lty = c(1,2,1), bty = "n")
grDevices::dev.off()
cat(sprintf("\nsaved comparison plot: %s\n", normalizePath(out_png)))
