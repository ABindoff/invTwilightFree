# =====================================================================
# H2: hierarchical partial-pooling sampler for a PANEL of tracks.
# Each individual's track is updated with the validated single-track block/DA
# kernel (mode="block") or a gold always-exact single-site kernel (mode="gold").
# The per-individual movement variance sig2_i and the population level beta are
# conjugate Gibbs draws:
#     sig2_i | x_i ~ InvGamma(a_pop + (T-1), beta + 0.5 * SS_i)
#     beta   | .   ~ Gamma(g0 + N*a_pop, h0 + sum(1/sig2_i))
# Gold check: block and gold must give the same posterior for the pooled
# parameters (population scale, per-individual sig2_i).
# =====================================================================
suppressPackageStartupMessages(library(invTwilightFree))
source("sim_panel.R")   # simulate_panel(); guarded, no side effects

if (!exists("CFG")) CFG <- list()
cfg <- utils::modifyList(list(
  N = 6, days = 14, step_hours = 12.0,
  mode = "block",            # "block" | "gold"
  block_len = 5L, polish = TRUE,   # block move + single-site wiggle polish (both exact)
  coarse_res = 1.0, mesh_pad = 12.0, inflate = 1.5,
  surrogate_diffusion = 45.0,   # km/sqrt(day) used ONLY to build the surrogate
  a_pop = 3.0, g0 = 1e-3, h0 = 1e-3,   # vague hyperprior (sig2 is in km^2 ~ O(1000))
  a0 = 1e-3, b0 = 1e-3,                 # vague per-track prior for the INDEPENDENT (unpooled) fit
  sweeps = 3000L, burn = 1000L, thin = 3L,
  seed = 1L
), CFG)
set.seed(cfg$seed)

CTR <- new.env(); CTR$calls <- 0L; CTR$emit <- 0L

# ---- per-individual setup: knots, calibration, obs bins (as TwilightFreeGrid) ----
setup_track <- function(df, start_lon, start_lat, step_hours) {
  light <- df$light; date_time <- df$time
  min_l <- stats::quantile(light, 0.05, na.rm = TRUE); max_l <- stats::quantile(light, 0.95, na.rm = TRUE)
  lsh <- pmax(0, light - min_l); maxs <- as.numeric(max_l - min_l)
  ci <- date_time < (date_time[1] + 3*24*3600)
  cz <- solar_zenith(as.numeric(date_time[ci]), rep(start_lon, sum(ci)), rep(start_lat, sum(ci)))
  ti <- which(lsh[ci] > 0 & lsh[ci] < maxs*0.95 & cz > 85 & cz < 100)
  if (length(ti) > 10) {
    fitc <- stats::lm(lsh[ci][ti] ~ cz[ti]); icpt <- stats::coef(fitc)[1]; slp <- -stats::coef(fitc)[2]
    if (is.na(slp) || slp <= 0) slp <- maxs/(96-85)
  } else { slp <- maxs/(96-85); icpt <- slp*96 }
  calibration <- as.numeric(c(icpt, slp)); likpar <- as.numeric(c(1/(maxs*0.5), maxs, 0.10))
  ut <- as.numeric(date_time); t0 <- ut[1]; t1 <- ut[length(ut)]
  ks <- ceiling((t1-t0)/(step_hours*3600)) + 1; tstep <- if (ks>1) (t1-t0)/(ks-1) else 0
  kt <- t0 + (0:(ks-1))*tstep
  valid <- !is.na(lsh); ot <- ut[valid]; ol <- lsh[valid]
  obin <- vector("list", ks)
  for (k in seq_len(ks)) {
    tc <- kt[k]; tp <- if (k==1) tc-(kt[2]-kt[1]) else kt[k-1]
    obin[[k]] <- which(ot > tp & ot <= tc)
  }
  list(calibration=calibration, likpar=likpar, kt=kt, K=ks, tstep=tstep, ot=ot, ol=ol, obin=obin)
}

gauss_from_weights <- function(w, lon, lat, inflate) {
  w <- w/sum(w); mu <- c(sum(w*lon), sum(w*lat))
  dx <- lon-mu[1]; dy <- lat-mu[2]
  Sig <- matrix(c(sum(w*dx*dx), sum(w*dx*dy), sum(w*dx*dy), sum(w*dy*dy)),2,2)*inflate + diag(1e-6,2)
  list(mu=mu, S=Sig, P=solve(Sig))
}

# ---- construct one individual as a self-contained environment ----
make_individual <- function(df, start_lon, start_lat) {
  E <- new.env()
  E$S <- setup_track(df, start_lon, start_lat, cfg$step_hours)
  K <- E$S$K; E$K <- K; E$start <- c(start_lon, start_lat)
  cal <- E$S$calibration; lp <- E$S$likpar
  E$intercept <- cal[1]; E$slope <- cal[2]
  E$lambda <- lp[1]; E$max_light <- lp[2]; E$prob_slab <- lp[3]
  # local metric: Cinv = diag(km_lon^2, km_lat^2); movement cov = sig2 * Cinv^{-1}
  lat_ref <- mean(df$true_lat)
  km_lat <- 111.0; km_lon <- 111.0*cos(lat_ref*pi/180)
  E$Cinv <- diag(c(km_lon^2, km_lat^2))     # P_move = Cinv / sig2

  # ll_exact via eval_logpk_grid (gold path)
  E$ll_exact <- function(k, lon, lat) {
    j <- E$S$obin[[k]]; if (length(j)==0L) return(rep(0, length(lon)))
    CTR$calls <- CTR$calls + 1L
    eval_logpk_grid(as.numeric(lon), as.numeric(lat), E$S$ot[j], E$S$ol[j], cal, lp)
  }
  # batched emitter (block path): all (knot,point) pairs in one solar_zenith call
  E$emit <- function(kv, lonv, latv) {
    m <- length(kv); counts <- lengths(E$S$obin[kv]); out <- numeric(m)
    if (sum(counts)==0L) return(out)
    pt <- rep(seq_len(m), counts); oj <- unlist(E$S$obin[kv], use.names=FALSE)
    z <- solar_zenith(E$S$ot[oj], lonv[pt], latv[pt]); CTR$emit <- CTR$emit + 1L
    expc <- pmin(pmax(E$intercept - E$slope*z, 0), E$max_light); ol <- E$S$ol[oj]
    raw <- ifelse(ol<=expc, E$lambda*exp(-E$lambda*(expc-ol)), E$lambda*exp(-E$lambda*2*(ol-expc)))
    norm <- pmax(1 - exp(-E$lambda*expc) + 0.5*(1 - exp(-2*E$lambda*(E$max_light-expc))), 1e-12)
    den <- (1-E$prob_slab)*(raw/norm) + E$prob_slab*(1/E$max_light)
    rs <- rowsum(log(den), pt); out[as.integer(rownames(rs))] <- rs[,1]; out
  }
  # coarse grid-HMM smoother surrogate (pure R), built once with a fixed diffusion
  build_surrogate <- function() {
    clon <- seq(min(df$true_lon)-cfg$mesh_pad, max(df$true_lon)+cfg$mesh_pad, by=cfg$coarse_res)
    clat <- seq(min(df$true_lat)-cfg$mesh_pad, max(df$true_lat)+cfg$mesh_pad, by=cfg$coarse_res)
    cg <- expand.grid(lon=clon, lat=clat); n <- nrow(cg)
    Emat <- matrix(1/n, K, n)
    for (k in seq_len(K)) { j <- E$S$obin[[k]]; if (!length(j)) next
      ll <- E$ll_exact(k, cg$lon, cg$lat); e <- exp(ll-max(ll)); Emat[k,] <- e/sum(e) }
    dtd <- E$S$tstep/86400; var_km <- cfg$surrogate_diffusion^2*dtd
    Pm <- E$Cinv/var_km
    A <- outer(cg$lon, cg$lon, function(a,b) b-a); B <- outer(cg$lat, cg$lat, function(a,b) b-a)
    quad <- Pm[1,1]*A*A + 2*Pm[1,2]*A*B + Pm[2,2]*B*B
    Tm <- exp(-0.5*(quad - apply(quad,1,min))); Tm <- Tm/rowSums(Tm); rm(A,B,quad)
    sc <- which.min((cg$lon-start_lon)^2 + (cg$lat-start_lat)^2)
    al <- matrix(0,K,n); a <- numeric(n); a[sc] <- 1; al[1,] <- a
    for (k in 2:K) { a <- as.numeric(crossprod(Tm,a))*Emat[k,]; a <- a/sum(a); al[k,] <- a }
    be <- matrix(0,K,n); be[K,] <- 1/n
    for (k in (K-1):1) { b <- as.numeric(Tm %*% (be[k+1,]*Emat[k+1,])); be[k,] <- b/sum(b) }
    mu <- matrix(NA_real_,K,2); Sl <- vector("list",K); Pl <- vector("list",K); info <- logical(K)
    for (k in seq_len(K)) { w <- al[k,]*be[k,]; sw <- sum(w)
      if (!is.finite(sw)||sw<=0) { Pl[[k]] <- matrix(0,2,2); next }
      g <- gauss_from_weights(w/sw, cg$lon, cg$lat, cfg$inflate)
      mu[k,] <- g$mu; Sl[[k]] <- g$S; Pl[[k]] <- g$P; info[k] <- TRUE }
    list(mu=mu, S=Sl, P=Pl, info=info)
  }
  E$sur <- build_surrogate()

  E$ll_cheap <- function(k, lon, lat) {
    if (!E$sur$info[k]) return(rep(0, length(lon)))
    P <- E$sur$P[[k]]; mu <- E$sur$mu[k,]; dx <- lon-mu[1]; dy <- lat-mu[2]
    -log(2*pi) - 0.5*log(det(E$sur$S[[k]])) - 0.5*(P[1,1]*dx*dx + 2*P[1,2]*dx*dy + P[2,2]*dy*dy)
  }
  # surrogate-mean init track
  E$init <- function() {
    x <- cbind(rep(start_lon,K), rep(start_lat,K))
    for (k in seq_len(K)) if (E$sur$info[k]) x[k,] <- E$sur$mu[k,]
    x[1,] <- c(start_lon, start_lat); x
  }
  E$draw2 <- function(mean, prec) as.numeric(mean + t(chol(solve(prec))) %*% rnorm(2))
  E
}

# ---- track update kernels (take current sig2 -> P_move) ----
gold_site <- function(E, x, t, Pm) {
  K <- E$K
  if (t==K) { m <- x[t-1,]; prec <- Pm } else { m <- (x[t-1,]+x[t+1,])/2; prec <- 2*Pm }
  prop <- E$draw2(m, prec); cur <- x[t,]
  le <- E$ll_exact(t, c(prop[1],cur[1]), c(prop[2],cur[2]))
  if (log(runif(1)) < le[1]-le[2]) x[t,] <- prop
  x
}
block_draw <- function(E, idx, xL, xR, Pm) {
  B <- length(idx); d <- 2*B; Q <- matrix(0,d,d); b <- numeric(d)
  for (a in seq_len(B)) { ia <- (2*a-1):(2*a)
    Q[ia,ia] <- Q[ia,ia] + 2*Pm
    if (a>1) { ip <- (2*a-3):(2*a-2); Q[ia,ip] <- Q[ia,ip]-Pm; Q[ip,ia] <- Q[ip,ia]-Pm }
    Pk <- E$sur$P[[idx[a]]]; Q[ia,ia] <- Q[ia,ia]+Pk
    if (E$sur$info[idx[a]]) b[ia] <- b[ia] + Pk %*% E$sur$mu[idx[a],]
  }
  b[1:2] <- b[1:2] + Pm %*% xL; b[(d-1):d] <- b[(d-1):d] + Pm %*% xR
  R <- chol(Q); mu <- solve(Q,b); matrix(as.numeric(mu + backsolve(R, rnorm(d))), B, 2, byrow=TRUE)
}
update_gold <- function(E, x, sig2) {
  Pm <- E$Cinv/sig2; for (t in 2:E$K) x <- gold_site(E, x, t, Pm); x
}
# Batched red-black single-site polish: interior knots of one parity are
# conditionally independent given the other parity, so each half-sweep proposes
# from the RW full-conditional and evaluates ALL of them in ONE emit call. Exact,
# and mixes the high-frequency wiggle that the smooth block proposal misses.
rb_polish <- function(E, x, Pm) {
  K <- E$K; L <- t(chol(solve(2*Pm)))
  for (par in c(0L, 1L)) {
    ts <- (2:(K-1L))[(2:(K-1L)) %% 2L == par]
    if (!length(ts)) next
    prop <- matrix(0, length(ts), 2)
    for (a in seq_along(ts)) { t <- ts[a]; prop[a, ] <- (x[t-1,]+x[t+1,])/2 + L %*% rnorm(2) }
    cur <- x[ts, , drop = FALSE]; B2 <- length(ts)
    le <- E$emit(c(ts, ts), c(prop[,1], cur[,1]), c(prop[,2], cur[,2]))
    acc <- log(runif(B2)) < le[seq_len(B2)] - le[B2 + seq_len(B2)]
    if (any(acc)) x[ts[acc], ] <- prop[acc, ]
  }
  x
}
update_block <- function(E, x, sig2, block_len) {
  K <- E$K; Pm <- E$Cinv/sig2
  phase <- sample(0:(block_len-1L), 1L); blk_max <- max(1L, block_len-phase)
  l <- 2L
  while (l <= K-1L) {
    if (!E$sur$info[l]) { x <- gold_site(E, x, l, Pm); l <- l+1L; next }
    r <- l; while (r < min(l+blk_max-1L, K-1L) && E$sur$info[r+1L]) r <- r+1L
    blk_max <- block_len; idx <- l:r; B2 <- length(idx)
    prop <- block_draw(E, idx, x[l-1L,], x[r+1L,], Pm); cur <- x[idx,,drop=FALSE]
    le <- E$emit(c(idx,idx), c(prop[,1],cur[,1]), c(prop[,2],cur[,2]))
    le_p <- le[seq_len(B2)]; le_c <- le[B2+seq_len(B2)]
    lc_p <- vapply(seq_len(B2), function(a) E$ll_cheap(idx[a], prop[a,1], prop[a,2]), numeric(1))
    lc_c <- vapply(seq_len(B2), function(a) E$ll_cheap(idx[a], cur[a,1], cur[a,2]), numeric(1))
    if (log(runif(1)) < sum(le_p-lc_p) - sum(le_c-lc_c)) x[idx,] <- prop
    l <- r+1L
  }
  x <- gold_site(E, x, K, Pm)
  # block mixes global position; the batched red-black polish mixes the
  # high-frequency wiggle that drives the movement-variance estimate. Both exact.
  if (isTRUE(cfg$polish)) x <- rb_polish(E, x, Pm)
  x
}

# SS_i = sum_t d_t' Cinv d_t over increments (for the conjugate sig2 draw)
ss_track <- function(E, x) {
  D <- x[2:E$K,,drop=FALSE] - x[1:(E$K-1),,drop=FALSE]
  sum(rowSums((D %*% E$Cinv) * D))
}

# ---------------------------------------------------------------------
# Build the panel and the individuals
# ---------------------------------------------------------------------
PAN <- simulate_panel(N=cfg$N, days=cfg$days, step_min=10, pop_km_per_day=45)
N <- cfg$N
message(sprintf("building %d individuals (surrogates) ...", N))
IND <- lapply(seq_len(N), function(i)
  make_individual(PAN$ind[[i]], PAN$start$lon[i], PAN$start$lat[i]))
Kv <- sapply(IND, function(E) E$K)
cat(sprintf("individuals built. knots per ind: %s\n", paste(Kv, collapse=",")))

# ---------------------------------------------------------------------
# Hierarchical Gibbs
# ---------------------------------------------------------------------
run_hier <- function(mode = cfg$mode, sweeps = cfg$sweeps, burn = cfg$burn, thin = cfg$thin,
                     pool = TRUE) {
  X <- lapply(IND, function(E) E$init())
  sig2 <- sapply(IND, function(E) ss_track(E, E$init())/(E$K-1))   # data-driven start
  beta <- 1
  keep_pop <- numeric(0); keep_sig <- matrix(0, 0, N)
  for (s in seq_len(sweeps)) {
    for (i in seq_len(N)) {
      E <- IND[[i]]
      if (mode == "block") X[[i]] <- update_block(E, X[[i]], sig2[i], cfg$block_len)
      else                 X[[i]] <- update_gold(E, X[[i]], sig2[i])
      SS <- ss_track(E, X[[i]])
      # conjugate per-individual movement variance. pool=TRUE shares the population
      # scale beta (partial pooling); pool=FALSE uses a vague per-track prior (no
      # information shared between individuals).
      if (pool) sig2[i] <- 1/rgamma(1, shape = cfg$a_pop + (E$K-1), rate = beta + 0.5*SS)
      else      sig2[i] <- 1/rgamma(1, shape = cfg$a0    + (E$K-1), rate = cfg$b0 + 0.5*SS)
    }
    if (pool) beta <- rgamma(1, shape = cfg$g0 + N*cfg$a_pop, rate = cfg$h0 + sum(1/sig2))
    if (s > burn && (s-burn) %% thin == 0) {
      keep_pop <- c(keep_pop, if (pool) sqrt(beta/(cfg$a_pop-1)) else NA_real_)
      keep_sig <- rbind(keep_sig, sqrt(sig2))              # per-ind per-step km sd
    }
  }
  list(pop = keep_pop, sig = keep_sig)
}

dt_day <- IND[[1]]$S$tstep/86400
summ_run <- function(R) list(pop = R$pop/sqrt(dt_day), sig = R$sig/sqrt(dt_day))

if (isTRUE(cfg$benefit)) {
  message("running POOLED (hierarchical) fit ...")
  set.seed(cfg$seed)
  Rp <- summ_run(run_hier("block", cfg$sweeps, cfg$burn, cfg$thin, pool = TRUE))
  message("running INDEPENDENT (unpooled) fit ...")
  set.seed(cfg$seed)
  Ri <- summ_run(run_hier("block", cfg$sweeps, cfg$burn, cfg$thin, pool = FALSE))
  truth <- PAN$sd_day_true; days <- PAN$days
  est_p <- colMeans(Rp$sig); est_i <- colMeans(Ri$sig)
  rmse <- function(e, s = seq_len(N)) sqrt(mean((e[s] - truth[s])^2))
  cat("\n=== POOLING BENEFIT (heterogeneous track lengths): sigma_i (km/day) ===\n")
  cat(sprintf("%-4s %6s %8s %12s %8s\n", "ind", "days", "truth", "independent", "pooled"))
  o <- order(days)
  for (k in o)
    cat(sprintf("%-4d %6.0f %8.1f %12.1f %8.1f\n", k, days[k], truth[k], est_i[k], est_p[k]))
  shortS <- which(days <= stats::median(days)); longS <- which(days > stats::median(days))
  cat(sprintf("\nRMSE vs truth        indep %.2f | pooled %.2f  (%.0f%% lower)\n",
      rmse(est_i), rmse(est_p), 100*(1 - rmse(est_p)/rmse(est_i))))
  cat(sprintf("  short tracks       indep %.2f | pooled %.2f  (%.0f%% lower)\n",
      rmse(est_i, shortS), rmse(est_p, shortS), 100*(1 - rmse(est_p,shortS)/rmse(est_i,shortS))))
  cat(sprintf("  long tracks        indep %.2f | pooled %.2f  (%.0f%% lower)\n",
      rmse(est_i, longS), rmse(est_p, longS), 100*(1 - rmse(est_p,longS)/rmse(est_i,longS))))
  shrink <- abs(est_i - est_p)
  cat(sprintf("mean |shrinkage|:    short %.1f | long %.1f km/day (partial pooling: short shrink more)\n",
      mean(shrink[shortS]), mean(shrink[longS])))

  out_png <- if (exists("OUT_PNG")) OUT_PNG else "pooling_benefit.png"
  grDevices::png(out_png, width = 1150, height = 500)
  graphics::par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  pop_lvl <- mean(Rp$pop, na.rm = TRUE)
  # colour ramp by track length (short = red, long = blue)
  drank <- (days - min(days)) / (max(days) - min(days) + 1e-9)
  col_i <- grDevices::rgb(1 - drank, 0, drank)
  # (a) estimate vs truth, indep -> pooled shrinkage arrows, coloured by length
  rng <- range(c(truth, est_i, est_p))
  plot(truth, est_i, pch = 1, col = col_i, cex = 1.4, lwd = 2, xlim = rng, ylim = rng,
       xlab = "true sigma_i (km/day)", ylab = "estimate (km/day)",
       main = "estimate vs truth (colour = track length)")
  graphics::abline(0, 1, col = "grey60", lwd = 2)
  graphics::abline(h = pop_lvl, col = "grey70", lty = 3)
  graphics::segments(truth, est_i, truth, est_p, col = "grey80")
  graphics::points(truth, est_p, pch = 19, col = col_i, cex = 1.2)
  graphics::legend("topleft", c("independent (open)", "pooled (filled)", "1:1", "pop level"),
                   col = c("black","black","grey60","grey70"), pch = c(1,19,NA,NA),
                   lty = c(NA,NA,1,3), bty = "n")
  # (b) partial-pooling signature: shrinkage vs track length
  plot(days, shrink, pch = 19, col = col_i, cex = 1.3,
       xlab = "track length (days)", ylab = "|shrinkage| = |indep - pooled| (km/day)",
       main = "shrinkage decreases with track length")
  if (length(unique(days)) > 2) {
    fit_s <- stats::lowess(days, shrink); graphics::lines(fit_s, col = "grey40", lwd = 2)
  }
  grDevices::dev.off()
  cat(sprintf("\nsaved %s\n", normalizePath(out_png)))
} else if (isTRUE(cfg$compare)) {
  gs <- if (!is.null(cfg$gold_sweeps)) cfg$gold_sweeps else cfg$sweeps
  message(sprintf("running GOLD hierarchy (%d sweeps) ...", gs))
  set.seed(cfg$seed); CTR$calls <- 0L; CTR$emit <- 0L; t0 <- Sys.time()
  G <- summ_run(run_hier("gold", gs, cfg$burn, cfg$thin)); gtime <- as.numeric(Sys.time()-t0, units="secs")
  g_eval <- CTR$calls
  message(sprintf("running BLOCK hierarchy (%d sweeps) ...", cfg$sweeps))
  set.seed(cfg$seed + 1L); CTR$calls <- 0L; CTR$emit <- 0L; t0 <- Sys.time()
  B <- summ_run(run_hier("block", cfg$sweeps, cfg$burn, cfg$thin))
  b_emit <- CTR$emit; b_eval <- CTR$calls; btime <- as.numeric(Sys.time()-t0, units="secs")

  cat("\n=== POOLED PARAMS: block vs gold (should coincide) vs truth ===\n")
  ksp <- suppressWarnings(ks.test(G$pop, B$pop)$p.value)
  cat(sprintf("pop scale km/day: gold %.1f [%.1f,%.1f] | block %.1f [%.1f,%.1f] | KS p=%.2f | truth %.0f\n",
      mean(G$pop), quantile(G$pop,.025), quantile(G$pop,.975),
      mean(B$pop), quantile(B$pop,.025), quantile(B$pop,.975), ksp, PAN$pop_km_per_day))
  for (i in seq_len(N)) {
    ks_i <- suppressWarnings(ks.test(G$sig[,i], B$sig[,i])$p.value)
    dz <- abs(mean(G$sig[,i])-mean(B$sig[,i]))/((sd(G$sig[,i])+sd(B$sig[,i]))/2+1e-9)
    cat(sprintf("sig ind %d: gold %.1f [%.1f,%.1f] | block %.1f [%.1f,%.1f] | dmu/sd %.2f KSp %.2f | truth %.1f\n",
        i, mean(G$sig[,i]), quantile(G$sig[,i],.025), quantile(G$sig[,i],.975),
        mean(B$sig[,i]), quantile(B$sig[,i],.025), quantile(B$sig[,i],.975), dz, ks_i, PAN$sd_day_true[i]))
  }
  cat(sprintf("\ncost: gold %.0fs (%d sw, %d eval FFI) | block %.0fs (%d sw, %d emit + %d eval FFI)\n",
      gtime, gs, g_eval, btime, cfg$sweeps, b_emit, b_eval))
} else {
  message(sprintf("running hierarchical sampler (mode=%s) ...", cfg$mode))
  t0 <- Sys.time(); R <- summ_run(run_hier())
  cat(sprintf("done in %.1fs ; kept %d draws ; FFI: %d eval + %d emit\n",
              as.numeric(Sys.time()-t0, units="secs"), length(R$pop), CTR$calls, CTR$emit))
  cat("\n=== per-individual movement (km/day): mean [95% CI] vs truth ===\n")
  for (i in seq_len(N))
    cat(sprintf("ind %d: %.1f [%.1f, %.1f]  (true %.1f)\n", i,
        mean(R$sig[,i]), quantile(R$sig[,i],.025), quantile(R$sig[,i],.975), PAN$sd_day_true[i]))
  cat(sprintf("\npopulation scale: %.1f [%.1f, %.1f] km/day (true %.0f)\n",
      mean(R$pop), quantile(R$pop,.025), quantile(R$pop,.975), PAN$pop_km_per_day))
}
