#' Hierarchical twilight-free geolocation over a panel of tags
#'
#' Fits a partial-pooling hierarchical model to a panel of light-logger tracks:
#' each tag's track is reconstructed with a surrogate-posterior block sampler
#' (delayed acceptance against the exact twilight-free light likelihood), while
#' the per-tag movement scale is pooled through a conjugate inverse-gamma /
#' gamma population model so short tracks borrow strength from the panel.
#'
#' The heavy sampling runs in the native Rust kernel; the coarse-HMM proposal
#' surrogate is built in R. Tracks may be of different lengths. Deployment and
#' (optional) retrieval locations are fixed at the first and last knot.
#'
#' @param data Either a list of data frames (one per tag) or a single long data
#'   frame with an id column (see \code{id}). Each frame holds a light time
#'   series; other columns (e.g. temperature) are ignored in this version.
#' @param locations A data frame with one row per tag: columns \code{id},
#'   \code{deploy_lon}, \code{deploy_lat}, and optionally \code{retrieve_lon},
#'   \code{retrieve_lat} (use \code{NA} to leave the retrieval end free).
#' @param columns Optional named character vector overriding column auto-detection,
#'   e.g. \code{c(time = "Date", light = "Light")}. Missing entries are auto-detected
#'   (a POSIXct/Date column for \code{time}; a column named light/lux/lig/irrad, else
#'   the first numeric non-time column, for \code{light}).
#' @param id Name of the id column when \code{data} is a single long data frame
#'   (default "id"). Ignored when \code{data} is a list.
#' @param step_hours Knot spacing in hours (uniform across tags so the pooled
#'   movement scale is in common units). Default 12.
#' @param calibration Optional \code{c(intercept, slope)} applied to all tags; if
#'   \code{NULL} (default) each tag is auto-calibrated from its first days at the
#'   deployment location.
#' @param likelihood_params Optional \code{c(lambda, max_light, prob_slab)}; if
#'   \code{NULL} it is derived per tag from the light range.
#' @param a_pop Fixed inverse-gamma shape for the per-tag movement variance (default 3).
#' @param hyperprior \code{c(g0, h0)} gamma hyperprior on the population scale beta
#'   (default vague \code{c(1e-3, 1e-3)}; movement variance is in km^2).
#' @param surrogate_diffusion Movement scale (km/sqrt(day)) used only to build the
#'   proposal surrogate (exactness-neutral). Default 100.
#' @param mesh_pad Degrees of padding around the deploy/retrieve box for the
#'   surrogate mesh (default 20; increase for wide-ranging animals).
#' @param coarse_res Surrogate mesh resolution in degrees (default 1.5).
#' @param inflate Surrogate variance inflation (>= 1; default 1.5).
#' @param block_len Maximum block length in knots (default 5).
#' @param sweeps,burn,thin MCMC length, burn-in and thinning (defaults 4000/1500/5).
#' @param polish Run the red-black single-site polish each sweep (default TRUE;
#'   needed for unbiased movement-variance estimation).
#' @param seed Integer RNG seed; \code{NULL} is non-deterministic.
#' @return An object of class \code{TwilightFreeHier} with elements \code{population}
#'   (population movement scale posterior, km/day), \code{movement} (per-tag
#'   movement posterior data frame), \code{tracks} (named list of per-knot location
#'   posteriors), \code{tags} (per-tag metadata) and \code{draws} (raw sig2/beta).
#' @importFrom stats quantile lm coef rgamma sd
#' @export
TwilightFreeHier <- function(data, locations,
                             columns = NULL, id = "id",
                             step_hours = 12,
                             calibration = NULL, likelihood_params = NULL,
                             a_pop = 3, hyperprior = c(1e-3, 1e-3),
                             surrogate_diffusion = 100, mesh_pad = 20, coarse_res = 1.5,
                             inflate = 1.5, block_len = 5L,
                             sweeps = 4000L, burn = 1500L, thin = 5L,
                             polish = TRUE, seed = NULL) {

  panel <- .tfh_as_list(data, id)
  tags  <- .tfh_resolve_ids(panel, locations)
  n <- length(tags$data)
  if (n < 1) stop("no tags in `data`")

  # ---- per-tag setup + surrogate ----
  ind <- vector("list", n)
  for (i in seq_len(n)) {
    df <- tags$data[[i]]
    cm <- .tfh_cols(df, columns)
    tm <- df[[cm$time]]; if (!inherits(tm, "POSIXct")) tm <- as.POSIXct(tm, tz = "UTC")
    lt <- as.numeric(df[[cm$light]])
    ok <- !is.na(tm) & !is.na(lt); tm <- tm[ok]; lt <- lt[ok]
    ord <- order(tm); tm <- tm[ord]; lt <- lt[ord]
    loc <- tags$loc[i, ]
    ind[[i]] <- .tfh_build_individual(tm, lt, loc, step_hours, calibration,
                                      likelihood_params, surrogate_diffusion,
                                      mesh_pad, coarse_res, inflate)
  }

  # ---- flatten to the global arrays the Rust kernel expects ----
  fl <- .tfh_flatten(ind)

  fit <- run_block_hier(
    as.integer(n), fl$knots_per_ind, fl$kob_start, fl$kob_len, fl$obs_times, fl$obs_light,
    fl$cal, fl$lp, fl$smu, fl$sp, fl$cinv, fl$slon, fl$slat, fl$elon, fl$elat,
    a_pop, hyperprior[1], hyperprior[2], as.integer(block_len),
    as.integer(sweeps), as.integer(burn), as.integer(thin), isTRUE(polish),
    as.numeric(if (is.null(seed)) 0 else seed))

  .tfh_assemble(fit, ind, tags, a_pop, step_hours)
}

# ---- input normalization -----------------------------------------------------
.tfh_as_list <- function(data, id) {
  if (is.data.frame(data)) {
    if (!id %in% names(data)) stop(sprintf("long-format `data` needs an id column '%s'", id))
    ids <- as.character(data[[id]])
    split(data[, setdiff(names(data), id), drop = FALSE], factor(ids, levels = unique(ids)))
  } else if (is.list(data)) {
    lapply(data, function(d) if (is.data.frame(d)) d else as.data.frame(d))
  } else stop("`data` must be a list of data frames or a single long data frame")
}

.tfh_resolve_ids <- function(panel, locations) {
  if (!is.data.frame(locations) || !all(c("deploy_lon", "deploy_lat") %in% names(locations)))
    stop("`locations` must be a data frame with deploy_lon, deploy_lat (and optional retrieve_lon/lat)")
  ids <- names(panel)
  if (is.null(ids) || any(ids == "")) {
    ids <- if ("id" %in% names(locations) && length(locations$id) == length(panel))
      as.character(locations$id) else paste0("tag", seq_along(panel))
  }
  names(panel) <- ids
  if ("id" %in% names(locations)) {
    loc <- locations[match(ids, as.character(locations$id)), , drop = FALSE]
    if (anyNA(loc$deploy_lon)) stop("some tag ids have no matching row in `locations`")
  } else {
    if (nrow(locations) != length(panel))
      stop("`locations` has no id column and its row count does not match the number of tags")
    loc <- locations
  }
  if (is.null(loc$retrieve_lon)) loc$retrieve_lon <- NA_real_
  if (is.null(loc$retrieve_lat)) loc$retrieve_lat <- NA_real_
  list(data = panel, loc = loc, ids = ids)
}

# ---- column detection --------------------------------------------------------
.tfh_cols <- function(df, columns) {
  tc <- if (!is.null(columns) && "time" %in% names(columns)) columns[["time"]] else {
    cand <- names(df)[vapply(df, function(x) inherits(x, "POSIXct") || inherits(x, "Date"), logical(1))]
    if (!length(cand)) stop("could not auto-detect a time column (POSIXct/Date); set columns=c(time=...)")
    cand[1]
  }
  lc <- if (!is.null(columns) && "light" %in% names(columns)) columns[["light"]] else {
    num <- names(df)[vapply(df, is.numeric, logical(1))]
    hit <- num[grepl("light|lux|^lig$|irrad|led", num, ignore.case = TRUE)]
    if (length(hit)) hit[1] else {
      num <- setdiff(num, tc)
      if (!length(num)) stop("could not auto-detect a light column; set columns=c(light=...)")
      num[1]
    }
  }
  list(time = tc, light = lc)
}

# ---- per-tag setup + coarse-HMM surrogate ------------------------------------
.tfh_build_individual <- function(tm, lt, loc, step_hours, calibration, likelihood_params,
                                  surrogate_diffusion, mesh_pad, coarse_res, inflate) {
  ut <- as.numeric(tm)
  dlon <- loc$deploy_lon; dlat <- loc$deploy_lat
  # calibration / likelihood params (fixed or auto)
  if (!is.null(calibration) && !is.null(likelihood_params)) {
    cal <- as.numeric(calibration); lp <- as.numeric(likelihood_params); lsh <- lt
  } else {
    min_l <- stats::quantile(lt, 0.05, na.rm = TRUE); max_l <- stats::quantile(lt, 0.95, na.rm = TRUE)
    lsh <- pmax(0, lt - min_l); maxs <- as.numeric(max_l - min_l)
    ci <- ut < (ut[1] + 3 * 24 * 3600)
    cz <- solar_zenith(ut[ci], rep(dlon, sum(ci)), rep(dlat, sum(ci)))
    ti <- which(lsh[ci] > 0 & lsh[ci] < maxs * 0.95 & cz > 85 & cz < 100)
    if (length(ti) > 10) {
      fc <- stats::lm(lsh[ci][ti] ~ cz[ti]); icpt <- stats::coef(fc)[1]; slp <- -stats::coef(fc)[2]
      if (is.na(slp) || slp <= 0) slp <- maxs / (96 - 85)
    } else { slp <- maxs / (96 - 85); icpt <- slp * 96 }
    cal <- as.numeric(c(icpt, slp)); lp <- as.numeric(c(1 / (maxs * 0.5), maxs, 0.10))
  }
  # uniform step_hours knots so per-knot sig2 is comparable across tags
  t0 <- ut[1]; t1 <- ut[length(ut)]; hstep <- step_hours * 3600
  K <- max(2L, as.integer(floor((t1 - t0) / hstep)) + 1L)
  kt <- t0 + (0:(K - 1)) * hstep
  # obs -> knot: (knot[k-1], knot[k]]; trailing obs fall in the last knot
  obin <- vector("list", K)
  for (k in seq_len(K)) {
    tp <- if (k == 1) kt[1] - hstep else kt[k - 1]
    tc <- if (k == K) t1 else kt[k]
    obin[[k]] <- which(ut > tp & ut <= tc)
  }
  # movement metric (local km/deg at the deployment latitude)
  km_lat <- 111.0; km_lon <- 111.0 * cos(dlat * pi / 180)
  Cinv <- c(km_lon^2, 0, 0, km_lat^2)
  # surrogate
  sur <- .tfh_surrogate(kt, obin, ut, lsh, cal, lp, loc, mesh_pad, coarse_res,
                        surrogate_diffusion, inflate, Cinv, hstep)
  list(K = K, kt = kt, ot = ut, ol = lsh, obin = obin, cal = cal, lp = lp,
       deploy = c(dlon, dlat), retrieve = c(loc$retrieve_lon, loc$retrieve_lat),
       Cinv = Cinv, sur = sur, tstep = hstep)
}

.tfh_surrogate <- function(kt, obin, ut, lsh, cal, lp, loc, mesh_pad, coarse_res,
                           surrogate_diffusion, inflate, Cinv, hstep) {
  K <- length(kt)
  lons <- c(loc$deploy_lon, loc$retrieve_lon); lats <- c(loc$deploy_lat, loc$retrieve_lat)
  lons <- lons[is.finite(lons)]; lats <- lats[is.finite(lats)]
  clon <- seq(min(lons) - mesh_pad, max(lons) + mesh_pad, by = coarse_res)
  clat <- seq(min(lats) - mesh_pad, max(lats) + mesh_pad, by = coarse_res)
  cg <- expand.grid(lon = clon, lat = clat); n <- nrow(cg)
  # emissions
  E <- matrix(1 / n, K, n)
  for (k in seq_len(K)) {
    j <- obin[[k]]; if (!length(j)) next
    ll <- eval_logpk_grid(cg$lon, cg$lat, ut[j], lsh[j], cal, lp)
    e <- exp(ll - max(ll)); E[k, ] <- e / sum(e)
  }
  # RW transition from surrogate_diffusion (proposal only)
  dtd <- hstep / 86400; var_km <- surrogate_diffusion^2 * dtd
  Pm <- matrix(Cinv, 2, 2, byrow = TRUE) / var_km
  A <- outer(cg$lon, cg$lon, function(a, b) b - a)
  B <- outer(cg$lat, cg$lat, function(a, b) b - a)
  quad <- Pm[1, 1] * A * A + 2 * Pm[1, 2] * A * B + Pm[2, 2] * B * B
  Tm <- exp(-0.5 * (quad - apply(quad, 1, min))); Tm <- Tm / rowSums(Tm); rm(A, B, quad)
  # forward-backward, pinned at deploy (and retrieve if given)
  sc <- which.min((cg$lon - loc$deploy_lon)^2 + (cg$lat - loc$deploy_lat)^2)
  al <- matrix(0, K, n); a <- numeric(n); a[sc] <- 1; al[1, ] <- a
  for (k in 2:K) { a <- as.numeric(crossprod(Tm, a)) * E[k, ]; a <- a / sum(a); al[k, ] <- a }
  be <- matrix(0, K, n)
  if (is.finite(loc$retrieve_lon)) {
    rc <- which.min((cg$lon - loc$retrieve_lon)^2 + (cg$lat - loc$retrieve_lat)^2)
    be[K, ] <- 0; be[K, rc] <- 1
  } else be[K, ] <- 1 / n
  for (k in (K - 1):1) { b <- as.numeric(Tm %*% (be[k + 1, ] * E[k + 1, ])); be[k, ] <- b / sum(b) }
  mu <- matrix(0, K, 2); P <- numeric(4 * K); info <- logical(K)
  for (k in seq_len(K)) {
    w <- al[k, ] * be[k, ]; sw <- sum(w)
    if (!is.finite(sw) || sw <= 0) next
    w <- w / sw
    m <- c(sum(w * cg$lon), sum(w * cg$lat))
    dx <- cg$lon - m[1]; dy <- cg$lat - m[2]
    S <- matrix(c(sum(w * dx * dx), sum(w * dx * dy), sum(w * dx * dy), sum(w * dy * dy)), 2, 2) *
      inflate + diag(1e-6, 2)
    mu[k, ] <- m; P[(4 * k - 3):(4 * k)] <- as.numeric(t(solve(S))); info[k] <- TRUE
  }
  list(mu = mu, P = P, info = info)
}

# ---- flatten to global arrays ------------------------------------------------
.tfh_flatten <- function(ind) {
  n <- length(ind)
  knots_per_ind <- as.integer(vapply(ind, function(E) E$K, integer(1)))
  obs_times <- numeric(0); obs_light <- numeric(0)
  kob_start <- integer(0); kob_len <- integer(0); smu <- numeric(0); sp <- numeric(0)
  cal <- numeric(0); lp <- numeric(0); cinv <- numeric(0)
  slon <- numeric(0); slat <- numeric(0); elon <- numeric(0); elat <- numeric(0)
  obs_off <- 0L
  for (i in seq_len(n)) {
    E <- ind[[i]]
    obs_times <- c(obs_times, E$ot); obs_light <- c(obs_light, E$ol)
    for (k in seq_len(E$K)) {
      j <- E$obin[[k]]
      if (length(j)) { kob_start <- c(kob_start, obs_off + min(j) - 1L); kob_len <- c(kob_len, length(j)) }
      else           { kob_start <- c(kob_start, obs_off); kob_len <- c(kob_len, 0L) }
      if (E$sur$info[k]) { smu <- c(smu, E$sur$mu[k, ]); sp <- c(sp, E$sur$P[(4 * k - 3):(4 * k)]) }
      else               { smu <- c(smu, c(0, 0)); sp <- c(sp, c(0, 0, 0, 0)) }
    }
    cal <- c(cal, E$cal); lp <- c(lp, E$lp); cinv <- c(cinv, E$Cinv)
    slon <- c(slon, E$deploy[1]); slat <- c(slat, E$deploy[2])
    elon <- c(elon, E$retrieve[1]); elat <- c(elat, E$retrieve[2])
    obs_off <- obs_off + length(E$ot)
  }
  # NA retrieval -> NaN so Rust treats the last knot as free
  elon[is.na(elon)] <- NaN; elat[is.na(elat)] <- NaN
  list(knots_per_ind = knots_per_ind, kob_start = as.integer(kob_start), kob_len = as.integer(kob_len),
       obs_times = obs_times, obs_light = obs_light, smu = smu, sp = sp, cal = cal, lp = lp,
       cinv = cinv, slon = slon, slat = slat, elon = elon, elat = elat)
}

# ---- assemble the result object ----------------------------------------------
.tfh_assemble <- function(fit, ind, tags, a_pop, step_hours) {
  n <- length(ind); ids <- tags$ids
  dt_day <- step_hours / 24
  kmday <- function(v) sqrt(pmax(v, 0)) / sqrt(dt_day)
  sig <- matrix(fit$sig2, nrow = fit$n_kept, ncol = n, byrow = TRUE)
  pop <- sqrt(pmax(fit$beta, 0) / (a_pop - 1)) / sqrt(dt_day)
  ci <- function(x) c(mean = mean(x), lower = stats::quantile(x, 0.025, names = FALSE),
                      upper = stats::quantile(x, 0.975, names = FALSE))
  movement <- data.frame(id = ids,
    sigma_kmday = apply(sig, 2, function(c) mean(kmday(c))),
    lower = apply(sig, 2, function(c) stats::quantile(kmday(c), 0.025, names = FALSE)),
    upper = apply(sig, 2, function(c) stats::quantile(kmday(c), 0.975, names = FALSE)),
    row.names = NULL, stringsAsFactors = FALSE)
  # per-tag track posteriors (global knot slices)
  koff <- cumsum(c(0, vapply(ind, function(E) E$K, integer(1))))
  tracks <- setNames(vector("list", n), ids)
  for (i in seq_len(n)) {
    gk <- (koff[i] + 1):koff[i + 1]
    tracks[[i]] <- data.frame(
      time = as.POSIXct(ind[[i]]$kt, origin = "1970-01-01", tz = "UTC"),
      lon = fit$mean_lon[gk], lat = fit$mean_lat[gk],
      lon_sd = fit$sd_lon[gk], lat_sd = fit$sd_lat[gk],
      row.names = NULL)
  }
  meta <- data.frame(id = ids,
    n_obs = vapply(ind, function(E) length(E$ot), integer(1)),
    n_knots = vapply(ind, function(E) E$K, integer(1)),
    retrieval = vapply(ind, function(E) all(is.finite(E$retrieve)), logical(1)),
    row.names = NULL, stringsAsFactors = FALSE)
  out <- list(population = ci(pop), movement = movement, tracks = tracks,
              tags = meta, draws = list(beta = pop, sig2_kmday = kmday(sig)))
  class(out) <- "TwilightFreeHier"
  out
}

#' @method print TwilightFreeHier
#' @export
print.TwilightFreeHier <- function(x, ...) {
  cat("\nHierarchical TwilightFree fit\n=============================================\n")
  cat(sprintf("Tags: %d   Knots: %d..%d   Obs: %d\n", nrow(x$tags),
              min(x$tags$n_knots), max(x$tags$n_knots), sum(x$tags$n_obs)))
  cat(sprintf("Population movement scale: %.1f km/day [%.1f, %.1f]\n",
              x$population["mean"], x$population["lower"], x$population["upper"]))
  cat("---------------------------------------------\nPer-tag movement (km/day):\n")
  m <- x$movement
  for (i in seq_len(nrow(m)))
    cat(sprintf("  %-10s %5.1f [%4.1f, %5.1f]\n", m$id[i], m$sigma_kmday[i], m$lower[i], m$upper[i]))
  cat("=============================================\n")
  invisible(x)
}

#' @method plot TwilightFreeHier
#' @param x A \code{TwilightFreeHier} object.
#' @param bands Draw +/- 2 sd latitude/longitude bands around each track (default TRUE).
#' @param ... Passed to the initial \code{plot}.
#' @importFrom grDevices rgb
#' @importFrom graphics lines points legend
#' @export
plot.TwilightFreeHier <- function(x, bands = TRUE, ...) {
  tr <- x$tracks; n <- length(tr)
  cols <- grDevices::hcl.colors(max(n, 2), "Dark 3")
  xr <- range(unlist(lapply(tr, function(d) d$lon)))
  yr <- range(unlist(lapply(tr, function(d) d$lat)))
  plot(NA, xlim = xr, ylim = yr, xlab = "longitude", ylab = "latitude",
       main = "TwilightFreeHier: reconstructed tracks", ...)
  for (i in seq_len(n)) {
    d <- tr[[i]]
    if (isTRUE(bands)) {
      graphics::polygon(c(d$lon - 2*d$lon_sd, rev(d$lon + 2*d$lon_sd)),
                        c(d$lat, rev(d$lat)), col = grDevices::adjustcolor(cols[i], 0.12), border = NA)
    }
    graphics::lines(d$lon, d$lat, col = cols[i], lwd = 1.5)
    graphics::points(d$lon[1], d$lat[1], pch = 19, col = cols[i])           # deploy
    if (x$tags$retrieval[i]) graphics::points(d$lon[nrow(d)], d$lat[nrow(d)], pch = 4, col = cols[i], lwd = 2)  # retrieve
  }
  graphics::legend("topright", legend = x$tags$id, col = cols[seq_len(n)], lwd = 1.5, bty = "n")
  invisible(x)
}
