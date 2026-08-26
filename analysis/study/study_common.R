# =============================================================================
# SHARED SCAFFOLDING FOR THE REPORTING RUN
#
# Panel assembly, calibration, scoring and instrumentation. Every arm of the
# study calls these; no arm reimplements one. The campaign's most expensive
# mistake was three verbatim copies of the panel block drifting apart, and a
# second scorer defined locally in four scripts -- which is why arm A_shipped of
# fit_final_29.R "reproduced" its reference on the mean while missing it on every
# single tag. One definition, sourced everywhere.
#
# Run from the package root.
# =============================================================================

suppressMessages({
  library(data.table)
  library(terra)
  library(invTwilightFree)
})
Sys.setlocale("LC_TIME", "C")

source("analysis/study/study_config.R")
# nes_common.R supplies the readers this study does not redefine: read_xlsx_grid,
# nes_meta, nes_argos, argos_at, departure_time, gc_km, lon360, dlon.
source("analysis/diagnostics/nes_common.R")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. MACHINE STATE
# =============================================================================

# A deterministic numeric workload, timed. Recorded at the start and end of every
# phase so the run certifies its own machine conditions rather than asserting
# them in prose. The engine is single-threaded (no rayon in Cargo.toml, no
# parallel:: in R), so what this detects is contention for cache and memory
# bandwidth, which was measured at 1.33-1.56x on a saturated 24-core box.
machine_probe <- function(reps = 3L) {
  f <- function() {
    set.seed(7)
    x <- 0
    for (i in 1:40) x <- x + sum(sqrt(abs(sin(seq_len(2e5) * 1e-4))))
    x
  }
  f()                                   # warm up, discarded
  t <- vapply(seq_len(reps), function(i) system.time(f())[["elapsed"]], 0)
  list(median = median(t), min = min(t), max = max(t),
       cores = parallel::detectCores())
}

log_line <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), sprintf(...)))
  flush.console()
}

# =============================================================================
# 2. PANEL
# =============================================================================

tag_family <- function(serial) {
  fifelse(is.na(serial), "unknown",
    fifelse(grepl("^219", serial), "Mk9_219",
      fifelse(grepl("^18A", serial), "F18A", "other")))
}

# A marine animal's ground truth should not contain positions on land.
#
# The speed filter alone does not catch every one: a single 2023041 fix sits at
# 255.66 E, 44.06 N -- inland in the western United States, 1,699 km from the
# colony -- and survives because it follows a long Argos gap, so the implied speed
# is not impossible. It is outside the study domain, and widening the domain to
# accommodate a bad fix would be chasing noise with the search space.
#
# Haul-out fixes at Ano Nuevo are legitimately on land, so the gate exempts
# anything within LAND_EXEMPT_KM of the colony. Bathymetry is the only input.
LAND_EXEMPT_KM <- 50

range_filter <- function(d, sea_rast) {
  if (!nrow(d)) return(d)
  in_domain <- d$lon >= DOMAIN[["xmin"]] & d$lon <= DOMAIN[["xmax"]] &
               d$lat >= DOMAIN[["ymin"]] & d$lat <= DOMAIN[["ymax"]]
  # The bathymetry raster covers 164-243 E, narrower than the search domain, so
  # extract() returns NA west of 164 E. NA means "no coverage", not "land", and
  # such a fix is kept -- the domain test above is what removes the inland one.
  at_sea  <- as.numeric(terra::extract(sea_rast, cbind(d$lon, d$lat))[, 1])
  on_land <- !is.na(at_sea) & at_sea <= 0
  near    <- gc_km(d$lon, d$lat, COLONY[["lon"]], COLONY[["lat"]]) <= LAND_EXEMPT_KM
  d[in_domain & (near | !on_land)]
}

# Forward speed filter, applied to BOTH deliveries so the ground truth is
# processed identically across seasons.
#
# This was not uniform in the campaign. `nes_argos()` returns the 2021 Locations
# fixes unfiltered, while the 2022/23 RawArgos fixes were filtered during ingest.
# Scoring 2021 against unfiltered fixes inflates its error asymmetrically against
# the other seasons, which is one of the three confounds that made the season
# contrast uninterpretable.
#
# Anchored on the first good-class fix: a class-0 solution at the head of a record
# otherwise becomes the reference every later fix is judged against.
speed_filter <- function(d, vmax = ARGOS_VMAX_KMH) {
  setorder(d, time)
  if (nrow(d) < 3) return(d)
  keep <- rep(TRUE, nrow(d)); last <- 1L
  for (i in 2:nrow(d)) {
    dt <- as.numeric(difftime(d$time[i], d$time[last], units = "hours"))
    if (dt <= 0) { keep[i] <- FALSE; next }
    if (gc_km(d$lon[last], d$lat[last], d$lon[i], d$lat[i]) / dt > vmax) keep[i] <- FALSE
    else last <- i
  }
  d[keep]
}

serial_table <- function() {
  old <- basename(list.files("data/nes_untracked", pattern = "^[0-9]+_.*[.]csv$"))
  new <- list.files("fvilches/extracted/TDR raw", pattern = "Archive[.]csv$")
  unique(rbind(
    data.table(topp = sub("^([0-9]+)_.*$", "\\1", old),
               serial = sub("^[0-9]+_(.+)[.]csv$", "\\1", old)),
    data.table(topp = sub("^([0-9]+)_.*$", "\\1", new),
               serial = sub("^[0-9]+_(.+)-out-Archive[.]csv$", "\\1", new))),
    by = "topp")
}

# The panel, with ONE clipping rule and ONE endpoint rule for all three seasons.
#
# The campaign used different rules per delivery: 2021 clipped to the deployment
# date and derived its end anchor from the last 72 h of Argos; 2022/23 clipped to
# the Argos span and took the manifest recovery position, which sits over 100 km
# from the colony on 4 of 19 deployments. That made the season contrast
# uninterpretable and put truth-derived positions into the estimation. Both are
# now uniform: clip to the intersection of the archive and Argos spans, and pin
# both ends at the colony.
build_panel <- function() {
  old_arch  <- lapply(readRDS(file.path(CACHE_DIR, "archives_v3_30min.rds")), `[[`, "main")
  new_arch  <- lapply(readRDS(file.path(CACHE_DIR, "fvilches_archives_v1_30min.rds")), `[[`, "main")
  new_argos <- readRDS(file.path(CACHE_DIR, "fvilches_argos_v1.rds"))
  new_man   <- fread("scratch/nes_calibration/fvilches_manifest.csv",
                     colClasses = list(character = "topp"))
  meta <- nes_meta(); ar <- nes_argos()
  ser  <- serial_table()

  sea <- load_layers()$sea
  n_raw <- 0L; n_kept <- 0L
  mk <- function(id, season, light, argos, argos_deploy, argos_recover) {
    n_raw  <<- n_raw + nrow(argos)
    argos  <- range_filter(speed_filter(argos), sea)
    n_kept <<- n_kept + nrow(argos)
    t0 <- max(min(light$time), min(argos$time))
    t1 <- min(max(light$time), max(argos$time))
    d  <- light[time >= t0 & time <= t1]
    if (nrow(d) < 1000 || nrow(argos) < 50) return(NULL)
    list(id = id, season = season,
         serial = ser[topp == id]$serial[1],
         family = tag_family(ser[topp == id]$serial[1]),
         light = d, argos = argos[time >= t0 & time <= t1],
         t0 = t0, t1 = t1,
         days = as.numeric(difftime(t1, t0, units = "days")),
         # FITTED endpoints: the colony, both ends, for every deployment.
         p0 = unname(COLONY), p1 = unname(COLONY),
         # DIAGNOSTIC ONLY -- never passed to an engine.
         argos_deploy = argos_deploy, argos_recover = argos_recover)
  }

  tags <- list()
  for (tg in names(old_arch)) {
    m <- meta[match(tg, meta$id), ]
    if (is.na(m$ptt) || !nzchar(m$ptt)) next
    df <- ar$deploy[Ptt == m$ptt]; a <- ar$fixes[Ptt == m$ptt]
    if (!nrow(df) || nrow(a) < 50) next
    a <- a[, .(time, lon, lat)][order(time)]
    tf <- a[time >= max(time) - 72 * 3600]
    tags[[tg]] <- mk(tg, "2021", as.data.table(old_arch[[tg]]), a,
                     c(df$deploy_lon[1], df$deploy_lat[1]),
                     c(median(tf$lon), median(tf$lat)))
  }
  for (tg in names(new_arch)) {
    mm <- new_man[topp == tg]; if (!nrow(mm)) next
    a <- as.data.table(new_argos[[tg]])[, .(time, lon, lat)][order(time)]
    tags[[tg]] <- mk(tg, as.character(mm$season[1]), as.data.table(new_arch[[tg]]), a,
                     c(mm$deploy_lon[1], mm$deploy_lat[1]),
                     c(mm$recover_lon[1], mm$recover_lat[1]))
  }
  message(sprintf("Argos QC: %d fixes in, %d kept (%.3f%% dropped by the speed and land filters)",
                  n_raw, n_kept, 100 * (1 - n_kept / n_raw)))
  Filter(Negate(is.null), tags)
}

# =============================================================================
# 3. CALIBRATION
#
# Three modes, differing in ONE step each, so the contrast between adjacent modes
# identifies that step and nothing else. The intensity-scale window is identical
# in all three: only the geometry changes.
#
#   study        geometry on [t0, min(t0 + 7 d, departure)], pooled within family
#   pooled_untrunc  geometry on [t0, t0 + 7 d], pooled within family
#   pertag_untrunc  geometry on [t0, t0 + 7 d], each tag its own
#
#   study vs pooled_untrunc  -> what the departure truncation is worth
#   pooled_untrunc vs pertag_untrunc -> what pooling is worth
# =============================================================================

geometry_window <- function(g, truncate) {
  t0  <- as.numeric(min(g$light$time))
  end <- t0 + CAL_GEOM_MAX_DAYS * 86400
  if (truncate) {
    dep <- departure_time(g$light, dive_m = DEPARTURE_DIVE_M)
    if (is.finite(dep)) end <- min(end, dep)
  }
  as.numeric(g$light$time) < end
}

fit_calibration <- function(tags, mode = c("study", "pooled_untrunc", "pertag_untrunc")) {
  mode <- match.arg(mode)
  truncate <- identical(mode, "study")
  pooled   <- mode %in% c("study", "pooled_untrunc")

  # Intensity scale: each tag's own first CAL_SCALE_DAYS. Identical in every mode.
  scale_fits <- lapply(tags, function(g) {
    k <- g$light$time <= min(g$light$time) + CAL_SCALE_DAYS * 86400
    try(fit_light_response(g$light$time[k], g$light$light[k],
                           g$p0[1], g$p0[2]), silent = TRUE)
  })
  # Geometry: fitted at the colony, which is where the animal is for as long as
  # the window says it is.
  geom_fits <- lapply(tags, function(g) {
    k <- geometry_window(g, truncate)
    if (sum(k) < CAL_GEOM_MIN_OBS) return(NULL)
    r <- try(fit_light_response(g$light$time[k], g$light$light[k],
                                g$p0[1], g$p0[2]), silent = TRUE)
    if (inherits(r, "try-error")) NULL else r
  })
  bad <- vapply(scale_fits, function(x) inherits(x, "try-error"), NA)
  if (any(bad)) stop("intensity-scale fit failed for: ",
                     paste(names(scale_fits)[bad], collapse = ", "))

  responses <- list()
  diag <- list()
  if (pooled) {
    for (fam in unique(vapply(tags, function(g) g$family, ""))) {
      ids <- names(tags)[vapply(tags, function(g) g$family, "") == fam]
      gf <- geom_fits[ids]; sf <- scale_fits[ids]
      n_ok <- sum(!vapply(gf, is.null, logical(1)))
      if (!n_ok) { message("family ", fam, ": no usable geometry, skipped"); next }
      p <- pool_light_responses(gf, sf)
      a <- attr(p, "pooled")
      diag[[fam]] <- data.table(family = fam, n_deployments = length(ids),
                                n_contributing = a[["n"]],
                                z50 = a[["z50"]], width_deg = 4.394 * a[["scale"]])
      for (i in ids) responses[[i]] <- p[[i]]
    }
  } else {
    # No pooling: each tag stands or falls on its own window. Tags without a
    # usable geometry get nothing, which is the honest consequence of the recipe
    # and is exactly what this control arm exists to show.
    for (i in names(tags)) {
      responses[[i]] <- if (is.null(geom_fits[[i]])) NULL
                        else graft_response(geom_fits[[i]], scale_fits[[i]])
    }
  }
  list(responses = responses, geom_fits = geom_fits, scale_fits = scale_fits,
       pooled_diag = if (length(diag)) rbindlist(diag) else NULL,
       mode = mode,
       n_contributing = sum(!vapply(geom_fits, is.null, logical(1))),
       n_usable = sum(!vapply(responses, is.null, logical(1))))
}

likelihood_params <- function(r) c(1 / (r$max_light * 0.5), r$max_light, PROB_SLAB)
tag_light <- function(g, r) pmax(0, g$light$light - r$baseline)

# =============================================================================
# 4. SCORING -- the only scorer in the study
# =============================================================================

# Per-knot errors against Argos. argos_at() returns NA when the bracketing gap
# exceeds ARGOS_MAX_GAP_H, so knots inside a long gap are dropped rather than
# scored against an interpolation. Those gaps are not random -- a seal transmits
# when it surfaces -- so `frac_scored` is reported alongside every error and must
# be read with it.
score_knots <- function(fit_time, fit_lon, fit_lat, sd_lat, sd_lon, g) {
  tm <- as.POSIXct(fit_time, origin = "1970-01-01", tz = "UTC")
  tt <- argos_at(as.data.table(g$argos), tm)
  e_lat <- fit_lat - tt$lat
  e_lon <- dlon(lon360(fit_lon), tt$lon)
  km    <- gc_km(lon360(fit_lon), fit_lat, tt$lon, tt$lat)
  ok    <- is.finite(e_lat) & is.finite(sd_lat) & sd_lat > 0
  data.table(time = tm, lon = lon360(fit_lon), lat = fit_lat,
             truth_lon = tt$lon, truth_lat = tt$lat,
             err_km = km, err_lat = e_lat, err_lon = e_lon,
             sd_lat = sd_lat, sd_lon = sd_lon, scored = ok)
}

EMPTY_SUMMARY <- function(n_knots) data.table(
  n_knots = n_knots, n_scored = 0L, frac_scored = 0, median_km = NA_real_,
  mean_km = NA_real_, q90_km = NA_real_, bias_lat = NA_real_, rmse_lat = NA_real_,
  bias_lon = NA_real_, rmse_lon = NA_real_, lat_sd = NA_real_, lon_sd = NA_real_,
  cover_lat = NA_real_, cover_lon = NA_real_)

# Always returns the full schema. A fit that scores nothing must still produce a
# full-width row: a short row makes the checkpoint file unreadable, which is how
# a resumed run silently refits work it had already done.
summarise_knots <- function(K) {
  ok <- K[scored == TRUE]
  if (!nrow(ok)) return(EMPTY_SUMMARY(nrow(K)))
  data.table(
    n_knots    = nrow(K),
    n_scored   = nrow(ok),
    frac_scored= round(nrow(ok) / nrow(K), 3),
    median_km  = round(median(ok$err_km), 1),
    mean_km    = round(mean(ok$err_km), 1),
    q90_km     = round(quantile(ok$err_km, 0.9), 1),
    bias_lat   = round(mean(ok$err_lat), 3),
    rmse_lat   = round(sqrt(mean(ok$err_lat^2)), 3),
    bias_lon   = round(mean(ok$err_lon), 3),
    rmse_lon   = round(sqrt(mean(ok$err_lon^2)), 3),
    lat_sd     = round(mean(ok$sd_lat), 3),
    lon_sd     = round(mean(ok$sd_lon), 3),
    cover_lat  = round(mean(abs(ok$err_lat) <= 1.96 * ok$sd_lat), 3),
    cover_lon  = round(mean(abs(ok$err_lon) <= 1.96 * ok$sd_lon), 3))
}

# Posterior marginal standard deviations per knot, from the full cell posterior.
posterior_sd <- function(f) {
  gp <- grid_posterior(f)
  n <- nrow(gp$P)
  s_lat <- s_lon <- numeric(n)
  for (i in seq_len(n)) {
    w <- gp$P[i, ]; s <- sum(w)
    if (!is.finite(s) || s <= 0) { s_lat[i] <- s_lon[i] <- NA_real_; next }
    w <- w / s
    s_lon[i] <- sqrt(sum(w * (gp$lon - sum(w * gp$lon))^2))
    s_lat[i] <- sqrt(sum(w * (gp$lat - sum(w * gp$lat))^2))
  }
  list(lat = s_lat, lon = s_lon, n_cells = ncol(gp$P))
}

# =============================================================================
# 5. ENVIRONMENTAL LAYERS AND FUSION TERMS
# =============================================================================

load_layers <- function() {
  bp <- file.path(CACHE_DIR, "bathy.tif")
  if (!file.exists(bp)) stop("bathymetry raster absent: ", bp)
  # bathy.tif stores ETOPO ALTITUDE: negative below sea level. Both
  # bathy_source() and floor_rule() want seafloor depth POSITIVE DOWN with land
  # at zero, and neither can tell it has been handed the wrong sign.
  #
  # The smoke test found this the hard way. With altitude passed straight
  # through, floor_rule() computes deficit = required - (-4000) > 0 for every
  # ocean cell and returns -Inf, so the hard arm excluded the entire ocean
  # (log_z = -Inf, no scorable knot) while the soft arm penalised deep water
  # ~103 log-units more than the shore and drove the track 12.7 degrees south.
  # A manual reconstruction of the same constraint left 2,021-2,267 of 5,000
  # cells alive per knot, which is what made the discrepancy visible.
  elev       <- terra::rast(bp)
  depth_rast <- -1 * elev             # positive down
  depth_rast[depth_rast < 0] <- 0     # land is zero depth, not negative
  names(depth_rast) <- "depth_m"
  sea_rast   <- elev < 0              # TRUE over water

  # A hard land mask and a colony-pinned endpoint are incompatible on a coarse
  # grid, and the smoke test found it: Ano Nuevo is a beach at 237.67 E, 37.11 N,
  # which on a 1 degree grid falls in the cell centred at 237.5 E, 37.5 N --
  # inland California. Masking that cell with -Inf makes the pinned endpoint
  # impossible and the whole fit returns log_z = -Inf and no scorable knot.
  #
  # The fix is not to move the endpoint offshore, which would change the
  # light-only arm too and make the arms incomparable. It is to let the mask
  # admit the colony, which is the one place the animal is known to be at both
  # ends of the record. Cells within LAND_EXEMPT_KM of the colony are treated as
  # navigable; everywhere else the hard mask is unchanged.
  # mask_rule() reads nonzero as navigable, zero as masked and NA as "no
  # coverage, contribute nothing", so the layer must be numeric 0/1.
  xy <- terra::crds(sea_rast, na.rm = FALSE)
  near_colony <- gc_km(xy[, 1], xy[, 2], COLONY[["lon"]], COLONY[["lat"]]) <= LAND_EXEMPT_KM
  sea_vals <- as.numeric(terra::values(sea_rast))
  sea_vals[near_colony] <- 1
  sea_num <- elev
  terra::values(sea_num) <- sea_vals
  names(sea_num) <- "sea"

  list(depth = depth_rast, sea = sea_num, n_exempt = sum(near_colony, na.rm = TRUE))
}

memo_static <- function(fun) {
  cached <- NULL
  src <- function(lon, lat, date = NULL) {
    if (is.null(cached)) cached <<- fun(lon, lat, NULL)
    cached
  }
  class(src) <- c("tf_source", "function")
  src
}

# `hard` declares a cell impossible, so the engine skips the light likelihood
# there entirely -- that is where the speed-up comes from, and a soft constraint
# carrying the same information buys none of it. The knot's deepest dive is
# reduced by BATHY_MARGIN_M first so a 7 km bathymetry grid cannot exclude truth.
make_terms <- function(g, layers, hard) {
  out <- list()

  # The pinned endpoint knots are exempt from the depth constraint.
  #
  # Ano Nuevo's grid cell holds 16 m of water, so a knot in which the animal
  # dives 500 m excludes the very cell the fit is pinned to, and the whole track
  # collapses to log_z = -Inf. A seal swimming in and hauling out can easily dive
  # deep inside its final 12 h knot, so this is not a rare case -- and it would
  # kill a tag silently, hours into a long run.
  #
  # No auxiliary sensor can add information at a knot whose position is already
  # known, and one that contradicts it is simply wrong. Setting the tag value to
  # NA over those knots uses floor_rule()'s own missing-data path, which returns
  # zero rather than a penalty. The coastline mask is handled separately, by
  # exempting cells near the colony in load_layers().
  depth_val <- g$light$depth_max
  edge <- as.numeric(g$light$time) < as.numeric(min(g$light$time)) + STEP_HOURS * 3600 |
          as.numeric(g$light$time) > as.numeric(max(g$light$time)) - STEP_HOURS * 3600
  depth_val[edge] <- NA_real_

  out$bathy <- location_term(
    name = "bathy",
    tag  = data.frame(time = g$light$time, value = depth_val),
    reduce = if (hard) function(x) {
               if (all(is.na(x))) NA_real_ else max(0, max(x, na.rm = TRUE) - BATHY_MARGIN_M)
             } else "max",
    source = memo_static(bathy_source(rast = layers$depth)),
    rule   = if (hard) floor_rule() else floor_rule(sd = BATHY_SOFT_SD))
  out$sea <- location_term(
    name = "sea",
    source = memo_static(raster_source(layers$sea)),
    rule = mask_rule(penalty = if (hard) -Inf else SOFT_MASK_PENALTY))
  unname(out)
}

# =============================================================================
# 6. CHECKPOINTED OUTPUT
# =============================================================================

append_rows <- function(dt, path) {
  fwrite(dt, path, append = file.exists(path))
}

done_keys <- function(path, cols) {
  if (!file.exists(path)) return(character(0))
  d <- fread(path, colClasses = list(character = cols))
  do.call(paste, c(lapply(cols, function(cc) as.character(d[[cc]])), sep = "|"))
}
