# =============================================================================
# THE REPORTING RUN -- 29 northern elephant seal deployments, start to finish.
#
#   Rscript analysis/study/preflight.R      # first, and read it
#   Rscript analysis/study/run_study.R      # then this
#
# Resumable: every fit is checkpointed on (arm, id) and appended to disk as it
# finishes, so an interrupted run continues rather than restarting. Persist to
# disk, never to the terminal -- a `Select-Object -Last 5` in a launcher once
# truncated a five-hour run's entire output.
#
# PHASES=G,H,F selects phases (grid / hierarchical / FFBS). Default all three.
#
# ---------------------------------------------------------------------------
# THE ARMS, and what each contrast identifies
#
#   grid_light          reference. Light only, study calibration.
#   grid_fusion         + HARD bathymetry and coastline. Against grid_light this
#                       measures what the tag's own depth channel is worth, in
#                       accuracy AND in wall clock: a hard constraint lets the
#                       engine skip the light likelihood on excluded cells, a
#                       soft one carrying the same information cannot.
#   grid_cal_untrunc    calibration geometry fitted on a fixed 7-day window with
#                       NO departure truncation, still pooled. Against
#                       grid_light this isolates the truncation -- the animals
#                       leave between day 0 and day 17, so a fixed window puts
#                       at-sea light into 16 of 29 calibrations.
#   grid_cal_pertag     as above, and NOT pooled. Against grid_cal_untrunc this
#                       isolates pooling within tag family, which is also the
#                       only thing that calibrates the 9 deployments with too
#                       little haul-out to calibrate themselves.
#   hier_light/fusion   the panel fitted jointly, movement scale pooled.
#   ffbs_light/fusion   continuous space, no grid and no coarse mesh, on the
#                       three tags with the densest Argos coverage.
#
# Together grid_light -> grid_cal_untrunc -> grid_cal_pertag is the calibration
# decomposition the manuscript needs, on the same 29 tags and one scorer.
# =============================================================================

source("analysis/study/study_common.R")

PHASES <- strsplit(Sys.getenv("PHASES", "G,H,F"), ",")[[1]]

# Smoke-test controls. SMOKE_TAGS limits the panel and SMOKE_DAYS truncates each
# light record, so the whole path -- every arm, every engine, the scorer and the
# output files -- can be exercised in minutes before committing to the real run.
# Both default to off. When either is set the outputs go to a separate directory
# so a smoke test can never be mistaken for the study.
SMOKE_TAGS <- as.integer(Sys.getenv("SMOKE_TAGS", "0"))
SMOKE_DAYS <- as.integer(Sys.getenv("SMOKE_DAYS", "0"))
SMOKE <- SMOKE_TAGS > 0 || SMOKE_DAYS > 0
if (SMOKE) {
  OUT_DIR <- file.path("analysis/output/study", paste0(STUDY_VERSION, "_smoke"))
  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
}

F_TAG   <- file.path(OUT_DIR, "study_tags.csv")
F_KNOT  <- file.path(OUT_DIR, "study_knots.csv")
F_TIME  <- file.path(OUT_DIR, "study_timing.csv")
F_PROV  <- file.path(OUT_DIR, "study_provenance.rds")

started <- Sys.time()
log_line("study %s starting | phases %s", STUDY_VERSION, paste(PHASES, collapse = "+"))

# ---- fixed inputs ------------------------------------------------------------
tags   <- build_panel()
if (SMOKE_TAGS > 0) tags <- tags[seq_len(min(SMOKE_TAGS, length(tags)))]
if (SMOKE_DAYS > 0) tags <- lapply(tags, function(g) {
  g$light <- g$light[time <= min(time) + SMOKE_DAYS * 86400]
  g$argos <- g$argos[time <= min(g$light$time) + SMOKE_DAYS * 86400]
  g$t1    <- max(g$light$time)
  g$days  <- as.numeric(difftime(g$t1, g$t0, units = "days"))
  g
})
if (SMOKE) log_line("SMOKE TEST: %d tags, %s days -- results go to %s and are NOT the study",
                    length(tags), if (SMOKE_DAYS > 0) SMOKE_DAYS else "all", OUT_DIR)
layers <- load_layers()
grid   <- terra::rast(xmin = DOMAIN[["xmin"]], xmax = DOMAIN[["xmax"]],
                      ymin = DOMAIN[["ymin"]], ymax = DOMAIN[["ymax"]],
                      resolution = GRID_CELL_DEG, crs = "EPSG:4326")
terra::values(grid) <- 1
N_CELLS <- terra::ncell(grid)

cal <- list(
  study          = fit_calibration(tags, "study"),
  pooled_untrunc = fit_calibration(tags, "pooled_untrunc"),
  pertag_untrunc = fit_calibration(tags, "pertag_untrunc"))

log_line("panel %d deployments | grid %d cells | calibration modes built",
         length(tags), N_CELLS)

probe_start <- machine_probe()
log_line("machine probe %.2f s on %d cores", probe_start$median, probe_start$cores)

# =============================================================================
# PHASE G -- grid HMM arms
# =============================================================================

GRID_ARMS <- data.table(
  arm  = c("grid_light", "grid_fusion", "grid_cal_untrunc", "grid_cal_pertag"),
  cal  = c("study",      "study",       "pooled_untrunc",   "pertag_untrunc"),
  hard = c(FALSE,        TRUE,          FALSE,              FALSE),
  fuse = c(FALSE,        TRUE,          FALSE,              FALSE))

run_grid <- function() {
  done <- done_keys(F_TAG, c("arm", "id"))
  todo <- nrow(GRID_ARMS) * length(tags)
  log_line("PHASE G: %d fits, %d already on disk", todo, sum(done != ""))

  for (ai in seq_len(nrow(GRID_ARMS))) {
    A <- GRID_ARMS[ai]
    responses <- cal[[A$cal]]$responses
    for (tg in names(tags)) {
      if (paste(A$arm, tg, sep = "|") %in% done) next
      g <- tags[[tg]]
      r <- responses[[tg]]
      if (is.null(r)) {
        # A control arm is allowed to fail on a tag; that failure IS the result.
        log_line("%s / %s : no usable calibration, recorded as a failure", A$arm, tg)
        append_rows(data.table(arm = A$arm, id = tg, season = g$season,
                               family = g$family, serial = g$serial %||% NA_character_,
                               calibrated = FALSE, n_knots = NA_integer_,
                               n_scored = NA_integer_, frac_scored = NA_real_,
                               median_km = NA_real_, mean_km = NA_real_, q90_km = NA_real_,
                               bias_lat = NA_real_, rmse_lat = NA_real_,
                               bias_lon = NA_real_, rmse_lon = NA_real_,
                               lat_sd = NA_real_, lon_sd = NA_real_,
                               cover_lat = NA_real_, cover_lon = NA_real_,
                               log_z = NA_real_), F_TAG)
        next
      }
      terms <- if (A$fuse) make_terms(g, layers, hard = A$hard) else list()
      log_line("%s / %s (%d obs, %.0f d)", A$arm, tg, nrow(g$light), g$days)
      t <- system.time(invisible(capture.output(
        f <- TwilightFreeGrid(g$light$time, tag_light(g, r), grid = grid,
          start_lon = g$p0[1], start_lat = g$p0[2],
          end_lon   = g$p1[1], end_lat   = g$p1[2],
          step_hours = STEP_HOURS, diffusion = DIFFUSION,
          calibration = r$calibration, likelihood_params = likelihood_params(r),
          shade_ratio = SHADE_RATIO, area_correction = AREA_CORRECTION,
          terms = terms))))[["elapsed"]]

      sd <- posterior_sd(f)
      K  <- score_knots(f$fit$time, f$fit$lon, f$fit$lat, sd$lat, sd$lon, g)
      S  <- summarise_knots(K)
      if (is.null(S)) { log_line("  no scorable knot -- recorded"); S <- data.table(n_knots = nrow(K)) }

      append_rows(cbind(data.table(arm = A$arm, id = tg, season = g$season,
                                   family = g$family, serial = g$serial %||% NA_character_,
                                   calibrated = TRUE),
                        S, data.table(log_z = round(f$log_z, 2))), F_TAG)
      append_rows(cbind(data.table(arm = A$arm, id = tg), K), F_KNOT)
      append_rows(data.table(arm = A$arm, id = tg, engine = "TwilightFreeGrid",
                             seconds = round(t, 2), n_knots = nrow(K),
                             n_cells = N_CELLS, n_obs = nrow(g$light),
                             sec_per_knot = round(t / nrow(K), 4),
                             sec_per_knot_cell = signif(t / (nrow(K) * N_CELLS), 4),
                             sec_per_tag_day = round(t / g$days, 4)), F_TIME)
      rm(f, sd, K); gc(verbose = FALSE)
      done <- c(done, paste(A$arm, tg, sep = "|"))
    }
  }
}

# =============================================================================
# PHASE H -- hierarchical
#
# TwilightFreeHier applies ONE calibration across the panel, so the per-tag
# responses are pooled by their median. Since the geometry is already pooled
# within family, this pools only what is left, the per-tag intensity scale.
# =============================================================================

run_hier <- function() {
  responses <- cal$study$responses
  ok <- Filter(Negate(is.null), responses)
  cal_panel <- c(median(vapply(ok, function(r) r$calibration[1], 0)),
                 median(vapply(ok, function(r) r$calibration[2], 0)))
  m <- median(vapply(ok, function(r) r$max_light, 0))
  lp_panel <- c(1 / (m * 0.5), m, PROB_SLAB)

  P <- list(
    data = setNames(lapply(tags, function(g)
      data.frame(Date = g$light$time, Light = tag_light(g, responses[[g$id]]))),
      names(tags)),
    locations = data.frame(id = names(tags),
      deploy_lon = COLONY[["lon"]], deploy_lat = COLONY[["lat"]],
      retrieve_lon = COLONY[["lon"]], retrieve_lat = COLONY[["lat"]],
      row.names = NULL))

  pad <- c(lon = 74, lat = 26)
  done <- done_keys(F_TAG, c("arm", "id"))

  for (arm in c("hier_light", "hier_fusion")) {
    if (any(grepl(paste0("^", arm, "\\|"), done))) { log_line("%s already on disk", arm); next }
    terms <- if (identical(arm, "hier_fusion"))
      function(id, df) make_terms(tags[[id]], layers, hard = FALSE) else NULL
    log_line("PHASE H: %s over %d deployments, %d sweeps", arm, length(tags), HIER_SWEEPS)
    t <- system.time(
      H <- TwilightFreeHier(data = P$data, locations = P$locations, terms = terms,
        calibration = cal_panel, likelihood_params = lp_panel,
        step_hours = STEP_HOURS, coarse_res = HIER_COARSE_RES, mesh_pad = pad,
        surrogate_diffusion = HIER_SURROGATE_DIFFUSION,
        sweeps = HIER_SWEEPS, burn = HIER_BURN, thin = HIER_THIN,
        metric = "spherical", seed = SEED))[["elapsed"]]

    for (tg in names(H$tracks)) {
      tr <- H$tracks[[tg]]; g <- tags[[tg]]
      K <- score_knots(as.numeric(tr$time), tr$lon, tr$lat, tr$lat_sd, tr$lon_sd, g)
      S <- summarise_knots(K)
      if (is.null(S)) next
      append_rows(cbind(data.table(arm = arm, id = tg, season = g$season,
                                   family = g$family, serial = g$serial %||% NA_character_,
                                   calibrated = TRUE),
                        S, data.table(log_z = NA_real_)), F_TAG)
      append_rows(cbind(data.table(arm = arm, id = tg), K), F_KNOT)
    }
    append_rows(data.table(arm = arm, id = "ALL", engine = "TwilightFreeHier",
                           seconds = round(t, 2), n_knots = NA_integer_,
                           n_cells = NA_integer_,
                           n_obs = sum(vapply(tags, function(g) nrow(g$light), 0)),
                           sec_per_knot = NA_real_, sec_per_knot_cell = NA_real_,
                           sec_per_tag_day = round(t / sum(vapply(tags, function(g) g$days, 0)), 4)),
                F_TIME)
    saveRDS(list(population = H$population, movement = H$movement,
                 persistence = H$persistence),
            file.path(OUT_DIR, paste0(arm, "_population.rds")))
    rm(H); gc(verbose = FALSE)
  }
}

# =============================================================================
# PHASE F -- FFBS on the densest-Argos subset
#
# The SMC engine rasterises auxiliary terms onto a background grid that defaults
# to 2 degrees spanning -180..180. This study works in 0-360, so without an
# explicit mask every particle falls outside that extent, is clamped to the edge
# column and receives a constant auxiliary value -- the terms would silently do
# nothing. The mask below is all-ones over the study domain, so it is the same
# mild "stay in the box" constraint in BOTH arms and cannot confound the
# light-versus-fusion contrast.
# =============================================================================

run_ffbs <- function() {
  responses <- cal$study$responses
  dens <- vapply(tags, function(g) nrow(g$argos) / max(g$days, 1), 0)
  ids  <- names(sort(dens, decreasing = TRUE))[seq_len(min(FFBS_N, length(tags)))]
  log_line("PHASE F: FFBS subset %s", paste(ids, collapse = ", "))

  mask <- terra::rast(xmin = DOMAIN[["xmin"]], xmax = DOMAIN[["xmax"]],
                      ymin = DOMAIN[["ymin"]], ymax = DOMAIN[["ymax"]],
                      resolution = 1, crs = "EPSG:4326")
  terra::values(mask) <- 1
  done <- done_keys(F_TAG, c("arm", "id"))

  for (arm in c("ffbs_light", "ffbs_fusion")) {
    for (tg in ids) {
      if (paste(arm, tg, sep = "|") %in% done) next
      g <- tags[[tg]]; r <- responses[[tg]]
      terms <- if (identical(arm, "ffbs_fusion")) make_terms(g, layers, hard = FALSE) else list()
      log_line("%s / %s", arm, tg)
      t <- system.time(
        f <- TwilightFreeSMC(date_time = g$light$time, light = tag_light(g, r),
          calibration = r$calibration, likelihood_params = likelihood_params(r),
          n_particles = FFBS_PARTICLES,
          start_lon = g$p0[1], start_lat = g$p0[2],
          end_lon = g$p1[1], end_lat = g$p1[2],
          method = "ffbs", step_hours = STEP_HOURS, diffusion = DIFFUSION,
          spatial_mask = mask, terms = terms, seed = SEED))[["elapsed"]]

      K <- score_knots(as.numeric(f$knot_times), f$lon, f$lat, f$lat_sd, f$lon_sd, g)
      S <- summarise_knots(K)
      if (!is.null(S)) {
        append_rows(cbind(data.table(arm = arm, id = tg, season = g$season,
                                     family = g$family, serial = g$serial %||% NA_character_,
                                     calibrated = TRUE),
                          S, data.table(log_z = NA_real_)), F_TAG)
        append_rows(cbind(data.table(arm = arm, id = tg), K), F_KNOT)
      }
      append_rows(data.table(arm = arm, id = tg, engine = "TwilightFreeSMC(ffbs)",
                             seconds = round(t, 2), n_knots = nrow(K),
                             n_cells = FFBS_PARTICLES, n_obs = nrow(g$light),
                             sec_per_knot = round(t / nrow(K), 4),
                             sec_per_knot_cell = signif(t / (nrow(K) * FFBS_PARTICLES), 4),
                             sec_per_tag_day = round(t / g$days, 4)), F_TIME)
      rm(f, K); gc(verbose = FALSE)
    }
  }
}

# =============================================================================
if ("G" %in% PHASES) run_grid()
if ("H" %in% PHASES) run_hier()
if ("F" %in% PHASES) run_ffbs()

probe_end <- machine_probe()
saveRDS(list(
  study_version = STUDY_VERSION,
  started = started, finished = Sys.time(),
  config = mget(c("COLONY", "DOMAIN", "GRID_CELL_DEG", "CAL_GEOM_MAX_DAYS",
                  "CAL_GEOM_MIN_OBS", "CAL_SCALE_DAYS", "ARGOS_MAX_GAP_H",
                  "ENDPOINT_RULE", "STEP_HOURS", "DIFFUSION", "SEED", "PROB_SLAB",
                  "SHADE_RATIO", "AREA_CORRECTION", "HIER_SWEEPS", "FFBS_N",
                  "FFBS_PARTICLES", "USE_SST"), envir = globalenv()),
  arms = GRID_ARMS,
  panel = rbindlist(lapply(tags, function(g) data.table(
    id = g$id, season = g$season, family = g$family,
    serial = g$serial %||% NA_character_, days = round(g$days, 1),
    n_obs = nrow(g$light), n_argos = nrow(g$argos)))),
  calibration = lapply(cal, function(x) list(mode = x$mode, pooled = x$pooled_diag,
                                             n_contributing = x$n_contributing,
                                             n_usable = x$n_usable)),
  machine = list(start = probe_start, end = probe_end),
  git = tryCatch(system("git rev-parse HEAD", intern = TRUE), error = function(e) NA),
  session = utils::sessionInfo()), F_PROV)

log_line("machine probe end %.2f s (start %.2f, drift %+.0f%%)",
         probe_end$median, probe_start$median,
         100 * (probe_end$median / probe_start$median - 1))
log_line("study finished in %.2f h", as.numeric(difftime(Sys.time(), started, units = "hours")))
cat(sprintf("\nwritten to %s\n  study_tags.csv  study_knots.csv  study_timing.csv  study_provenance.rds\n",
            OUT_DIR))
