# =============================================================================
# REFIT grid_light, KEEPING THE CELL POSTERIORS.
#
# The reporting run kept per-knot summaries and threw the posteriors away, which
# is right for a 12-hour run — 29 tags x ~470 knots x 5000 cells is 545 MB of
# doubles — but it means the exact `occupancy_map(method = "posterior")` cannot
# be built from its output. The Gaussian-kernel approximation from the reported
# marginal standard deviations stands in for it, and for a figure in the paper
# the exact surface is worth 2-3 hours.
#
# WHAT IS SAVED. Not the raw posteriors: the per-knot cell probabilities summed
# within CALENDAR MONTH, one matrix of (month x cell) per deployment. About
# 320 kB a tag instead of 19 MB. Any coarser period — season, whole track —
# is a sum of months. A period boundary that falls mid-month is the one thing
# this cannot serve, and no seasonal occupancy figure needs one.
#
# REPRODUCTION CONTROL. TwilightFreeGrid is deterministic, so this refit must
# return the reporting run's grid_light track bit-for-bit. It is asserted PER
# TAG on the full knot vectors, not on a summary. That distinction is the whole
# lesson of fit_final_29.R, whose declared control matched its reference on two
# means while missing every one of 29 tags.
#
#   Rscript analysis/study/refit_posteriors.R
#
# Resumable: checkpointed per deployment.
# =============================================================================

source("analysis/study/study_common.R")

# The package's own accumulator is internal; a local copy keeps this script
# runnable against an installed build without reaching into the namespace.
accumulate_cells <- function(cells, mass, n) {
  out <- numeric(n)
  o <- order(cells); cs <- cells[o]; ms <- mass[o]
  r <- rle(cs); ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1L
  out[r$values] <- vapply(seq_along(r$values),
                          function(i) sum(ms[starts[i]:ends[i]]), 0)
  out
}

OUT   <- file.path(OUT_DIR, "posteriors")
F_CHK <- file.path(OUT_DIR, "refit_check.csv")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

REF <- fread(file.path(OUT_DIR, "study_knots.csv"), colClasses = list(character = "id"))
REF <- REF[arm == "grid_light"]
# fwrite wrote these as ISO strings with no zone; force UTC rather than let the
# session locale decide, or the timestamp assertion below compares two clocks.
REF[, time := as.POSIXct(time, tz = "UTC")]

REFIT_TAGS <- as.integer(Sys.getenv("REFIT_TAGS", "0"))   # 0 = all; smoke-test control
if (!nrow(REF)) stop("no grid_light knots in ", OUT_DIR, " -- run run_study.R first")

started <- Sys.time()
log_line("refit starting | %d deployments to reproduce", uniqueN(REF$id))

tags   <- build_panel()
cal    <- fit_calibration(tags, "study")
grid   <- terra::rast(xmin = DOMAIN[["xmin"]], xmax = DOMAIN[["xmax"]],
                      ymin = DOMAIN[["ymin"]], ymax = DOMAIN[["ymax"]],
                      resolution = GRID_CELL_DEG, crs = "EPSG:4326")
terra::values(grid) <- 1
N_CELLS <- terra::ncell(grid)
probe_start <- machine_probe()
log_line("grid %d cells | machine probe %.2f s", N_CELLS, probe_start$median)

done <- if (file.exists(F_CHK)) as.character(fread(F_CHK)$id) else character(0)
log_line("%d of %d already on disk", length(done), length(tags))

todo <- setdiff(names(tags), done)
if (REFIT_TAGS > 0) {
  todo <- head(todo, REFIT_TAGS)
  log_line("SMOKE: limiting to %d deployment(s)", length(todo))
}

for (tg in todo) {
  if (tg %in% done) next
  g <- tags[[tg]]; r <- cal$responses[[tg]]
  if (is.null(r)) { log_line("%s: no calibration, skipped", tg); next }
  log_line("refit %s (%d obs, %.0f d)", tg, nrow(g$light), g$days)

  t <- system.time(invisible(capture.output(
    f <- TwilightFreeGrid(g$light$time, tag_light(g, r), grid = grid,
      start_lon = g$p0[1], start_lat = g$p0[2],
      end_lon   = g$p1[1], end_lat   = g$p1[2],
      step_hours = STEP_HOURS, diffusion = DIFFUSION,
      calibration = r$calibration, likelihood_params = likelihood_params(r),
      shade_ratio = SHADE_RATIO, area_correction = AREA_CORRECTION,
      terms = list()))))[["elapsed"]]

  # ---- the control, per tag, on the full vectors ----------------------------
  ref <- REF[id == tg][order(time)]
  tm  <- as.POSIXct(f$fit$time, origin = "1970-01-01", tz = "UTC")
  ok_n   <- nrow(ref) == length(tm)
  d_lon  <- if (ok_n) max(abs(lon360(f$fit$lon) - lon360(ref$lon))) else NA_real_
  d_lat  <- if (ok_n) max(abs(f$fit$lat - ref$lat)) else NA_real_
  d_time <- if (ok_n) max(abs(as.numeric(tm) - as.numeric(ref$time))) else NA_real_
  # POSITIONS must be bit-exact: the engine is deterministic and this is the same
  # fit. TIMESTAMPS are compared to within a second, because the reference came
  # back through a CSV. Knots are spaced 43158.29 s apart here, not a round 12 h
  # -- the engine divides the record span into equal intervals -- and ISO output
  # does not carry that sub-second part. The tolerance is on the serialisation,
  # not on the fit, and the observed maximum is recorded per tag so a real drift
  # cannot hide behind it.
  exact  <- isTRUE(ok_n && d_lon == 0 && d_lat == 0 && d_time < 1)
  if (!exact)
    stop("refit does not reproduce grid_light for ", tg,
         ": n ", nrow(ref), " vs ", length(tm),
         ", max |dlon| ", format(d_lon), ", max |dlat| ", format(d_lat),
         ", max |dt| ", format(d_time), " s",
         ". The refit and the reporting run are not the same fit; stop and find out why.")

  # ---- accumulate the posterior within calendar month -----------------------
  gp <- grid_posterior(f)
  # grid_posterior returns cells in the engine's own order; map them once onto
  # the raster's cell index so the saved matrix is addressable by cell number.
  cell <- terra::cellFromXY(grid, cbind(gp$lon, gp$lat))
  mon  <- format(tm, "%Y-%m")
  months <- sort(unique(mon))
  P <- matrix(0, nrow = length(months), ncol = N_CELLS,
              dimnames = list(months, NULL))
  nk <- integer(length(months)); names(nk) <- months
  for (k in seq_len(nrow(gp$P))) {
    w <- gp$P[k, ]; s <- sum(w)
    if (!is.finite(s) || s <= 0) next
    w <- w / s
    i <- match(mon[k], months)
    keep <- !is.na(cell)
    P[i, ] <- P[i, ] + accumulate_cells(cell[keep], w[keep], N_CELLS)
    nk[i] <- nk[i] + 1L
  }

  saveRDS(list(id = tg, season = g$season, family = g$family, serial = g$serial,
               months = months, knots_per_month = nk, P = P,
               track = data.frame(time = tm, lon = lon360(f$fit$lon), lat = f$fit$lat),
               grid = list(ext = as.vector(terra::ext(grid)),
                           res = terra::res(grid), ncell = N_CELLS)),
          file.path(OUT, paste0(tg, ".rds")), compress = "xz")

  append_rows(data.table(id = tg, season = g$season, n_knots = length(tm),
                         n_months = length(months), exact = exact,
                         max_dlon = d_lon, max_dlat = d_lat, max_dt_s = d_time,
                         seconds = round(t, 2)), F_CHK)
  rm(f, gp, P); gc(verbose = FALSE)
}

chk <- fread(F_CHK)
probe_end <- machine_probe()
log_line("machine probe end %.2f s (start %.2f, drift %+.0f%%)",
         probe_end$median, probe_start$median,
         100 * (probe_end$median / probe_start$median - 1))
cat("\n=== REPRODUCTION CONTROL ===\n")
cat(sprintf("deployments refitted            : %d\n", nrow(chk)))
cat(sprintf("bit-exact against grid_light    : %d\n", sum(chk$exact)))
cat(sprintf("max |dlon| over all deployments : %s\n", format(max(chk$max_dlon))))
cat(sprintf("max |dlat| over all deployments : %s\n", format(max(chk$max_dlat))))
if (all(chk$exact)) {
  cat("PASS -- the saved posteriors belong to exactly the tracks the study reported.\n")
} else {
  cat("FAIL -- do not use these posteriors.\n")
}
log_line("refit finished in %.2f h",
         as.numeric(difftime(Sys.time(), started, units = "hours")))
