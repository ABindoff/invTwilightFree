# =============================================================================
# PREFLIGHT -- run this, and read it, BEFORE launching the study.
#
# Every check here is an assertion that has already failed once in this project.
# The run is long; a failure discovered at hour nine costs a day.
#
#   Rscript analysis/study/preflight.R
#
# Set COST_PROBE=0 to skip the single timed fit at the end (it takes ~7 min and
# is what turns the cost estimate from a guess into a measurement).
# =============================================================================

source("analysis/study/study_common.R")

fails <- 0L
warns <- 0L
ok   <- function(...) cat(sprintf("  PASS  %s\n", sprintf(...)))
bad  <- function(...) { fails <<- fails + 1L; cat(sprintf("  FAIL  %s\n", sprintf(...))) }
warn <- function(...) { warns <<- warns + 1L; cat(sprintf("  WARN  %s\n", sprintf(...))) }

cat("\n=========================================================\n")
cat(sprintf("PREFLIGHT  study %s  %s\n", STUDY_VERSION,
            format(Sys.time(), "%Y-%m-%d %H:%M")))
cat("=========================================================\n")

# ---- 1. the build that will actually run ------------------------------------
# The trap this prevents cost sixteen hours in August: R CMD INSTALL was run,
# THEN the source was edited, then pkgload::load_all() passed (it reads source)
# and a long job was launched against the stale INSTALLED binary. This study uses
# the installed package deliberately, and asserts its features here.
cat("\n1. INSTALLED BUILD\n")
cat(sprintf("     invTwilightFree %s | %s\n",
            as.character(packageVersion("invTwilightFree")), R.version.string))
need_grid <- c("area_correction", "shade_ratio", "lambda_scale", "terms", "calibration")
missing <- setdiff(need_grid, names(formals(TwilightFreeGrid)))
if (length(missing)) { bad("TwilightFreeGrid lacks: %s", paste(missing, collapse = ", "))
} else { ok("TwilightFreeGrid exposes every argument the study passes") }
for (fn in c("pool_light_responses", "fit_light_response", "grid_posterior",
             "location_term", "floor_rule", "mask_rule", "bathy_source",
             "raster_source", "TwilightFreeHier", "TwilightFreeSMC")) {
  if (!exists(fn, asNamespace("invTwilightFree"))) bad("%s() not exported", fn)
}
if (!fails) ok("every entry point the study calls is present")

gs <- system("git rev-parse --short HEAD", intern = TRUE, ignore.stderr = TRUE)
dirty <- length(system("git status --porcelain -- R src DESCRIPTION NAMESPACE",
                       intern = TRUE, ignore.stderr = TRUE))
if (dirty) { warn("package sources are modified relative to HEAD (%s); the installed binary may not match the tree", gs[1])
} else { ok("package tree clean at %s", gs[1]) }

# ---- 2. inputs ---------------------------------------------------------------
cat("\n2. INPUTS\n")
need_files <- c(file.path(CACHE_DIR, "archives_v3_30min.rds"),
                file.path(CACHE_DIR, "fvilches_archives_v1_30min.rds"),
                file.path(CACHE_DIR, "fvilches_argos_v1.rds"),
                file.path(CACHE_DIR, "bathy.tif"),
                "scratch/nes_calibration/fvilches_manifest.csv",
                "data/nes_untracked")
for (f in need_files) if (!file.exists(f)) bad("missing: %s", f)
if (!fails) ok("every input present")

# ---- 3. the panel ------------------------------------------------------------
cat("\n3. PANEL\n")
tags <- build_panel()
cat(sprintf("     %d deployments assembled\n", length(tags)))
if (length(tags) != 29) { bad("expected 29 deployments, assembled %d", length(tags))
} else { ok("29 deployments") }
tb <- table(vapply(tags, function(g) g$season, ""),
            vapply(tags, function(g) g$family, ""))
cat("     season x family:\n")
print(tb)
cat(sprintf("     distinct loggers: %d (deployments are NOT independent at logger level)\n",
            length(unique(na.omit(vapply(tags, function(g) g$serial, ""))))))

# ---- 4. the search domain does not touch the animals -------------------------
cat("\n4. SEARCH DOMAIN\n")
allfix <- rbindlist(lapply(tags, function(g) g$argos[, .(lon, lat)]))
m <- c(lon_lo = min(allfix$lon) - DOMAIN[["xmin"]],
       lon_hi = DOMAIN[["xmax"]] - max(allfix$lon),
       lat_lo = min(allfix$lat) - DOMAIN[["ymin"]],
       lat_hi = DOMAIN[["ymax"]] - max(allfix$lat))
cat(sprintf("     margins (deg): lon %.1f / %.1f, lat %.1f / %.1f\n",
            m[1], m[2], m[3], m[4]))
if (any(m < 0)) { bad("the domain does not contain every Argos fix")
} else if (any(m < 1)) { warn("a domain edge is within 1 degree of an animal")
} else { ok("domain contains every fix with at least 1 degree of margin") }

# ---- 5. endpoints carry no ground truth --------------------------------------
cat("\n5. ENDPOINT RULE\n")
dev <- rbindlist(lapply(tags, function(g) data.table(
  id = g$id, season = g$season,
  dep_km = gc_km(g$argos_deploy[1],  g$argos_deploy[2],  COLONY[["lon"]], COLONY[["lat"]]),
  rec_km = gc_km(g$argos_recover[1], g$argos_recover[2], COLONY[["lon"]], COLONY[["lat"]]))))
cat(sprintf("     both ends pinned at the colony (%.2f E, %.2f N) for all %d deployments\n",
            COLONY[["lon"]], COLONY[["lat"]], length(tags)))
cat("     deviation of the Argos-derived positions from the colony (diagnostic only):\n")
cat(sprintf("       deploy : median %.1f km, q90 %.1f, max %.1f\n",
            median(dev$dep_km), quantile(dev$dep_km, .9), max(dev$dep_km)))
cat(sprintf("       recover: median %.1f km, q90 %.1f, max %.1f\n",
            median(dev$rec_km), quantile(dev$rec_km, .9), max(dev$rec_km)))
far <- dev[dep_km > ENDPOINT_MAX_DEVIATION_KM | rec_km > ENDPOINT_MAX_DEVIATION_KM]
if (nrow(far)) {
  warn("%d deployment(s) have an Argos endpoint beyond %d km of the colony; the colony pin may be wrong for these and they should be reported separately",
       nrow(far), ENDPOINT_MAX_DEVIATION_KM)
  cat("       affected: ", paste(far$id, collapse = ", "), "\n")
} else ok("no deployment's Argos endpoints contradict the colony pin")
fwrite(dev, file.path(OUT_DIR, "preflight_endpoints.csv"))

# ---- 6. the scorer -----------------------------------------------------------
# The failure this catches: two definitions of ground truth were live in the
# campaign at once. argos_at() returns NA across a long gap; a local
# approx(rule = 2) variant never does, and extrapolates past both track ends.
# Numbers from scripts using different scorers are not comparable, which is how a
# reproduction control passed on the mean while missing every tag.
cat("\n6. SCORER\n")
syn <- data.table(
  time = as.POSIXct("2021-06-01", tz = "UTC") + c(0, 3600, 30 * 3600, 31 * 3600),
  lon = c(237, 237.5, 240, 240.5), lat = c(37, 37.2, 38, 38.2))
probe <- as.POSIXct("2021-06-01", tz = "UTC") + c(1800, 15 * 3600, 100 * 3600)
got <- argos_at(syn, probe)
if (!is.finite(got$lat[1])) { bad("argos_at returned NA inside a 1 h gap")
} else { ok("interpolates inside a short gap") }
if (is.finite(got$lat[2])) { bad("argos_at interpolated across a 29 h gap; the >%d h rule is not live", ARGOS_MAX_GAP_H)
} else { ok("returns NA across a %d h gap", ARGOS_MAX_GAP_H) }
if (is.finite(got$lat[3])) { bad("argos_at extrapolated past the last fix")
} else { ok("returns NA past the end of the record") }

# ---- 7. the engine is deterministic ------------------------------------------
cat("\n7. DETERMINISM\n")
g <- tags[[1]]
short <- g$light[time <= min(time) + 20 * 86400]
cal0 <- fit_light_response(short$time, short$light, COLONY[["lon"]], COLONY[["lat"]])
gsmall <- terra::rast(xmin = 225, xmax = 245, ymin = 30, ymax = 45,
                      resolution = 1, crs = "EPSG:4326")
terra::values(gsmall) <- 1
run1 <- function() {
  invisible(capture.output(
    f <- TwilightFreeGrid(short$time, pmax(0, short$light - cal0$baseline), grid = gsmall,
                     start_lon = COLONY[["lon"]], start_lat = COLONY[["lat"]],
                     end_lon = COLONY[["lon"]], end_lat = COLONY[["lat"]],
                     step_hours = STEP_HOURS, diffusion = DIFFUSION,
                     calibration = cal0$calibration,
                     likelihood_params = likelihood_params(cal0),
                     area_correction = AREA_CORRECTION)))
  f
}
a <- run1(); b <- run1()
if (identical(a$fit$lat, b$fit$lat) && identical(a$fit$lon, b$fit$lon) &&
    identical(a$log_z, b$log_z)) { ok("grid HMM is bit-identical across repeat calls")
} else { bad("grid HMM is NOT deterministic") }

# ---- 8. the area double-count guard is live ----------------------------------
cat("\n8. AREA CORRECTION GUARD\n")
guard <- try(suppressWarnings(TwilightFreeGrid(
  short$time, pmax(0, short$light - cal0$baseline), grid = gsmall,
  start_lon = COLONY[["lon"]], start_lat = COLONY[["lat"]],
  end_lon = COLONY[["lon"]], end_lat = COLONY[["lat"]],
  step_hours = STEP_HOURS, diffusion = DIFFUSION,
  calibration = cal0$calibration, likelihood_params = likelihood_params(cal0),
  area_correction = TRUE, terms = list(area_prior()))), silent = TRUE)
if (inherits(guard, "try-error")) { ok("passing area_prior() alongside area_correction errors, as it must")
} else { bad("the cos-lat factor can still be applied twice silently") }

# ---- 9. the three calibration modes ------------------------------------------
cat("\n9. CALIBRATION\n")
for (mode in c("study", "pooled_untrunc", "pertag_untrunc")) {
  cal <- fit_calibration(tags, mode)
  cat(sprintf("     %-15s contributing %2d/%d | usable responses %2d/%d\n",
              mode, cal$n_contributing, length(tags), cal$n_usable, length(tags)))
  if (!is.null(cal$pooled_diag)) {
    for (i in seq_len(nrow(cal$pooled_diag))) {
      r <- cal$pooled_diag[i]
      cat(sprintf("       family %-8s n=%2d from %2d  z50 %.2f deg  width %.1f deg\n",
                  r$family, r$n_deployments, r$n_contributing, r$z50, r$width_deg))
      if (r$width_deg < 5 || r$width_deg > 60)
        bad("pooled transition width %.1f deg is physically implausible (%s, %s)",
            r$width_deg, mode, r$family)
    }
  }
  if (identical(mode, "study") && cal$n_usable != length(tags))
    bad("the study arm must give every deployment a response; %d have none", length(tags) - cal$n_usable)
}
ok("calibration modes built")

# ---- 10. machine state -------------------------------------------------------
cat("\n10. MACHINE\n")
p <- machine_probe()
cat(sprintf("     %d logical cores | probe %.2f s (min %.2f, max %.2f)\n",
            p$cores, p$median, p$min, p$max))
cat("     The engine is single-threaded: no rayon in Cargo.toml, no parallel:: in R.\n")
cat("     Free cores do not change the fit; contention does, by 1.33-1.56x measured.\n")
if (p$max / p$min > 1.15) {
  warn("probe repeats vary by %.0f%%; something else is using this machine",
       100 * (p$max / p$min - 1))
} else { ok("probe is stable; the machine looks quiescent") }

# ---- 11. cost ----------------------------------------------------------------
if (!identical(Sys.getenv("COST_PROBE"), "0")) {
  cat("\n11. COST (one real fit, timed)\n")
  cal <- fit_calibration(tags, "study")
  lens <- vapply(tags, function(g) nrow(g$light), 0)
  mid  <- names(tags)[order(lens)][ceiling(length(tags) / 2)]
  gm   <- tags[[mid]]; rm_ <- cal$responses[[mid]]
  grid <- terra::rast(xmin = DOMAIN[["xmin"]], xmax = DOMAIN[["xmax"]],
                      ymin = DOMAIN[["ymin"]], ymax = DOMAIN[["ymax"]],
                      resolution = GRID_CELL_DEG, crs = "EPSG:4326")
  terra::values(grid) <- 1
  cat(sprintf("     probe tag %s (%d obs, %.0f days), grid %d cells\n",
              mid, nrow(gm$light), gm$days, terra::ncell(grid)))
  t <- system.time(invisible(capture.output(
    f <- TwilightFreeGrid(gm$light$time, tag_light(gm, rm_), grid = grid,
      start_lon = gm$p0[1], start_lat = gm$p0[2],
      end_lon = gm$p1[1], end_lat = gm$p1[2],
      step_hours = STEP_HOURS, diffusion = DIFFUSION,
      calibration = rm_$calibration, likelihood_params = likelihood_params(rm_),
      area_correction = AREA_CORRECTION))))[["elapsed"]]
  n_knots <- length(f$fit$lat)
  cat(sprintf("     %.0f s for %d knots = %.3f s/knot\n", t, n_knots, t / n_knots))
  grid_arm <- t * length(tags) / 3600
  cat(sprintf("\n     ESTIMATED WALL CLOCK, quiescent machine\n"))
  cat(sprintf("       grid arm (29 fits)          %5.1f h  x 3 light arms = %.1f h\n",
              grid_arm, 3 * grid_arm))
  cat(sprintf("       grid fusion arm (~18%% less) %5.1f h\n", 0.82 * grid_arm))
  cat(sprintf("       hierarchical, 2 batches      ~0.5 h\n"))
  cat(sprintf("       FFBS, 3 tags x 2 arms        ~0.3 h\n"))
  cat(sprintf("       ------------------------------------\n"))
  cat(sprintf("       TOTAL                       %5.1f h\n", 3 * grid_arm + 0.82 * grid_arm + 0.8))
}

cat("\n=========================================================\n")
cat(sprintf("PREFLIGHT: %d failure(s), %d warning(s)\n", fails, warns))
if (fails) {
  cat("DO NOT LAUNCH. Fix the failures above.\n")
} else {
  cat("Clear to launch:  Rscript analysis/study/run_study.R\n")
}
cat("=========================================================\n")
quit(status = if (fails) 1L else 0L)
