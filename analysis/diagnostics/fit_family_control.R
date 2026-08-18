# REPRODUCTION CONTROL for response-family comparisons.
#
# WHY THIS EXISTS. A first attempt to compare response families on real light
# returned 706 km median against the pipeline's 227 on the same six tags -- 2 to 7x
# off on five of them. Both arms were badly calibrated, so the contrast between them
# was meaningless. The failure was not subtle and it was entirely mine: I rebuilt the
# harness from notes instead of from the pipeline, and got six things wrong at once
#   * per-tag battery geometry instead of POOLED-WITHIN-FAMILY haul-out geometry
#   * q05 of the whole record as baseline instead of the fitted response's own
#   * a per-tag grid around the truth instead of the fixed 150-250 E, 20-70 N domain
#   * Argos-at-first-observation as endpoints instead of the deployment/recovery
#     positions from metadata
#   * `approx` on pre-interpolated positions instead of `argos_at()`
#   * light resampled onto foreign timestamps instead of the tag's own series
#
# So this script does NOT rebuild anything. The panel assembly and calibration below
# are copied VERBATIM from `fit_all29.R`, which is the validated pipeline and which
# produced `all29_tags.csv`. The only thing that will vary, once arms are added, is
# the response handed to the engine.
#
# THE GUARD AGAINST COPY DRIFT is the control itself: arm `tangent` must reproduce
# `all29_tags.csv` per tag. If it does not, the copy is unfaithful and nothing built
# on it is readable. That check runs at the end and is the whole point of this file.
#
# ARMS ARE DELIBERATELY NOT ADDED YET. Build the control, verify it, then extend.
# That ordering is the lesson of the day it cost to learn.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))
OUT <- "scratch/nes_calibration/family_arms.csv"
# THE CORRECT REFERENCE. `all29_tags.csv` is the PRE-area-correction baseline --
# SESSION_STATE is explicit that it is "the archived tangent on the BIASED kernel,
# 242 km", and that "242 was two bugs agreeing". Comparing against it made a faithful
# harness look broken (ratios up to 1.58). `rescore29_results.csv` arm
# `tangent_areaON` is the corrected-kernel run and is what a current harness must
# reproduce.
REF <- "scratch/nes_calibration/rescore29_results.csv"
REF_ARM <- "tangent_areaON"

STEP_HOURS <- 12; DIFFUSION <- 110; CELL <- 1; CAL_DAYS <- 15
# fit only these, for turnaround; calibration still pools over the WHOLE panel,
# exactly as the pipeline does
SUBSET <- c("2021033", "2021027", "2023037", "2023032", "2021032", "2022041")

# ================= copied verbatim from fit_all29.R: panel ====================
old_arch <- lapply(readRDS("analysis/cache/nes/archives_v3_30min.rds"), `[[`, "main")
new_arch <- lapply(readRDS("analysis/cache/nes/fvilches_archives_v1_30min.rds"), `[[`, "main")
new_argos <- readRDS("analysis/cache/nes/fvilches_argos_v1.rds")
new_man <- fread("scratch/nes_calibration/fvilches_manifest.csv",
                 colClasses = list(character = "topp"))
meta <- nes_meta(); ar <- nes_argos()

ser <- rbind(
  data.table(topp = names(old_arch),
             serial = sub("^[0-9]+_(.+)\\.csv$", "\\1",
                          basename(list.files("data/nes_untracked", pattern = "^[0-9]+_.*\\.csv$")))[
                            match(names(old_arch),
                                  sub("^([0-9]+)_.*$", "\\1",
                                      basename(list.files("data/nes_untracked",
                                                          pattern = "^[0-9]+_.*\\.csv$"))))]),
  data.table(topp = sub("^([0-9]+)_.*$", "\\1",
                        list.files("fvilches/extracted/TDR raw", pattern = "Archive[.]csv$")),
             serial = sub("^[0-9]+_(.+)-out-Archive[.]csv$", "\\1",
                          list.files("fvilches/extracted/TDR raw", pattern = "Archive[.]csv$"))),
  fill = TRUE)
ser <- unique(ser, by = "topp")
family <- function(s) fifelse(is.na(s), "unknown",
                       fifelse(grepl("^219", s), "Mk9_219",
                        fifelse(grepl("^18A", s), "F18A", "other")))

tags <- list()
for (tg in names(old_arch)) {
  m <- meta[match(tg, meta$id), ]
  if (is.na(m$ptt) || !nzchar(m$ptt)) next
  df <- ar$deploy[Ptt == m$ptt]; a <- ar$fixes[Ptt == m$ptt]
  if (!nrow(df) || nrow(a) < 50) next
  x <- as.data.table(old_arch[[tg]])
  t0 <- max(min(x$time), as.POSIXct(max(m$deploy_date, df$user_date), tz = "UTC"))
  t1 <- min(max(x$time), max(a$time))
  d <- x[time >= t0 & time <= t1]
  if (nrow(d) < 1000) next
  tf <- a[time >= max(time) - 72*3600]
  tags[[tg]] <- list(id = tg, season = "2021", light = d, argos = a,
                     p0 = c(df$deploy_lon, df$deploy_lat),
                     p1 = c(median(tf$lon), median(tf$lat)))
}
for (tg in names(new_arch)) {
  mm <- new_man[topp == tg]
  if (!nrow(mm)) next
  a <- as.data.table(new_argos[[tg]])
  d <- as.data.table(new_arch[[tg]])
  d <- d[time >= min(a$time) & time <= max(a$time)]
  if (nrow(d) < 1000 || nrow(a) < 50) next
  tags[[tg]] <- list(id = tg, season = as.character(mm$season[1]), light = d, argos = a,
                     p0 = c(mm$deploy_lon[1], mm$deploy_lat[1]),
                     p1 = c(mm$recover_lon[1], mm$recover_lat[1]))
}
for (tg in names(tags)) tags[[tg]]$family <- family(ser[topp == tg]$serial[1])
cat(sprintf("panel: %d deployments\n", length(tags)))

# ================= copied verbatim from fit_all29.R: calibration ==============
scale_fits <- lapply(tags, function(g) {
  k <- g$light$time <= min(g$light$time) + CAL_DAYS * 86400
  fit_light_response(g$light$time[k], g$light$light[k], g$p0[1], g$p0[2])
})
geom_fits <- lapply(tags, function(g) {
  dep <- departure_time(g$light)
  k <- is.finite(dep) & as.numeric(g$light$time) < dep
  if (sum(k) < 200) return(NULL)
  fit_light_response(g$light$time[k], g$light$light[k], g$p0[1], g$p0[2])
})
cat(sprintf("haul-out geometry available for %d of %d deployments\n",
            sum(!vapply(geom_fits, is.null, logical(1))), length(tags)))

responses <- list()
for (fam in unique(vapply(tags, function(g) g$family, ""))) {
  ids <- names(tags)[vapply(tags, function(g) g$family, "") == fam]
  gf <- geom_fits[ids]; sf <- scale_fits[ids]
  if (!sum(!vapply(gf, is.null, logical(1)))) next
  pooled <- pool_light_responses(gf, sf)
  p <- attr(pooled, "pooled")
  cat(sprintf("family %-8s: %d deployments, pooled from %d -> z50 %.2f, width %.1f\n",
              fam, length(ids), p[["n"]], p[["z50"]], 4.394 * p[["scale"]]))
  for (i in ids) responses[[i]] <- pooled[[i]]
}

grid <- rast(xmin = 150, xmax = 250, ymin = 20, ymax = 70,
             resolution = CELL, crs = "EPSG:4326")
values(grid) <- 1

# ================= response arms ==============================================
# GOMPERTZ, pooled the same way the shipped recipe pools the logistic: geometry
# (z50, scale) fitted per tag to its HAUL-OUT envelope then pooled within family,
# intensity (floor, amp) kept per tag from its own 15-day scale fit. Anything else
# would confound the family comparison with a change of pooling.
#
# Handed to the engine as a LOOKUP TABLE (leading -1 sentinel), which
# `expected_light` already supports. No Rust change: a 4-element Gompertz would
# collide with the logistic branch, and the table represents it exactly to
# interpolation error on a 1-degree grid.
gomp <- function(z, th) th[1] + th[2] * exp(-exp((z - th[3]) / th[4]))
fit_gomp <- function(env, base, ml) {
  z <- env$zenith; y <- (env$light - base) / ml
  k <- is.finite(z) & is.finite(y)
  if (sum(k) < 20) return(NULL)
  z <- z[k]; y <- pmin(pmax(y[k], 0), 1)
  obj <- function(th) sum((gomp(z, th) - y)^2)
  best <- NULL
  for (j in 1:6) {
    st <- c(0.1, 0.8, 95, 6) * (if (j == 1) 1 else runif(4, 0.7, 1.4))
    o <- try(optim(st, obj, control = list(maxit = 8000, reltol = 1e-12)), silent = TRUE)
    if (!inherits(o, "try-error") && (is.null(best) || o$value < best$value)) best <- o
  }
  best$par
}

set.seed(1)
gpar <- list()
for (id in names(geom_fits)) {
  gfit <- geom_fits[[id]]
  if (is.null(gfit)) next
  p <- try(fit_gomp(gfit$envelope, gfit$baseline, gfit$max_light), silent = TRUE)
  if (!inherits(p, "try-error") && !is.null(p)) gpar[[id]] <- p
}
gpool <- list()
for (fam in unique(vapply(tags, function(g) g$family, ""))) {
  ids <- intersect(names(tags)[vapply(tags, function(g) g$family, "") == fam], names(gpar))
  if (!length(ids)) next
  M <- do.call(rbind, gpar[ids])
  gpool[[fam]] <- c(z50 = median(M[, 3]), scale = median(M[, 4]))
  cat(sprintf("gompertz %-8s: pooled from %d -> z50 %.2f, scale %.2f\n",
              fam, length(ids), gpool[[fam]]["z50"], gpool[[fam]]["scale"]))
}

ZGRID <- seq(30, 140, by = 1)
make_table <- function(r, fam) {
  gp <- gpool[[fam]]
  if (is.null(gp)) return(NULL)
  # per-tag intensity from the shipped pooled response, pooled geometry from above
  th <- c(r$calibration_logistic[1] / r$max_light,
          r$calibration_logistic[2] / r$max_light, gp["z50"], gp["scale"])
  y <- pmin(pmax(gomp(ZGRID, th), 0), 1) * r$max_light
  c(-1, ZGRID[1], ZGRID[2] - ZGRID[1], y)
}

ARMS <- c("tangent", "gompertz", "logistic")

# ================= the fit, tangent arm only ==================================
done <- if (file.exists(OUT)) {
  d0 <- fread(OUT, colClasses = list(character = "id")); paste(d0$id, d0$arm)
} else character(0)
todo <- intersect(SUBSET, names(tags))
cat(sprintf("\nfitting %d of %d (id, arm) pairs\n\n",
            length(ARMS) * length(todo) - length(done), length(ARMS) * length(todo)))

for (arm in ARMS) for (tg in todo) {
  if (paste(tg, arm) %in% done) next
  g <- tags[[tg]]; r <- responses[[tg]]
  if (is.null(r)) { message(tg, ": no response, skipped"); next }
  cal <- switch(arm,
                tangent  = r$calibration,
                logistic = r$calibration_logistic,
                gompertz = make_table(r, g$family))
  if (is.null(cal)) { message(tg, "/", arm, ": no calibration, skipped"); next }
  message("fit: ", tg, " / ", arm)
  invisible(capture.output(
    f <- TwilightFreeGrid(g$light$time, pmax(0, g$light$light - r$baseline), grid = grid,
      start_lon = unname(g$p0[1]), start_lat = unname(g$p0[2]),
      end_lon = unname(g$p1[1]), end_lat = unname(g$p1[2]),
      step_hours = STEP_HOURS, diffusion = DIFFUSION, calibration = cal,
      likelihood_params = c(1/(r$max_light*0.5), r$max_light, 0.10))))
  gp <- grid_posterior(f)
  sdp <- numeric(nrow(gp$P))
  for (i in seq_len(nrow(gp$P))) {
    w <- gp$P[i, ]; s <- sum(w)
    if (!is.finite(s) || s <= 0) { sdp[i] <- NA; next }
    w <- w / s; sdp[i] <- sqrt(sum(w * (gp$lat - sum(w * gp$lat))^2))
  }
  tm <- as.POSIXct(f$fit$time, origin = "1970-01-01", tz = "UTC")
  tr <- argos_at(g$argos, tm)
  ep <- f$fit$lat - tr$lat
  e  <- gc_km(lon360(f$fit$lon), f$fit$lat, tr$lon, tr$lat)
  # mask copied from fit_rescore29.R:216 -- an earlier version used a different one
  # and the scored-knot sets diverged, which alone moved median km by a few percent
  ok <- is.finite(ep) & is.finite(sdp) & sdp > 0
  fwrite(data.table(id = tg, arm = arm, n = sum(ok),
                    median_km = median(e[ok]),
                    bias_lat = mean(ep[ok]),
                    rmse_lat = sqrt(mean(ep[ok]^2)),
                    cover_lat = mean(abs(ep[ok]) <= 1.96 * sdp[ok]),
                    log_z = f$log_z),
         OUT, append = file.exists(OUT))
}

# ================= THE CONTROL ================================================
D <- fread(OUT, colClasses = list(character = "id"))
R <- fread(REF, colClasses = list(character = "id"))[arm == REF_ARM]
m <- merge(D[arm == "tangent", .(id, mine_km = median_km, mine_bias = bias_lat)],
           R[, .(id, ref_km = median_km, ref_bias = bias_lat)], by = "id")
m[, ratio := mine_km / ref_km]
cat("\n=== REPRODUCTION CONTROL: tangent arm vs rescore29 tangent_areaON ===\n")
print(as.data.frame(m[, .(id, ref_km = round(ref_km), mine_km = round(mine_km),
                          ratio = round(ratio, 4),
                          d_bias = round(mine_bias - ref_bias, 3))]), row.names = FALSE)
cat(sprintf("\n  max |ratio - 1| = %.4f | max |d bias| = %.4f deg\n",
            max(abs(m$ratio - 1)), max(abs(m$mine_bias - m$ref_bias))))
ok_ctl <- max(abs(m$ratio - 1)) < 0.02
cat(sprintf("  VERDICT: %s\n", if (ok_ctl)
      "REPRODUCED -- the arm comparison below is readable." else
      "*** NOT REPRODUCED -- treat the comparison below with suspicion."))

cat("\n=== RESPONSE FAMILY COMPARISON (real light, 6 tags) ===\n")
S <- D[, .(n = .N, median_km = round(median(median_km)),
           bias = round(mean(bias_lat), 3),
           abs_bias = round(mean(abs(bias_lat)), 3),
           rmse_lat = round(mean(rmse_lat), 3),
           cover = round(mean(cover_lat), 3),
           log_z = round(mean(log_z))), by = arm]
print(as.data.frame(S[match(ARMS, arm)]), row.names = FALSE)

for (a in setdiff(ARMS, "tangent")) {
  w <- merge(D[arm == "tangent", .(id, k0 = median_km, b0 = bias_lat, z0 = log_z)],
             D[arm == a, .(id, k1 = median_km, b1 = bias_lat, z1 = log_z)], by = "id")
  if (!nrow(w)) next
  cat(sprintf("\n  %s vs tangent, paired (n = %d):\n", a, nrow(w)))
  cat(sprintf("    d km      %+0.0f   better on %d/%d\n",
              mean(w$k1 - w$k0), sum(w$k1 < w$k0), nrow(w)))
  cat(sprintf("    d |bias|  %+0.3f  better on %d/%d\n",
              mean(abs(w$b1)) - mean(abs(w$b0)),
              sum(abs(w$b1) < abs(w$b0)), nrow(w)))
  cat(sprintf("    d log_z   %+0.0f  (evidence prefers %s)\n",
              mean(w$z1 - w$z0), if (mean(w$z1) > mean(w$z0)) toupper(a) else "tangent"))
}
cat("\n  Envelope-fit prediction was gompertz 6.8x better than the tangent in the\n")
cat("  twilight band, at no cost in gauge inflation (1.01 vs 1.39). This is the\n")
cat("  fit-level test of that prediction.\n")
