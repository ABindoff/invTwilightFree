# =============================================================================
# EXACT POSTERIOR OCCUPANCY, against Argos and against the approximation.
#
# refit_posteriors.R reproduced grid_light bit-for-bit and kept the per-knot cell
# posteriors summed within calendar month. This builds the occupancy map from
# those surfaces and asks three things:
#
#   1. how close is the Gaussian-kernel approximation, built from the reported
#      marginal standard deviations, to the exact surface? If it is close, the
#      approximation is fine for engines that report only a mean and an sd.
#   2. how well does either recover the Argos occupancy distribution, at cell
#      sizes from 1 to 5 degrees?
#   3. where does the residual sit? A displaced mode and a noisy one are
#      different failures and only one of them is fixable by more tags.
#
# ONE MISMATCH, STATED. The saved posteriors cover EVERY knot; the Argos
# reference exists only at SCORED knots (94% overall, but under half on five
# deployments, and the gaps concentrate offshore where a seal surfaces least).
# The exact map therefore carries mass in places the reference cannot. Section 2
# measures how much that is worth by running the approximation both ways.
#
#   Rscript analysis/study/occupancy_exact.R
# =============================================================================

suppressMessages({library(data.table); library(terra); library(invTwilightFree)})
`%||%` <- function(a, b) if (is.null(a)) b else a
source("analysis/study/study_config.R")

POST <- file.path(OUT_DIR, "posteriors")
REP  <- file.path(OUT_DIR, "report")
dir.create(REP, recursive = TRUE, showWarnings = FALSE)
files <- list.files(POST, pattern = "[.]rds$", full.names = TRUE)
if (!length(files)) stop("no posteriors in ", POST, " -- run refit_posteriors.R first")
cat(sprintf("%d deployments with saved posteriors\n", length(files)))

base <- rast(xmin = DOMAIN[["xmin"]], xmax = DOMAIN[["xmax"]],
             ymin = DOMAIN[["ymin"]], ymax = DOMAIN[["ymax"]],
             resolution = GRID_CELL_DEG, crs = "EPSG:4326")
values(base) <- 1

season_of <- function(month_label) {
  m <- as.integer(substr(month_label, 6, 7))
  c("winter","winter","spring","spring","spring","summer",
    "summer","summer","autumn","autumn","autumn","winter")[m]
}

# ---- 1. the exact map, one tag at a time, equal weight per tag ---------------
# Same accumulation occupancy_map() performs for a grid fit: normalise each
# tag's contribution before folding it in, so a 313-day deployment cannot
# outvote a 137-day one.
exact_layers <- list()
for (f in files) {
  o <- readRDS(f)
  sea <- season_of(o$months)
  for (s in unique(sea)) {
    v <- colSums(o$P[sea == s, , drop = FALSE])
    tot <- sum(v)
    if (!is.finite(tot) || tot <= 0) next
    exact_layers[[s]] <- (exact_layers[[s]] %||% numeric(ncell(base))) + v / tot
  }
  exact_layers[["all"]] <- (exact_layers[["all"]] %||% numeric(ncell(base))) +
    colSums(o$P) / sum(o$P)
}
to_rast <- function(v) { r <- base; values(r) <- v / sum(v); r }
exact_all <- to_rast(exact_layers[["all"]])

# ---- 2. the reference and the approximation, on matched knots ---------------
K <- fread(file.path(OUT_DIR, "study_knots.csv"), colClasses = list(character = "id"))
K <- K[arm == "grid_light"]
K[, time := as.POSIXct(time, tz = "UTC")]
KS <- K[scored == TRUE]

tracks <- function(D, lon, lat, sd = FALSE) {
  setNames(lapply(unique(D$id), function(i) {
    d <- D[id == i]
    out <- data.frame(time = d$time, lon = d[[lon]], lat = d[[lat]])
    if (sd) { out$lon_sd <- d$sd_lon; out$lat_sd <- d$sd_lat }
    out
  }), unique(D$id))
}

cat("\n=== 1. HOW MUCH DOES RESTRICTING TO SCORED KNOTS MATTER? ===\n")
cat("The approximation run over ALL knots and over SCORED knots only, both\n")
cat("against the Argos map. The difference bounds the caveat on the exact map.\n\n")
argos <- occupancy_map(tracks(KS, "truth_lon", "truth_lat"), grid = base,
                       method = "point", weight = "tag")
gauss_all    <- occupancy_map(tracks(K,  "lon", "lat", sd = TRUE), grid = base,
                              method = "posterior", weight = "tag")
gauss_scored <- occupancy_map(tracks(KS, "lon", "lat", sd = TRUE), grid = base,
                              method = "posterior", weight = "tag")
cat(sprintf("  approximation, ALL knots    vs Argos: overlap %.3f\n",
            occupancy_overlap(argos, gauss_all)$overlap))
cat(sprintf("  approximation, SCORED knots vs Argos: overlap %.3f\n",
            occupancy_overlap(argos, gauss_scored)$overlap))

cat("\n=== 2. IS THE GAUSSIAN APPROXIMATION FAITHFUL TO THE EXACT SURFACE? ===\n")
cat("Both over ALL knots, so this isolates the approximation alone.\n\n")
o <- occupancy_overlap(exact_all, gauss_all)
cat(sprintf("  exact vs approximation: overlap %.3f  bhattacharyya %.3f  cor %.3f\n",
            o$overlap, o$bhattacharyya, o$correlation))

cat("\n=== 3. RECOVERY OF THE ARGOS OCCUPANCY DISTRIBUTION ===\n")
point <- occupancy_map(tracks(KS, "lon", "lat"), grid = base,
                       method = "point", weight = "tag")
res <- rbindlist(lapply(c(1, 2, 5), function(cs) {
  agg <- function(r) if (cs == GRID_CELL_DEG) r else
    aggregate(r, fact = cs / GRID_CELL_DEG, fun = "sum", na.rm = TRUE)
  a <- agg(argos)
  rbindlist(list(
    data.table(cell_deg = cs, form = "point",       occupancy_overlap(a, agg(point))[, -1]),
    data.table(cell_deg = cs, form = "approximation", occupancy_overlap(a, agg(gauss_scored))[, -1]),
    data.table(cell_deg = cs, form = "exact",       occupancy_overlap(a, agg(exact_all))[, -1])))
}))
print(as.data.frame(res[, .(cell_deg, form, overlap = round(overlap, 3),
                            bhattacharyya = round(bhattacharyya, 3),
                            correlation = round(correlation, 3))]), row.names = FALSE)
fwrite(res, file.path(REP, "occupancy_exact.csv"))

cat("\n=== 4. BY SEASON (2 degree cells) ===\n")
seasons <- setdiff(names(exact_layers), "all")
a_s <- occupancy_map(tracks(KS, "truth_lon", "truth_lat"), grid = base,
                     method = "point", weight = "tag",
                     by = function(t) season_of(format(t, "%Y-%m")))
# Both forms, season by season, and BOTH ON THE SCORED KNOTS the Argos reference
# uses. This matters more per season than it does overall. The exact map
# integrates every knot, and the scored fraction is 0.83, 0.83 and 0.77 in
# summer, autumn and winter but only 0.60 in spring, on 250 knots. Comparing an
# all-knot model map with a scored-knot reference made the posterior look far
# better than the point form in spring (0.382 against 0.248); on matched knots
# that reversal disappears (0.239 against 0.248). Quote the matched columns.
p_s <- occupancy_map(tracks(KS, "lon", "lat"), grid = base,
                     method = "point", weight = "tag",
                     by = function(t) season_of(format(t, "%Y-%m")))
gs <- occupancy_map(tracks(KS, "lon", "lat", sd = TRUE), grid = base,
                    method = "posterior", weight = "tag",
                    by = function(t) season_of(format(t, "%Y-%m")))
nk <- K[, .(all = .N, scored = sum(scored)),
        by = .(sname = season_of(format(time, "%Y-%m")))]
sres <- rbindlist(lapply(seasons, function(s) {
  if (!s %in% names(a_s)) return(NULL)
  agg <- function(r) aggregate(r, fact = 2 / GRID_CELL_DEG, fun = "sum", na.rm = TRUE)
  a <- agg(a_s[[s]])
  ex <- occupancy_overlap(a, agg(to_rast(exact_layers[[s]])))
  pt <- if (s %in% names(p_s)) occupancy_overlap(a, agg(p_s[[s]])) else NULL
  gg <- if (s %in% names(gs)) occupancy_overlap(a, agg(gs[[s]])) else NULL
  n  <- nk[sname == s]
  data.table(season = s,
             n_knots  = if (nrow(n)) n$all[1] else NA_integer_,
             n_scored = if (nrow(n)) n$scored[1] else NA_integer_,
             overlap_point = if (is.null(pt)) NA_real_ else pt$overlap,
             overlap_gauss = if (is.null(gg)) NA_real_ else gg$overlap,
             overlap_exact_unmatched = ex$overlap,
             cor_exact = ex$correlation)
}))
print(as.data.frame(sres[, .(season, n_knots, n_scored,
                             point = round(overlap_point, 3),
                             gaussian = round(overlap_gauss, 3),
                             exact_unmatched = round(overlap_exact_unmatched, 3))]),
      row.names = FALSE)
cat("point and gaussian share the Argos knot set; exact_unmatched does not.
")
cat("Spring rests on 250 knots and supports no inference.
")
fwrite(sres, file.path(REP, "occupancy_exact_season.csv"))

cat("\n=== 5. WHERE IS THE RESIDUAL? latitude profile, 5 degree bands ===\n")
prof <- function(r, nm) {
  v <- values(r)[, 1]; xy <- crds(r, na.rm = FALSE)
  b <- tapply(v, cut(xy[, 2], seq(20, 70, by = 5)), sum, na.rm = TRUE)
  setNames(data.table(band = names(b), round(as.numeric(b), 3)), c("band", nm))
}
P <- Reduce(function(a, b) merge(a, b, by = "band", sort = FALSE),
            list(prof(argos, "argos"), prof(point, "point"), prof(exact_all, "exact")))
P[, `:=`(d_point = round(point - argos, 3), d_exact = round(exact - argos, 3))]
print(as.data.frame(P), row.names = FALSE)
fwrite(P, file.path(REP, "occupancy_latitude_profile.csv"))

saveRDS(list(exact = exact_layers, grid = list(ext = as.vector(ext(base)),
             res = res(base))), file.path(OUT_DIR, "occupancy_exact.rds"))
cat(sprintf("\nwritten to %s\n", REP))
