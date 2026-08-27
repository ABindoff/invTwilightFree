# =============================================================================
# CONTROL: is the pooling gain in occupancy overlap real, or is it sample size?
#
# Pooling 29 tracks raises overlap against Argos from a per-tag median of 0.346
# to 0.615. Two things changed at once and they are not the same claim:
#
#   (a) 29 animals' errors partly cancel        <- a real benefit of a panel
#   (b) 10,909 knots fill a 2-D map that 376    <- an artefact of sample size,
#       knots cannot                                and it would happen even if
#                                                   every error were identical
#
# The control holds sample size fixed and varies only how many animals the knots
# come from. At matched n:
#
#   ONE      all n knots from a single deployment
#   SPREAD   n knots drawn evenly across all 29 deployments
#
# Both compare a light map against an Argos map built on the SAME knots, so
# sparsity penalises both arms equally. If SPREAD beats ONE at the same n, the
# panel is buying error cancellation. If they match, the gain was sample size.
#
# A sparsity curve for SPREAD shows where overlap saturates, which is the number
# a reader needs to know how many tags this argument requires.
#
#   Rscript analysis/study/occupancy_subsample.R
# =============================================================================

suppressMessages({library(data.table); library(terra); library(invTwilightFree)})
source("analysis/study/study_config.R")

SEED <- 42
REPS <- 200
set.seed(SEED)

K <- fread(file.path(OUT_DIR, "study_knots.csv"), colClasses = list(character = "id"))
K <- K[arm == "grid_light" & scored == TRUE]
K[, time := as.POSIXct(time, tz = "UTC")]
ids <- unique(K$id)
cat(sprintf("%d scored knots, %d deployments\n", nrow(K), length(ids)))

CELL <- 2
g <- rast(xmin = DOMAIN[["xmin"]], xmax = DOMAIN[["xmax"]],
          ymin = DOMAIN[["ymin"]], ymax = DOMAIN[["ymax"]],
          resolution = CELL, crs = "EPSG:4326")
values(g) <- 1
cat(sprintf("grid %d x %d at %g degrees\n\n", ncol(g), nrow(g), CELL))

# Overlap of the light map against the Argos map, on one set of knots. Both maps
# are built the same way from the same rows, so any sparsity penalty applies to
# both. weight = "knot" throughout: a subsample is not a panel of whole tracks
# and re-weighting per tag would reintroduce the thing being controlled for.
ov <- function(D) {
  if (nrow(D) < 2) return(NA_real_)
  a <- occupancy_map(data.frame(time = D$time, lon = D$truth_lon, lat = D$truth_lat),
                     grid = g, method = "point", weight = "knot")
  b <- occupancy_map(data.frame(time = D$time, lon = D$lon, lat = D$lat),
                     grid = g, method = "point", weight = "knot")
  occupancy_overlap(a, b)$overlap
}

per_tag_n <- K[, .N, by = id]
n_typical <- as.integer(median(per_tag_n$N))
cat(sprintf("median knots per deployment: %d\n", n_typical))

# ---- ONE: a single deployment, subsampled to n ------------------------------
one_at <- function(n, reps = REPS) {
  elig <- per_tag_n[N >= n]$id
  if (!length(elig)) return(NA_real_)
  v <- vapply(seq_len(reps), function(i) {
    d <- K[id == sample(elig, 1)]
    ov(d[sample(.N, n)])
  }, 0)
  median(v, na.rm = TRUE)
}

# ---- SPREAD: n knots drawn evenly across all deployments --------------------
spread_at <- function(n, reps = REPS) {
  per <- max(1L, floor(n / length(ids)))
  v <- vapply(seq_len(reps), function(i) {
    d <- K[, .SD[sample(.N, min(.N, per))], by = id]
    if (nrow(d) > n) d <- d[sample(.N, n)]
    ov(d)
  }, 0)
  median(v, na.rm = TRUE)
}

cat("\n=== MATCHED SAMPLE SIZE: one deployment vs the whole panel ===\n")
ladder <- c(50, 100, 200, n_typical, 800, 2000, 5000)
ladder <- sort(unique(ladder[ladder <= nrow(K)]))
res <- rbindlist(lapply(ladder, function(n) {
  data.table(n_knots = n,
             one = round(one_at(n), 3),
             spread = round(spread_at(n), 3),
             n_eligible_tags = per_tag_n[N >= n, .N])
}))
res[, gain := round(spread - one, 3)]
print(as.data.frame(res), row.names = FALSE)

cat(sprintf("\nfull pool, all %d knots: %.3f\n", nrow(K), ov(K)))
cat(sprintf("observed per-tag median (each tag's own knots, no subsampling): %.3f\n",
            median(vapply(ids, function(i) ov(K[id == i]), 0))))

cat("\n=== READING ===\n")
at_typ <- res[n_knots == n_typical]
cat(sprintf("At the typical single-deployment sample size (%d knots):\n", n_typical))
cat(sprintf("  one deployment      %.3f\n", at_typ$one))
cat(sprintf("  spread over 29 tags %.3f\n", at_typ$spread))
cat(sprintf("  difference          %+.3f\n", at_typ$gain))
cat(sprintf("\nfull-pool overlap minus one-deployment-at-matched-n = %+.3f\n",
            ov(K) - at_typ$one))
cat(sprintf("of which spreading across animals accounts for %+.3f and the extra\n",
            at_typ$gain))
cat(sprintf("%d knots for %+.3f\n", nrow(K) - n_typical, ov(K) - at_typ$spread))

saveRDS(list(seed = SEED, reps = REPS, cell = CELL, ladder = res,
             full_pool = ov(K), n_typical = n_typical),
        file.path(OUT_DIR, "occupancy_subsample.rds"))
fwrite(res, file.path(OUT_DIR, "report", "occupancy_subsample.csv"))
cat(sprintf("\nwritten to %s\n", file.path(OUT_DIR, "occupancy_subsample.rds")))
