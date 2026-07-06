# M3 gate: a pure-R per-obs light log-likelihood evaluator that batches a whole
# block into ONE solar_zenith call, validated bit-close against eval_logpk_grid.
# If this does not match to ~machine precision, batching would break exactness.
suppressPackageStartupMessages(library(invTwilightFree))
set.seed(1)

calibration <- c(525.97, 5.458)
likpar      <- c(0.0312, 64.0, 0.10)
intercept <- calibration[1]; slope <- calibration[2]
lambda <- likpar[1]; max_light <- likpar[2]; prob_slab <- likpar[3]

# Rust replica: spike-and-slab log-density summed over a set of obs, for one point.
emit_point <- function(lon, lat, otimes, olight) {
  z <- solar_zenith(otimes, rep(lon, length(otimes)), rep(lat, length(otimes)))
  expc <- pmin(pmax(intercept - slope * z, 0), max_light)
  raw <- ifelse(olight <= expc,
                lambda * exp(-lambda * (expc - olight)),
                lambda * exp(-lambda * 2 * (olight - expc)))
  norm <- pmax(1 - exp(-lambda * expc) + 0.5 * (1 - exp(-2 * lambda * (max_light - expc))), 1e-12)
  spike <- raw / norm
  den <- (1 - prob_slab) * spike + prob_slab * (1 / max_light)
  sum(log(den))
}

# random obs series and random points
t0 <- as.numeric(as.POSIXct("2024-06-15", tz = "UTC"))
otimes <- t0 + sort(runif(60, 0, 12 * 3600))
olight <- pmin(pmax(rnorm(60, 30, 15), 0), 64)
pts <- cbind(runif(12, -40, 0), runif(12, -60, 60))  # lon, lat

# reference: eval_logpk_grid over all points in one call
ref <- eval_logpk_grid(pts[, 1], pts[, 2], otimes, olight, calibration, likpar)
mine <- vapply(seq_len(nrow(pts)), function(i) emit_point(pts[i, 1], pts[i, 2], otimes, olight), numeric(1))

cat(sprintf("max abs diff (R emitter vs eval_logpk_grid): %.3e\n", max(abs(ref - mine))))
cat(sprintf("max rel diff: %.3e\n", max(abs((ref - mine) / (abs(ref) + 1e-9)))))
cat("sample ref: ", paste(round(head(ref, 4), 6), collapse = ", "), "\n")
cat("sample mine:", paste(round(head(mine, 4), 6), collapse = ", "), "\n")
