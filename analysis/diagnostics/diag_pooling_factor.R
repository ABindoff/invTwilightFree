# Pooling factors, per tag, in the two places this analysis pools.
#
# DEFINITION USED (Gelman & Pardoe 2006, Technometrics). For group-level errors
# eps_j = theta_j - mu, with V_j the finite-population variance across groups
# within a draw and E the expectation over draws,
#
#     lambda = 1 - V_j( E(eps_j) ) / E( V_j(eps_j) )
#
# lambda -> 1 is complete pooling (the group estimates are pulled to the
# population and the between-group spread of posterior MEANS is small relative to
# the posterior uncertainty); lambda -> 0 is no pooling. The "pooling fraction"
# in the sense of weight left on a group's own data is 1 - lambda.
#
# This may not coincide with the definition in Bindoff (2026); it is stated
# explicitly so the two can be mapped.
#
# TWO PLACES POOLING HAPPENS HERE, in very different states:
#
#   1. MOVEMENT SCALE. Genuinely hierarchical: sig2_i has an inverse-gamma prior
#      whose parameters are estimated across the panel. lambda is computable from
#      the posterior draws, which the fit retains.
#
#   2. LIGHT-RESPONSE GEOMETRY. Currently pooled COMPLETELY by construction
#      (the panel median replaces each tag's own estimate), i.e. an imposed
#      lambda = 1. The interesting question is what lambda SHOULD be, which the
#      normal-normal algebra answers from two quantities we have measured:
#      the between-tag spread of the true response, and the error with which any
#      one tag's response can be measured.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))

# ---- 1. movement scale, from the posterior draws ----------------------------
f <- NULL
for (nm in c("analysis/cache/nes/fit_hier_light_v2.rds",
             "analysis/cache/nes/fit_hier_light.rds")) {
  if (file.exists(nm)) { f <- readRDS(nm); cat("using", nm, "\n"); break }
}
if (!is.null(f)) {
  fit <- if (!is.null(f$value)) f$value else f
  S <- fit$draws$sig2_kmday          # draws x tags, movement scale in km/day
  if (!is.null(S)) {
    S <- as.matrix(S)
    ids <- if (!is.null(fit$movement$id)) fit$movement$id else colnames(S)
    mu  <- rowMeans(S)                       # population level, per draw
    eps <- S - mu                            # group errors, per draw
    num <- var(colMeans(eps))                # V_j of E(eps_j)
    den <- mean(apply(eps, 1, function(r) var(r)))   # E of V_j(eps_j)
    lambda <- 1 - num / den
    cat(sprintf("\n=== movement scale: panel pooling factor lambda = %.3f ===\n", lambda))
    cat(sprintf("    (pooling fraction, weight on own data = %.3f)\n", 1 - lambda))

    # a per-tag reading: how far each tag's posterior mean sits from the
    # population, relative to its own posterior sd
    per <- data.frame(id = ids,
      post_mean = round(colMeans(S), 1),
      post_sd   = round(apply(S, 2, sd), 1),
      z_from_pop = round((colMeans(S) - mean(mu)) / apply(S, 2, sd), 2),
      row.names = NULL)
    per$shrunk_by <- round(1 - abs(per$post_mean - mean(mu)) /
                             pmax(abs(per$post_mean - mean(mu)), 1e-9), 3)
    cat("\nper tag (z_from_pop large => that tag's data overrode the population):\n")
    print(per[order(-abs(per$z_from_pop)), c("id","post_mean","post_sd","z_from_pop")],
          row.names = FALSE)
    cat(sprintf("\npopulation mean %.1f km/day; between-tag sd of posterior means %.1f;\n",
                mean(mu), sd(colMeans(S))))
    cat(sprintf("mean within-tag posterior sd %.1f\n", mean(apply(S, 2, sd))))
  } else cat("no sig2 draws retained in this fit\n")
} else cat("no hierarchical fit on disk yet (re-render in progress)\n")

# ---- 2. light-response geometry: what SHOULD the pooling factor be? ---------
# Between-tag spread of the TRUE at-sea response, and the error with which the
# haul-out measures it. Both were measured in diag_truegeom.R.
true_z50 <- c(90.66, 91.07, 91.51, 91.60, 91.60, 92.51, 91.70, 91.74, 92.33, 92.26)
mean_abs_diff <- 0.87        # mean |haul-out z50 - at-sea z50| over the 7 with both

tau  <- sd(true_z50)                     # between-tag sd of the true response
# for d ~ N(0, s^2), E|d| = s * sqrt(2/pi)
s_meas <- mean_abs_diff / sqrt(2/pi)     # measurement sd of a single tag's fit

lambda_geom <- s_meas^2 / (s_meas^2 + tau^2)
cat(sprintf("\n=== light-response geometry ===\n"))
cat(sprintf("between-tag sd of the TRUE response (tau)      : %.3f deg\n", tau))
cat(sprintf("measurement sd of one tag's haul-out fit (s)   : %.3f deg\n", s_meas))
cat(sprintf("implied shrinkage lambda = s^2/(s^2+tau^2)     : %.3f\n", lambda_geom))
cat(sprintf("implied POOLING FRACTION (weight on own data)  : %.3f\n", 1 - lambda_geom))
cat("\nThe analysis currently imposes lambda = 1 (complete pooling). The algebra\n")
cat("says the optimal value is close to that but not equal to it, so a small\n")
cat("gain is available from partial pooling -- and it confirms why per-tag\n")
cat(sprintf("fitting (lambda = 0) was worse: it puts %.0f%% of the weight on an\n",
            100 * 1))
cat(sprintf("estimate that deserves %.0f%%.\n", 100 * (1 - lambda_geom)))

# what the partial-pooling adjustment would actually be, per tag
haul <- c("2021023"=91.84, "2021025"=92.00, "2021027"=92.09, "2021028"=92.99,
          "2021032"=93.02, "2021033"=93.06, "2021035"=91.92)
pooled <- 92.09
adj <- round((1 - lambda_geom) * (haul - pooled), 3)
cat("\nper-tag z50 adjustment implied by partial pooling (deg):\n")
print(adj)
cat(sprintf("largest adjustment %.3f deg of zenith\n", max(abs(adj))))
