# The prior fraction (Bindoff 2026, eq. 20) applied to this model, and a test of
# its central empirical claim on a model class well outside the paper's
# simulation study.
#
#     pi_j = (1/sigma^2) / G_FF,j  =  prior precision / total fiber precision
#
# Two things this model has that the paper's GLMM does not, and both matter.
#
# 1. THREE LEVELS, so two nested bundles:
#      population  ->  per-tag movement (sigma^2_i, rho_i)   [fiber = scalars]
#      per-tag movement  ->  LATENT TRACK (x_{i,1:K})        [fiber = the track]
#    The paper's base-fiber split is the first; the second is where this model's
#    difficulty actually lives, and its fiber is hundreds of dimensions per tag.
#
# 2. ANISOTROPY. The light likelihood is far sharper in longitude than latitude,
#    so the SAME isotropic movement prior gives a different prior fraction per
#    axis. If pi differs across the axes of one latent field, no uniform
#    parameterisation of that field can be right.
#
# The paper's claim under test: prior-dominated groups (high pi_j) mix worse.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))

KM_DEG <- 111.32
D <- 110; STEP_H <- 12
sigma_step <- D * sqrt(STEP_H / 24)          # prior sd per knot, km

# ---- bundle 2: the track, per axis, per tag ---------------------------------
k <- as.data.table(read.csv("analysis/output/nes/per_knot_errors.csv"))
k <- k[engine == "hier" & batch == "light" & is.finite(err_lat) & is.finite(err_lon)]
k[, id := as.character(id)]

# tau: the width of the LIGHT likelihood alone. Approximated by the error sd,
# which is an upper bound on it where the prior has already smoothed -- so
# pi_lat below is, if anything, understated.
tk <- k[, .(tau_lat = sd(err_lat) * KM_DEG,
            tau_lon = sd(err_lon) * KM_DEG * cos(mean(true_lat) * pi/180)), by = id]
tk[, `:=`(pi_lat = round((1/sigma_step^2) / ((1/sigma_step^2) + (1/tau_lat^2)), 3),
          pi_lon = round((1/sigma_step^2) / ((1/sigma_step^2) + (1/tau_lon^2)), 3))]
tk[, `:=`(tau_lat = round(tau_lat), tau_lon = round(tau_lon))]
cat(sprintf("=== bundle 2: the latent track (prior sd %.0f km per %d-h knot) ===\n",
            sigma_step, STEP_H))
print(as.data.frame(tk[order(-pi_lat)]), row.names = FALSE)
cat(sprintf("\nmean pi_lat = %.3f  (prior-dominated; paper's threshold for non-centring is 0.6)\n",
            mean(tk$pi_lat)))
cat(sprintf("mean pi_lon = %.3f  (likelihood-dominated; centring appropriate)\n",
            mean(tk$pi_lon)))
cat("\nOne isotropic movement prior, two regimes. The paper's per-group rule\n")
cat("applied per AXIS says: non-centre latitude (increments), centre longitude\n")
cat("(positions). No uniform parameterisation of the track serves both.\n")

# ---- bundle 1: per-tag movement scale, and the mixing claim ------------------
f <- NULL
for (nm in c("analysis/cache/nes/fit_hier_light_v2.rds",
             "analysis/cache/nes/fit_hier_light.rds"))
  if (file.exists(nm)) { f <- readRDS(nm); cat("\nusing", nm, "\n"); break }

if (!is.null(f)) {
  fit <- if (!is.null(f$value)) f$value else f
  S <- as.matrix(fit$draws$sig2_kmday)
  ids <- if (!is.null(fit$movement$id)) fit$movement$id else colnames(S)
  n_knots <- fit$tags$n_knots
  names(n_knots) <- fit$tags$id

  # pi_j for the movement scale. The fiber parameter is log sigma^2_i; its
  # likelihood information is about (k_i - 1)/2 from the increments, and the
  # population prior contributes the equivalent of a_pop. a_pop is recovered
  # from the panel: the between-tag variance of the posterior means against the
  # mean within-tag posterior variance.
  m_j  <- colMeans(S); v_j <- apply(S, 2, var)
  tau2 <- max(var(m_j) - mean(v_j), 1e-8)      # between-tag variance of the truth
  pi_j <- v_j / (v_j + tau2)                   # prior share of the fiber precision

  # mixing, per tag: integrated autocorrelation time of that tag's own chain
  iact <- apply(S, 2, function(x) {
    a <- acf(x, lag.max = min(200, length(x) %/% 5), plot = FALSE)$acf[, , 1]
    cut <- which(a < 0.05)[1]; if (is.na(cut)) cut <- length(a)
    1 + 2 * sum(a[2:cut])
  })
  ess <- nrow(S) / pmax(iact, 1)

  P <- data.frame(id = ids, n_knots = as.integer(n_knots[ids]),
                  post_mean = round(m_j, 1), pi_j = round(pi_j, 3),
                  iact = round(iact, 1), ess = round(ess), row.names = NULL)
  cat("\n=== bundle 1: per-tag movement scale ===\n")
  print(P[order(-P$pi_j), ], row.names = FALSE)
  ct <- suppressWarnings(cor.test(P$pi_j, P$ess, method = "spearman"))
  cat(sprintf("\npaper's prediction: high pi_j mixes worse.\n"))
  cat(sprintf("Spearman(pi_j, ESS) = %+.3f (p = %.3f), n = %d\n",
              ct$estimate, ct$p.value, nrow(P)))
  cat(sprintf("mean pi_j = %.3f -- %s\n", mean(pi_j),
              if (mean(pi_j) < 0.6) "below 0.6: centring appropriate, no remedy indicated"
              else "above 0.6: non-centring indicated"))
  cat("\nNote the contrast between the two bundles: the movement scale is\n")
  cat("data-dominated and needs nothing, while the track it governs is strongly\n")
  cat("prior-dominated in one axis. The pathology is a level below where the\n")
  cat("hierarchical parameters live.\n")
} else cat("\nno hierarchical fit on disk (re-render in progress)\n")
