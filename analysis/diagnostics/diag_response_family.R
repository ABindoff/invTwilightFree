# WHICH RESPONSE FAMILY? Fit gained, against degeneracy added.
#
# THE MEASURED CURVE (relative light, pooled over tags at Argos-known positions):
# 0.90 at zenith 71 declining GENTLY to 0.81 at 86 -- a sloping daytime, not a
# plateau -- then a collapse to 0.12 by 104, then a FLOOR near 0.10 rather than zero.
# The engine's own comment says as much: "a clamped line has two parameters and
# cannot place a floor at all; a logistic can, but is still committed to a symmetric
# transition."
#
# TWO CRITERIA, WHICH PULL IN OPPOSITE DIRECTIONS.
#
# 1. FIT WHERE IT MATTERS. Latitude is read from the twilight band, so a family that
#    wins in daylight and loses at twilight is worse than useless. Reported
#    separately for 85-100 deg.
#
# 2. DEGENERACY ADDED. The emission is already sloppy -- condition number 3.8e4 with
#    latitude at weight 0.881 in the softest eigenvector. A transition-asymmetry
#    parameter is exactly the sort of thing that trades against day length, hence
#    against latitude. So each family is also scored on the GAUGE INFLATION it
#    induces: se(latitude offset) with the response unknown, divided by the same with
#    the response known. 1.0 = latitude untouched by response uncertainty.
#
# The right family maximises twilight fit per unit of inflation, not fit alone.
#
# FAMILIES
#   linear     clamp(a - b z)                       2 par, current engine default
#   logistic   floor + amp/(1+exp((z-z50)/s))       4 par, symmetric transition
#   gompertz   floor + amp*exp(-exp((z-z50)/s))     4 par, ASYMMETRIC
#   richards   floor + amp/(1+nu*exp((z-z50)/s))^(1/nu)  5 par, asymmetry as a free
#                                                   parameter; logistic at nu = 1
#   logslope   logistic + a linear daytime term     5 par, holds the sloping day
#
# Fitted to each tag's CLEAR-SKY ENVELOPE (upper quantile of baselined light per
# zenith bin), which is what `fit_light_response` targets: shading only ever pushes
# light down, so the envelope is the signal and the spread below it is the noise.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))
suppressMessages(devtools::load_all(".", quiet = TRUE))
BAT <- readRDS("scratch/nes_calibration/light_battery.rds")
REC <- Filter(function(z) z$offset == 0, BAT)
arch <- lapply(readRDS("analysis/cache/nes/archives_v3_30min.rds"), `[[`, "main")
fv <- try(lapply(readRDS("analysis/cache/nes/fvilches_archives_v1_30min.rds"),
                 `[[`, "main"), silent = TRUE)
if (!inherits(fv, "try-error")) arch <- c(arch, fv)

# ---- build each tag's clear-sky envelope on the baselined scale ---------------
ENV <- list()
for (b in REC) {
  g <- arch[[b$id]]; if (is.null(g) || !nrow(g)) next
  ML <- b$max_light
  gt <- as.numeric(g$time); bt <- as.numeric(b$time)
  base <- as.numeric(quantile(as.numeric(g$light), 0.05, na.rm = TRUE))
  y <- pmax(0, approx(gt, as.numeric(g$light), bt, rule = 2)$y - base)
  ok <- b$supported & is.finite(y)
  z <- solar_zenith(bt[ok], lon360(b$lon)[ok], b$lat[ok])
  d <- data.table(z = z, y = y[ok])
  d[, zb := round(z)]
  e <- d[, .(n = .N, y = quantile(y, 0.90)), by = zb][n >= 20][order(zb)]
  if (nrow(e) < 25) next
  ENV[[b$id]] <- data.table(id = b$id, z = as.numeric(e$zb), y = e$y / ML)  # units of ML
}
E <- rbindlist(ENV)
cat(sprintf("%d tags, %d envelope points, zenith %.0f-%.0f\n\n",
            uniqueN(E$id), nrow(E), min(E$z), max(E$z)))

# ---- the families ------------------------------------------------------------
fam <- list(
  linear = list(p = 2, f = function(z, th) pmin(pmax(th[1] - th[2] * z, 0), 1),
                s = function(d) c(3, 0.03)),
  logistic = list(p = 4,
    f = function(z, th) pmin(pmax(th[1] + th[2] / (1 + exp((z - th[3]) / th[4])), 0), 1),
    s = function(d) c(0.1, 0.8, 95, 6)),
  gompertz = list(p = 4,
    f = function(z, th) pmin(pmax(th[1] + th[2] * exp(-exp((z - th[3]) / th[4])), 0), 1),
    s = function(d) c(0.1, 0.8, 95, 6)),
  richards = list(p = 5,
    f = function(z, th) {
      nu <- max(th[5], 0.05)
      pmin(pmax(th[1] + th[2] / (1 + nu * exp((z - th[3]) / th[4]))^(1 / nu), 0), 1)
    }, s = function(d) c(0.1, 0.8, 95, 6, 1)),
  logslope = list(p = 5,
    f = function(z, th) pmin(pmax(th[1] + th[2] / (1 + exp((z - th[3]) / th[4])) -
                                    th[5] * (z - 70), 0), 1),
    s = function(d) c(0.1, 0.8, 95, 6, 0.002)))

fit_fam <- function(d, fm) {
  obj <- function(th) sum((fm$f(d$z, th) - d$y)^2)
  best <- NULL
  for (jit in 1:6) {
    st <- fm$s(d) * (if (jit == 1) 1 else runif(fm$p, 0.7, 1.4))
    o <- try(optim(st, obj, method = "Nelder-Mead",
                   control = list(maxit = 8000, reltol = 1e-12)), silent = TRUE)
    if (!inherits(o, "try-error") && (is.null(best) || o$value < best$value)) best <- o
  }
  best
}

set.seed(1)
res <- list()
for (tg in unique(E$id)) {
  d <- E[id == tg]
  tw <- d$z >= 85 & d$z <= 100
  for (nm in names(fam)) {
    o <- fit_fam(d, fam[[nm]])
    if (is.null(o)) next
    r <- fam[[nm]]$f(d$z, o$par) - d$y
    res[[length(res) + 1]] <- data.table(
      id = tg, family = nm, p = fam[[nm]]$p,
      rmse_all = sqrt(mean(r^2)),
      rmse_twi = sqrt(mean(r[tw]^2)),
      aic = nrow(d) * log(mean(r^2)) + 2 * fam[[nm]]$p)
  }
}
R <- rbindlist(res)
cat("=== fit to the clear-sky envelope (units of max_light) ===\n")
S <- R[, .(p = p[1],
           rmse_all = round(mean(rmse_all), 4),
           rmse_twilight = round(mean(rmse_twi), 4),
           aic = round(mean(aic), 1)), by = family][order(rmse_twilight)]
print(as.data.frame(S), row.names = FALSE)
cat("\n  ranked by TWILIGHT rmse (85-100 deg), which is where latitude is read.\n")
cat(sprintf("  best twilight fit: %s (%.4f) vs the current linear (%.4f) -- %.1fx better\n",
            S$family[1], S$rmse_twilight[1],
            S[family == "linear"]$rmse_twilight,
            S[family == "linear"]$rmse_twilight / S$rmse_twilight[1]))

# ---- criterion 2: what does each family do to the gauge? ---------------------
cat("\n=== gauge inflation induced by each family ===\n")
cat("  se(latitude offset) with the response UNKNOWN / with it KNOWN, on a full\n")
cat("  240-day track. 1.0 means latitude is untouched by response uncertainty.\n\n")
b <- REC[[1]]
ML <- b$max_light; lam <- 1 / (ML * 0.5); LH <- lam * 2; PS <- 0.10
tn <- as.numeric(b$time)
sel <- which(b$supported)[seq(1, sum(b$supported), by = 8)]
tt <- tn[sel]; yy <- b$noisy[sel]; la <- b$lat[sel]; lo <- lon360(b$lon)[sel]

# UNCLAMPED forms for the information matrix. The fitted families clamp to [0,1],
# and a numerical Hessian taken across a clamp boundary is meaningless -- the first
# attempt at this returned an indefinite matrix for every multi-parameter family for
# exactly that reason. The clamp is irrelevant to identifiability anyway: it is a
# range restriction on the output, not a constraint linking parameters.
raw_f <- list(
  linear   = function(z, th) th[1] - th[2] * z,
  logistic = function(z, th) th[1] + th[2] / (1 + exp((z - th[3]) / th[4])),
  gompertz = function(z, th) th[1] + th[2] * exp(-exp((z - th[3]) / th[4])),
  richards = function(z, th) {
    nu <- max(th[5], 0.05); th[1] + th[2] / (1 + nu * exp((z - th[3]) / th[4]))^(1 / nu)
  },
  logslope = function(z, th) th[1] + th[2] / (1 + exp((z - th[3]) / th[4])) -
                             th[5] * (z - 70))

# PER-OBSERVATION log-likelihood (a vector, not a sum) -- needed for BHHH below.
logLi <- function(th, np, nm) {
  z <- solar_zenith(tt, lo, la + th[1])
  mu <- pmin(pmax(raw_f[[nm]](z, th[2:(np + 1)]), 1e-6), 1 - 1e-6) * ML
  raw <- ifelse(yy <= mu, lam * exp(-lam * (mu - yy)), lam * exp(-LH * (yy - mu)))
  lo_ <- (1 - exp(-lam * mu)) / lam
  hi_ <- (1 - exp(-LH * (ML - mu))) / LH
  log((1 - PS) * (raw / pmax((lo_ + hi_) * lam, 1e-12)) + PS / ML)
}

# INFORMATION BY OUTER PRODUCT OF GRADIENTS (BHHH), not by a numerical Hessian.
#
# The asymmetric-Laplace spike has a kink at y = mu, so the log-likelihood is C0 but
# not C2 in the parameters. A numerical second derivative therefore picks up kink
# crossings and gets NOISIER as the step shrinks -- which is exactly what happened
# here: every family came back indefinite, and shrinking the steps made the one that
# had worked (linear) fail too. Earlier Fisher work in this directory escaped it only
# because it used noiseless data, where every observation sits coherently at the kink.
#
# First derivatives exist almost everywhere and are well behaved, so BHHH --
# I = sum_i s_i s_i' with s_i the per-observation score -- is the right estimator
# here, and is guaranteed positive semi-definite by construction.
bhhh <- function(np, nm, th, h) {
  p <- length(th)
  S <- matrix(0, length(yy), p)
  for (i in seq_len(p)) {
    e <- numeric(p); e[i] <- h[i]
    S[, i] <- (logLi(th + e, np, nm) - logLi(th - e, np, nm)) / (2 * h[i])
  }
  crossprod(S)
}
d1 <- E[id == b$id]
for (nm in names(fam)) {
  fm <- fam[[nm]]; np <- fm$p
  o <- fit_fam(d1, fm); if (is.null(o)) next
  th <- c(0, o$par)
  # Steps must reflect each parameter's MEANINGFUL scale, not its magnitude. A fixed
  # fraction fails badly here: 5% of z50 (~95) is a 4.75 degree finite difference,
  # far wider than the curve's own structure, which is what made the first two
  # attempts return indefinite matrices. floor/amp live in units of max_light (0-1),
  # z50 and scale in degrees, nu is dimensionless.
  hstep <- switch(nm,
    linear   = c(0.02, 2e-4),
    logistic = ,
    gompertz = c(0.005, 0.01, 0.25, 0.10),
    richards = c(0.005, 0.01, 0.25, 0.10, 0.05),
    logslope = c(0.005, 0.01, 0.25, 0.10, 2e-4))
  h <- c(0.25, hstep)
  F <- bhhh(np, nm, th, h)
  ev <- eigen(F, symmetric = TRUE)$values
  if (min(ev) <= 0 || !all(is.finite(ev))) {
    cat(sprintf("  %-9s  (information matrix singular)\n", nm)); next
  }
  infl <- sqrt(solve(F)[1, 1]) * sqrt(F[1, 1])
  cat(sprintf("  %-9s p=%d  cond %8.2g  se(lat) %.4f  inflation %.2f\n",
              nm, np, max(ev) / min(ev), sqrt(solve(F)[1, 1]), infl))
}

cat("\nREADING\n")
cat("  Take the family with the best TWILIGHT fit whose inflation stays near 1.\n")
cat("  A family that fits better but inflates the gauge has bought precision in the\n")
cat("  response at the cost of identifiability in latitude, which is the trade this\n")
cat("  model can least afford.\n")
