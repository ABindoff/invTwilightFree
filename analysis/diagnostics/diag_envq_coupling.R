# env_q IS COUPLED TO THE DECIMATION STATISTIC.
#
# fit_light_response() takes the env_q (default 0.95) quantile of light in each
# 1-degree zenith bin as the clear-sky envelope. Calibration bins hold only 3-9
# observations, so a 0.95 quantile there is effectively the per-bin MAXIMUM.
#
# The envelope's PURPOSE is to reject shaded observations. That made sense when
# each observation was a 30-minute max (already near-surface but still
# contaminated). Under surface conditioning the rejection has already happened, so
# a 0.95 quantile no longer removes contamination -- it selects the upper noise
# tail, and a monotone fit to a noisy envelope comes out shallow and wide. That is
# the likely source of width 31.0 -> 43.0 deg and z50 91.92 -> 89.56.
#
# If that account is right, LOWERING env_q on the surface decimation should recover
# a physically sensible width, and the recovered z50 should then differ from the
# shipped 91.92 by something near the +1.42 deg of measured max bias rather than
# by 2.36 deg.
#
# Reported per env_q: pooled z50 and width, plus how well the fitted response
# predicts the OBSERVED light at known zenith over the calibration window. The
# residual is the model-selection criterion; width alone is only a plausibility
# check.
suppressMessages({ library(data.table) })
source("analysis/diagnostics/nes_common.R")
CAL_DAYS <- 15
O <- readRDS("analysis/cache/nes/panel_v1.rds")
N <- readRDS("analysis/cache/nes/panel_surface_v1.rds")

# light predicted by a fitted response at zenith z: the package parameterises the
# transition by calibration = c(intercept, slope) on light, so invert that.
resid_of <- function(P, resp, tg) {
  g <- P$tags[[tg]]; r <- resp[[tg]]
  if (is.null(r)) return(NA_real_)
  L <- as.data.table(g$light)
  k <- as.numeric(L$time) < min(as.numeric(L$time)) + CAL_DAYS * 86400
  L <- L[k]; if (nrow(L) < 50) return(NA_real_)
  z <- solar_zenith(as.numeric(L$time), rep(unname(g$p0[1]), nrow(L)),
                    rep(unname(g$p0[2]), nrow(L)))
  # fitted mean light at that zenith, on the same scale as the observations
  pred <- r$baseline + (r$max_light - r$baseline) /
          (1 + exp((z - r$z50) / max(1e-6, r$scale)))
  obs <- as.numeric(L$light)
  ok <- is.finite(pred) & is.finite(obs)
  if (!any(ok)) return(NA_real_)
  sqrt(mean((obs[ok] - pred[ok])^2))
}

fit_panel <- function(P, env_q) {
  sf <- lapply(P$tags, function(g) {
    k <- g$light$time <= min(g$light$time) + CAL_DAYS * 86400
    fit_light_response(g$light$time[k], g$light$light[k], g$p0[1], g$p0[2], env_q = env_q)
  })
  gf <- lapply(P$tags, function(g) {
    dep <- departure_time(g$light)
    k <- is.finite(dep) & as.numeric(g$light$time) < dep
    if (sum(k) < 200) return(NULL)
    fit_light_response(g$light$time[k], g$light$light[k], g$p0[1], g$p0[2], env_q = env_q)
  })
  fam <- vapply(P$tags, function(g) g$family, "")
  out <- list(); z50 <- c(); wid <- c()
  for (f in unique(fam)) {
    ids <- names(P$tags)[fam == f]
    if (!sum(!vapply(gf[ids], is.null, logical(1)))) next
    pooled <- pool_light_responses(gf[ids], sf[ids])
    p <- attr(pooled, "pooled")
    z50[f] <- p[["z50"]]; wid[f] <- 4.394 * p[["scale"]]
    for (i in ids) out[[i]] <- pooled[[i]]
  }
  list(resp = out, z50 = z50, width = wid,
       ngeom = sum(!vapply(gf, is.null, logical(1))))
}

cat("=== SHIPPED decimation, env_q 0.95 (the published configuration) ===\n")
r0 <- fit_panel(O, 0.95)
cat(sprintf("  Mk9_219 z50 %6.2f  width %5.1f | geom %d/29 | mean RMSE vs observed light %5.1f units\n",
            r0$z50["Mk9_219"], r0$width["Mk9_219"], r0$ngeom,
            mean(vapply(names(O$tags), function(t) resid_of(O, r0$resp, t), 0), na.rm = TRUE)))

cat("\n=== SURFACE decimation, sweeping env_q ===\n")
for (q in c(0.95, 0.85, 0.70, 0.50)) {
  r <- fit_panel(N, q)
  rms <- mean(vapply(names(N$tags), function(t) resid_of(N, r$resp, t), 0), na.rm = TRUE)
  cat(sprintf("  env_q %.2f : Mk9_219 z50 %6.2f  width %5.1f | geom %d/29 | RMSE %5.1f units | z50 shift vs shipped %+6.2f deg\n",
              q, r$z50["Mk9_219"], r$width["Mk9_219"], r$ngeom, rms,
              r$z50["Mk9_219"] - r0$z50["Mk9_219"]))
}
cat("\nmeasured max-selection bias to be recovered: +1.42 deg of zenith\n")
cat("shipped width for reference: 31.0 deg\n")
