# WHY DID THE FITTED RESPONSE WIDEN 40%?
#
# Refitting on the surface-conditioned decimation moved pooled z50 from 91.92 to
# 89.56 and the width from 31.0 to 43.0 deg. z50 and width are not independent, so
# the z50 change cannot be read as "the max bias coming out" while the shape is
# also moving, and a transition spanning 68-111 deg of zenith is not physical:
# light should still be saturated at 22 deg solar elevation.
#
# `fit_light_response()` fits to the EMPIRICAL ENVELOPE -- the env_q quantile of
# light in each 1-degree zenith bin, forced non-increasing. So the fitted shape is
# downstream of that envelope, and the envelope is directly comparable between the
# two decimations. This prints it, per tag, over the calibration window, at the
# deployment position where zenith is known.
#
# The question it answers: did the surface statistic genuinely reveal a wider
# response that the max was artificially sharpening, or has the envelope become
# noisy/flat because each zenith bin now holds fewer and dimmer observations?
suppressMessages({ library(data.table) })
source("analysis/diagnostics/nes_common.R")

ENV_Q <- 0.95     # matches fit_light_response's default
CAL_DAYS <- 15
O <- readRDS("analysis/cache/nes/panel_v1.rds")
N <- readRDS("analysis/cache/nes/panel_surface_v1.rds")

env_of <- function(P, tg) {
  g <- P$tags[[tg]]
  L <- as.data.table(g$light)
  dep <- departure_time(L)
  k <- if (is.finite(dep)) as.numeric(L$time) < dep else
       as.numeric(L$time) < min(as.numeric(L$time)) + CAL_DAYS * 86400
  L <- L[k]
  if (nrow(L) < 50) return(NULL)
  # solar_zenith does not recycle scalars the way fit_light_response does
  z <- solar_zenith(as.numeric(L$time),
                    rep(unname(g$p0[1]), nrow(L)), rep(unname(g$p0[2]), nrow(L)))
  d <- data.table(zb = round(z), li = as.numeric(L$light))
  d[, .(n = .N, env = as.numeric(quantile(li, ENV_Q)),
        med = as.numeric(median(li))), by = zb][n >= 3][order(zb)]
}

tags <- intersect(names(O$tags), names(N$tags))
tags <- tags[vapply(tags, function(t) !is.null(O$geom_fits[[t]]), TRUE)][1:3]

for (tg in tags) {
  eo <- env_of(O, tg); en <- env_of(N, tg)
  if (is.null(eo) || is.null(en)) next
  m <- merge(eo[, .(zb, n_old = n, env_old = env, med_old = med)],
             en[, .(zb, n_new = n, env_new = env, med_new = med)], by = "zb", all = TRUE)
  m <- m[zb >= 60 & zb <= 115]
  cat(sprintf("\n=== tag %s : empirical envelope (%.2f quantile of light per zenith degree) ===\n",
              tg, ENV_Q))
  cat("  zenith  n_old n_new | env_old env_new  d_env | med_old med_new\n")
  s <- m[seq(1, nrow(m), by = 4)]
  for (i in seq_len(nrow(s))) cat(sprintf("   %5.0f  %5s %5s | %7s %7s %6s | %7s %7s\n",
    s$zb[i], s$n_old[i], s$n_new[i],
    round(s$env_old[i], 1), round(s$env_new[i], 1),
    round(s$env_new[i] - s$env_old[i], 1),
    round(s$med_old[i], 1), round(s$med_new[i], 1)))
  # where does the envelope leave saturation, and where does it reach the floor?
  sat <- function(e) { if (all(is.na(e$env))) return(NA_real_)
    top <- max(e$env, na.rm = TRUE); e$zb[which(e$env < 0.95 * top)[1]] }
  flo <- function(e) { if (all(is.na(e$env))) return(NA_real_)
    top <- max(e$env, na.rm = TRUE); bot <- min(e$env, na.rm = TRUE)
    e$zb[which(e$env < bot + 0.05 * (top - bot))[1]] }
  cat(sprintf("  leaves saturation at zenith: old %5.0f  new %5.0f\n", sat(eo), sat(en)))
  cat(sprintf("  reaches floor at zenith    : old %5.0f  new %5.0f\n", flo(eo), flo(en)))
  cat(sprintf("  implied transition span    : old %5.0f  new %5.0f deg\n",
              flo(eo) - sat(eo), flo(en) - sat(en)))
}
cat("\nIf the NEW span is genuinely wider, the max was sharpening the response and\n")
cat("43 deg may be right. If the new ENVELOPE is merely noisier at the same span,\n")
cat("the width change is a fitting artefact of fewer/dimmer observations per bin.\n")
