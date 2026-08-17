# DOES SEASONAL SPAN, OR POOLING ACROSS TAGS, BREAK THE GAUGE?
#
# Established (diag_fisher_sloppy.R, 553 one-day windows): the emission is SLOPPY --
# condition number 3.8e4 -- and latitude is essentially THE soft direction within a
# day (weight 0.881 in the softest eigenvector, 0.005 in the stiffest). The exchange
# rate against the threshold z50 CHANGES SIGN across the year (+1.69 in NH winter,
# -4.43 at equinox), so the soft direction ROTATES seasonally.
#
# That rotation is the hope: a single (z50, amp) error cannot satisfy both exchange
# rates at once, so a track or a panel spanning seasons should IDENTIFY the
# calibration that no single day can. This measures whether it actually does.
#
# THE PARAMETERISATION, chosen so the question is well posed. Over a multi-day window
# the animal's latitude changes, so a single "latitude" parameter is meaningless.
# Instead the free parameter is a UNIFORM LATITUDE OFFSET dlat applied to the known
# path, alongside shared calibration offsets. That asks exactly the right question:
# can a global latitude displacement be absorbed by a change in calibration?
#
#     theta = (dlat, dz50, dscale, damp, dfloor)
#
# THE DECISION-RELEVANT NUMBER is the GAUGE INFLATION FACTOR:
#     se(dlat) with calibration UNKNOWN (profiled out, from the inverse FIM)
#   / se(dlat) with calibration KNOWN   (from the FIM diagonal)
# = how much not knowing the calibration costs you in latitude. 1.0 means no gauge
# problem; large means latitude is hostage to calibration.
#
# PART 1: inflation vs window span within one track (1 day to the full 237).
# PART 2: inflation vs number of tags pooled, each with its OWN dlat but a SHARED
#         calibration offset -- the partial-pooling structure, as a joint estimation
#         rather than a pre-fit.
#
# PRE-REGISTERED: if the seasonal contrast identifies the calibration, inflation
# should fall steeply once the window spans a solstice-to-equinox range (~90 days),
# and fall further with pooled tags. If it stays high at full span AND with 6 tags,
# the gauge survives, pooling will not fix it, and the honest conclusion is that
# latitude is only weakly identified in this model class.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))
suppressMessages(devtools::load_all(".", quiet = TRUE))
BAT <- readRDS("scratch/nes_calibration/light_battery.rds")
REC <- Filter(function(z) z$offset == 0, BAT)
PSLAB <- 0.10; RATIO <- 2

# log-likelihood for a list of tag-windows sharing calibration offsets.
# theta = (dlat_1..dlat_n, dz50, dscale, damp, dfloor)
mk_logL <- function(W) {
  n <- length(W)
  function(th) {
    dz50 <- th[n + 1]; dsc <- th[n + 2]; damp <- th[n + 3]; dflo <- th[n + 4]
    tot <- 0
    for (i in seq_len(n)) {
      w <- W[[i]]; ML <- w$ML; lam <- w$lam; lam_hi <- lam * RATIO; r <- w$r
      z <- solar_zenith(w$t, w$lon, w$lat + th[i])
      mu <- pmin(pmax(r[1] + dflo * ML + r[2] * (1 + damp) /
                        (1 + exp((z - (r[3] + dz50)) / (r[4] * (1 + dsc)))), 0), ML)
      raw <- ifelse(w$y <= mu, lam * exp(-lam * (mu - w$y)),
                    lam * exp(-lam_hi * (w$y - mu)))
      lo <- (1 - exp(-lam * mu)) / lam
      hi <- (1 - exp(-lam_hi * (ML - mu))) / lam_hi
      tot <- tot + sum(log((1 - PSLAB) * (raw / pmax((lo + hi) * lam, 1e-12)) +
                             PSLAB / ML))
    }
    tot
  }
}

hess <- function(f, th, h) {
  p <- length(th); H <- matrix(0, p, p)
  for (i in seq_len(p)) for (j in i:p) {
    ei <- ej <- numeric(p); ei[i] <- h[i]; ej[j] <- h[j]
    H[i, j] <- H[j, i] <-
      (f(th + ei + ej) - f(th + ei - ej) - f(th - ei + ej) + f(th - ei - ej)) /
      (4 * h[i] * h[j])
  }
  H
}

win <- function(b, days, thin = 4) {
  tn <- as.numeric(b$time)
  sel <- which(tn <= min(tn) + days * 86400 & b$supported)
  sel <- sel[seq(1, length(sel), by = thin)]
  list(t = tn[sel], y = b$perfect[sel],
       lon = lon360(b$lon)[sel], lat = b$lat[sel],
       ML = b$max_light, lam = 1 / (b$max_light * 0.5), r = b$response,
       span = diff(range(solar_declination(b$time[sel]))))
}

gauge <- function(W) {
  n <- length(W); p <- n + 4
  f <- mk_logL(W)
  th <- numeric(p)
  h <- c(rep(0.25, n), 0.25, 0.02, 0.01, 0.01)
  F <- -hess(f, th, h)
  F <- (F + t(F)) / 2
  ev <- eigen(F, symmetric = TRUE)$values
  if (min(ev) <= 0) return(NULL)
  Fi <- solve(F)
  se_unknown <- sqrt(diag(Fi)[seq_len(n)])
  se_known   <- 1 / sqrt(diag(F)[seq_len(n)])
  list(cond = max(ev) / min(ev),
       infl = mean(se_unknown / se_known),
       se_lat = mean(se_unknown),
       se_z50 = sqrt(Fi[n + 1, n + 1]))
}

cat("=== PART 1: one track, increasing window span ===\n")
b <- REC[[1]]
cat(sprintf("  track %s\n", b$id))
cat(sprintf("%8s %10s %12s %12s %12s\n",
            "days", "decl span", "cond", "se(dlat)", "INFLATION"))
for (d in c(1, 2, 5, 10, 20, 45, 90, 160, 240)) {
  W <- list(win(b, d))
  if (length(W[[1]]$t) < 50) next
  g <- gauge(W)
  if (is.null(g)) { cat(sprintf("%8.0f  (indefinite)\n", d)); next }
  cat(sprintf("%8.0f %10.1f %12.3g %12.4f %12.2f\n",
              d, W[[1]]$span, g$cond, g$se_lat, g$infl))
}

cat("\n=== PART 2: full tracks, pooling tags with a SHARED calibration ===\n")
cat(sprintf("%8s %12s %12s %12s %12s\n",
            "n tags", "cond", "se(dlat)", "se(dz50)", "INFLATION"))
for (n in 1:length(REC)) {
  W <- lapply(REC[seq_len(n)], win, days = 240, thin = 8)
  g <- gauge(W)
  if (is.null(g)) { cat(sprintf("%8d  (indefinite)\n", n)); next }
  cat(sprintf("%8d %12.3g %12.4f %12.4f %12.2f\n",
              n, g$cond, g$se_lat, g$se_z50, g$infl))
}

cat("\nREADING\n")
cat("  INFLATION = se(latitude offset) with calibration unknown, divided by the same\n")
cat("  with calibration known. 1.0 = no gauge problem. Large = latitude is hostage\n")
cat("  to calibration, and no within-track experiment can identify it.\n")
cat("  Falling with span or with pooled tags => the seasonal contrast breaks the\n")
cat("  gauge, and estimating calibration as a shared latent is well founded.\n")
