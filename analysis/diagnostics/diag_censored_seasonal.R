# Censoring halves the tilt. Does it remove the SEASONAL part -- the sign flip?
#
# The headline pathology is not the mean tilt, it is that the bias points
# equatorward in northern winter (-2.156 deg at declination -7) and POLEWARD in
# northern summer (+0.101 at +22). A fix that halves the average but leaves the
# seasonal swing intact would still produce season-dependent errors, which is the
# thing that would be picked apart in a habitat or climate application.
#
# So: per-knot tilt under both likelihoods, keyed on SIGNED declination, over all
# four date offsets and two records. Pure numerics on the validated R twin.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))
suppressMessages(devtools::load_all(".", quiet = TRUE))
BAT <- readRDS("scratch/nes_calibration/light_battery.rds")
PSLAB <- 0.10; RATIO <- 2

tilt_record <- function(b) {
  ML <- b$max_light; LAM <- 1 / (ML * 0.5); LAM_HI <- LAM * RATIO
  CC <- LAM * LAM_HI / (LAM + LAM_HI)
  r <- b$response
  exl <- function(z) pmin(pmax(r[1] + r[2] / (1 + exp((z - r[3]) / r[4])), 0), ML)

  ltrunc <- function(y, mu) {
    raw <- ifelse(y <= mu, LAM * exp(-LAM * (mu - y)), LAM * exp(-LAM_HI * (y - mu)))
    lo <- (1 - exp(-LAM * mu)) / LAM
    hi <- (1 - exp(-LAM_HI * (ML - mu))) / LAM_HI
    log((1 - PSLAB) * (raw / pmax((lo + hi) * LAM, 1e-12)) + PSLAB / ML)
  }
  lcens <- function(y, mu) {
    dens <- ifelse(y <= mu, CC * exp(-LAM * (mu - y)), CC * exp(-LAM_HI * (y - mu)))
    a_lo <- (LAM_HI / (LAM + LAM_HI)) * exp(-LAM * mu)
    a_hi <- (LAM  / (LAM + LAM_HI)) * exp(-LAM_HI * (ML - mu))
    sp <- ifelse(y >= ML - 1e-9, a_hi, ifelse(y <= 1e-9, a_lo, dens))
    log((1 - PSLAB) * sp + PSLAB / ML)
  }

  tn <- as.numeric(b$time)
  kb <- seq(min(tn), max(tn), by = 12 * 3600)
  dl <- 0.5
  out <- list()
  for (k in seq_len(length(kb) - 1)) {
    sel <- tn > kb[k] & tn <= kb[k + 1]
    if (sum(sel) < 6) next
    t_k <- tn[sel]; y_k <- b$perfect[sel]; m <- sum(sel)
    lat0 <- approx(tn, b$lat, mean(t_k), rule = 2)$y
    lon0 <- approx(tn, b$lon, mean(t_k), rule = 2)$y
    mu_up <- exl(solar_zenith(t_k, rep(lon0, m), rep(lat0 + dl, m)))
    mu_dn <- exl(solar_zenith(t_k, rep(lon0, m), rep(lat0 - dl, m)))
    out[[length(out) + 1]] <- data.table(
      id = b$id, offset = b$offset,
      decl = solar_declination(as.POSIXct(mean(t_k), origin = "1970-01-01", tz = "UTC")),
      lat = lat0,
      trunc = (sum(ltrunc(y_k, mu_up)) - sum(ltrunc(y_k, mu_dn))) / (2 * dl),
      cens  = (sum(lcens(y_k, mu_up))  - sum(lcens(y_k, mu_dn)))  / (2 * dl))
  }
  rbindlist(out)
}

RECS <- Filter(function(z) z$id %in% c("2021033", "2023032"), BAT)
D <- rbindlist(lapply(RECS, tilt_record))
cat(sprintf("%d knots, %d records x %d offsets\n\n",
            nrow(D), uniqueN(D$id), uniqueN(D$offset)))

D[, sband := cut(decl, c(-24, -15, -7, 7, 15, 24), include.lowest = TRUE,
                 labels = c("<-15 (NH winter)", "-15..-7", "-7..+7 (equinox)",
                            "+7..+15", ">+15 (NH summer)"))]
cat("=== per-knot latitude tilt (nats/deg) by SIGNED declination ===\n")
print(as.data.frame(D[, .(n = .N,
                          truncated = round(mean(trunc), 4),
                          censored  = round(mean(cens), 4),
                          reduction = sprintf("%.0f%%", 100 * (1 - abs(mean(cens)) / abs(mean(trunc)))))
                      , by = sband][order(sband)]), row.names = FALSE)

sw <- function(v) diff(range(D[, .(m = mean(get(v))), by = sband]$m))
cat(sprintf("\n  SEASONAL SWING (max - min across bands):\n"))
cat(sprintf("    truncated %.4f nats/deg\n    censored  %.4f nats/deg   (%.0f%% smaller)\n",
            sw("trunc"), sw("cens"), 100 * (1 - sw("cens") / sw("trunc"))))

cat(sprintf("\n  overall mean: truncated %+0.4f | censored %+0.4f\n",
            mean(D$trunc), mean(D$cens)))
cat(sprintf("  sign flips across season? truncated %s | censored %s\n",
            if (any(D[, mean(trunc), by = sband]$V1 > 0) &&
                any(D[, mean(trunc), by = sband]$V1 < 0)) "YES" else "no",
            if (any(D[, mean(cens), by = sband]$V1 > 0) &&
                any(D[, mean(cens), by = sband]$V1 < 0)) "YES" else "no"))
