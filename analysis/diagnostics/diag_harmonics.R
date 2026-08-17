# A position-free measure of how much LATITUDE information a stretch of light
# actually carries, derived from the transform rather than guessed at.
#
# THE ARGUMENT. Geolocation from light is phase-and-duty demodulation of a
# one-cycle-per-day carrier: longitude is the phase, latitude is the duty cycle
# (day length). For a duty-cycle-D waveform the Fourier coefficients go as
# sin(n pi D)/(n pi), so at D = 0.5 -- a twelve-hour day -- EVERY EVEN HARMONIC
# IS IDENTICALLY ZERO. The day-length information lives in the even harmonics
# and nulls exactly at the equinox, which is the same degeneracy as the
# vanishing Jacobian dH/dphi = -tan(delta) sec^2(phi)/sin(H), but expressed as
# something measurable in the data.
#
# WHY THIS MATTERS. Every previous attempt at a position-free reliability
# measure was a hand-built statistic and reached only r ~ 0.54, p = 0.11. This
# one is not guessed: it is the amplitude of the very harmonics that carry the
# latitude signal, computed by FFT from the light series with NO position, no
# track and no ground truth.
#
# It cannot fix bias. The claim under test is narrower and is the actual
# coverage problem: the model reports near-identical uncertainty regardless of
# how much latitude information it has, and this would tell it.
#
# TEST, per tag per 20-day window:
#   E2 = power at 2 cycles/day relative to 1 cycle/day  (the even harmonic)
#   does E2 predict the realised latitude error in that window?
# Prediction: MORE even-harmonic power -> more day-length information -> SMALLER
# latitude error. A negative correlation is the result being sought.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))
WIN_DAYS <- 20

declination <- function(t) {
  doy <- as.numeric(format(t, "%j")); g <- 2*pi/365 * (doy - 1)
  (0.006918 - 0.399912*cos(g) + 0.070257*sin(g) - 0.006758*cos(2*g) +
     0.000907*sin(2*g) - 0.002697*cos(3*g) + 0.00148*sin(3*g)) * 180/pi
}

# Power at k cycles/day, by direct projection onto the complex exponential.
# Robust to gaps and uneven spacing, unlike an FFT on a padded series.
harm <- function(tt, y, k) {
  y <- y - mean(y)
  w <- 2*pi*k*(as.numeric(tt)/86400)
  sqrt(sum(y*cos(w))^2 + sum(y*sin(w))^2) / length(y)
}

arch <- lapply(readRDS("analysis/cache/nes/archives_v3_30min.rds"), `[[`, "main")
k <- as.data.table(read.csv("analysis/output/nes/per_knot_errors.csv"))
k <- k[engine == "hier" & batch == "light" & is.finite(err_lat)]
k[, `:=`(id = as.character(id), time = as.POSIXct(time, tz = "UTC"))]

rows <- list()
for (tg in names(arch)) {
  x <- as.data.table(arch[[tg]]); kk <- k[id == tg]
  if (!nrow(kk)) next
  x[, w := floor(as.numeric(time) / (WIN_DAYS*86400))]
  for (ww in unique(x$w)) {
    s <- x[w == ww]
    if (nrow(s) < WIN_DAYS*24) next                 # need most of the window
    rng <- range(s$time)
    e <- kk[time >= rng[1] & time <= rng[2]]
    if (nrow(e) < 10) next
    h1 <- harm(s$time, s$light, 1)
    h2 <- harm(s$time, s$light, 2)
    h3 <- harm(s$time, s$light, 3)
    rows[[length(rows)+1]] <- data.frame(
      id = tg, win = ww,
      decl = round(mean(declination(s$time)), 1),
      h1 = round(h1, 2), h2 = round(h2, 2),
      E2 = round(h2 / h1, 4),                       # even-harmonic share
      E2_abs = round(h2, 3),
      odd_even = round(h3 / pmax(h2, 1e-9), 3),
      abs_err_lat = round(mean(abs(e$err_lat)), 3),
      sd_err_lat = round(sd(e$err_lat), 3),
      n = nrow(e), row.names = NULL)
  }
}
D <- as.data.table(rbindlist(rows))
cat(sprintf("%d windows of %d days across %d tags\n\n", nrow(D), WIN_DAYS, uniqueN(D$id)))

cat("=== does the even harmonic null at the equinox, as the algebra says? ===\n")
D[, dbin := cut(abs(decl), breaks = c(0, 5, 10, 15, 20, 24))]
print(as.data.frame(D[!is.na(dbin), .(windows = .N, mean_decl = round(mean(abs(decl)), 1),
  E2 = round(mean(E2), 4), abs_err_lat = round(mean(abs_err_lat), 2)),
  by = dbin][order(dbin)]), row.names = FALSE)

cat("\n=== does E2 predict the latitude error? ===\n")
for (v in c("E2", "E2_abs", "h2")) {
  ct <- cor.test(D[[v]], D$abs_err_lat, method = "spearman")
  cat(sprintf("  %-7s vs |latitude error|: rho = %+.3f (p = %.4f)\n",
              v, ct$estimate, ct$p.value))
}
# within tag, so that between-animal differences cannot produce it
D[, `:=`(E2_c = E2 - mean(E2), err_c = abs_err_lat - mean(abs_err_lat)), by = id]
ct <- cor.test(D$E2_c, D$err_c, method = "spearman")
cat(sprintf("\n  WITHIN TAG (each animal centred): rho = %+.3f (p = %.4f, n = %d)\n",
            ct$estimate, ct$p.value, nrow(D)))

cat("\n=== is it just declination in disguise? ===\n")
m1 <- lm(abs_err_lat ~ abs(decl), data = D)
m2 <- lm(abs_err_lat ~ abs(decl) + E2, data = D)
cat(sprintf("|declination| alone      R2 = %.3f\n", summary(m1)$r.squared))
cat(sprintf("|declination| + E2       R2 = %.3f  (increment %.3f, p = %.4g)\n",
            summary(m2)$r.squared, summary(m2)$r.squared - summary(m1)$r.squared,
            anova(m1, m2)$`Pr(>F)`[2]))
cat("\nE2 is worth having only if it adds to what the DATE already tells you.\n")
cat("Declination is known exactly with no data at all; a harmonic measure has to\n")
cat("beat that to be interesting, by responding to the tag and the animal.\n")
saveRDS(D, file.path(SCRATCH, "harmonics.rds"))
