# The even-harmonic latitude-information measure, re-tested on all 29
# deployments across three seasons.
#
# At n = 10 it gave rho = -0.34 (p = 0.0002) and added to declination alone
# (R2 0.024 -> 0.081, p = 0.0098). Everything else in this project that reached
# p ~ 0.1 at n = 10 has since either failed or weakened with more data, so this
# is the honest test: does it survive tripling the panel and adding two seasons
# it was never fitted on?
#
# The measure: for a duty-cycle-D daily light waveform the Fourier coefficients
# go as sin(n pi D)/(n pi), so the even harmonics null at D = 0.5 -- a twelve
# hour day. Day-length information lives in those harmonics, and E2 = |H2|/|H1|
# measures how far from the null this stretch of light sits. Computed from the
# light series alone: no position, no track, no ground truth.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
suppressMessages({ library(data.table) })
WIN_DAYS <- 20

declination <- function(t) {
  doy <- as.numeric(format(t, "%j")); g <- 2*pi/365 * (doy - 1)
  (0.006918 - 0.399912*cos(g) + 0.070257*sin(g) - 0.006758*cos(2*g) +
     0.000907*sin(2*g) - 0.002697*cos(3*g) + 0.00148*sin(3*g)) * 180/pi
}
harm <- function(tt, y, k) {
  y <- y - mean(y); w <- 2*pi*k*(as.numeric(tt)/86400)
  sqrt(sum(y*cos(w))^2 + sum(y*sin(w))^2) / length(y)
}

old <- lapply(readRDS("analysis/cache/nes/archives_v3_30min.rds"), `[[`, "main")
new <- lapply(readRDS("analysis/cache/nes/fvilches_archives_v1_30min.rds"), `[[`, "main")
arch <- c(old, new)
k <- fread("scratch/nes_calibration/all29_knots.csv")
k <- k[is.finite(err_lat)]
k[, `:=`(id = as.character(id), time = as.POSIXct(time, tz = "UTC"))]

rows <- list()
for (tg in names(arch)) {
  x <- as.data.table(arch[[tg]]); kk <- k[id == tg]
  if (!nrow(kk)) next
  x[, w := floor(as.numeric(time) / (WIN_DAYS*86400))]
  for (ww in unique(x$w)) {
    s <- x[w == ww]
    if (nrow(s) < WIN_DAYS*24) next
    rng <- range(s$time)
    e <- kk[time >= rng[1] & time <= rng[2]]
    if (nrow(e) < 10) next
    h1 <- harm(s$time, s$light, 1); h2 <- harm(s$time, s$light, 2)
    rows[[length(rows)+1]] <- data.frame(
      id = tg, season = substr(tg, 1, 4), win = ww,
      decl = round(mean(declination(s$time)), 1),
      E2 = round(h2 / h1, 4),
      abs_err_lat = round(mean(abs(e$err_lat)), 3),
      lat_sd = round(mean(e$lat_sd, na.rm = TRUE), 3),
      n = nrow(e), row.names = NULL)
  }
}
D <- as.data.table(rbindlist(rows))
cat(sprintf("%d windows across %d deployments, %d seasons\n\n",
            nrow(D), uniqueN(D$id), uniqueN(D$season)))

cat("=== the null, by |declination| ===\n")
D[, dbin := cut(abs(decl), breaks = c(0, 5, 10, 15, 20, 24))]
print(as.data.frame(D[!is.na(dbin), .(windows = .N, E2 = round(mean(E2), 4),
  abs_err_lat = round(mean(abs_err_lat), 2)), by = dbin][order(dbin)]), row.names = FALSE)

cat("\n=== does E2 predict the latitude error? ===\n")
ct <- cor.test(D$E2, D$abs_err_lat, method = "spearman")
cat(sprintf("  pooled          rho = %+.3f (p = %.5f, n = %d)\n",
            ct$estimate, ct$p.value, nrow(D)))
D[, `:=`(E2_c = E2 - mean(E2), err_c = abs_err_lat - mean(abs_err_lat)), by = id]
ct2 <- cor.test(D$E2_c, D$err_c, method = "spearman")
cat(sprintf("  within tag      rho = %+.3f (p = %.5f)\n", ct2$estimate, ct2$p.value))
for (s in sort(unique(D$season))) {
  ds <- D[season == s]
  cs <- cor.test(ds$E2, ds$abs_err_lat, method = "spearman")
  cat(sprintf("  season %s     rho = %+.3f (p = %.4f, n = %d)\n",
              s, cs$estimate, cs$p.value, nrow(ds)))
}

cat("\n=== does it add to the calendar? ===\n")
m1 <- lm(abs_err_lat ~ abs(decl), data = D)
m2 <- lm(abs_err_lat ~ abs(decl) + E2, data = D)
cat(sprintf("|decl| alone     R2 = %.3f\n|decl| + E2      R2 = %.3f  (increment %.3f, p = %.4g)\n",
            summary(m1)$r.squared, summary(m2)$r.squared,
            summary(m2)$r.squared - summary(m1)$r.squared, anova(m1, m2)$`Pr(>F)`[2]))

cat("\n=== the practical question: does the model KNOW when it is uninformed? ===\n")
cat("If the reported interval already tracked the information, lat_sd would\n")
cat("correlate with E2 as strongly as the error does. Any gap is what a\n")
cat("harmonic-scaled interval would recover.\n")
c1 <- cor.test(D$E2, D$abs_err_lat, method = "spearman")$estimate
c2 <- cor.test(D$E2, D$lat_sd, method = "spearman")$estimate
cat(sprintf("  E2 vs actual |error| : rho = %+.3f\n", c1))
cat(sprintf("  E2 vs reported sd    : rho = %+.3f\n", c2))
saveRDS(D, file.path(SCRATCH, "harmonics29.rds"))
