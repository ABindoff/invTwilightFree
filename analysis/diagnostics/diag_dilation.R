# Why does a looped track come back DILATED rather than displaced?
#
# Two hypotheses, and they are separable because they are organised by different
# variables.
#
#   OBSERVATION MODEL. Latitude is read from day length. A fixed error in the
#   assumed threshold zenith enters as a fixed perturbation to cos(z), and the
#   latitude error is that divided by
#       f'(phi) = cos(phi) sin(delta) - sin(phi) cos(delta) cos(H)
#   which changes SIGN with declination: positive in northern summer (delta > 0,
#   H > 90) and negative in northern winter. So one calibration error produces
#   opposite-signed latitude bias on the two sides of the equinox. A post-moult
#   trip leaves in summer and returns in winter, so the loop inflates.
#   PREDICTION: bias is organised by DATE, and flips sign at the equinox,
#   independently of where the animal is or how fast it is moving.
#
#   MOVEMENT MODEL. A single diffusion cannot serve both transit and area-
#   restricted search, so the prior over-smooths one and under-smooths the other.
#   PREDICTION: bias is organised by SPEED, independently of date.
#
# The clean discriminator is the crossing itself: pairs of knots close in SPACE
# but far apart in TIME. Same place, different season, and different behavioural
# mode too. If the observation model is responsible, the latitude errors of such
# a pair differ in proportion to their declination difference.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
suppressMessages({ library(data.table) })
gc_km <- function(l1, p1, l2, p2, R = 6371.0088) { d <- pi/180
  2*R*asin(pmin(1, sqrt(sin((p2-p1)*d/2)^2 + cos(p1*d)*cos(p2*d)*sin((l2-l1)*d/2)^2))) }

# Spencer's series for solar declination, accurate to about 0.1 degrees.
declination <- function(t) {
  doy <- as.numeric(format(t, "%j"))
  hr  <- as.numeric(format(t, "%H")) + as.numeric(format(t, "%M")) / 60
  g <- 2*pi/365 * (doy - 1 + (hr - 12)/24)
  (0.006918 - 0.399912*cos(g) + 0.070257*sin(g) - 0.006758*cos(2*g) +
     0.000907*sin(2*g) - 0.002697*cos(3*g) + 0.00148*sin(3*g)) * 180/pi
}

k <- as.data.table(read.csv("analysis/output/nes/per_knot_errors.csv"))
k[, time := as.POSIXct(time, tz = "UTC")]
k <- k[engine == "hier" & batch == "light" & is.finite(err_lat) & is.finite(true_lat)]
k[, decl := declination(time)]
setorder(k, id, time)

cat("=== 1. latitude bias against solar declination ===\n")
k[, dbin := cut(decl, breaks = c(-24, -18, -12, -6, 0, 6, 12, 18, 24))]
print(as.data.frame(k[!is.na(dbin), .(knots = .N,
  mean_decl = round(mean(decl), 1), mean_true_lat = round(mean(true_lat), 1),
  bias_lat = round(mean(err_lat), 2), sd_lat = round(sd(err_lat), 2)),
  by = dbin][order(dbin)]), row.names = FALSE)
cat(sprintf("\ncorrelation of latitude error with declination: %.3f (n = %d)\n",
            cor(k$decl, k$err_lat), nrow(k)))

# Is it really declination, or just that the animal is further north in summer?
cat("\n=== 2. is it date or is it place? (bias by declination WITHIN latitude bands)\n")
k[, lbin := cut(true_lat, breaks = c(30, 40, 45, 50, 60))]
tab <- k[!is.na(lbin) & !is.na(dbin), .(n = .N, bias = round(mean(err_lat), 2)),
         by = .(lbin, sum_win = ifelse(decl > 0, "summer", "winter"))]
print(as.data.frame(dcast(tab, lbin ~ sum_win, value.var = c("n", "bias"))),
      row.names = FALSE)

# Partial correlations: which variable survives controlling for the other?
m1 <- lm(err_lat ~ decl + true_lat, data = k)
cat("\nerr_lat ~ declination + true latitude:\n")
print(round(summary(m1)$coefficients, 4))

cat("\n=== 3. the crossing test: same place, different season ===\n")
# pairs within 150 km in space but more than 60 days apart in time
pairs <- rbindlist(lapply(unique(k$id), function(tg) {
  d <- k[id == tg]
  n <- nrow(d); if (n < 20) return(NULL)
  ii <- CJ(i = seq_len(n), j = seq_len(n))[i < j]
  ii[, km := gc_km(d$true_lon[i], d$true_lat[i], d$true_lon[j], d$true_lat[j])]
  ii[, dt := abs(as.numeric(difftime(d$time[j], d$time[i], units = "days")))]
  ii <- ii[km < 150 & dt > 60]
  if (!nrow(ii)) return(NULL)
  data.table(id = tg, km = ii$km, dt = ii$dt,
             d_decl = d$decl[ii$j] - d$decl[ii$i],
             d_err  = d$err_lat[ii$j] - d$err_lat[ii$i])
}))
if (nrow(pairs)) {
  cat(sprintf("%d self-crossing pairs across %d tags (within 150 km, over 60 days apart)\n",
              nrow(pairs), length(unique(pairs$id))))
  cat(sprintf("correlation of the latitude-error difference with the declination difference: %.3f\n",
              cor(pairs$d_decl, pairs$d_err)))
  fit <- lm(d_err ~ d_decl, data = pairs)
  cat(sprintf("slope %.4f deg latitude error per degree of declination (p = %.3g)\n",
              coef(fit)[2], summary(fit)$coefficients[2, 4]))
  pairs[, db := cut(d_decl, breaks = c(-50, -30, -15, 0, 15, 30, 50))]
  print(as.data.frame(pairs[!is.na(db), .(pairs = .N,
    mean_d_decl = round(mean(d_decl), 1), mean_d_err = round(mean(d_err), 2)),
    by = db][order(db)]), row.names = FALSE)
} else cat("no qualifying pairs\n")

cat("\n=== 4. the movement hypothesis: is bias organised by SPEED? ===\n")
k[, spd := {
  v <- c(NA_real_, gc_km(head(true_lon, -1), head(true_lat, -1),
                         tail(true_lon, -1), tail(true_lat, -1)) /
           as.numeric(difftime(tail(time, -1), head(time, -1), units = "hours")))
  v }, by = id]
k[, sbin := cut(spd, breaks = c(0, 1, 2, 3, 4, 20))]
print(as.data.frame(k[!is.na(sbin), .(knots = .N, mean_kmh = round(mean(spd), 2),
  bias_lat = round(mean(err_lat), 2), abs_err = round(mean(abs(err_lat)), 2)),
  by = sbin][order(sbin)]), row.names = FALSE)
m2 <- lm(err_lat ~ decl + spd, data = k[is.finite(spd)])
cat("\nerr_lat ~ declination + speed:\n")
print(round(summary(m2)$coefficients, 4))
saveRDS(list(k = k, pairs = pairs), file.path(SCRATCH, "dilation.rds"))
