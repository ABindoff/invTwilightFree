# Day-length diagnostic, done properly.
#
# The first pass thresholded each day at half of THAT day's amplitude, which is
# wrong: a clouded or shaded day never reaches clear-sky maximum, so its own
# half sits below the response's half, the crossing falls at a deeper zenith,
# and a positive bias appears out of nothing. Threshold on the response's
# absolute half level instead, and keep only days bright enough to cross it.
#
# Two response fits are compared: one at the colony, as the analysis does, and
# one along the Argos track over the same window. If the animal has already left
# when the response is fitted, the colony fit is measuring the curve at the wrong
# zenith angles, and that is a calibration-window defect rather than a model one.
suppressMessages({ library(data.table); library(invTwilightFree) })
Sys.setlocale("LC_TIME", "C")
lon360 <- function(x) (x %% 360 + 360) %% 360
dlon <- function(a, b) ((a - b + 180) %% 360) - 180
gc_km <- function(l1, p1, l2, p2, R = 6371.0088) { d <- pi/180
  2*R*asin(pmin(1, sqrt(sin((p2-p1)*d/2)^2 + cos(p1*d)*cos(p2*d)*sin((l2-l1)*d/2)^2))) }

arch <- lapply(readRDS("analysis/cache/nes/archives_v3_30min.rds"), `[[`, "main")
locs <- list.files("data/nes_untracked/argos", pattern = "-Locations[.]csv$",
                   recursive = TRUE, full.names = TRUE)
A <- rbindlist(lapply(locs, fread, showProgress = FALSE,
                      colClasses = list(character = c("Ptt","Quality"))), fill = TRUE)
A[, time := as.POSIXct(Date, format = "%H:%M:%S %d-%b-%Y", tz = "UTC")]
A[, `:=`(lon = lon360(Longitude), lat = Latitude)]
Aq <- A[Type == "Argos" & Quality %in% c("3","2","1","0","A","B") & !is.na(time)]
xw <- c("2021025"="214393","2021027"="214396","2021033"="214399")
interp <- function(a, tm) {
  setorder(a, time); at <- as.numeric(a$time); qt <- as.numeric(tm)
  j <- findInterval(qt, at); ok <- j >= 1 & j < length(at)
  lo <- la <- rep(NA_real_, length(qt)); jj <- j[ok]
  good <- (at[jj+1]-at[jj])/3600 <= 24
  w <- (qt[ok]-at[jj])/(at[jj+1]-at[jj]); idx <- which(ok)[good]
  lo[idx] <- lon360(a$lon[jj[good]] + w[good]*dlon(a$lon[jj[good]+1], a$lon[jj[good]]))
  la[idx] <- a$lat[jj[good]] + w[good]*(a$lat[jj[good]+1]-a$lat[jj[good]])
  list(lon = lo, lat = la)
}
COL_LON <- 237.67; COL_LAT <- 37.11
T0 <- as.POSIXct("2021-07-01", tz = "UTC"); T1 <- as.POSIXct("2021-09-15", tz = "UTC")
LATS <- seq(10, 75, by = 0.1)

cross_time <- function(tt, y, thr, rising) {
  n <- length(y); if (n < 3) return(NA_real_)
  s <- if (rising) which(y[-n] < thr & y[-1] >= thr) else which(y[-n] >= thr & y[-1] < thr)
  if (!length(s)) return(NA_real_)
  i <- if (rising) s[1] else s[length(s)]
  d <- y[i+1] - y[i]; if (!is.finite(d) || d == 0) return(as.numeric(tt[i]))
  as.numeric(tt[i]) + (thr - y[i]) / d * (as.numeric(tt[i+1]) - as.numeric(tt[i]))
}

daylen_bias <- function(d, a, r, tag) {
  floor_l <- r$calibration_logistic[1]; amp <- r$calibration_logistic[2]
  z50 <- r$calibration_logistic[3]
  thr <- floor_l + amp / 2
  d <- copy(d)[, day := as.integer(floor((as.numeric(time) - 8*3600) / 86400))]
  res <- list()
  for (dy in unique(d$day)) {
    w <- d[day == dy]
    if (nrow(w) < 24) next
    if (max(w$light) < floor_l + 0.8 * amp) next     # too dim to cross cleanly
    tr <- interp(a, mean(w$time)); if (!is.finite(tr$lat)) next
    t_up <- cross_time(w$time, w$light, thr, TRUE)
    t_dn <- cross_time(w$time, w$light, thr, FALSE)
    if (!is.finite(t_up) || !is.finite(t_dn) || t_dn <= t_up) next
    z_up <- solar_zenith(t_up, tr$lon, tr$lat)
    z_dn <- solar_zenith(t_dn, tr$lon, tr$lat)
    zu <- solar_zenith(rep(t_up, length(LATS)), rep(tr$lon, length(LATS)), LATS)
    zd <- solar_zenith(rep(t_dn, length(LATS)), rep(tr$lon, length(LATS)), LATS)
    lat_hat <- LATS[which.min((zu - z50)^2 + (zd - z50)^2)]
    res[[length(res)+1]] <- data.frame(id = tag, z50 = z50, z_obs = (z_up+z_dn)/2,
      lat_true = tr$lat, lat_implied = lat_hat, daylen_h = (t_dn-t_up)/3600,
      row.names = NULL)
  }
  if (!length(res)) return(NULL)
  as.data.table(do.call(rbind, res))
}

cal_days <- 15
out <- list(); dist_tab <- list()
for (id in names(xw)) {
  x <- as.data.table(arch[[id]])
  a <- Aq[Ptt == xw[[id]]]
  k <- x$time <= min(x$time) + cal_days*86400
  tc <- interp(a, x$time[k])
  dkm <- gc_km(tc$lon, tc$lat, COL_LON, COL_LAT)
  dist_tab[[id]] <- data.frame(id = id, n_cal = sum(k),
    argos_cover = round(mean(is.finite(dkm)), 2),
    med_km_from_colony = round(median(dkm, na.rm = TRUE)),
    max_km_from_colony = round(max(dkm, na.rm = TRUE)), row.names = NULL)

  r_col <- fit_light_response(x$time[k], x$light[k], COL_LON, COL_LAT)
  gg <- is.finite(tc$lon) & is.finite(tc$lat)
  r_trk <- if (sum(gg) >= 200)
    fit_light_response(x$time[k][gg], x$light[k][gg], tc$lon[gg], tc$lat[gg]) else NULL

  d <- x[time >= T0 & time <= T1]
  for (nm in c("colony", "argos_track")) {
    r <- if (nm == "colony") r_col else r_trk
    if (is.null(r)) next
    o <- daylen_bias(d, a, r, id)
    if (is.null(o)) next
    o[, source := nm]
    o[, `:=`(dz = z_obs - z50, dlat = lat_implied - lat_true)]
    out[[length(out)+1]] <- o
  }
}
cat("=== is the animal still at the colony during the calibration window? ===\n")
print(do.call(rbind, dist_tab), row.names = FALSE)

O <- rbindlist(out)
cat("\n=== day-length bias, thresholding on the response's own half level ===\n")
print(as.data.frame(O[, .(days = .N, z50 = round(mean(z50), 1),
      z_obs = round(mean(z_obs), 1), dz = round(mean(dz), 2),
      lat_true = round(mean(lat_true), 1), lat_implied = round(mean(lat_implied), 1),
      dlat = round(mean(dlat), 2), dlat_sd = round(sd(dlat), 2)),
      by = .(source, id)][order(source, id)]), row.names = FALSE)
cat("\n=== pooled by response fit ===\n")
print(as.data.frame(O[, .(days = .N, dz = round(mean(dz), 2),
      dlat = round(mean(dlat), 2), dlat_sd = round(sd(dlat), 2)), by = source]),
      row.names = FALSE)
