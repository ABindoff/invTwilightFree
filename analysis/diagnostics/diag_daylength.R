# Where does the residual latitude bias come from? Latitude is read off day
# length, so isolate day length and drop the HMM entirely.
#
# For each day: find the zenith at which the observed light actually crosses
# half of its daily amplitude, at the KNOWN Argos position. The response model
# says that crossing happens at z50. Any systematic gap between the two is a
# response-shape error, and it converts straight into a latitude bias -- which
# is measured here by asking which latitude WOULD have put the crossings at z50.
# No sampler, no prior, no grid, so nothing else can be blamed.
suppressMessages({ library(data.table); library(invTwilightFree) })
Sys.setlocale("LC_TIME", "C")
lon360 <- function(x) (x %% 360 + 360) %% 360
dlon <- function(a, b) ((a - b + 180) %% 360) - 180

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
LATS <- seq(15, 70, by = 0.1)

# Linear interpolation of the time at which a series crosses `thr`.
cross_time <- function(tt, y, thr, rising) {
  n <- length(y); if (n < 3) return(NA_real_)
  s <- if (rising) which(y[-n] < thr & y[-1] >= thr) else which(y[-n] >= thr & y[-1] < thr)
  if (!length(s)) return(NA_real_)
  i <- if (rising) s[1] else s[length(s)]
  d <- y[i+1] - y[i]; if (!is.finite(d) || d == 0) return(as.numeric(tt[i]))
  as.numeric(tt[i]) + (thr - y[i]) / d * (as.numeric(tt[i+1]) - as.numeric(tt[i]))
}

res <- list()
for (id in names(xw)) {
  x <- as.data.table(arch[[id]])
  k <- x$time <= min(x$time) + 15*86400
  r <- fit_light_response(x$time[k], x$light[k], COL_LON, COL_LAT)
  z50 <- r$calibration_logistic[3]
  d <- x[time >= T0 & time <= T1]
  a <- Aq[Ptt == xw[[id]]]
  # Days are centred on the local solar day so a whole daylight arc sits in one
  # window; at ~150-220 W the offset from UTC is roughly 10-15 hours.
  d[, day := as.integer(floor((as.numeric(time) - 8*3600) / 86400))]
  for (dy in unique(d$day)) {
    w <- d[day == dy]
    if (nrow(w) < 24) next
    tr <- interp(a, mean(w$time))
    if (!is.finite(tr$lat)) next
    lo <- min(w$light); hi <- max(w$light)
    if (hi - lo < 20) next                        # no usable diurnal contrast
    thr <- lo + 0.5 * (hi - lo)
    t_up <- cross_time(w$time, w$light, thr, TRUE)
    t_dn <- cross_time(w$time, w$light, thr, FALSE)
    if (!is.finite(t_up) || !is.finite(t_dn) || t_dn <= t_up) next
    # The zenith the model expects at a crossing is z50; the zenith actually
    # obtaining at the known position is this:
    z_up <- solar_zenith(t_up, tr$lon, tr$lat)
    z_dn <- solar_zenith(t_dn, tr$lon, tr$lat)
    # Which latitude would have placed both crossings at z50?
    zu <- solar_zenith(rep(t_up, length(LATS)), rep(tr$lon, length(LATS)), LATS)
    zd <- solar_zenith(rep(t_dn, length(LATS)), rep(tr$lon, length(LATS)), LATS)
    ss <- (zu - z50)^2 + (zd - z50)^2
    lat_hat <- LATS[which.min(ss)]
    res[[length(res)+1]] <- data.frame(
      id = id, z50 = z50, lat_true = tr$lat, lat_implied = lat_hat,
      z_obs = (z_up + z_dn)/2, daylen_h = (t_dn - t_up)/3600, row.names = NULL)
  }
}
o <- as.data.table(do.call(rbind, res))
o[, `:=`(dz = z_obs - z50, dlat = lat_implied - lat_true)]
cat("=== half-amplitude crossing: model z50 vs the zenith actually obtaining ===\n")
print(as.data.frame(o[, .(days = .N, z50 = round(mean(z50), 1),
      z_obs = round(mean(z_obs), 1), dz = round(mean(dz), 2),
      lat_true = round(mean(lat_true), 1), lat_implied = round(mean(lat_implied), 1),
      dlat = round(mean(dlat), 2), dlat_sd = round(sd(dlat), 2)), by = id]),
      row.names = FALSE)
cat("\npooled: dz = ", round(mean(o$dz), 2), " deg zenith, dlat = ",
    round(mean(o$dlat), 2), " deg latitude (n = ", nrow(o), " days)\n", sep = "")
cat("\n=== does the bias track latitude? (binned by true latitude) ===\n")
o[, bin := cut(lat_true, breaks = seq(25, 60, by = 5))]
print(as.data.frame(o[!is.na(bin), .(days = .N, dz = round(mean(dz), 2),
      dlat = round(mean(dlat), 2), daylen_h = round(mean(daylen_h), 1)),
      by = bin][order(bin)]), row.names = FALSE)
saveRDS(o, file.path(tempdir(), "daylen.rds"))
