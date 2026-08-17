# The light scale is set from quantiles of the calibration window:
# baseline = q05, max_light = q95 - q05. That makes it depend on how much
# DIVING happened to fall in the window, which is why the haul-out window (no
# dives, so a high q05) wrecked tag 2021027 while the 15-day window was fine.
# It is an accident of behaviour standing in for a property of the sensor.
#
# The envelope already measures the same two quantities physically: `floor` is
# the dark reading at high zenith and `amp` the clear-sky amplitude above it.
# The code prefers the quantiles, with a comment recording that the envelope
# roughly doubled the error -- but that comparison was made with the colony
# calibration, which we now know is broken for two of these three animals. It is
# worth re-running against a calibration that is right.
suppressMessages({ library(data.table); library(invTwilightFree); library(terra) })
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
grid <- rast(xmin = 164, xmax = 244, ymin = 20, ymax = 66, resolution = 2, crs = "EPSG:4326")
values(grid) <- 1
T0 <- as.POSIXct("2021-07-01", tz = "UTC"); T1 <- as.POSIXct("2021-09-15", tz = "UTC")

departure_day <- function(d, dive_m = 50) {
  dd <- copy(d)[, day := as.integer(floor(as.numeric(time)/86400))]
  s <- dd[, .(deep = mean(depth_max > dive_m, na.rm = TRUE)), by = day][order(day)]
  at_sea <- s$deep > 0.5
  if (!any(at_sea)) return(NA_real_)
  last_ashore <- if (any(!at_sea)) max(which(!at_sea)) else 0L
  i <- if (last_ashore >= length(at_sea)) which(at_sea)[1] else last_ashore + 1L
  s$day[i] * 86400
}

# Take the light scale from the fitted envelope rather than window quantiles.
envelope_scale <- function(r) {
  if (is.null(r)) return(NULL)
  r$baseline <- r$floor
  r$max_light <- r$amp
  slope <- r$amp / (4 * r$scale); zero <- r$z50 + 2 * r$scale
  r$calibration <- c(slope * zero, slope)
  r
}

score <- function(id, r, label) {
  if (is.null(r)) return(NULL)
  x <- as.data.table(arch[[id]]); a <- Aq[Ptt == xw[[id]]]
  d <- x[time >= T0 & time <= T1]
  t0 <- interp(a, d$time[1]); t1 <- interp(a, d$time[nrow(d)])
  invisible(capture.output(
    f <- TwilightFreeGrid(d$time, pmax(0, d$light - r$baseline), grid = grid,
      start_lon = t0$lon, start_lat = t0$lat, end_lon = t1$lon, end_lat = t1$lat,
      step_hours = 12, diffusion = 250, calibration = r$calibration,
      likelihood_params = c(1/(r$max_light*0.5), r$max_light, 0.10))))
  gp <- grid_posterior(f)
  sdp <- vapply(seq_len(nrow(gp$P)), function(i) { w <- gp$P[i, ]; s <- sum(w)
    if (!is.finite(s) || s <= 0) return(NA_real_); w <- w/s
    sqrt(sum(w*(gp$lat - sum(w*gp$lat))^2)) }, 0)
  tm <- as.POSIXct(f$fit$time, origin = "1970-01-01", tz = "UTC")
  tr <- interp(a, tm)
  ep <- f$fit$lat - tr$lat
  e <- gc_km(lon360(f$fit$lon), f$fit$lat, tr$lon, tr$lat)
  g <- is.finite(ep) & is.finite(sdp) & sdp > 0
  data.frame(id = id, recipe = label,
    baseline = round(r$baseline, 1), max_light = round(r$max_light, 1),
    median_km = round(median(e[g])), bias_lat = round(mean(ep[g]), 2),
    rmse_lat = round(sqrt(mean(ep[g]^2)), 2), post_sd_lat = round(median(sdp[g]), 2),
    cover_lat = round(mean(abs(ep[g]) <= 1.96*sdp[g]), 2), row.names = NULL)
}

res <- list()
for (id in names(xw)) {
  x <- as.data.table(arch[[id]]); a <- Aq[Ptt == xw[[id]]]
  k15 <- x$time <= min(x$time) + 15*86400
  kpre <- as.numeric(x$time) < departure_day(x)
  tc <- interp(a, x$time[k15]); gg <- is.finite(tc$lon)

  r15  <- fit_light_response(x$time[k15], x$light[k15], COL_LON, COL_LAT)
  rpre <- if (sum(kpre) >= 200)
    fit_light_response(x$time[kpre], x$light[kpre], COL_LON, COL_LAT) else NULL
  rarg <- if (sum(gg) >= 200)
    fit_light_response(x$time[k15][gg], x$light[k15][gg], tc$lon[gg], tc$lat[gg]) else NULL

  res[[length(res)+1]] <- score(id, r15,  "A 15d colony, quantile scale")
  res[[length(res)+1]] <- score(id, envelope_scale(r15),  "J 15d colony, envelope scale")
  res[[length(res)+1]] <- score(id, rpre, "B haulout, quantile scale")
  res[[length(res)+1]] <- score(id, envelope_scale(rpre), "K haulout, envelope scale")
  res[[length(res)+1]] <- score(id, rarg, "E argos, quantile scale")
  res[[length(res)+1]] <- score(id, envelope_scale(rarg), "I argos, envelope scale")
}
o <- as.data.table(do.call(rbind, res))
cat("=== grid fits, July-September ===\n")
print(as.data.frame(o[order(recipe, id)]), row.names = FALSE)
cat("\n=== pooled ===\n")
print(as.data.frame(o[, .(median_km = round(mean(median_km)), worst_km = max(median_km),
      bias_lat = round(mean(bias_lat), 2), rmse_lat = round(mean(rmse_lat), 2),
      post_sd_lat = round(mean(post_sd_lat), 2), cover_lat = round(mean(cover_lat), 2)),
      by = recipe][order(recipe)]), row.names = FALSE)
