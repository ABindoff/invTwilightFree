# Splitting geometry from intensity failed: every recipe that took the scale
# from the whole record drifted north (deep dives drag the 5th percentile down,
# the light looks brighter, the days look longer). Both halves of the response
# have to come from the same window. The only thing wrong with that window is
# the ASSUMED POSITION -- the animal is not at the colony.
#
# Argos supplies the true positions and takes the fit from 321 km to 221 km with
# coverage 0.99, but that is ground truth and cannot ship. Every deployment does
# know one thing though: where the tag went on the animal, at t = 0. So run a
# first pass anchored there, read the positions it gives over the calibration
# window, and refit the response against those.
#
# This is refine_light_response() confined to the calibration window. Refining
# over the WHOLE record already failed (421 -> 1178 km) because latitude error
# feeds back into the envelope. Over the first few weeks the animal is still
# near a known anchor, so the feedback has much less room to run.
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

fit_track <- function(tm, lig, r, s_lon, s_lat, e_lon = NA_real_, e_lat = NA_real_) {
  invisible(capture.output(
    f <- TwilightFreeGrid(tm, pmax(0, lig - r$baseline), grid = grid,
      start_lon = s_lon, start_lat = s_lat, end_lon = e_lon, end_lat = e_lat,
      step_hours = 12, diffusion = 250, calibration = r$calibration,
      likelihood_params = c(1/(r$max_light*0.5), r$max_light, 0.10))))
  f
}

score <- function(id, r, label) {
  if (is.null(r)) return(NULL)
  x <- as.data.table(arch[[id]]); a <- Aq[Ptt == xw[[id]]]
  d <- x[time >= T0 & time <= T1]
  t0 <- interp(a, d$time[1]); t1 <- interp(a, d$time[nrow(d)])
  f <- fit_track(d$time, d$light, r, t0$lon, t0$lat, t1$lon, t1$lat)
  gp <- grid_posterior(f)
  sdp <- vapply(seq_len(nrow(gp$P)), function(i) { w <- gp$P[i, ]; s <- sum(w)
    if (!is.finite(s) || s <= 0) return(NA_real_); w <- w/s
    sqrt(sum(w*(gp$lat - sum(w*gp$lat))^2)) }, 0)
  tm <- as.POSIXct(f$fit$time, origin = "1970-01-01", tz = "UTC")
  tr <- interp(a, tm)
  ep <- f$fit$lat - tr$lat
  e <- gc_km(lon360(f$fit$lon), f$fit$lat, tr$lon, tr$lat)
  g <- is.finite(ep) & is.finite(sdp) & sdp > 0
  data.frame(id = id, recipe = label, z50 = round(r$z50, 1),
    width = round(4.394*r$scale, 1), median_km = round(median(e[g])),
    bias_lat = round(mean(ep[g]), 2), rmse_lat = round(sqrt(mean(ep[g]^2)), 2),
    post_sd_lat = round(median(sdp[g]), 2),
    cover_lat = round(mean(abs(ep[g]) <= 1.96*sdp[g]), 2), row.names = NULL)
}

res <- list(); cal <- list()
for (id in names(xw)) {
  x <- as.data.table(arch[[id]]); a <- Aq[Ptt == xw[[id]]]
  k15 <- x$time <= min(x$time) + 15*86400
  tc <- interp(a, x$time[k15]); gg <- is.finite(tc$lon)

  r15  <- fit_light_response(x$time[k15], x$light[k15], COL_LON, COL_LAT)
  rarg <- if (sum(gg) >= 200)
    fit_light_response(x$time[k15][gg], x$light[k15][gg], tc$lon[gg], tc$lat[gg]) else NULL

  # first pass over the opening 45 days, anchored only at the release site
  k45 <- x$time <= min(x$time) + 45*86400
  fp <- fit_track(x$time[k45], x$light[k45], r15, COL_LON, COL_LAT)
  trk <- data.frame(time = as.POSIXct(fp$fit$time, origin = "1970-01-01", tz = "UTC"),
                    lon = lon360(fp$fit$lon), lat = fp$fit$lat)
  rboot <- refine_light_response(x$time[k15], x$light[k15], trk)
  # a second turn of the same crank, to see whether it settles or runs away
  fp2 <- fit_track(x$time[k45], x$light[k45], rboot, COL_LON, COL_LAT)
  trk2 <- data.frame(time = as.POSIXct(fp2$fit$time, origin = "1970-01-01", tz = "UTC"),
                     lon = lon360(fp2$fit$lon), lat = fp2$fit$lat)
  rboot2 <- refine_light_response(x$time[k15], x$light[k15], trk2)

  # how far off were the first-pass positions over the calibration window?
  ti <- interp(a, trk$time)
  kk <- trk$time <= min(x$time) + 15*86400
  cal[[length(cal)+1]] <- data.frame(id = id,
    firstpass_km = round(median(gc_km(trk$lon[kk], trk$lat[kk], ti$lon[kk], ti$lat[kk]),
                                na.rm = TRUE)),
    colony_km = round(median(gc_km(COL_LON, COL_LAT, ti$lon[kk], ti$lat[kk]), na.rm = TRUE)),
    z50_colony = round(r15$z50, 1), z50_boot = round(rboot$z50, 1),
    z50_boot2 = round(rboot2$z50, 1),
    z50_argos = if (is.null(rarg)) NA else round(rarg$z50, 1),
    w_colony = round(r15$width_deg, 1), w_boot = round(rboot$width_deg, 1),
    w_argos = if (is.null(rarg)) NA else round(rarg$width_deg, 1), row.names = NULL)

  res[[length(res)+1]] <- score(id, r15,    "A 15d colony (current)")
  res[[length(res)+1]] <- score(id, rboot,  "G first-pass refit x1")
  res[[length(res)+1]] <- score(id, rboot2, "H first-pass refit x2")
  res[[length(res)+1]] <- score(id, rarg,   "E argos (ceiling)")
}
cat("=== how wrong is the assumed position over the calibration window? ===\n")
print(do.call(rbind, cal), row.names = FALSE)
o <- as.data.table(do.call(rbind, res))
cat("\n=== grid fits, July-September ===\n")
print(as.data.frame(o[order(recipe, id)]), row.names = FALSE)
cat("\n=== pooled ===\n")
print(as.data.frame(o[, .(median_km = round(mean(median_km)), worst_km = max(median_km),
      bias_lat = round(mean(bias_lat), 2), rmse_lat = round(mean(rmse_lat), 2),
      post_sd_lat = round(mean(post_sd_lat), 2), cover_lat = round(mean(cover_lat), 2)),
      by = recipe][order(recipe)]), row.names = FALSE)
