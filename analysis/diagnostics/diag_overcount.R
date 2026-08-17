# Composite-likelihood test. The prediction is sharp and falsifiable:
#   posterior sd scales as sqrt(overcount)
#   the point estimate does NOT move
#   coverage rises toward 0.95 near overcount = 5
# If the point estimate moves, the adjustment is doing something other than
# correcting for over-counted information and should be rejected.
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

res <- list()
for (id in names(xw)) {
  x <- as.data.table(arch[[id]])
  k <- x$time <= min(x$time) + 15*86400
  r <- fit_light_response(x$time[k], x$light[k], COL_LON, COL_LAT)
  d <- x[time >= T0 & time <= T1]
  a <- Aq[Ptt == xw[[id]]]
  t0 <- interp(a, d$time[1]); t1 <- interp(a, d$time[nrow(d)])
  for (oc in c(1, 2, 4, 6, 9)) {
    invisible(capture.output(
      f <- TwilightFreeGrid(d$time, pmax(0, d$light - r$baseline), grid = grid,
        start_lon = t0$lon, start_lat = t0$lat, end_lon = t1$lon, end_lat = t1$lat,
        step_hours = 12, diffusion = 250, calibration = r$calibration,
        likelihood_params = c(1/(r$max_light*0.5), r$max_light, 0.10),
        overcount = oc)))
    gp <- grid_posterior(f)
    sdl <- sdp <- numeric(nrow(gp$P))
    for (i in seq_len(nrow(gp$P))) {
      w <- gp$P[i, ]; sm <- sum(w)
      if (!is.finite(sm) || sm <= 0) { sdl[i] <- sdp[i] <- NA; next }
      w <- w/sm
      ml <- sum(w*gp$lon); mp <- sum(w*gp$lat)
      sdl[i] <- sqrt(sum(w*(gp$lon-ml)^2)); sdp[i] <- sqrt(sum(w*(gp$lat-mp)^2))
    }
    tm <- as.POSIXct(f$fit$time, origin = "1970-01-01", tz = "UTC")
    tr <- interp(a, tm)
    el <- dlon(lon360(f$fit$lon), tr$lon); ep <- f$fit$lat - tr$lat
    e <- gc_km(lon360(f$fit$lon), f$fit$lat, tr$lon, tr$lat)
    g <- is.finite(ep) & is.finite(sdp) & sdp > 0
    res[[length(res)+1]] <- data.frame(
      id = id, overcount = oc, median_km = round(median(e[g])),
      post_sd_lat = round(median(sdp[g]), 2), err_sd_lat = round(sd(ep[g]), 2),
      cover_lat = round(mean(abs(ep[g]) <= 1.96*sdp[g]), 2),
      cover_lon = round(mean(abs(el[g]) <= 1.96*sdl[g]), 2),
      bias_lat = round(mean(ep[g]), 2), row.names = NULL)
  }
}
out <- do.call(rbind, res)
cat("=== composite-likelihood weight, grid HMM exact posterior ===\n")
print(as.data.frame(as.data.table(out)[, .(
  median_km = round(mean(median_km)), bias_lat = round(mean(bias_lat), 2),
  post_sd_lat = round(mean(post_sd_lat), 2), err_sd_lat = round(mean(err_sd_lat), 2),
  cover_lat = round(mean(cover_lat), 2), cover_lon = round(mean(cover_lon), 2)),
  by = overcount][order(overcount)]), row.names = FALSE)
b <- as.data.table(out)[, .(s = mean(post_sd_lat)), by = overcount][order(overcount)]
cat("\nposterior sd relative to overcount = 1, against the sqrt(overcount) prediction:\n")
for (i in seq_len(nrow(b)))
  cat(sprintf("  overcount %2.0f : observed %.2fx, predicted %.2fx\n",
              b$overcount[i], b$s[i]/b$s[1], sqrt(b$overcount[i])))
