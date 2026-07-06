# Probe the exact per-knot light likelihood as a function of latitude
# (lon fixed at the true track), for sim_equinox. Tells us whether the
# per-knot latitude posterior is genuinely bimodal (=> mixture / fallback)
# or unimodal-but-broad (=> Laplace fit to the mode).
suppressPackageStartupMessages(library(invTwilightFree))

e <- new.env(); utils::data("sim_equinox", package = "invTwilightFree", envir = e)
df <- e$sim_equinox; df <- df[order(df$time), ]
start_lon <- df$true_lon[1]; start_lat <- df$true_lat[1]

# replicate calibration + knots (same as main script)
light <- df$light; date_time <- df$time
min_l <- quantile(light,.05,na.rm=TRUE); max_l <- quantile(light,.95,na.rm=TRUE)
lsh <- pmax(0, light-min_l); maxs <- as.numeric(max_l-min_l)
ci <- date_time < (date_time[1]+3*24*3600)
cz <- solar_zenith(as.numeric(date_time[ci]), rep(start_lon,sum(ci)), rep(start_lat,sum(ci)))
ti <- which(lsh[ci]>0 & lsh[ci]<maxs*.95 & cz>85 & cz<100)
fit <- lm(lsh[ci][ti] ~ cz[ti]); icpt <- coef(fit)[1]; slp <- -coef(fit)[2]
if (is.na(slp)||slp<=0) { slp <- maxs/(96-85); icpt <- slp*96 }
calibration <- as.numeric(c(icpt,slp))
likpar <- as.numeric(c(1/(maxs*.5), maxs, .10))
ut <- as.numeric(date_time)
t0 <- ut[1]; t1 <- ut[length(ut)]
ks <- ceiling((t1-t0)/(12*3600))+1; tstep <- (t1-t0)/(ks-1)
knot_times <- t0 + (0:(ks-1))*tstep

lat_grid <- seq(-70, 60, by = 0.5)
probe_knot <- function(k) {
  tc <- knot_times[k]; tp <- if (k==1) tc-(knot_times[2]-knot_times[1]) else knot_times[k-1]
  j <- which(ut>tp & ut<=tc & !is.na(lsh))
  if (!length(j)) return(NULL)
  # lon fixed at the true lon nearest this knot
  lon_k <- df$true_lon[which.min(abs(ut - tc))]
  lat_true <- df$true_lat[which.min(abs(ut - tc))]
  ll <- eval_logpk_grid(rep(lon_k,length(lat_grid)), lat_grid, ut[j], lsh[j], calibration, likpar)
  list(ll = ll, lat_true = lat_true, lon = lon_k)
}

sel <- round(seq(2, ks-1, length.out = 6))
grDevices::png("equinox_ll_profiles.png", width = 1400, height = 460)
par(mfrow = c(2,3), mar = c(4,4,2,1))
for (k in sel) {
  pr <- probe_knot(k); if (is.null(pr)) next
  ll <- pr$ll - max(pr$ll)
  plot(lat_grid, exp(ll), type = "l", lwd = 2,
       xlab = "lat", ylab = "rel. likelihood",
       main = sprintf("knot %d/%d (true lat=%.1f)", k, ks, pr$lat_true))
  abline(v = pr$lat_true, col = "grey50", lwd = 2)
  abline(v = -pr$lat_true, col = "blue", lty = 3)  # mirror hemisphere
  # moment-matched sd over this grid
  w <- exp(ll); w <- w/sum(w); mu <- sum(w*lat_grid); sdv <- sqrt(sum(w*(lat_grid-mu)^2))
  legend("topright", c(sprintf("mm mu=%.1f sd=%.1f", mu, sdv)), bty="n", cex=0.9)
}
grDevices::dev.off()
cat("saved equinox_ll_profiles.png\n")
cat("knots probed:", sel, "of", ks, "\n")
