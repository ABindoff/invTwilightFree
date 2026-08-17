# Figures for eye-judging the battery: real track, perfect light, degraded light.
# Written to a local folder and assembled into a local HTML file -- NOT published,
# because these are tracks and light records from the restricted deliveries.
suppressMessages({ library(data.table) })
# run from the package root (setwd removed for portability)
suppressMessages(devtools::load_all(".", quiet = TRUE))
BAT <- readRDS("scratch/nes_calibration/light_battery.rds")
DIR <- "scratch/nes_calibration/battery_fig"
dir.create(DIR, showWarnings = FALSE, recursive = TRUE)

ids <- unique(vapply(BAT, function(b) b$id, ""))
offs <- sort(unique(vapply(BAT, function(b) b$offset, 0)))
get <- function(id, off) BAT[[which(vapply(BAT, function(b) b$id == id && b$offset == off, TRUE))[1]]]
COL_P <- "#1f6feb"; COL_D <- "#d1242f"

# ---- 1. all tracks on one map ----------------------------------------------
png(file.path(DIR, "tracks.png"), width = 1500, height = 820, res = 130)
par(mar = c(4, 4, 2, 1))
b1 <- lapply(ids, function(i) get(i, 0))
xr <- range(unlist(lapply(b1, function(b) b$lon)))
yr <- range(unlist(lapply(b1, function(b) b$lat)))
plot(NA, xlim = xr, ylim = yr, xlab = "longitude (deg E)", ylab = "latitude (deg N)",
     main = "Argos tracks used for the battery (real movement, subsampled 1:20)")
grid(col = "grey88")
pal <- hcl.colors(length(ids), "Dark 3")
for (k in seq_along(b1)) {
  b <- b1[[k]]; s <- seq(1, length(b$lon), by = 20)
  lines(b$lon[s], b$lat[s], col = pal[k], lwd = 1.6)
  points(b$lon[1], b$lat[1], pch = 19, col = pal[k], cex = 1.1)
}
legend("topright", legend = ids, col = pal, lwd = 2, bty = "n", cex = 0.8)
dev.off()

# ---- 2. per track: perfect vs degraded, full record, all offsets ------------
for (id in ids) {
  png(file.path(DIR, sprintf("light_%s.png", id)), width = 1500, height = 1100, res = 130)
  par(mfrow = c(length(offs), 1), mar = c(2.4, 4, 2, 1), oma = c(2, 0, 2, 0))
  for (o in offs) {
    b <- get(id, o)
    s <- seq(1, length(b$time), by = 3)     # thin for legibility only
    plot(b$time[s], b$perfect[s], type = "l", col = COL_P, lwd = 0.5,
         ylim = c(0, b$max_light), xlab = "", ylab = "light",
         main = sprintf("%s  |  date offset %+d days  |  %s to %s", id, o,
                        format(min(b$time), "%d %b %Y"), format(max(b$time), "%d %b %Y")),
         cex.main = 0.9)
    lines(b$time[s], b$degraded[s], col = adjustcolor(COL_D, 0.6), lwd = 0.4)
    dec <- abs(solar_declination(b$time[s]))
    eq <- which(diff(sign(diff(dec))) > 0) + 1     # local minima of |declination|
    abline(v = b$time[s][eq], col = "grey45", lty = 3)
  }
  mtext(sprintf("%s: clear-sky (blue) vs shaded (red). Dotted = equinox.", id),
        outer = TRUE, cex = 0.95)
  dev.off()
}

# ---- 3. diel zoom: equinox vs solstice, one track --------------------------
for (id in ids[1:2]) {
  b <- get(id, 0)
  dec <- abs(solar_declination(b$time))
  i_eq <- which.min(dec); i_so <- which.max(dec)
  png(file.path(DIR, sprintf("zoom_%s.png", id)), width = 1500, height = 700, res = 130)
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  for (nm in c("near equinox", "near solstice")) {
    i0 <- if (nm == "near equinox") i_eq else i_so
    w <- which(abs(as.numeric(b$time) - as.numeric(b$time[i0])) < 3 * 86400)
    plot(b$time[w], b$perfect[w], type = "l", col = COL_P, lwd = 1.4,
         ylim = c(0, b$max_light), xlab = "", ylab = "light",
         main = sprintf("%s, %s (|dec| %.1f, lat %.1f)", id, nm,
                        dec[i0], b$lat[i0]), cex.main = 0.95)
    lines(b$time[w], b$degraded[w], col = COL_D, lwd = 1.2)
    legend("topleft", c("clear sky", "shaded"), col = c(COL_P, COL_D),
           lwd = 2, bty = "n", cex = 0.8)
  }
  dev.off()
}

# ---- 4. the measured attenuation itself -------------------------------------
png(file.path(DIR, "attenuation.png"), width = 1500, height = 620, res = 130)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
sl <- lapply(ids, function(i) { a <- get(i, 0)$atten; a[is.finite(a)] })
hist(unlist(sl), breaks = 60, col = "grey75", border = "white",
     xlab = "attenuation factor (observed sky / clear sky)", main = "Measured attenuation, pooled")
abline(v = 1, col = COL_P, lwd = 2)
boxplot(sl, names = ids, las = 2, col = "grey85", ylab = "attenuation factor",
        main = "Per deployment", cex.axis = 0.7)
abline(h = 1, col = COL_P, lwd = 2)
dev.off()

cat("figures written to", DIR, "\n")
print(list.files(DIR))

# summary table for the document
S <- rbindlist(lapply(BAT, function(b) data.table(
  id = b$id, offset = b$offset, n = length(b$time),
  start = format(min(b$time), "%Y-%m-%d"), days = round(as.numeric(diff(range(b$time)), units = "days")),
  lat_min = round(min(b$lat), 1), lat_max = round(max(b$lat), 1),
  lon_min = round(min(b$lon), 1), lon_max = round(max(b$lon), 1),
  atten_med = round(median(b$atten, na.rm = TRUE), 2),
  atten_zero = round(mean(b$atten < 0.02, na.rm = TRUE), 3),
  perfect_med = round(median(b$perfect), 1), degraded_med = round(median(b$degraded), 1))))
fwrite(S, "scratch/nes_calibration/battery_summary.csv")
print(as.data.frame(S), row.names = FALSE)
