# =============================================================================
# FIGURE 1: reconstructed tracks for the three seal noise scenarios of Table 2.
#
# Reads the first-repetition fits cached by inst/paper/rerun_benchmarks.R, so
# the figure and the table are drawn from the same run rather than from
# separately cached fits made months apart.
#
#   R-4.6.0/bin/Rscript inst/paper/fig_benchmark_tracks.R
#
# Longitudes are drawn on 0-360. SGAT's excursions cross the dateline, and on
# -180..180 the wrap draws a false line straight across the map. The panels
# share one extent: SGAT leaves the region entirely on the two noisy scenarios,
# and a per-panel extent would quietly rescale that away.
# =============================================================================

FITS <- "inst/paper/fig1_fits.rds"
CSV  <- "inst/paper/benchmark_rerun.csv"
OUT  <- "inst/paper/fig_benchmark_tracks.png"
if (!file.exists(FITS)) stop("run inst/paper/rerun_benchmarks.R first: ", FITS, " missing")
fits <- readRDS(FITS)

SCEN  <- c("Cloudy", "Shaded (ARS diving)", "ALAN near colony")
PANEL <- c("(a) cloud", "(b) deep shading (ARS diving)", "(c) ALAN near colony")
METH  <- list(truth   = list(lab = "true path",       col = "grey55",  lwd = 4.0),
              invTF   = list(lab = "invTwilightFree", col = "#0072B2", lwd = 1.8),
              SGAT    = list(lab = "SGAT",            col = "#D55E00", lwd = 1.3),
              FLightR = list(lab = "FLightR",         col = "#009E73", lwd = 1.3))

wrap360 <- function(x) x %% 360
have <- function(sc, m) {
  d <- fits[[paste0(sc, "|", m)]]
  if (is.null(d) || !nrow(d)) return(NULL)
  d$lon <- wrap360(d$lon)
  d[order(as.numeric(d$time)), ]
}

# median RMSE per scenario/method, for the panel annotation
rms <- NULL
if (file.exists(CSV)) {
  z <- read.csv(CSV); z <- z[z$table == 2, ]
  rms <- tapply(z$rmse_km, list(z$scenario, z$method), median, na.rm = TRUE)
}

xs <- ys <- numeric(0)
for (sc in SCEN) for (m in names(METH)) {
  d <- have(sc, m); if (is.null(d)) next
  xs <- c(xs, range(d$lon, na.rm = TRUE)); ys <- c(ys, range(d$lat, na.rm = TRUE))
}
padx <- diff(range(xs)) * 0.05; pady <- diff(range(ys)) * 0.05
xlim <- c(min(xs) - padx, max(xs) + padx)
ylim <- c(min(ys) - pady, max(ys) + pady)

land <- tryCatch({
  w <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
  sf::st_geometry(sf::st_shift_longitude(w))
}, error = function(e) NULL)

png(OUT, width = 3300, height = 1500, res = 300)
op <- par(mfrow = c(1, 3), mar = c(3.2, 3.3, 2.5, 0.8), mgp = c(2.0, 0.6, 0),
          cex.axis = 0.85, cex.lab = 0.95)

for (i in seq_along(SCEN)) {
  sc <- SCEN[i]
  plot(NA, xlim = xlim, ylim = ylim, xlab = "longitude (°E)",
       ylab = if (i == 1) "latitude (°)" else "", axes = FALSE)
  if (!is.null(land)) plot(land, add = TRUE, col = "grey92", border = "grey78", lwd = 0.4)
  abline(h = pretty(ylim), v = pretty(xlim), col = "grey93", lwd = 0.5)
  axis(1); axis(2, las = 1); box(col = "grey45")

  for (m in c("truth", "invTF", "SGAT", "FLightR")) {
    d <- have(sc, m); if (is.null(d)) next
    lines(d$lon, d$lat, col = METH[[m]]$col, lwd = METH[[m]]$lwd)
  }
  tr <- have(sc, "truth")
  if (!is.null(tr)) points(tr$lon[1], tr$lat[1], pch = 21, bg = "white", col = "grey20", cex = 1.2)
  mtext(PANEL[i], side = 3, line = 0.6, adj = 0, cex = 0.85, font = 2)

  if (!is.null(rms) && sc %in% rownames(rms)) {
    lab <- sapply(c("invTF", "SGAT", "FLightR"), function(m)
      if (!is.na(rms[sc, m])) sprintf("%s %.0f km", METH[[m]]$lab, rms[sc, m]) else "")
    legend("topleft", bty = "n", cex = 0.72, text.col = c("#0072B2", "#D55E00", "#009E73"),
           legend = lab[lab != ""])
  }
  if (i == 1)
    legend("bottomright", bty = "n", cex = 0.78, seg.len = 1.8,
           legend = sapply(METH, `[[`, "lab"),
           col    = sapply(METH, `[[`, "col"),
           lwd    = sapply(METH, `[[`, "lwd"))
}
par(op); invisible(dev.off())
cat("wrote", OUT, sprintf("(lon %.0f-%.0f, lat %.0f-%.0f)\n", xlim[1], xlim[2], ylim[1], ylim[2]))
