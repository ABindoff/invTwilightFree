# FIGURE: habitat occupancy, invTwilightFree hierarchical against Argos, by season.
#
# The argument the figure has to carry: per-fix accuracy is the wrong yardstick for
# what these tags are deployed to answer. What matters is which water the animals
# used and when. So the figure compares OCCUPANCY DISTRIBUTIONS, not tracks, and
# splits by season because the "when" is half the claim.
#
# Form (following the data-viz procedure):
#   maps      magnitude over space -> ONE sequential hue, light to dark, shared
#             scale across all six panels so they are directly comparable.
#             Identity (Argos vs model) is carried by the column header, not by hue.
#   marginals two distributions compared -> categorical slots 1 and 2, validated:
#             #2a78d6 / #eb6834, CVD dE 24.7, normal-vision 33.6, both >= 3:1 on
#             the surface. Line TYPE also differs, so identity is never colour
#             alone in print or for CVD readers.
#
# Zero-occupancy cells are left at the surface colour rather than given the
# lightest ramp step, so "no data" and "a little" are not confusable.
suppressMessages({ library(data.table) })

D <- as.data.table(readRDS("scratch/nes_calibration/hier_occupancy.rds"))
D[, season := fifelse(month %in% 6:8, "Summer (Jun-Aug)",
             fifelse(month %in% 9:11, "Autumn (Sep-Nov)",
             fifelse(month %in% c(12, 1, 2), "Winter (Dec-Feb)", NA_character_)))]
D <- D[!is.na(season)]
SEAS <- c("Summer (Jun-Aug)", "Autumn (Sep-Nov)", "Winter (Dec-Feb)")

SURFACE <- "#fcfcfb"; INK <- "#0b0b0b"; INK2 <- "#52514e"; GRID <- "#e6e5e2"
ARGOS <- "#2a78d6"; MODEL <- "#eb6834"
RAMP <- colorRampPalette(c("#cde2fb", "#9ec5f4", "#5598e7", "#2a78d6",
                           "#1c5cab", "#104281", "#0d366b"))(64)

LONB <- seq(170, 245, by = 5)
LATB <- seq(32, 60, by = 2)
occ <- function(lon, lat) {
  t <- table(cut(lon, LONB), cut(lat, LATB))
  as.matrix(t) / sum(t)
}
ovl <- function(a, b) sum(pmin(a, b))
latden <- function(x) { d <- density(x, from = 32, to = 60, bw = 1.1, n = 256)
  d$y <- d$y / sum(d$y); d }

# shared colour scale across every map panel
mx <- max(vapply(SEAS, function(s) {
  d <- D[season == s]; max(c(occ(d$true_lon, d$true_lat), occ(d$est_lon, d$est_lat)))
}, 0))

png("inst/paper/occupancy_by_season.png", width = 2700, height = 2500, res = 300)
op <- par(bg = SURFACE, family = "sans")
layout(matrix(c(1:9, 10, 10, 11), nrow = 4, byrow = TRUE),
       heights = c(1, 1, 1, 0.30))

for (si in seq_along(SEAS)) {
  s <- SEAS[si]; d <- D[season == s]
  oa <- occ(d$true_lon, d$true_lat); om <- occ(d$est_lon, d$est_lat)
  for (which in c("argos", "model")) {
    m <- if (which == "argos") oa else om
    par(mar = c(2.6, if (which == "argos") 3.4 else 1.4, 2.4, 0.6),
        mgp = c(1.7, 0.5, 0), tcl = -0.25, col.axis = INK2, col.lab = INK2)
    image(LONB[-length(LONB)] + 2.5, LATB[-length(LATB)] + 1,
          ifelse(m == 0, NA, m), zlim = c(0, mx), col = RAMP,
          xlab = "", ylab = "", axes = FALSE,
          xlim = range(LONB), ylim = range(LATB), useRaster = FALSE)
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],
         border = GRID, lwd = 0.8)
    axis(1, at = seq(180, 240, 20), labels = paste0(seq(180, 240, 20), "°E"),
         cex.axis = 0.8, lwd = 0, lwd.ticks = 0.8, col.ticks = GRID)
    if (which == "argos")
      axis(2, at = seq(35, 55, 5), labels = paste0(seq(35, 55, 5), "°N"),
           las = 1, cex.axis = 0.8, lwd = 0, lwd.ticks = 0.8, col.ticks = GRID)
    if (si == 1) mtext(if (which == "argos") "Argos" else "invTwilightFree (hierarchical)",
                       side = 3, line = 0.7, cex = 0.85, col = INK, font = 2)
    if (which == "argos")
      mtext(s, side = 2, line = 2.3, cex = 0.8, col = INK, font = 2)
    # the 2-D overlap belongs with the MAPS, which is what it measures. Annotating
    # it on the latitude panel invites reading it as describing that curve, and the
    # two differ substantially (0.49 vs 0.78 in autumn) because a 5x2-deg cell is a
    # far stricter test than a 5-deg latitude band.
    if (which == "model")
      legend("topleft", bty = "n", cex = 0.72, text.col = INK2,
             legend = sprintf("2-D overlap %.2f", ovl(oa, om)))
  }
  # latitude marginal
  par(mar = c(2.6, 3.0, 2.4, 1.0), mgp = c(1.7, 0.5, 0))
  da <- latden(d$true_lat); dm <- latden(d$est_lat)
  plot(NA, xlim = c(32, 60), ylim = c(0, max(da$y, dm$y) * 1.28),
       xlab = "", ylab = "", axes = FALSE)
  abline(v = seq(35, 55, 5), col = GRID, lwd = 0.7)
  lines(da$x, da$y, col = ARGOS, lwd = 2.2)
  lines(dm$x, dm$y, col = MODEL, lwd = 2.2, lty = 2)
  axis(1, at = seq(35, 55, 5), labels = paste0(seq(35, 55, 5), "°N"),
       cex.axis = 0.8, lwd = 0, lwd.ticks = 0.8, col.ticks = GRID)
  mtext("proportion of time", side = 2, line = 0.9, cex = 0.7, col = INK2)
  if (si == 1) mtext("latitude occupancy", side = 3, line = 0.7, cex = 0.85,
                     col = INK, font = 2)
  lat5 <- seq(10, 80, by = 5)
  pa <- prop.table(table(cut(d$true_lat, lat5)))
  pb <- prop.table(table(cut(d$est_lat, lat5)))
  sh <- mean(d$est_lat) - mean(d$true_lat)
  legend("topleft", bty = "n", cex = 0.72, text.col = INK2,
         legend = c(sprintf("latitude overlap %.2f", sum(pmin(pa, pb))),
                    sprintf("mean shift %+.2f°", sh),
                    sprintf("n = %d", nrow(d))))
}

# shared colour bar
par(mar = c(2.6, 3.4, 1.4, 1.4))
z <- seq(0, mx, length.out = 200)
image(z, 1, matrix(z, ncol = 1), col = RAMP, axes = FALSE, xlab = "", ylab = "")
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], border = GRID, lwd = 0.8)
axis(1, at = pretty(c(0, mx), 4), labels = sprintf("%.1f%%", pretty(c(0, mx), 4) * 100),
     cex.axis = 0.78, lwd = 0, lwd.ticks = 0.8, col.ticks = GRID, col.axis = INK2)
mtext("proportion of time in cell (5° × 2°)", side = 3, line = 0.2,
      cex = 0.75, col = INK2, adj = 0)

par(mar = c(2.6, 1.0, 1.4, 1.0))
plot.new()
legend("center", bty = "n", cex = 0.95, seg.len = 2.4, text.col = INK,
       legend = c("Argos", "invTwilightFree"), col = c(ARGOS, MODEL),
       lwd = 2.2, lty = c(1, 2))
par(op); dev.off()
cat("written: inst/paper/occupancy_by_season.png\n\n")

cat("=== numbers behind the figure ===\n")
for (s in SEAS) {
  d <- D[season == s]
  oa <- occ(d$true_lon, d$true_lat); om <- occ(d$est_lon, d$est_lat)
  cat(sprintf("  %-18s n=%4d | 2-D occupancy overlap %.3f | lat shift %+.2f deg | median |lat err| %.2f deg\n",
              s, nrow(d), ovl(oa, om), mean(d$est_lat) - mean(d$true_lat),
              median(abs(d$est_lat - d$true_lat))))
}
oa <- occ(D$true_lon, D$true_lat); om <- occ(D$est_lon, D$est_lat)
cat(sprintf("  %-18s n=%4d | 2-D occupancy overlap %.3f | lat shift %+.2f deg\n",
            "ALL SEASONS", nrow(D), ovl(oa, om), mean(D$est_lat) - mean(D$true_lat)))
