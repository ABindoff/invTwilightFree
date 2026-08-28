# =============================================================================
# FIGURE: seasonal habitat occupancy, invTwilightFree against Argos.
#
# The argument the figure carries: per-knot distance is the estimand of a
# methods paper, not of the ecology these tags are deployed to do. What matters
# is which water the animals used and when, so the comparison is between
# OCCUPANCY DISTRIBUTIONS, and it is split by season because the "when" is half
# the claim.
#
# The model column is the EXACT integrated posterior (analysis/study/
# refit_posteriors.R), not a track of point estimates: every knot contributes
# its whole cell-probability surface, so a position the light does not determine
# spreads its mass instead of being placed confidently in one cell.
#
# Form, following the data-visualisation procedure:
#   maps       magnitude over space -> ONE sequential hue, light to dark, with a
#              single scale shared by every panel so they are directly
#              comparable. Identity (Argos or model) is carried by the column
#              header, never by hue.
#   marginals  two distributions compared -> categorical slots 1 and 2,
#              #2a78d6 / #eb6834, which separate under both normal vision and
#              the common CVD forms. Line TYPE differs as well, so identity is
#              never colour alone in print.
#   extent     an HMM posterior assigns SOME probability to every cell, so a
#              raw surface is a full-domain wash and "never went there" cannot be
#              read off it. Both columns are therefore masked to the smallest set
#              of cells holding 95% of the mass -- the 95% utilisation
#              distribution, the standard object in this literature -- and the
#              same rule is applied to Argos, so the columns stay comparable.
#              The overlap statistics annotated on the panels are computed on the
#              UNMASKED surfaces, so the display choice cannot flatter them.
#   scale      one shared scale across every panel, on a square-root transform.
#              A point-estimate reference concentrates an order of magnitude more
#              mass in its peak cell than an integrated posterior does; on a
#              linear shared scale the posterior column renders as blank paper,
#              which is a statement about the colour ramp rather than about the
#              animals.
#
#   Rscript inst/paper/fig_occupancy.R
# =============================================================================

suppressMessages({library(data.table); library(terra); library(invTwilightFree)})
source("analysis/study/study_config.R")

POST <- file.path(OUT_DIR, "posteriors")
OUT  <- "inst/paper/fig_occupancy.png"
CELL <- 2                      # degrees; the resolution the text reports
files <- list.files(POST, pattern = "[.]rds$", full.names = TRUE)
if (!length(files)) stop("no posteriors -- run analysis/study/refit_posteriors.R")

SURFACE <- "#fcfcfb"; INK <- "#0b0b0b"; INK2 <- "#52514e"; GRID <- "#e6e5e2"
LAND <- "#eceae5"; LANDB <- "#d6d3cc"
ARGOS <- "#2a78d6"; MODEL <- "#eb6834"
RAMP <- colorRampPalette(c("#cde2fb", "#9ec5f4", "#5598e7", "#2a78d6",
                          "#1c5cab", "#104281", "#0d366b"))(64)

season_of <- function(m) c("Winter","Winter","Spring","Spring","Spring","Summer",
                           "Summer","Summer","Autumn","Autumn","Autumn","Winter")[m]

base <- rast(xmin = DOMAIN[["xmin"]], xmax = DOMAIN[["xmax"]],
             ymin = DOMAIN[["ymin"]], ymax = DOMAIN[["ymax"]],
             resolution = GRID_CELL_DEG, crs = "EPSG:4326")
values(base) <- 1

# ---- model: exact posterior, equal weight per deployment --------------------
model_v <- list(); n_knot <- list()
for (f in files) {
  o <- readRDS(f)
  s <- season_of(as.integer(substr(o$months, 6, 7)))
  for (ss in unique(s)) {
    v <- colSums(o$P[s == ss, , drop = FALSE]); tot <- sum(v)
    if (!is.finite(tot) || tot <= 0) next
    model_v[[ss]] <- (if (is.null(model_v[[ss]])) 0 else model_v[[ss]]) + v / tot
    n_knot[[ss]] <- (if (is.null(n_knot[[ss]])) 0 else n_knot[[ss]]) +
      sum(o$knots_per_month[s == ss])
  }
}

# ---- reference: Argos at the scored knots ----------------------------------
K <- fread(file.path(OUT_DIR, "study_knots.csv"), colClasses = list(character = "id"))
K <- K[arm == "grid_light" & scored == TRUE]
K[, time := as.POSIXct(time, tz = "UTC")]
K[, season := season_of(as.integer(format(time, "%m")))]
argos_r <- occupancy_map(
  setNames(lapply(unique(K$id), function(i)
    data.frame(time = K[id == i]$time, lon = K[id == i]$truth_lon,
               lat = K[id == i]$truth_lat)), unique(K$id)),
  grid = base, method = "point", weight = "tag",
  by = function(t) season_of(as.integer(format(t, "%m"))))

# Seasons carried by enough of the panel to plot. Spring is the tail of the
# longest deployments only and is reported in the text rather than shown.
SEAS <- c("Summer", "Autumn", "Winter")
SEAS <- SEAS[SEAS %in% names(model_v) & SEAS %in% names(argos_r)]

agg <- function(r) aggregate(r, fact = CELL / GRID_CELL_DEG, fun = "sum", na.rm = TRUE)
as_r <- function(v) { r <- base; values(r) <- v / sum(v); agg(r) }
A <- lapply(SEAS, function(s) agg(argos_r[[s]]))
M <- lapply(SEAS, function(s) as_r(model_v[[s]]))
names(A) <- names(M) <- SEAS

# The smallest set of cells holding `p` of the mass: the 95% utilisation
# distribution. Everything outside is set NA and drawn as bare surface.
hdr <- function(r, p = 0.95) {
  v <- values(r)[, 1]
  ok <- is.finite(v) & v > 0
  o <- order(v[ok], decreasing = TRUE)
  cum <- cumsum(v[ok][o]) / sum(v[ok])
  cut <- v[ok][o][min(which(cum >= p))]
  out <- r; vv <- v; vv[!is.finite(vv) | vv < cut] <- NA
  values(out) <- vv
  out
}
Am <- lapply(A, hdr); Mm <- lapply(M, hdr)
names(Am) <- names(Mm) <- SEAS

# Plot extent: the union of the two 95% distributions, not the search domain.
occupied <- Reduce(function(a, b) {
  x <- a; values(x) <- pmax(values(a)[, 1], values(b)[, 1], na.rm = TRUE); x
}, c(Am, Mm))
xy <- crds(occupied, na.rm = FALSE); vv <- values(occupied)[, 1]
keep <- is.finite(vv)
XL <- range(xy[keep, 1]) + c(-CELL, CELL)
YL <- range(xy[keep, 2]) + c(-CELL, CELL)

# Square-root colour scale, shared. See the header for why linear fails here.
mx <- max(vapply(c(Am, Mm), function(r) max(values(r), na.rm = TRUE), 0))
sq <- function(v) sqrt(pmax(v, 0) / mx)
ovl <- function(a, b) {
  va <- values(a)[, 1]; vb <- values(b)[, 1]
  ok <- is.finite(va) & is.finite(vb)
  sum(pmin(va[ok] / sum(va[ok]), vb[ok] / sum(vb[ok])))
}
land <- try(suppressWarnings(
  terra::project(terra::vect(rnaturalearth::ne_countries(scale = "medium",
                                                         returnclass = "sf")),
                 "EPSG:4326")), silent = TRUE)
draw_land <- function() {
  if (inherits(land, "try-error")) return(invisible(NULL))
  for (sh in c(0, 360))
    try(plot(land, add = TRUE, col = LAND, border = LANDB, lwd = 0.4,
             ext = c(XL[1] - sh, XL[2] - sh, YL)), silent = TRUE)
}

lat_marg <- function(r) {
  v <- values(r)[, 1]; xy <- crds(r, na.rm = FALSE)
  b <- tapply(v, xy[, 2], sum, na.rm = TRUE)
  list(lat = as.numeric(names(b)), p = as.numeric(b) / sum(b, na.rm = TRUE))
}

png(OUT, width = 2700, height = 2500, res = 300)
op <- par(bg = SURFACE, family = "sans")
layout(matrix(c(1:(3 * length(SEAS)),
                rep(3 * length(SEAS) + 1, 2), 3 * length(SEAS) + 2),
              nrow = length(SEAS) + 1, byrow = TRUE),
       heights = c(rep(1, length(SEAS)), 0.30))

for (si in seq_along(SEAS)) {
  s <- SEAS[si]
  for (which in c("argos", "model")) {
    r <- if (which == "argos") Am[[s]] else Mm[[s]]
    par(mar = c(2.6, if (which == "argos") 3.4 else 1.4, 2.4, 0.6),
        mgp = c(1.7, 0.5, 0), tcl = -0.25, col.axis = INK2)
    plot(NA, xlim = XL, ylim = YL, xlab = "", ylab = "", axes = FALSE,
         xaxs = "i", yaxs = "i")
    draw_land()
    rr <- r; values(rr) <- sq(values(r)[, 1])
    image(rr, add = TRUE, zlim = c(0, 1), col = RAMP, useRaster = FALSE)
    # Ano Nuevo. Both ends of every track are pinned here, so it is the one
    # position in the figure that carries no estimation error, and a reader
    # needs it to see that the winter distribution is a return to the colony.
    points(COLONY[["lon"]], COLONY[["lat"]], pch = 21, cex = 0.9,
           col = INK, bg = SURFACE, lwd = 1.1)
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],
         border = GRID, lwd = 0.8)
    at_x <- pretty(XL, 4)
    axis(1, at = at_x, labels = paste0(at_x, "°E"), cex.axis = 0.8,
         lwd = 0, lwd.ticks = 0.8, col.ticks = GRID)
    if (which == "argos") {
      at_y <- pretty(YL, 4)
      axis(2, at = at_y, labels = paste0(at_y, "°N"), las = 1,
           cex.axis = 0.8, lwd = 0, lwd.ticks = 0.8, col.ticks = GRID)
      mtext(s, side = 2, line = 2.3, cex = 0.8, col = INK, font = 2)
    }
    if (si == 1) {
      mtext(if (which == "argos") "Argos" else "invTwilightFree (posterior)",
            side = 3, line = 0.7, cex = 0.85, col = INK, font = 2)
      if (which == "argos")
        text(COLONY[["lon"]] - 1.5, COLONY[["lat"]] - 2.5, "Ano Nuevo",
             cex = 0.62, col = INK2, adj = 1)
    }
    # The 2-D overlap belongs with the MAPS, which is what it measures. Putting
    # it on the latitude panel invites reading it as describing that curve, and
    # the two differ substantially because a 2x2-degree cell is a far stricter
    # test than a latitude band.
    #
    # "(all knots)" is not decoration. The model surface integrates every knot
    # while the Argos reference exists only where a satellite fix falls within
    # 24 h, so this number is not on a matched knot set. It is the right number
    # for the maps drawn here; the matched-knot figures the text quotes are in
    # analysis/study/occupancy_exact.R, and they differ by up to 0.02 in the
    # seasons shown.
    if (which == "model")
      legend("topleft", bty = "n", cex = 0.72, text.col = INK2,
             bg = SURFACE,
             legend = sprintf("2-D overlap %.2f (all knots)", ovl(A[[s]], M[[s]])))
  }
  par(mar = c(2.6, 3.0, 2.4, 1.0), mgp = c(1.7, 0.5, 0))
  ma <- lat_marg(A[[s]]); mm <- lat_marg(M[[s]])
  plot(NA, xlim = YL, ylim = c(0, max(ma$p, mm$p) * 1.30),
       xlab = "", ylab = "", axes = FALSE)
  abline(v = pretty(YL, 4), col = GRID, lwd = 0.7)
  lines(ma$lat, ma$p, col = ARGOS, lwd = 2.2)
  lines(mm$lat, mm$p, col = MODEL, lwd = 2.2, lty = 2)
  at_y <- pretty(YL, 4)
  axis(1, at = at_y, labels = paste0(at_y, "°N"), cex.axis = 0.8,
       lwd = 0, lwd.ticks = 0.8, col.ticks = GRID)
  mtext("proportion of time", side = 2, line = 0.9, cex = 0.7, col = INK2)
  if (si == 1) mtext("latitude occupancy", side = 3, line = 0.7, cex = 0.85,
                     col = INK, font = 2)
  shift <- sum(mm$lat * mm$p) - sum(ma$lat * ma$p)
  legend("topleft", bty = "n", cex = 0.72, text.col = INK2,
         legend = c(sprintf("latitude overlap %.2f", sum(pmin(ma$p, mm$p))),
                    sprintf("mean shift %+.2f°", shift),
                    sprintf("n = %s knots", format(n_knot[[s]], big.mark = ","))))
}

par(mar = c(2.6, 3.4, 1.4, 1.4))
z <- seq(0, 1, length.out = 200)
image(z, 1, matrix(z, ncol = 1), col = RAMP, axes = FALSE, xlab = "", ylab = "")
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], border = GRID, lwd = 0.8)
tk <- c(0, 0.01, 0.05, 0.10, 0.15)
tk <- tk[tk <= mx]
axis(1, at = sqrt(tk / mx), labels = sprintf("%.0f%%", tk * 100), cex.axis = 0.78,
     lwd = 0, lwd.ticks = 0.8, col.ticks = GRID, col.axis = INK2)
mtext(sprintf("proportion of time in cell (%d° × %d°, square-root scale)", CELL, CELL),
      side = 3, line = 0.2, cex = 0.75, col = INK2, adj = 0)

par(mar = c(2.6, 1.0, 1.4, 1.0))
plot.new()
legend("center", bty = "n", cex = 0.95, seg.len = 2.4, text.col = INK,
       legend = c("Argos", "invTwilightFree"), col = c(ARGOS, MODEL),
       lwd = 2.2, lty = c(1, 2))
par(op); dev.off()

cat(sprintf("wrote %s\n", OUT))
for (s in SEAS)
  cat(sprintf("  %-7s 2-D overlap %.3f | knots %s\n", s, ovl(A[[s]], M[[s]]),
              format(n_knot[[s]], big.mark = ",")))
