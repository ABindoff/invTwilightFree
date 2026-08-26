#' Population occupancy from reconstructed tracks
#'
#' Per-knot position error is the wrong yardstick for most of the science
#' archival tags are deployed to do. Few questions in movement ecology need to
#' know where an animal was on a Tuesday; they need to know which water it used
#' and when. `occupancy_map()` reduces one or many reconstructed tracks to a
#' utilisation distribution over a grid, optionally split by period, so that the
#' estimand is the one the ecology is organised around.
#'
#' Two forms are provided and they answer slightly different questions.
#'
#' \describe{
#'   \item{`method = "posterior"`}{Integrates the **whole per-knot posterior**.
#'     For a [TwilightFreeGrid()] fit this is exact: the forward-backward
#'     recursions already give the marginal probability of every cell at every
#'     knot, and the occupancy map is their weighted sum. A location the model
#'     is unsure about spreads its mass, which is the honest representation of
#'     what a light record supports.}
#'   \item{`method = "point"`}{Bins the reported point estimate, by default the
#'     per-knot maximum a posteriori cell. This is what a reader gets by
#'     plotting the track, and it is deliberately over-confident: every knot
#'     contributes its whole weight to one cell however uncertain it was.}
#' }
#'
#' Reporting both is informative. Where they agree, the occupancy estimate is
#' not an artefact of how the posterior was summarised; where the point map is
#' concentrated and the posterior map is not, the apparent structure is a
#' readout of the movement prior rather than of the light.
#'
#' @section Which errors average out, and which do not:
#' Pooling many tracks divides the *independent* part of position error by
#' roughly the square root of the number of tags, and the *shared* part not at
#' all. Calibration in this package is pooled across a tag family by
#' construction (see [pool_light_responses()]), so per-tag errors are correlated
#' rather than independent and a systematic displacement survives any amount of
#' pooling. An occupancy map that matches a reference is therefore evidence that
#' the residual bias is small **against the scale of the habitat structure**,
#' not evidence that the sample was large. Those are different claims and are
#' worth separating in print. Comparing the per-tag overlap distribution with
#' the pooled overlap ([occupancy_overlap()]) distinguishes them: if the error
#' were independent, the pooled map would beat the typical per-tag map.
#'
#' @param x Tracks to summarise. Accepts a single [TwilightFreeGrid()],
#'   [TwilightFreeSMC()] or [TwilightFreeHier()] fit, a named list of fits, or a
#'   `data.frame` with columns `time`, `lon` and `lat` (and optionally `lon_sd`
#'   and `lat_sd`). A plain data frame is how a reference track — Argos, GPS —
#'   is passed in for comparison.
#' @param grid A `SpatRaster` giving the target cells. Values are ignored; only
#'   the geometry is used. Cells masked `NA` are dropped from the output.
#' @param method `"posterior"` to integrate the per-knot posterior,
#'   `"point"` to bin the point estimate. See Details.
#' @param by `NULL` for a single map, or a function taking a `POSIXct` vector
#'   and returning a period label per knot, for example
#'   `function(t) format(t, "%b")` for calendar months. One layer is produced
#'   per distinct label.
#' @param weight `"tag"` gives every track equal total weight, so a 300-day
#'   deployment does not outvote a 130-day one; `"knot"` weights by time
#'   recorded, so the map reflects animal-days. The choice is ecological, not
#'   technical, and should be stated.
#' @param per_area Divide each cell's mass by its area, giving occupancy per
#'   square kilometre rather than per cell. Over a wide latitude span this
#'   matters: a one-degree cell at 65 degrees north holds under half the area of
#'   one at 20 degrees. Defaults to `FALSE`, which keeps layers summing to one.
#' @param sd_floor_km Smallest per-knot standard deviation, in km, used when
#'   building a Gaussian kernel for a track fit. Guards against a degenerate
#'   posterior contributing a spike to a single cell. Ignored for grid fits and
#'   for `method = "point"`.
#'
#' @return A `SpatRaster` with one layer per period, named for it. Each layer
#'   sums to one over its non-`NA` cells (before `per_area`), so layers are
#'   directly comparable and can be passed to [occupancy_overlap()].
#'
#' @seealso [occupancy_overlap()] to compare two maps, [grid_posterior()] for
#'   the underlying per-knot surfaces.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("terra", quietly = TRUE)) {
#'   g <- makeGrid(lon = c(140, 180), lat = c(-70, -40), cell.size = 2)
#'   # a reference track summarised the same way as a fit
#'   ref <- data.frame(
#'     time = seq(as.POSIXct("2021-06-01", tz = "UTC"), by = "12 hours",
#'                length.out = 60),
#'     lon = seq(150, 170, length.out = 60),
#'     lat = seq(-60, -50, length.out = 60))
#'   m <- occupancy_map(ref, grid = g, method = "point")
#'   terra::global(m, "sum", na.rm = TRUE)
#' }
#' }
#' @export
occupancy_map <- function(x, grid, method = c("posterior", "point"),
                          by = NULL, weight = c("tag", "knot"),
                          per_area = FALSE, sd_floor_km = 25) {
  method <- match.arg(method)
  weight <- match.arg(weight)
  if (!inherits(grid, "SpatRaster"))
    stop("`grid` must be a SpatRaster giving the target cells")
  if (!is.null(by) && !is.function(by))
    stop("`by` must be NULL or a function of a POSIXct vector")

  tracks <- as_track_list(x)
  if (!length(tracks)) stop("`x` contains no tracks")

  xy    <- terra::crds(grid, na.rm = FALSE)
  keep  <- !is.na(terra::values(grid)[, 1])
  ncell <- terra::ncell(grid)
  xmin  <- terra::xmin(grid); xmax <- terra::xmax(grid)

  # Longitude conventions differ between a fit and a grid: the SMC engine
  # reports -180..180 while a Pacific analysis works in 0..360. Silently
  # indexing across that mismatch is a bug this package has already shipped
  # once, so wrap explicitly here.
  wrap_lon <- function(lon) {
    lon <- ((lon - xmin) %% 360) + xmin
    lon[lon > xmax] <- lon[lon > xmax] - 360
    lon
  }

  labels_of <- function(times) {
    if (is.null(by)) rep("all", length(times)) else as.character(by(times))
  }

  # Accumulate one tag at a time into its own per-period vectors, normalise
  # those, and only then fold them into the panel total. Doing it in that order
  # is what makes `weight` mean anything: a tag's contribution is rescaled while
  # it is still separable. A full posterior is knots x cells, so a panel of them
  # held at once would be gigabytes for no reason.
  acc <- new.env(parent = emptyenv())

  for (nm in names(tracks)) {
    tr  <- tracks[[nm]]
    per <- labels_of(tr$time)
    mine <- new.env(parent = emptyenv())   # period -> numeric(ncell), this tag
    put <- function(period, cells, mass) {
      ok <- !is.na(cells) & is.finite(mass) & mass > 0
      if (!any(ok)) return(invisible(NULL))
      v <- mine[[period]]
      if (is.null(v)) v <- numeric(ncell)
      mine[[period]] <- v + tapply_sum(cells[ok], mass[ok], ncell)
      invisible(NULL)
    }

    if (method == "posterior" && !is.null(tr$posterior)) {
      # Exact: the grid HMM's own marginal posterior over cells, per knot.
      P <- tr$posterior
      src <- terra::cellFromXY(grid, cbind(wrap_lon(P$lon), P$lat))
      for (k in seq_len(nrow(P$P))) {
        w <- P$P[k, ]; s <- sum(w)
        if (!is.finite(s) || s <= 0) next
        put(per[k], src, w / s)
      }
    } else if (method == "posterior") {
      # Approximate: a Gaussian kernel per knot from the reported marginal
      # spreads. Honest for an engine that reports only a mean and an sd, but it
      # is a normal approximation to a posterior that need not be normal.
      if (is.null(tr$lat_sd) || is.null(tr$lon_sd))
        stop("`", nm, "` carries neither a cell posterior nor per-knot standard ",
             "deviations, so method = \"posterior\" is not available for it")
      for (k in seq_along(tr$lon)) {
        m <- knot_kernel(xy, wrap_lon(tr$lon[k]), tr$lat[k],
                         tr$lon_sd[k], tr$lat_sd[k], sd_floor_km)
        if (is.null(m)) next
        put(per[k], m$cell, m$mass)
      }
    } else {
      cells <- terra::cellFromXY(grid, cbind(wrap_lon(tr$lon), tr$lat))
      for (p in unique(per)) {
        i <- per == p
        put(p, cells[i], rep(1, sum(i)))
      }
    }

    for (p in ls(mine)) {
      v <- mine[[p]]
      tot <- sum(v)
      if (!is.finite(tot) || tot <= 0) next
      # "tag": every track contributes one unit to every period it appears in,
      # so a 313-day deployment cannot outvote a 137-day one. "knot": leave the
      # raw mass, so the map reflects animal-days.
      if (weight == "tag") v <- v / tot
      acc[[p]] <- (if (is.null(acc[[p]])) numeric(ncell) else acc[[p]]) + v
    }
  }

  periods <- sort(ls(acc))
  if (!length(periods)) stop("no knot fell inside `grid`")

  out <- terra::rast(grid, nlyrs = length(periods))
  names(out) <- periods
  for (i in seq_along(periods)) {
    v <- acc[[periods[i]]]
    v[!keep] <- NA_real_
    s <- sum(v, na.rm = TRUE)
    if (is.finite(s) && s > 0) v <- v / s
    terra::values(out[[i]]) <- v
  }
  if (per_area) out <- out / cell_area_km2(grid)

  attr(out, "occupancy") <- list(method = method, weight = weight,
                                 n_tracks = length(tracks), per_area = per_area)
  out
}

# tapply-free accumulation: fast and allocation-light for long cell vectors.
tapply_sum <- function(cells, mass, n) {
  out <- numeric(n)
  o <- order(cells)
  cs <- cells[o]; ms <- mass[o]
  r <- rle(cs)
  ends <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1L
  sums <- vapply(seq_along(r$values), function(i) sum(ms[starts[i]:ends[i]]), 0)
  out[r$values] <- sums
  out
}

# A Gaussian kernel for one knot, evaluated only on cells within 4 sd so a
# panel of long tracks does not become an O(knots x cells) dense operation.
knot_kernel <- function(xy, lon, lat, lon_sd, lat_sd, sd_floor_km) {
  if (!is.finite(lon) || !is.finite(lat)) return(NULL)
  km_per_deg_lat <- 111.32
  km_per_deg_lon <- 111.32 * cos(lat * pi / 180)
  s_lat <- max(if (is.finite(lat_sd)) lat_sd else 0, sd_floor_km / km_per_deg_lat)
  s_lon <- max(if (is.finite(lon_sd)) lon_sd else 0,
               sd_floor_km / max(km_per_deg_lon, 1e-6))
  near <- abs(xy[, 2] - lat) <= 4 * s_lat &
          abs(((xy[, 1] - lon + 180) %% 360) - 180) <= 4 * s_lon
  if (!any(near)) return(NULL)
  d_lat <- (xy[near, 2] - lat) / s_lat
  d_lon <- (((xy[near, 1] - lon + 180) %% 360) - 180) / s_lon
  m <- exp(-0.5 * (d_lat^2 + d_lon^2))
  s <- sum(m)
  if (!is.finite(s) || s <= 0) return(NULL)
  list(cell = which(near), mass = m / s)
}

cell_area_km2 <- function(grid) {
  a <- terra::cellSize(grid, unit = "km")
  a
}

# Normalise the many shapes a fit can arrive in to time / lon / lat (/ sds /
# posterior). Keeping this in one place is what lets occupancy_map() accept a
# grid fit, an SMC fit, a hierarchical panel and a reference data frame without
# any of them being a special case downstream.
as_track_list <- function(x) {
  one <- function(f, nm) {
    if (inherits(f, "TwilightFreeGrid")) {
      gp <- try(grid_posterior(f), silent = TRUE)
      list(time = as.POSIXct(f$fit$time, origin = "1970-01-01", tz = "UTC"),
           lon = f$fit$lon, lat = f$fit$lat,
           lon_sd = NULL, lat_sd = NULL,
           posterior = if (inherits(gp, "try-error")) NULL else gp)
    } else if (inherits(f, "TwilightFreeTrack") ||
               (is.list(f) && !is.null(f$knot_times))) {
      list(time = as.POSIXct(f$knot_times, origin = "1970-01-01", tz = "UTC"),
           lon = f$lon, lat = f$lat, lon_sd = f$lon_sd, lat_sd = f$lat_sd,
           posterior = NULL)
    } else if (is.data.frame(f)) {
      need <- c("time", "lon", "lat")
      if (!all(need %in% names(f)))
        stop("a data.frame track needs columns time, lon and lat")
      list(time = as.POSIXct(f$time, tz = "UTC"), lon = f$lon, lat = f$lat,
           lon_sd = f$lon_sd, lat_sd = f$lat_sd, posterior = NULL)
    } else if (is.list(f) && !is.null(f$lon) && !is.null(f$lat) && !is.null(f$time)) {
      list(time = as.POSIXct(f$time, origin = "1970-01-01", tz = "UTC"),
           lon = f$lon, lat = f$lat, lon_sd = f$lon_sd, lat_sd = f$lat_sd,
           posterior = NULL)
    } else {
      stop("cannot read a track out of element ", nm)
    }
  }
  if (inherits(x, "TwilightFreeHier")) {
    return(stats::setNames(lapply(names(x$tracks), function(i) one(x$tracks[[i]], i)),
                           names(x$tracks)))
  }
  if (inherits(x, c("TwilightFreeGrid", "TwilightFreeTrack")) || is.data.frame(x))
    return(list(track = one(x, "track")))
  if (is.list(x)) {
    nms <- names(x)
    if (is.null(nms)) nms <- paste0("track", seq_along(x))
    return(stats::setNames(lapply(seq_along(x), function(i) one(x[[i]], nms[i])), nms))
  }
  stop("`x` must be a fit, a list of fits, or a data.frame of time/lon/lat")
}

#' Overlap between two occupancy maps
#'
#' Quantifies how much of one utilisation distribution falls where another one
#' does. This is the number that turns "the occupancy map looks about right"
#' into a result, and it is the natural way to report a light-geolocation
#' reconstruction against an independent reference: per-fix error can be
#' hundreds of kilometres while the occupancy distributions agree closely, and
#' that combination is the substantive claim.
#'
#' Two statistics are returned because they weight disagreement differently.
#' `overlap` is the shared mass, \eqn{\sum_i \min(a_i, b_i)}, equal to one minus
#' the total-variation distance; it is interpretable as "the proportion of time
#' the two distributions place in the same cells". `bhattacharyya` is
#' \eqn{\sum_i \sqrt{a_i b_i}}, which is more forgiving of a distribution that
#' is diffuse rather than displaced, so a reconstruction that is honestly
#' uncertain scores better on it than on `overlap`. Reporting both separates
#' "in the wrong place" from "not sure where".
#'
#' @param a,b `SpatRaster` occupancy maps on the same geometry, as returned by
#'   [occupancy_map()]. Multi-layer maps are compared layer by layer, matched by
#'   name where both are named.
#'
#' @return A `data.frame` with one row per compared layer and columns `layer`,
#'   `overlap`, `bhattacharyya` and `correlation` (Pearson, over non-`NA`
#'   cells).
#'
#' @seealso [occupancy_map()]
#' @export
occupancy_overlap <- function(a, b) {
  if (!inherits(a, "SpatRaster") || !inherits(b, "SpatRaster"))
    stop("`a` and `b` must both be SpatRaster occupancy maps")
  if (!terra::compareGeom(a, b, stopOnError = FALSE))
    stop("`a` and `b` must share the same grid geometry")
  la <- names(a); lb <- names(b)
  common <- if (!is.null(la) && !is.null(lb) && length(intersect(la, lb)))
    intersect(la, lb) else NULL
  idx <- if (!is.null(common)) lapply(common, function(n) c(match(n, la), match(n, lb)))
         else lapply(seq_len(min(terra::nlyr(a), terra::nlyr(b))), function(i) c(i, i))
  nms <- if (!is.null(common)) common else la[vapply(idx, `[`, 1L, 1)]

  do.call(rbind, lapply(seq_along(idx), function(i) {
    va <- terra::values(a[[idx[[i]][1]]])[, 1]
    vb <- terra::values(b[[idx[[i]][2]]])[, 1]
    ok <- is.finite(va) & is.finite(vb)
    va <- va[ok]; vb <- vb[ok]
    if (!length(va)) return(NULL)
    va <- va / sum(va); vb <- vb / sum(vb)
    data.frame(layer = nms[i],
               overlap = sum(pmin(va, vb)),
               bhattacharyya = sum(sqrt(va * vb)),
               correlation = if (stats::sd(va) > 0 && stats::sd(vb) > 0)
                 stats::cor(va, vb) else NA_real_,
               row.names = NULL)
  }))
}
