#' Posterior extent of a fitted grid, for refining the search domain
#'
#' Returns the smallest bounding box containing `1 - epsilon` of the posterior mass
#' at every knot, expanded by a movement margin. Use it to shrink a coarse search
#' domain before re-fitting at the resolution you actually want.
#'
#' @details
#' Longitude is a carrier phase and is identified far more sharply than latitude: on
#' 2184 knots of elephant-seal track over a 100-degree domain, the smallest
#' contiguous longitude band holding all but 1e-6 of the marginal posterior was 15
#' columns wide at the median and never exceeded 21. Most of a generous starting
#' domain is therefore longitudes that never hold any mass, and evaluating the
#' emission and the transition there is wasted work.
#'
#' **Why a margin is not optional.** The transition kernel is evaluated out to five
#' movement standard deviations and is *not* renormalised per source cell, so
#' probability directed outside the domain is silently dropped. Shrink the domain to
#' the posterior alone and every knot near the new edge loses part of its outgoing
#' mass, which is a position-dependent tilt rather than a harmless approximation. The
#' returned box is therefore padded by `n_sigma` movement standard deviations so that
#' the mass crossing the new boundary is negligible by construction. Pad in time as
#' well by widening `epsilon` if the track is long.
#'
#' **Take the envelope, never the mode.** Near an equinox the latitude posterior is
#' genuinely bimodal, and a box drawn around the MAP path would discard the competing
#' hemisphere. This function works from the posterior mass at every knot, so a bimodal
#' knot contributes both modes.
#'
#' @param fit A `TwilightFreeGrid` object from a coarse fit.
#' @param epsilon Posterior mass allowed to fall outside the box, per knot per axis
#'   (default 1e-6). Smaller is safer and slower. At 1e-9 the box was within one cell
#'   of the full-domain answer on every track tested.
#' @param diffusion,step_hours The movement scale and knot spacing the FINE fit will
#'   use, needed to size the margin. Read from the fit when not supplied.
#' @param n_sigma Movement standard deviations of padding (default 5, matching the
#'   engine's own transition cutoff).
#' @param cell_margin Extra padding in degrees, for the resolution of the fit the box
#'   was drawn from. A coarse fit cannot localise the posterior better than its own
#'   cell size, so a box taken from it must be padded by at least one cell or it will
#'   clip. [TwilightFreeGridRefine()] sets this to `coarse.size` automatically.
#' @param clock_margin_deg Extra longitude padding, degrees. Longitude is set by the
#'   tag clock, so a drifting clock displaces it; if the clock is uncalibrated,
#'   1 degree per 4 minutes of expected drift is the conversion.
#' @return Numeric vector `c(lon_min, lon_max, lat_min, lat_max)` with attributes
#'   `retained` (fraction of the original cells the box covers) and `epsilon`.
#' @seealso [makeGrid()], [TwilightFreeGrid()]
#' @export
posterior_extent <- function(fit, epsilon = 1e-6, diffusion = NULL, step_hours = NULL,
                             n_sigma = 5, clock_margin_deg = 0, cell_margin = 0) {
  stopifnot(inherits(fit, "TwilightFreeGrid"))
  if (epsilon <= 0 || epsilon >= 1) stop("`epsilon` must be in (0, 1)")
  gp <- grid_posterior(fit)
  ulon <- sort(unique(gp$lon)); ulat <- sort(unique(gp$lat))

  # Per knot, the smallest CONTIGUOUS interval on each axis holding 1 - epsilon of
  # that axis's marginal. Contiguous because the fine grid is a box: a rule that
  # returned a scattered set could not be expressed as one, and the union over knots
  # is what the box must cover anyway.
  span <- function(vals, uniq) {
    lo <- rep(NA_real_, nrow(gp$P)); hi <- lo
    for (k in seq_len(nrow(gp$P))) {
      w <- gp$P[k, ]; s <- sum(w)
      if (!is.finite(s) || s <= 0) next
      m <- vapply(uniq, function(v) sum(w[vals == v]), 0) / s
      ord <- order(m, decreasing = TRUE)
      keep <- uniq[ord[seq_len(which(cumsum(m[ord]) >= 1 - epsilon)[1])]]
      lo[k] <- min(keep); hi[k] <- max(keep)
    }
    c(min(lo, na.rm = TRUE), max(hi, na.rm = TRUE))
  }
  lon_r <- span(gp$lon, ulon)
  lat_r <- span(gp$lat, ulat)

  if (is.null(diffusion))  diffusion  <- fit$diffusion
  if (is.null(step_hours)) step_hours <- fit$step_hours
  if (is.null(diffusion) || is.null(step_hours))
    stop("`diffusion` and `step_hours` are needed to size the movement margin; ",
         "pass them explicitly if the fit does not carry them")
  sig_km <- max(diffusion) * sqrt(max(step_hours) / 24)
  pad_km <- n_sigma * sig_km
  pad_lat <- pad_km / 111.32
  # widen longitude at the poleward edge of the box, where a degree is shortest
  worst_lat <- max(abs(lat_r + c(-pad_lat, pad_lat)))
  pad_lon <- pad_km / (111.32 * max(cos(worst_lat * pi / 180), 0.05)) + clock_margin_deg
  # a fit cannot localise the posterior better than its own cell size
  pad_lat <- pad_lat + cell_margin
  pad_lon <- pad_lon + cell_margin

  out <- c(lon_r[1] - pad_lon, lon_r[2] + pad_lon,
           lat_r[1] - pad_lat, lat_r[2] + pad_lat)
  # never expand beyond the domain that was actually searched
  out <- c(max(out[1], min(ulon) - 0.5), min(out[2], max(ulon) + 0.5),
           max(out[3], min(ulat) - 0.5), min(out[4], max(ulat) + 0.5))
  names(out) <- c("lon_min", "lon_max", "lat_min", "lat_max")
  attr(out, "retained") <- ((out[2] - out[1]) * (out[4] - out[3])) /
    ((diff(range(ulon)) + 1) * (diff(range(ulat)) + 1))
  attr(out, "epsilon") <- epsilon
  attr(out, "pad_km") <- pad_km
  attr(out, "pad_deg") <- c(lon = pad_lon, lat = pad_lat)
  out
}


#' Fit a grid HMM by coarse search then refinement
#'
#' Runs [TwilightFreeGrid()] twice: once on a generous domain at a coarse cell size
#' to find where the posterior lives, then again at the cell size you want, on a box
#' drawn around it. The domain is a property of the data and is rarely known in
#' advance, so it is cheaper to measure it than to guess a wide one and pay for it at
#' full resolution.
#'
#' @details
#' The saving comes from both halves of the engine: the emission costs
#' `cells x observations` per knot and the forward-backward costs
#' `cells x neighbours`, so removing cells helps both. On elephant-seal tracks over a
#' 100 x 50 degree domain the refined box is typically a third to a quarter of the
#' area, for a three- to four-fold saving.
#'
#' **`log_z` is not comparable across different domains.** The evidence is a sum over
#' cells; change which cells exist and it changes. Compare `log_z` only between fits
#' sharing a grid, exactly as with `area_correction`.
#'
#' @param ... Arguments for [TwilightFreeGrid()], including the data. `grid` is
#'   supplied by this function and must not be passed.
#' @param lon,lat Length-2 search bounds for the COARSE pass, `c(min, max)`.
#' @param cell.size Cell size for the FINE fit.
#' @param coarse.size Cell size for the coarse pass (default 4x `cell.size`). The
#'   coarse fit only has to locate the posterior, not resolve it.
#' @param coarse_inflate Factor by which the coarse pass's `diffusion` is inflated
#'   (default 3). This is not a tuning knob, it is what makes the coarse pass safe.
#'   A coarse grid does NOT blur the posterior -- it sharpens it, because widely
#'   spaced cells differ so much in day length that one dominates its neighbours
#'   outright. Measured on a real track, a 4-degree coarse fit put its entire
#'   latitude posterior in a SINGLE cell (38.0 to 38.0) while the 1-degree fit spread
#'   over 24.5 to 50.5, so the coarse box clipped the answer by 2.5 degrees. Inflating
#'   the movement scale broadens the coarse posterior into an honest envelope; the
#'   coarse pass has to bound the region, not resolve it.
#' @param epsilon,n_sigma,clock_margin_deg Passed to [posterior_extent()].
#' @param mask Passed to [makeGrid()] for both passes.
#' @param verbose Report the box and the saving (default `TRUE`).
#' @return The FINE `TwilightFreeGrid` fit, with attributes `extent` (the refined box)
#'   and `coarse` (the coarse fit, for inspection).
#' @seealso [posterior_extent()], [TwilightFreeGrid()], [makeGrid()]
#' @export
TwilightFreeGridRefine <- function(..., lon, lat, cell.size, coarse.size = 4 * cell.size,
                                   coarse_inflate = 3, epsilon = 1e-6, n_sigma = 5,
                                   clock_margin_deg = 0,
                                   mask = c("none", "sea", "land"), verbose = TRUE) {
  mask <- match.arg(mask)
  dots <- list(...)
  if (!is.null(dots$grid)) stop("`grid` is built by TwilightFreeGridRefine(); pass `lon`, `lat` and `cell.size` instead")
  if (coarse.size < cell.size)
    stop("`coarse.size` must be at least `cell.size`; the coarse pass exists to be cheap")

  if (coarse_inflate < 1) stop("`coarse_inflate` must be at least 1")
  g0 <- makeGrid(lon = lon, lat = lat, cell.size = coarse.size, mask = mask)
  d0 <- dots
  d0$diffusion <- (if (is.null(dots$diffusion)) formals(TwilightFreeGrid)$diffusion
                   else dots$diffusion) * coarse_inflate
  if (verbose) message(sprintf("coarse pass at %g degrees, diffusion inflated %gx",
                               coarse.size, coarse_inflate))
  f0 <- do.call(TwilightFreeGrid, c(d0, list(grid = g0)))

  # size the margin from what the COARSE fit actually used, not from `dots`, so a
  # caller relying on TwilightFreeGrid's own defaults still gets a correct margin
  # size the margin on the FINE fit's movement scale, not the inflated coarse one,
  # and pad by a full coarse cell for the coarse pass's own resolution
  ext <- posterior_extent(f0, epsilon = epsilon, n_sigma = n_sigma,
                          diffusion = if (is.null(dots$diffusion))
                              eval(formals(TwilightFreeGrid)$diffusion) else dots$diffusion,
                          step_hours = if (is.null(dots$step_hours))
                              eval(formals(TwilightFreeGrid)$step_hours) else dots$step_hours,
                          clock_margin_deg = clock_margin_deg,
                          cell_margin = coarse.size)
  # a refined box smaller than a few fine cells means the coarse pass collapsed;
  # fall back to the original domain rather than fit in a sliver
  if ((ext[2] - ext[1]) < 3 * cell.size || (ext[4] - ext[3]) < 3 * cell.size) {
    warning("refined extent is smaller than 3 fine cells; falling back to the full domain")
    ext <- c(lon[1], lon[2], lat[1], lat[2])
  }
  if (verbose)
    message(sprintf("refined box: lon %.1f-%.1f, lat %.1f-%.1f (%.0f%% of the area, margin %.0f km)",
                    ext[1], ext[2], ext[3], ext[4], 100 * attr(ext, "retained"),
                    attr(ext, "pad_km")))

  g1 <- makeGrid(lon = c(ext[1], ext[2]), lat = c(ext[3], ext[4]),
                 cell.size = cell.size, mask = mask)
  f1 <- do.call(TwilightFreeGrid, c(dots, list(grid = g1)))
  attr(f1, "extent") <- ext
  attr(f1, "coarse") <- f0
  f1
}
