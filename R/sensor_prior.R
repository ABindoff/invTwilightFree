#' Hemisphere (seasonal latitude) prior
#'
#' Builds a location-prior source that softly (or hard) confines a track to one
#' hemisphere as a function of date. This is the cheapest fix for the equinox
#' latitude ambiguity ("hemisphere swapping"), where day length alone cannot
#' distinguish north from south and the light likelihood is genuinely bimodal.
#'
#' The returned object is a source closure consumed by `location_term()` with
#' `identity_rule()`. For each knot it returns, per candidate grid cell, a
#' log-prior density of `0` (i.e. `log(1)`) for cells in the allowed hemisphere
#' and `log(softness)` for cells in the disallowed hemisphere. These values are
#' added to the light log-likelihood of each cell inside the grid HMM.
#'
#' @param hemisphere_by_date A function taking a single date (`Date` or
#'   `POSIXct`) and returning one of `"N"`, `"S"`, or `"both"`. `"N"` favours
#'   the northern side (latitude >= `boundary`), `"S"` the southern side
#'   (latitude <= `boundary`), and `"both"` applies no constraint for that knot
#'   (e.g. during a migration window when the animal may be crossing).
#' @param softness Relative prior weight on the disallowed hemisphere, in
#'   `[0, 1]`. `1` means no penalty; `0` is a hard cut (`-Inf`). The default
#'   `1e-3` downweights the wrong hemisphere by about 6.9 log-units: strong
#'   enough to break an equinox tie, soft enough to be overridden by decisive
#'   light data.
#' @param boundary Latitude (degrees) of the dividing line, default `0` (the
#'   equator). Use e.g. `boundary = 10` with `"S"` to allow a little slack north
#'   of the equator.
#'
#' @return A source closure of class `tf_source`: `function(lon, lat, date)`
#'   returning a numeric vector of per-cell log-prior densities, where `lon` and
#'   `lat` are equal-length vectors of candidate cell coordinates and `date` is
#'   a single knot date.
#'
#' @seealso [identity_rule()]
#' @examples
#' # Austral breeder: confined south of the equator outside Nov-Feb,
#' # unconstrained during the summer breeding window.
#' season <- function(d) {
#'   if (as.integer(format(d, "%m")) %in% c(11, 12, 1, 2)) "both" else "S"
#' }
#' src <- hemisphere_prior(season, softness = 1e-3)
#'
#' # Evaluate on a tiny lon/lat grid for one winter date.
#' g <- expand.grid(lon = c(140, 150), lat = c(-45, -10, 30))
#' src(g$lon, g$lat, as.Date("2024-06-21"))
#' @export
hemisphere_prior <- function(hemisphere_by_date, softness = 1e-3, boundary = 0) {
  stopifnot(is.function(hemisphere_by_date))
  if (length(softness) != 1L || is.na(softness) || softness < 0 || softness > 1) {
    stop("`softness` must be a single value in [0, 1]")
  }
  log_pen <- if (softness <= 0) -Inf else log(softness)

  src <- function(lon, lat, date) {
    side <- hemisphere_by_date(date)
    if (length(side) != 1L || !side %in% c("N", "S", "both")) {
      stop("`hemisphere_by_date` must return one of 'N', 'S', or 'both'")
    }
    out <- numeric(length(lat))                 # 0 = log(1), the allowed value
    if (identical(side, "both")) return(out)
    allowed <- if (side == "N") lat >= boundary else lat <= boundary
    out[!allowed] <- log_pen
    out
  }
  class(src) <- c("tf_source", "function")
  src
}


#' Spherical area (stationary) prior over grid cells (deprecated)
#'
#' Deprecated. Use `TwilightFreeGrid(area_correction = TRUE)`, the default,
#' instead. This function is retained only so that existing scripts keep running,
#' and warns on use.
#'
#' Returns the log of each cell's relative area, `log(cos(latitude))`, as a
#' location-prior source: a lon/lat grid has equal angular cells but unequal area
#' (cells shrink poleward by `cos(lat)`), so this favours equatorward cells in
#' proportion to the area they represent.
#'
#' Why it is deprecated: routed through [identity_rule()], this is *mathematically
#' identical* to the engine's `area_correction`, which applies the same
#' `log(cos(lat))` to the same destination cell in the initial, forward and
#' backward passes. Verified bit-for-bit: the two routes return identical
#' latitudes and identical `log_z`. Two names for one operator invited
#' double-application, and because `area_correction` now defaults to `TRUE`,
#' supplying this term as well applies the factor **twice**. Passing both is an
#' error; see [TwilightFreeGrid()].
#'
#' Correction to earlier documentation: this term was previously described here as
#' "far too weak to bridge the planar-vs-spherical movement gap", citing an SBC
#' experiment in `notes/topology/sbc_design.md`. That conclusion does not
#' generalise. It was measured on a design that makes latitude easy -- 41 daily
#' knots with **both** endpoints anchored, at 50 degrees south in strong austral
#' winter light -- where a per-knot log-prior tilt moves the posterior mean very
#' little. On a 12-hourly track of 273-494 knots with a free retrieval endpoint
#' and near-equinox segments, the identical term is worth 0.4 to 2.2 degrees of
#' latitude bias. The effect scales with knot count, endpoint anchoring and how
#' flat the per-knot likelihood is, so a null result from a sharply identified
#' design says nothing about a weakly identified one.
#'
#' What remains true is the narrower original point: this is the STATIONARY part
#' of the spherical area measure only. Isotropic Brownian motion on a sphere also
#' carries an equatorward drift in its transition (a `tan(lat)` term) that no
#' per-cell prior can reproduce, and supplying the area factor does not by itself
#' make the SBC rank shape for the spherical generator uniform.
#'
#' @return A source closure of class `tf_source`: `function(lon, lat, date)`
#'   returning `log(cos(lat))` per cell (`date` is ignored). Carries the attribute
#'   `tf_area_prior = TRUE` so [TwilightFreeGrid()] can detect the double-count.
#' @seealso [TwilightFreeGrid()] and its `area_correction` argument.
#' @examples
#' suppressWarnings(area_prior()(lon = c(0, 0), lat = c(-60, -30)))
#' @export
area_prior <- function() {
  warning("`area_prior()` is deprecated: the engine now applies the cell-area ",
          "factor itself via `TwilightFreeGrid(area_correction = TRUE)` (the ",
          "default), and the two are the same operator. Passing both applies it ",
          "twice. Drop the term.", call. = FALSE)
  src <- function(lon, lat, date = NULL) {
    log(pmax(cos(lat * pi / 180), 1e-6))
  }
  class(src) <- c("tf_source", "function")
  attr(src, "tf_area_prior") <- TRUE
  src
}


#' Identity rule: treat the source field as a (log-)prior density
#'
#' A pass-through rule for `location_term()`s whose source already encodes the
#' contribution directly, such as a prior surface or [hemisphere_prior()]. No
#' tag observation is used, so `obs` is ignored.
#'
#' @param is_log If `TRUE` (default) the source values are already log-densities
#'   and are returned unchanged. If `FALSE` the source values are treated as
#'   non-negative densities and logged (zeros map to `-Inf`).
#'
#' @return A rule closure of class `tf_rule`: `function(obs, expected)`
#'   returning a numeric vector the same length as `expected`. Cells with `NA`
#'   in `expected` (for example outside a prior raster's extent) contribute `0`,
#'   so an incomplete prior is treated as uninformative there rather than
#'   impossible.
#'
#' @seealso [hemisphere_prior()]
#' @examples
#' r <- identity_rule()
#' r(obs = NULL, expected = c(0, log(1e-3), NA))      # 0.0, -6.9, 0.0
#'
#' rd <- identity_rule(is_log = FALSE)
#' rd(obs = NULL, expected = c(1, exp(1), 0))         # 0, 1, -Inf
#' @export
identity_rule <- function(is_log = TRUE) {
  rule <- function(obs, expected) {
    val <- if (isTRUE(is_log)) {
      as.numeric(expected)
    } else {
      e <- as.numeric(expected)
      ifelse(e > 0, log(e), -Inf)
    }
    val[is.na(val)] <- 0                          # undefined -> uninformative
    val
  }
  class(rule) <- c("tf_rule", "function")
  rule
}
