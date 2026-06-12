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
