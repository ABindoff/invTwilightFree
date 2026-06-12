#' Floor rule: a one-sided constraint from a bounding surface
#'
#' For sensors where the observed value is a *lower bound* on the expected field
#' at the true location. The motivating case is dive depth versus bathymetry: an
#' animal that dived to `obs` metres must have been somewhere the seabed is at
#' least that deep, so candidate cells shallower than the dive are penalised,
#' while cells that are deeper are not (the animal need not reach the bottom).
#'
#' Conventions: both `obs` and `expected` are positive-down depths in the same
#' units (metres). A cell is *too shallow* when `obs > expected`; the shortfall
#' is `deficit = max(0, obs - expected)`.
#'
#' Generalisation: this single rule expresses any "the solid surface must lie
#' beyond the animal's extreme excursion" constraint, not just diving. If a tag
#' records altitude, the same rule rules out terrain that pokes above the animal
#' (a tag at sea level cannot be over a mountain). Keep everything positive-down:
#' express the surface and the observation as depth below sea level, so a
#' +2000 m summit is `-2000`, a +50 m flight altitude is `-50`, and a 1200 m
#' seabed is `+1200`. Under that convention the rule is always "surface depth
#' `expected` is at least the animal's least depth `obs`". This is why there is
#' no separate "ceiling" rule: an upper-bound constraint is this rule with
#' negated inputs.
#'
#' Pair with [bathy_source()] (diving) or a terrain-elevation source expressed
#' positive-down (altimetry). For a clean land exclusion, also add a sea mask;
#' the floor rule only constrains cells via depth, so a flat-light day with no
#' deep dives leaves shallow cells unpenalised.
#'
#' @param sd Tolerance (same units as the depths) for the soft constraint. The
#'   penalty is `-0.5 * (deficit / sd)^2`, a half-normal on the shortfall. `sd`
#'   should absorb bathymetry-grid error and depth-sensor error. If `NULL`
#'   (default) the constraint is hard: `0` where the cell is deep enough and
#'   `-Inf` where it is too shallow.
#'
#' @return A rule closure of class `tf_rule`: `function(obs, expected)` where
#'   `obs` is the single per-knot reduced value (e.g. the knot's maximum dive
#'   depth) and `expected` is the length-`n` vector of field values at the
#'   candidate cells. Returns a numeric vector of per-cell log-likelihood
#'   contributions. A missing `obs` (`NA` or length 0) yields all zeros (no
#'   information), and cells with `NA` in `expected` (no bathymetry there)
#'   contribute `0` rather than being excluded.
#'
#' @seealso [bathy_source()], [identity_rule()]
#' @examples
#' # Hard constraint: a 300 m dive rules out cells shallower than 300 m.
#' fr <- floor_rule()
#' fr(obs = 300, expected = c(50, 300, 1200))      # -Inf, 0, 0
#'
#' # Soft constraint with 200 m tolerance.
#' fs <- floor_rule(sd = 200)
#' round(fs(obs = 300, expected = c(100, 250, 1000)), 3)
#'
#' # Altimetry, positive-down: a tag at sea level (obs = 0) over terrain whose
#' # elevation expressed as depth is c(summit -2000, plain -100, shore 0, sea 1200).
#' floor_rule()(obs = 0, expected = c(-2000, -100, 0, 1200))   # -Inf, -Inf, 0, 0
#' @export
floor_rule <- function(sd = NULL) {
  if (!is.null(sd) && (length(sd) != 1L || is.na(sd) || sd <= 0)) {
    stop("`sd` must be NULL (hard constraint) or a single positive value")
  }
  rule <- function(obs, expected) {
    expected <- as.numeric(expected)
    if (length(obs) == 0L || is.na(obs[1L])) {
      return(rep(0, length(expected)))            # no reading -> no information
    }
    deficit <- pmax(0, obs[1L] - expected)        # how much too shallow the cell is
    out <- if (is.null(sd)) {
      ifelse(deficit > 0, -Inf, 0)                # hard: cell must clear the dive
    } else {
      -0.5 * (deficit / sd)^2                     # soft: half-normal on shortfall
    }
    out[is.na(expected)] <- 0                     # no bathymetry data -> uninformative
    out
  }
  class(rule) <- c("tf_rule", "function")
  rule
}


#' Student-t matching rule: robust comparison of a tag reading to a field
#'
#' A heavy-tailed alternative to a Gaussian match, for sensors whose value
#' should equal the environmental field at the true location but is noisy. The
#' motivating case is tag-recorded temperature versus remotely sensed sea
#' surface temperature (SST): the residual `(obs - expected)` is modelled as a
#' scaled Student-t, so a few badly wrong readings (a deep dive logged as a
#' "surface" temperature, a frontal mismatch, a bad pixel) do not collapse the
#' likelihood the way squared error would. This matters because logged tag
#' temperature is noisy and sampled sparsely near the surface.
#'
#' The Gaussian match is the `df -> Inf` limit of this rule, so there is no
#' separate gaussian_rule; use a large `df` (say 30) if you want near-Gaussian
#' behaviour.
#'
#' @param sd Scale of the residual in the field's units (e.g. degrees C). Sets
#'   how far an observed temperature can sit from the cell's SST before the cell
#'   is downweighted. Should reflect combined tag-sensor noise and SST-product
#'   error, typically order 1 to 2 degrees C for surface temperature.
#' @param df Degrees of freedom of the Student-t. Lower is more robust (heavier
#'   tails); default `4` gives finite variance with markedly fatter tails than a
#'   normal.
#'
#' @return A rule closure of class `tf_rule`: `function(obs, expected)` where
#'   `obs` is the single per-knot reduced tag value and `expected` is the
#'   length-`n` vector of field values at the candidate cells. Returns per-cell
#'   log-likelihood contributions. A missing `obs` (`NA` or length 0) yields all
#'   zeros, and cells with `NA` in `expected` contribute `0`, so a knot with no
#'   surface reading or a cell with no SST simply adds no information.
#'
#' @seealso [sst_source()]
#' @examples
#' sr <- student_rule(sd = 1.5, df = 4)
#' # Tag reads 12 C; cells whose SST is 12, 13, 18 C.
#' round(sr(obs = 12, expected = c(12, 13, 18)), 3)
#' @export
student_rule <- function(sd, df = 4) {
  if (length(sd) != 1L || is.na(sd) || sd <= 0) {
    stop("`sd` must be a single positive value")
  }
  if (length(df) != 1L || is.na(df) || df <= 0) {
    stop("`df` must be a single positive value")
  }
  rule <- function(obs, expected) {
    expected <- as.numeric(expected)
    if (length(obs) == 0L || is.na(obs[1L])) {
      return(rep(0, length(expected)))            # no reading -> no information
    }
    resid <- (obs[1L] - expected) / sd
    out <- stats::dt(resid, df = df, log = TRUE) - log(sd)
    out[is.na(expected)] <- 0                     # no field value -> uninformative
    out
  }
  class(rule) <- c("tf_rule", "function")
  rule
}


#' Gaussian matching rule: symmetric comparison of a tag reading to a field
#'
#' A standard squared-error match, for sensors whose value should equal the
#' environmental field at the true location. This is the `df -> Inf` limit of
#' [student_rule()]; prefer `student_rule()` with a large `df` when in doubt,
#' since heavy tails are almost always safer for noisy tag sensors.
#'
#' @param sd Scale of the residual in the field's units. Sets how far an
#'   observed value can sit from the expected value before the cell is
#'   penalised.
#'
#' @return A rule closure of class `tf_rule`: `function(obs, expected)` returning
#'   per-cell log-likelihood contributions. Missing `obs` yields all zeros; cells
#'   with `NA` in `expected` contribute `0`.
#'
#' @seealso [student_rule()]
#' @examples
#' gr <- gaussian_rule(sd = 2)
#' round(gr(obs = 10, expected = c(10, 12, 6)), 3)
#' @export
gaussian_rule <- function(sd) {
  if (length(sd) != 1L || is.na(sd) || sd <= 0) {
    stop("`sd` must be a single positive value")
  }
  rule <- function(obs, expected) {
    expected <- as.numeric(expected)
    if (length(obs) == 0L || is.na(obs[1L])) {
      return(rep(0, length(expected)))
    }
    out <- -0.5 * ((obs[1L] - expected) / sd)^2 - log(sd)
    out[is.na(expected)] <- 0
    out
  }
  class(rule) <- c("tf_rule", "function")
  rule
}


#' Ceiling rule: a one-sided constraint from an upper bounding surface
#'
#' Mirror of [floor_rule()]: for sensors where the observed value is an
#' *upper bound* on the expected field at the true location. A cell is infeasible
#' when `expected > obs` (i.e. the field value pokes above the observation).
#' The deficit is `pmax(0, expected - obs)`.
#'
#' @param sd Tolerance for the soft constraint (`NULL` for a hard cut).
#'
#' @return A rule closure of class `tf_rule`. Missing `obs` or `NA` in `expected`
#'   contribute `0` (uninformative).
#'
#' @seealso [floor_rule()]
#' @examples
#' cr <- ceiling_rule()
#' cr(obs = 100, expected = c(50, 100, 200))  # 0, 0, -Inf
#' @export
ceiling_rule <- function(sd = NULL) {
  if (!is.null(sd) && (length(sd) != 1L || is.na(sd) || sd <= 0)) {
    stop("`sd` must be NULL (hard constraint) or a single positive value")
  }
  rule <- function(obs, expected) {
    expected <- as.numeric(expected)
    if (length(obs) == 0L || is.na(obs[1L])) {
      return(rep(0, length(expected)))
    }
    deficit <- pmax(0, expected - obs[1L])
    out <- if (is.null(sd)) {
      ifelse(deficit > 0, -Inf, 0)
    } else {
      -0.5 * (deficit / sd)^2
    }
    out[is.na(expected)] <- 0
    out
  }
  class(rule) <- c("tf_rule", "function")
  rule
}


#' Mask rule: exclude cells where the field is zero (or NA)
#'
#' Assigns `0` to cells where `expected != 0` and is not `NA` (i.e. the cell
#' is in the permitted region), and `penalty` (default `-Inf`) elsewhere. Pair
#' with [sea_mask_source()] to hard-exclude land cells, or any binary raster
#' where `1 = allowed`, `0 = forbidden`.
#'
#' @param penalty Log-likelihood penalty for forbidden cells. Default `-Inf`
#'   (hard exclusion). Use a large negative finite value (e.g. `-20`) for a
#'   soft mask that can still be overridden by very strong light data.
#'
#' @return A rule closure of class `tf_rule`. `obs` is ignored (masks need no
#'   tag channel). Cells with `NA` in `expected` are treated as allowed (`0`),
#'   so an incomplete mask raster is uninformative outside its extent.
#'
#' @seealso [sea_mask_source()]
#' @examples
#' mr <- mask_rule()
#' mr(obs = NULL, expected = c(1, 0, NA))   # 0, -Inf, 0
#'
#' mr_soft <- mask_rule(penalty = -20)
#' mr_soft(obs = NULL, expected = c(1, 0))  # 0, -20
#' @export
mask_rule <- function(penalty = -Inf) {
  if (length(penalty) != 1L || is.na(penalty)) {
    stop("`penalty` must be a single numeric value (use -Inf for a hard mask)")
  }
  rule <- function(obs, expected) {
    expected <- as.numeric(expected)
    out <- ifelse(expected != 0 & !is.na(expected), 0, penalty)
    out[is.na(expected)] <- 0
    out
  }
  class(rule) <- c("tf_rule", "function")
  rule
}


#' Custom rule: bring your own log-likelihood function
#'
#' An escape hatch for any sensor or constraint not covered by the built-in
#' rules. Supply a function `fn(obs, expected)` that accepts a scalar `obs`
#' (the per-knot reduced tag value, possibly `NA`) and a length-`n` numeric
#' vector `expected` (field values at candidate cells), and returns a length-`n`
#' numeric vector of log-likelihood contributions.
#'
#' @param fn A function `function(obs, expected)` returning a numeric vector of
#'   the same length as `expected`. It is called once per knot with a single
#'   scalar `obs` value.
#'
#' @return A rule closure of class `tf_rule`.
#'
#' @examples
#' # A simple absolute-deviation rule.
#' fn <- function(obs, expected) -abs(obs - expected)
#' cr <- custom_rule(fn)
#' cr(obs = 5, expected = c(4, 5, 8))
#' @export
custom_rule <- function(fn) {
  if (!is.function(fn)) stop("`fn` must be a function(obs, expected)")
  rule <- function(obs, expected) fn(obs, expected)
  class(rule) <- c("tf_rule", "function")
  rule
}
