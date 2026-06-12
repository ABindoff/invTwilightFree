# Sensor-fusion plumbing: assemble a list of location_term() objects into the
# additive log-likelihood matrix that TwilightFreeGrid() hands to the Rust HMM.
# A "term" is (tag channel -> per-knot value) x (environmental field) x (rule).
# Priors and masks have no tag channel. See HANDOFF_2_SENSOR_FUSION.md.

#' Define an auxiliary location term for sensor fusion
#'
#' Bundles a tag channel, an environmental source, and a rule into one additive
#' constraint on candidate locations, for use via the `terms` argument of
#' [TwilightFreeGrid()]. Each term contributes a log-likelihood at every grid
#' cell and knot, added to the light likelihood. Priors and masks omit `tag`.
#'
#' @param name Short identifier for the term (character).
#' @param tag Optional `data.frame` with columns `time` (POSIXct) and `value`
#'   (numeric) giving the tag sensor channel (temperature, depth, ...). `NULL`
#'   for priors and masks, whose source needs no observation.
#' @param reduce How to reduce the tag values falling in each knot interval to a
#'   single number: one of `"median"`, `"mean"`, `"max"`, `"min"`, `"last"`, or
#'   a function of a numeric vector. Ignored when `tag` is `NULL`.
#' @param source A `tf_source` (e.g. [sst_source()], [bathy_source()],
#'   [hemisphere_prior()]) or a bare `function(lon, lat, date)` returning the
#'   expected field value at each cell.
#' @param rule A `tf_rule` (e.g. [student_rule()], [floor_rule()],
#'   [identity_rule()]) or a bare `function(obs, expected)` returning a per-cell
#'   log-likelihood vector.
#' @param weight Multiplier on this term's contribution (default `1`). Use to
#'   up- or down-weight a sensor relative to the light likelihood.
#' @param missing Contribution for a knot whose reduced tag value is `NA` (no
#'   reading in that interval). Default `0`, i.e. the term adds no information
#'   there rather than penalising.
#'
#' @return An object of class `tf_term`.
#' @seealso [TwilightFreeGrid()], [hemisphere_prior()], [sst_source()],
#'   [bathy_source()]
#' @examples
#' season <- function(d) "S"
#' location_term("hemisphere", source = hemisphere_prior(season), rule = identity_rule())
#' @export
location_term <- function(name, tag = NULL, reduce = "median", source, rule,
                          weight = 1, missing = 0) {
  if (missing(source) || !is.function(source)) {
    stop("`source` must be a tf_source or a function(lon, lat, date)")
  }
  if (missing(rule) || !is.function(rule)) {
    stop("`rule` must be a tf_rule or a function(obs, expected)")
  }
  if (!is.null(tag)) {
    if (!is.data.frame(tag) || !all(c("time", "value") %in% names(tag))) {
      stop("`tag` must be a data.frame with columns 'time' and 'value'")
    }
  }
  structure(
    list(name = name, tag = tag, reduce = reduce, source = source,
         rule = rule, weight = weight, missing = missing),
    class = "tf_term"
  )
}


# Resolve a reduce specification to a function of a numeric vector.
get_reducer <- function(reduce) {
  if (is.function(reduce)) return(reduce)
  switch(
    reduce,
    median = function(x) stats::median(x, na.rm = TRUE),
    mean   = function(x) mean(x, na.rm = TRUE),
    max    = function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE),
    min    = function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE),
    last   = function(x) x[length(x)],
    stop("unknown `reduce`: ", reduce)
  )
}


# Lower bound of each knot interval, matching run_grid_hmm's convention:
# t_prev[k] = knot_times[k-1], and for the first knot the interval has the same
# width as the first step. Observations are assigned by (t_prev, t_curr].
knot_lower_bounds <- function(knot_times) {
  K <- length(knot_times)
  lower <- numeric(K)
  step1 <- if (K > 1) knot_times[2] - knot_times[1] else 1
  for (k in seq_len(K)) {
    lower[k] <- if (k == 1) knot_times[1] - step1 else knot_times[k - 1]
  }
  lower
}


# Reduce a term's tag channel to one value per knot (NA where no reading).
reduce_tag <- function(term, knot_times) {
  K <- length(knot_times)
  if (is.null(term$tag)) return(rep(NA_real_, K))
  tt <- as.numeric(term$tag$time)
  vv <- as.numeric(term$tag$value)
  lower <- knot_lower_bounds(knot_times)
  redfun <- get_reducer(term$reduce)
  out <- rep(NA_real_, K)
  for (k in seq_len(K)) {
    sel <- tt > lower[k] & tt <= knot_times[k]
    if (any(sel)) out[k] <- redfun(vv[sel])
  }
  out
}


#' Assemble the auxiliary log-likelihood matrix for a set of location terms
#'
#' Internal helper used by [TwilightFreeGrid()]. For each term it reduces the
#' tag channel per knot, evaluates the environmental source at every grid cell
#' (cached per calendar date to avoid refetching), applies the rule, and sums
#' the weighted contributions. Knots with a missing tag reading take the term's
#' `missing` value and skip the (potentially expensive) source evaluation.
#'
#' @param terms List of `tf_term` objects.
#' @param lon,lat Numeric vectors of grid-cell coordinates (length `n`).
#' @param knot_times Numeric vector of knot times (Unix seconds, length `K`).
#' @return A `K` by `n` numeric matrix of summed log-likelihood contributions.
#' @keywords internal
#' @export
build_aux_matrix <- function(terms, lon, lat, knot_times) {
  K <- length(knot_times)
  n <- length(lon)
  aux <- matrix(0.0, nrow = K, ncol = n)
  if (length(terms) == 0) return(aux)

  knot_dates <- as.Date(as.POSIXct(knot_times, origin = "1970-01-01", tz = "UTC"))

  for (term in terms) {
    if (!inherits(term, "tf_term")) stop("each element of `terms` must be a location_term()")
    obs_by_knot <- reduce_tag(term, knot_times)
    needs_obs <- !is.null(term$tag)
    cache <- list()                                  # date string -> expected field vector

    for (k in seq_len(K)) {
      obs_k <- obs_by_knot[k]
      if (needs_obs && is.na(obs_k)) {
        aux[k, ] <- aux[k, ] + term$weight * term$missing
        next
      }
      key <- as.character(knot_dates[k])
      expected <- cache[[key]]
      if (is.null(expected)) {
        expected <- term$source(lon, lat, knot_dates[k])
        cache[[key]] <- expected
      }
      ll <- term$rule(obs_k, expected)
      aux[k, ] <- aux[k, ] + term$weight * ll
    }
  }
  aux
}
