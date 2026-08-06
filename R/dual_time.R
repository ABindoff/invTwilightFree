# Dual (time-conditioned) primitives: P(time | location, light).
#
# The location grid HMM grids space and conditions on the known clock, giving
# P(location | time). These flip it: hold a location fixed and read the light
# likelihood as a function of a time offset applied to the clock. Because the
# light's temporal edges (twilight, noon) are sharp, this direction is far better
# conditioned than latitude, so it pins time/longitude precisely. Anchored at the
# known deployment and retrieval fixes (location AND time), it calibrates the tag
# clock, which removes the longitude gauge across the whole interior.

#' Dual light log-likelihood at a fixed location over a grid of time offsets
#'
#' For a single candidate location, evaluates the spike-and-slab light
#' log-likelihood as a function of a time offset `tau` added to the observation
#' clock. This is the transpose of [eval_logpk_grid()]: location fixed, time
#' varied. A sharp peak in `tau` locates the solar time the light implies at this
#' location; against the known clock it is the clock error (and, equivalently,
#' the longitude line of position expressed in time).
#'
#' @param lon,lat Scalar candidate location (decimal degrees).
#' @param times Observation times (POSIXct or numeric Unix seconds).
#' @param light Observed light values (same length as `times`).
#' @param offsets Time offsets to evaluate, in seconds. Default a grid of
#'   +/- 3 hours at 2-minute steps.
#' @param calibration `c(intercept, slope)` mapping zenith to expected light.
#' @param likelihood_params `c(lambda, max_light, prob_slab)`.
#' @param shade_ratio Ratio of the spike's upper-arm rate to its shading rate
#'   (default 2), matching the engines.
#'
#' @return A data frame with columns `offset` (seconds) and `loglik`.
#' @seealso [calibrate_clock()], [eval_logpk_grid()]
#' @examples
#' \dontrun{
#' d <- eval_logpt_loc(150, -50, track$time, track$light,
#'                     calibration = c(64, 64/90), likelihood_params = c(0.5, 64, 0.05))
#' d$offset[which.max(d$loglik)]   # clock error implied at (150, -50), seconds
#' }
#' @export
eval_logpt_loc <- function(lon, lat, times, light,
                           offsets = seq(-3 * 3600, 3 * 3600, by = 120),
                           calibration, likelihood_params, shade_ratio = 2) {
  tt <- as.numeric(times)
  intercept <- calibration[1]; slope <- calibration[2]
  lambda <- likelihood_params[1]; max_light <- likelihood_params[2]
  prob_slab <- if (length(likelihood_params) >= 3) likelihood_params[3] else 0.05
  slab <- 1 / max_light
  n <- length(tt)
  lonv <- rep(lon, n); latv <- rep(lat, n)

  # Normalised over [0, max_light], matching the engine's spike_density(): the
  # normaliser depends on the expected light and therefore on the offset, so
  # omitting it changes the shape of the profile, not just its level.
  lam_hi <- lambda * shade_ratio
  normaliser <- function(mu) {
    pmax(1e-12, (1 - exp(-lambda * mu)) +
                (lambda / lam_hi) * (1 - exp(-lam_hi * (max_light - mu))))
  }
  loglik <- vapply(offsets, function(tau) {
    z <- solar_zenith(tt + tau, lonv, latv)
    expected <- .tf_expected_light(z, calibration, max_light)
    raw <- ifelse(light <= expected,
                  lambda * exp(-lambda * (expected - light)),
                  lambda * exp(-lam_hi * (light - expected)))       # engine convention
    den <- (1 - prob_slab) * raw / normaliser(expected) + prob_slab * slab
    sum(log(den))
  }, numeric(1))

  data.frame(offset = offsets, loglik = loglik)
}

# Parabolic-refined argmax of a (offset, loglik) curve -> (offset_hat, curvature).
.peak_refine <- function(d) {
  i <- which.max(d$loglik)
  if (i == 1L || i == nrow(d)) {
    return(list(offset = d$offset[i], curvature = NA_real_))
  }
  x <- d$offset[(i - 1):(i + 1)]; y <- d$loglik[(i - 1):(i + 1)]
  # vertex of the quadratic through three points
  denom <- (y[1] - 2 * y[2] + y[3])
  off <- if (denom != 0) x[2] - 0.5 * (x[3] - x[1]) * (y[3] - y[1]) / (2 * denom) else x[2]
  h <- x[2] - x[1]
  curv <- denom / (h * h)                 # second derivative (negative at a max)
  list(offset = off, curvature = curv)
}

#' Calibrate the tag clock from the known deployment and retrieval fixes
#'
#' Exploits that `deployed.at` and `retrieved.at` are known in both location and
#' time. At each end the animal's location is known, so the dual likelihood
#' ([eval_logpt_loc()]) over a window of light there is sharply peaked at the
#' tag's clock error. Two fixes give a linear drift model `tau(t)`, whose removal
#' makes longitude essentially exact across the interior (it eliminates the
#' longitude holonomy globally rather than per knot).
#'
#' @param deploy,retrieve Lists with `lon`, `lat`, and `time` (POSIXct) for the
#'   known release and recovery fixes.
#' @param times,light Observation times (POSIXct) and light values.
#' @param calibration,likelihood_params As in [eval_logpt_loc()].
#' @param window_days Days of light around each fix used for calibration.
#' @param search Half-width of the offset search, seconds (default 2 h).
#' @param by Offset step, seconds (default 60).
#'
#' @return A list with the per-end offset estimates (`offset_deploy`,
#'   `offset_retrieve`, seconds), the linear drift coefficients, and `correct()`,
#'   a function mapping observed POSIXct times to clock-corrected times.
#' @seealso [eval_logpt_loc()]
#' @export
calibrate_clock <- function(deploy, retrieve, times, light,
                            calibration, likelihood_params,
                            window_days = 3, search = 2 * 3600, by = 60) {
  tt <- as.numeric(times)
  offsets <- seq(-search, search, by = by)
  td <- as.numeric(deploy$time); tr <- as.numeric(retrieve$time)

  end_offset <- function(fix, t_centre) {
    keep <- abs(tt - t_centre) <= window_days * 86400
    if (sum(keep) < 5) stop("too few observations within window_days of a fix")
    d <- eval_logpt_loc(fix$lon, fix$lat, tt[keep], light[keep],
                        offsets, calibration, likelihood_params)
    .peak_refine(d)
  }

  od <- end_offset(deploy, td)
  or <- end_offset(retrieve, tr)

  # Linear drift tau(t) through the two endpoint offsets.
  slope <- if (tr != td) (or$offset - od$offset) / (tr - td) else 0
  intercept <- od$offset                                  # at t = td
  tau <- function(t) intercept + slope * (as.numeric(t) - td)

  list(
    offset_deploy = od$offset,
    offset_retrieve = or$offset,
    curvature = c(deploy = od$curvature, retrieve = or$curvature),
    drift = c(intercept = intercept, slope_per_sec = slope),
    correct = function(t) {
      .POSIXct(as.numeric(t) + tau(t), tz = "UTC")
    }
  )
}

# Rough light calibration for clock estimation only, used when the caller
# supplies none.
#
# The clock offset is identified by the TIMING OF TWILIGHT, so the calibration
# must put the light's fall to zero where twilight actually is. Expected light
# therefore runs from the full observed range down to zero across civil twilight
# (zenith 85 to 96), matching the engine's own auto-calibration fallback, and the
# light is baselined to its lower quantile so that "dark" really is zero.
#
# An earlier version spread the same fall across the whole 0-96 degree zenith
# range (`calibration = c(hi, rng/96)`). That leaves a gentle ramp with no
# twilight edge, so the profile likelihood in the offset has no interior optimum
# and the search runs to its boundary. Tested on ten elephant-seal deployments,
# every tag returned exactly +/- the search half-width; with the twilight slope
# every tag returns an interior peak. Do not reintroduce a shallower slope here.
.quick_light_calibration <- function(light) {
  lo <- as.numeric(stats::quantile(light, 0.05, na.rm = TRUE))
  hi <- as.numeric(stats::quantile(light, 0.95, na.rm = TRUE))
  rng <- max(hi - lo, .Machine$double.eps)
  slope <- rng / (96 - 85)
  list(baseline = lo,
       calibration = c(slope * 96, slope),   # expected = slope*(96 - zenith), clamped
       likelihood_params = c(1 / (rng * 0.5), rng, 0.1))
}

#' Estimate the clock correction from the known deployment and retrieval fixes
#'
#' Engine-agnostic front end: clock calibration is a transform on the observation
#' timestamps, computed before any engine runs, so the same correction serves
#' [TwilightFreeSMC()] and [TwilightFreeGrid()] identically (both expose it via
#' `calibrate = TRUE`). Returns `NULL` with a warning if either endpoint location
#' is unknown. If `calibration`/`likelihood_params` are not supplied, a rough
#' light calibration is derived internally, which is adequate because the clock
#' estimate depends on twilight timing rather than absolute light scale.
#'
#' @inheritParams calibrate_clock
#' @param date_time Observation times (POSIXct).
#' @param light Observed light values.
#' @param start_lon,start_lat,end_lon,end_lat Known release and recovery
#'   locations (decimal degrees).
#' @return A clock object as returned by [calibrate_clock()], or `NULL`.
#' @seealso [calibrate_clock()], [eval_logpt_loc()]
#' @export
calibrate_clock_from_endpoints <- function(date_time, light,
                                           start_lon, start_lat, end_lon, end_lat,
                                           calibration = NULL, likelihood_params = NULL,
                                           window_days = 3, search = 2 * 3600, by = 60) {
  if (anyNA(c(start_lon, start_lat, end_lon, end_lat))) {
    warning("clock calibration needs known start AND end locations; skipping")
    return(NULL)
  }
  if (is.null(calibration) || is.null(likelihood_params)) {
    # Fit the tag's own response curve at the known release position rather than
    # assuming one. The clock offset is identified by the timing of twilight, so
    # a response with the wrong transition width puts the profile's peak in the
    # wrong place, or removes it altogether.
    resp <- tryCatch(fit_light_response(date_time, light, start_lon, start_lat),
                     error = function(e) NULL)
    if (!is.null(resp)) {
      if (is.null(calibration)) calibration <- resp$calibration
      if (is.null(likelihood_params))
        likelihood_params <- c(1 / (resp$amp * 0.5), resp$max_light, 0.10)
    } else {
      qc <- .quick_light_calibration(light)
      if (is.null(calibration)) {
        calibration <- qc$calibration
        # This fallback assumes baselined light, so baseline it. A user-supplied
        # calibration is left alone: it belongs to their light scale.
        light <- pmax(0, light - qc$baseline)
      }
      if (is.null(likelihood_params)) likelihood_params <- qc$likelihood_params
    }
  }
  calibrate_clock(
    deploy   = list(lon = start_lon, lat = start_lat, time = min(date_time)),
    retrieve = list(lon = end_lon,   lat = end_lat,   time = max(date_time)),
    times = date_time, light = light,
    calibration = calibration, likelihood_params = likelihood_params,
    window_days = window_days, search = search, by = by)
}

#' Format a fitted clock model for reporting
#'
#' Returns a short, human-readable summary of a clock object (from
#' [calibrate_clock()] / [calibrate_clock_from_endpoints()]): the offset at each
#' end and the drift rate. Used by the `print` methods and suitable for a methods
#' section. Reports a fast clock (events timestamped early) as a positive
#' correction, since the correction adds time.
#'
#' @param clock A clock object, or `NULL`.
#' @return A character vector of summary lines (empty if `clock` is `NULL`).
#' @export
format_clock <- function(clock) {
  if (is.null(clock)) return(character(0))
  drift_per_day <- clock$drift[["slope_per_sec"]] * 86400
  c(
    sprintf("Clock offset (deploy): %+.1f min", clock$offset_deploy / 60),
    sprintf("Clock offset (retrieve): %+.1f min", clock$offset_retrieve / 60),
    sprintf("Clock drift: %+.2f min/day", drift_per_day / 60)
  )
}
