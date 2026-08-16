# Fitting the tag's light response curve.
#
# The engine needs to know what the clear-sky light record looks like as a
# function of solar zenith angle. That is a property of the LIGHT CHANNEL, and
# it differs between tag models and between individual tags of the same model,
# so it has to be measured rather than assumed.
#
# The original clamped-linear response (light falls from saturation to zero
# across civil twilight, zenith 85 to 96) is right for a logger whose output is
# linear in irradiance. Many archival tags, including the Wildlife Computers
# Mk9, report a compressed, roughly logarithmic light level whose clear-sky
# curve declines smoothly over tens of degrees of zenith and settles on a
# non-zero dark reading. Fitting a clamped line to such a channel gets the
# twilight gradient wrong by a factor of three, and latitude is read from
# exactly that gradient.

#' Fit a tag's clear-sky light response against solar zenith angle
#'
#' Estimates the four-parameter logistic response used by the geolocation
#' engines, from a stretch of light for which the location is known (normally
#' the deployment period at the release site), or from a whole record once a
#' first-pass track is available (see [refine_light_response()]).
#'
#' @details
#' The response model is
#' \deqn{\mu(z) = \mathrm{floor} + \frac{\mathrm{amp}}{1 + \exp((z - z_{50})/s)}}
#' where \eqn{z} is the solar zenith angle in degrees, `floor` is the sensor's
#' reading in darkness, `amp` the clear-sky amplitude above it, `z50` the zenith
#' at half amplitude, and `s` the transition scale (the 90 to 10 per cent span
#' is about `4.39 * s`). It degenerates to the clamped-linear response as `s`
#' approaches zero, so it is a strict generalisation.
#'
#' Fitting uses the **upper envelope** of light against zenith rather than its
#' centre. Shading, cloud and (for a diving animal) depth can only ever reduce
#' the measured light, so the clear-sky curve is the top of the point cloud, not
#' its middle. The envelope is taken as a high quantile within each zenith bin
#' and forced to be non-increasing, since clear-sky light cannot rise as the sun
#' sets.
#'
#' The dark level is read from the envelope's asymptote at high zenith, not from
#' a low quantile of the whole record. For a continuously diving animal a low
#' quantile is a deep daytime dive, which is darker than night at the surface,
#' and using it makes the model expect zero where the sensor reads well above
#' zero.
#'
#' @section Two quantities, two windows:
#' The fit estimates two things that want different stretches of record, and
#' `scale_light` exists so they do not have to share one.
#'
#' The **geometry** (`z50` and `scale`) is where twilight sits in zenith and how
#' fast the channel falls through it. It needs a window over which the position
#' is genuinely known. An animal that has already left the release site smears
#' the envelope across zenith and flattens the fitted curve, and since latitude
#' is read from that gradient the error is first order. On ten northern elephant
#' seal deployments, assuming the release site for a fixed fifteen days while
#' two of the animals were already hundreds of kilometres away flattened one
#' tag's transition from 30 to 54 degrees and cost 543 km of median error
#' against 289 km fitted at the positions the animal was actually at.
#'
#' The **intensity scale** (`baseline` and `max_light`) is the range the light is
#' expressed on. It needs a window whose *shading* regime matches the record
#' being tracked, which is a quite different requirement. For a diving animal
#' the natural known-position window is a haul-out, and a hauled-out animal
#' never dives, so its fifth percentile is night at the surface rather than a
#' deep daytime dive. Carrying that scale into a record full of dives cost 447 km
#' against 104 km when only the geometry was taken from the haul-out.
#'
#' So for a diving animal the recipe is: geometry from the pre-departure
#' haul-out, where the animal is stationary, ashore and under open sky, and
#' `scale_light` from the record about to be tracked. This needs no independent
#' positions beyond the release site. Where several animals carry the same tag
#' model, pool the geometry across them with [pool_light_responses()]; that is
#' worth more again, and it is the only option for an animal that departs too
#' fast to be calibrated at all.
#'
#' Do not instead try to bootstrap the geometry from a first-pass track. Latitude
#' error biases the zenith angles the envelope is built from, which flattens the
#' response, which biases latitude further; the loop amplifies rather than
#' corrects. Measured over the calibration window alone it still ran away, from
#' 321 km to 391 km to 755 km on successive turns. See [refine_light_response()].
#'
#' @param time Observation times (POSIXct, or numeric Unix seconds).
#' @param light Observed light, on the tag's own scale. Do not baseline it: the
#'   fitted `floor` is the baseline.
#' @param lon,lat Known location for this stretch of record. Either scalars (a
#'   stationary calibration period) or vectors as long as `time` (a known or
#'   first-pass track).
#' @param scale_light Optional light values from which to take the intensity
#'   scale (`baseline` and `max_light`), when the window that pins the
#'   *geometry* is not representative of the record to be fitted. Give it the
#'   light from the record you are about to track. See the details.
#' @param env_q Quantile taken within each zenith bin to form the upper
#'   envelope. Default `0.95`. Lower it if the record is short.
#' @param bin Zenith bin width in degrees (default `1`).
#' @param night_zenith Zenith above which the record is treated as dark for the
#'   purpose of estimating `floor` (default `108`, i.e. past astronomical
#'   twilight).
#' @param min_bin Minimum observations in a zenith bin for it to be used
#'   (default `3`).
#'
#' @return An object of class `tf_light_response`. The fields to pass to an
#'   engine are `calibration` (`c(intercept, slope)`, the tangent to the fitted
#'   curve at half amplitude), `baseline` (subtract it from the light before
#'   fitting a track; the 5th percentile, not the dark level) and `max_light`.
#'   Those last two, and hence the slope, come from `scale_light` when it is
#'   supplied. Also returned
#'   are `calibration_logistic` (`c(floor, amp, z50, scale)`, for engines given
#'   the four-parameter form directly; its `floor` is relative to `baseline`, so
#'   the whole vector is on the scale of the baselined light the engine sees, and
#'   is not the raw `dark_level`), `width_deg` (the
#'   90 to 10 per cent zenith span), `saturate_at` and `zero_at`, and `envelope`,
#'   the binned envelope the fit was made to. `NULL` if the record is too short
#'   or too dark to support a fit.
#'
#' @seealso [pool_light_responses()], [refine_light_response()],
#'   [TwilightFreeGrid()]
#' @examples
#' # A simulated logger that is linear in irradiance: the fit finds a narrow
#' # transition, ending near the zenith the line really reaches zero at (96).
#' t <- seq(as.POSIXct("2021-06-01", tz = "UTC"), by = "5 min", length.out = 10 * 288)
#' z <- solar_zenith(as.numeric(t), rep(150, length(t)), rep(-45, length(t)))
#' lig <- pmin(pmax(558.5 - 5.818 * z, 0), 64)
#' fit <- fit_light_response(t, lig, 150, -45)
#' round(c(width = fit$width_deg, slope = fit$slope, zero = fit$zero_at), 2)
#' # The slope comes back steeper than the 5.818 the data were built with, and
#' # that is the construction rather than a failure: a logistic fitted to a
#' # clamped line has to bend at both ends, so its gradient at half amplitude
#' # exceeds the line's by about 1.37. On a channel that really is smooth, which
#' # is the case this function exists for, the tangent is the true local
#' # gradient. Read `width_deg` to tell a linear channel (about 11 degrees) from
#' # a log-scaled one (several tens); do not read `slope` as the line's own.
#' @importFrom stats quantile median optim
#' @export
fit_light_response <- function(time, light, lon, lat,
                               env_q = 0.95, bin = 1, night_zenith = 108,
                               min_bin = 3, scale_light = NULL) {
  tt <- as.numeric(time)
  n <- length(tt)
  lon <- if (length(lon) == 1L) rep(lon, n) else lon
  lat <- if (length(lat) == 1L) rep(lat, n) else lat
  ok <- is.finite(tt) & is.finite(light) & is.finite(lon) & is.finite(lat)
  if (sum(ok) < 200) return(NULL)
  tt <- tt[ok]; li <- as.numeric(light)[ok]
  sl <- if (is.null(scale_light)) li else as.numeric(scale_light)
  sl <- sl[is.finite(sl)]
  if (!length(sl)) return(NULL)
  q05 <- as.numeric(stats::quantile(sl, 0.05))
  rng <- as.numeric(stats::quantile(sl, 0.95)) - q05
  if (!is.finite(rng) || rng <= 0) return(NULL)

  z <- solar_zenith(tt, lon[ok], lat[ok])
  zb <- round(z / bin) * bin
  cnt <- table(zb)
  keep <- as.numeric(names(cnt))[cnt >= min_bin]
  if (length(keep) < 8) return(NULL)
  env <- vapply(keep, function(b) stats::quantile(li[zb == b], env_q, names = FALSE), 0)
  o <- order(keep); keep <- keep[o]; env <- env[o]

  # Work down from the brightest bin and force the envelope to be non-increasing:
  # clear-sky light cannot rise as the sun sets. Anchoring on the envelope's own
  # maximum keeps one dim bin at low zenith (fog, a shaded animal) from
  # collapsing everything below it.
  # Keep the envelope as measured, before the monotone constraint below. The
  # constraint is right for fitting the transition, but it is wrong at night:
  # the dark reading drifts slightly UPWARD with zenith on both tag families,
  # and `cummin` pins it to the running minimum, biasing the fitted floor low.
  # `light_response_table()` needs the unconstrained shape.
  imax <- which.max(env)
  keep <- keep[imax:length(keep)]
  env_raw <- env[imax:length(env)]
  env <- cummin(env_raw)

  night <- env[keep >= night_zenith]
  floor_l <- if (length(night)) stats::median(night) else min(env)
  top <- env[1]
  amp <- top - floor_l
  if (!is.finite(amp) || amp <= 0) return(NULL)

  # Initialise from where the envelope crosses half amplitude, then fit z50 and
  # the scale by least squares on the envelope.
  half <- keep[which.min(abs(env - (floor_l + amp / 2)))]
  q9 <- keep[which.min(abs(env - (floor_l + 0.9 * amp)))]
  q1 <- keep[which.min(abs(env - (floor_l + 0.1 * amp)))]
  s0 <- max((q1 - q9) / 4.394, 0.5)

  sse <- function(p) {
    mu <- floor_l + amp / (1 + exp((keep - p[1]) / max(p[2], 1e-3)))
    sum((env - mu)^2)
  }
  fit <- tryCatch(stats::optim(c(half, s0), sse, method = "Nelder-Mead",
                               control = list(reltol = 1e-10, maxit = 500)),
                  error = function(e) NULL)
  par <- if (is.null(fit)) c(half, s0) else fit$par
  z50 <- par[1]; scale <- max(par[2], 1e-3)

  # The engine's clamped-linear response is the TANGENT to this logistic at half
  # amplitude: slope amp/(4*scale), leaving saturation at z50 - 2*scale and
  # reaching zero at z50 + 2*scale. Validated against Argos on ten elephant-seal
  # deployments, that tangent line beats the logistic itself, because the engine's
  # spike normaliser is defined over [0, max_light] and behaves better when the
  # expected curve actually attains both ends of that interval. So the logistic is
  # used as a stable way to MEASURE the response, and the tangent is what the
  # engine is given. `calibration_logistic` is kept for engines given the
  # four-parameter form directly.
  # The tangent saturates at `4 * slope * scale`, so the slope is set from `rng`
  # (which is `max_light`) rather than from the envelope's `amp`. That makes the
  # expected curve reach exactly the top of the range the light is expressed on;
  # tie it to `amp` instead and the curve tops out somewhere else, leaving the
  # engine unable to explain the brightest observations as clear sky. The two
  # differ by up to a quarter even measured on the same window, and it is worth
  # 473 km against 363 km of median error on ten elephant seal deployments,
  # rising to 411 km against 259 km once `scale_light` puts them on different
  # windows.
  lin_slope <- rng / (4 * scale)
  lin_zero  <- z50 + 2 * scale
  # `baseline` is what to subtract from the light before fitting a track, and it
  # is the low quantile, NOT the fitted dark level. Subtracting the dark level
  # clips every deep-dive observation to zero and destroys the gradient the
  # engine reads: scored against Argos it roughly doubles the error (832 km
  # against 421 km on the same tags). `floor` is still reported, as the measured
  # dark reading, because it is a useful diagnostic of the channel.
  # `calibration_logistic` must be on the scale the ENGINE works in. The engine
  # is handed `pmax(0, light - baseline)` and evaluates
  # `floor + amp / (1 + exp((z - z50) / scale))`, so its floor has to be the dark
  # level RELATIVE to `baseline`, not the raw dark reading. Supplying the raw
  # `floor_l` puts the whole expected curve `floor_l - q05` too high at every
  # zenith: measured across ten elephant-seal deployments that is 18.2 units,
  # 11.4% of the response range, and it is why the logistic form scored 833 km
  # against the tangent's 354 and was set aside. That comparison was not a fair
  # one. Clamped at zero because the observations are, and because a negative
  # expectation falls outside the spike normaliser's [0, max_light] support.
  out <- list(calibration = c(lin_slope * lin_zero, lin_slope),
              calibration_logistic = c(max(floor_l - q05, 0), amp, z50, scale),
              baseline = q05, max_light = rng, dark_level = floor_l,
              slope = lin_slope, zero_at = lin_zero,
              saturate_at = z50 - 2 * scale,
              width_deg = 4.394 * scale,
              z50 = z50, scale = scale, floor = floor_l, amp = amp,
              envelope = data.frame(zenith = keep, light = env,
                                    measured = env_raw),
              n_used = sum(ok))
  class(out) <- "tf_light_response"
  out
}

#' Pool response geometry across tags of the same model
#'
#' Replaces each tag's fitted `z50` and `scale` with the median across a set of
#' fits, while every tag keeps its own intensity scale. Use it when a study
#' deploys several tags of the same model, which is the usual case.
#'
#' @details
#' The response geometry is a property of the **light channel**, so tags of one
#' model are repeat measurements of nearly the same curve, and any single
#' deployment measures it noisily. Pooling is worth more than it sounds. On ten
#' northern elephant seal deployments the per-animal haul-out fits agreed
#' closely (`z50` spanning 91.8 to 93.1 degrees, transition width 19.8 to 30.5),
#' and substituting their median improved even the animals that had contributed
#' a perfectly good fit of their own, taking median error from 259 km to 209 km.
#' Across all ten it went from 473 km to 237 km, the worst tag from 1255 km to
#' 483 km, and latitude interval coverage from 0.68 to 0.89.
#'
#' It also covers the animals that cannot be calibrated at all. Three of those
#' ten left within a day of being tagged, leaving no stretch of record over
#' which their position was known; pooled geometry is the only thing that gives
#' them a defensible response.
#'
#' The intensity scale is deliberately **not** pooled. It depends on how much
#' shading and diving the individual record contains, which is a property of the
#' animal rather than the tag, and transferring it between animals makes things
#' worse (see the details of [fit_light_response()]).
#'
#' @param geometry A list of `tf_light_response` objects to take the geometry
#'   from, normally fitted over each animal's known-position window. `NULL`
#'   entries are ignored, so a list with gaps can be passed straight in.
#' @param scale A list of `tf_light_response` objects supplying each tag's own
#'   `baseline` and `max_light`, normally fitted over the record about to be
#'   tracked. Defaults to `geometry`. The returned list follows this one, in its
#'   order and with its names.
#'
#' @return A list of `tf_light_response` objects, one per element of `scale`,
#'   each carrying the pooled geometry and its own intensity scale. Elements are
#'   `NULL` wherever `scale` is `NULL`. The pooled `z50` and `scale` are also
#'   attached as the attribute `"pooled"`.
#'
#' @seealso [fit_light_response()]
#' @examples
#' # Two tags of one model, one calibrated over a long stationary period and one
#' # over a short noisy one: pooling lends the second the first's precision.
#' t <- seq(as.POSIXct("2021-06-01", tz = "UTC"), by = "10 min", length.out = 3000)
#' z <- solar_zenith(as.numeric(t), rep(150, length(t)), rep(-45, length(t)))
#' mk <- function(seed) {
#'   set.seed(seed)
#'   fit_light_response(t, pmin(pmax(558.5 - 5.818 * z, 0), 64) *
#'                        runif(length(t), 0.8, 1), 150, -45)
#' }
#' pooled <- pool_light_responses(list(a = mk(1), b = mk(2)))
#' attr(pooled, "pooled")
#' @export
pool_light_responses <- function(geometry, scale = geometry) {
  if (!is.list(geometry) || !is.list(scale))
    stop("`geometry` and `scale` must both be lists of tf_light_response objects")
  g <- Filter(Negate(is.null), geometry)
  if (!length(g)) stop("`geometry` contains no fitted responses")
  bad <- !vapply(g, inherits, logical(1), "tf_light_response")
  if (any(bad)) stop("`geometry` must contain tf_light_response objects")
  z50_p <- stats::median(vapply(g, function(r) r$z50, 0))
  sc_p  <- stats::median(vapply(g, function(r) r$scale, 0))

  out <- lapply(scale, function(r) {
    if (is.null(r)) return(NULL)
    if (!inherits(r, "tf_light_response"))
      stop("`scale` must contain tf_light_response objects")
    # Same construction as fit_light_response(): the tangent at half amplitude,
    # rising to the top of the range the light is expressed on.
    slope <- r$max_light / (4 * sc_p)
    zero  <- z50_p + 2 * sc_p
    r$z50 <- z50_p; r$scale <- sc_p
    r$calibration <- c(slope * zero, slope)
    # As in fit_light_response(): relative to `baseline`, not the raw dark level.
    r$calibration_logistic <- c(max(r$floor - r$baseline, 0), r$amp, z50_p, sc_p)
    r$slope <- slope; r$zero_at <- zero
    r$saturate_at <- z50_p - 2 * sc_p
    r$width_deg <- 4.394 * sc_p
    r$pooled_from <- length(g)
    r
  })
  attr(out, "pooled") <- c(z50 = z50_p, scale = sc_p, n = length(g))
  out
}

#' Refine a light response against a first-pass track
#'
#' The deployment period is the natural place to fit a response curve, because
#' the location is known there, but it is often a poor one: a few days is a
#' small sample, and a release site can be persistently foggy or the animal
#' shaded while it remains ashore. Once a track has been estimated, the whole
#' record becomes usable, spanning far more zenith angles and far more weather.
#'
#' This is the refit half of an alternating scheme: fit a response, fit a track,
#' refit the response along that track, refit the track.
#'
#' **Use with care, and check it against ground truth before trusting it.** On
#' ten northern elephant seal deployments a single refit pass made the fit worse,
#' not better (median error 421 km rising to 1178 km): residual latitude error in
#' the first-pass track biases the zenith angles the envelope is built from,
#' which flattens the fitted response, which biases latitude further. The loop
#' amplifies an initial bias instead of correcting it. It is provided because it
#' is the right idea where the calibration period is genuinely unusable, but the
#' deployment-site fit is the safer default.
#'
#' @param time,light As in [fit_light_response()].
#' @param track A data frame of the first-pass track with columns `time`
#'   (POSIXct), `lon` and `lat`. Positions are interpolated to the observation
#'   times.
#' @param ... Passed to [fit_light_response()].
#'
#' @return A `tf_light_response`, or `NULL` if the fit fails.
#' @seealso [fit_light_response()]
#' @export
refine_light_response <- function(time, light, track, ...) {
  if (!all(c("time", "lon", "lat") %in% names(track)))
    stop("`track` needs columns time, lon and lat")
  tt <- as.numeric(time)
  kt <- as.numeric(track$time)
  o <- order(kt); kt <- kt[o]
  lon <- stats::approx(kt, track$lon[o], xout = tt, rule = 2)$y
  lat <- stats::approx(kt, track$lat[o], xout = tt, rule = 2)$y
  fit_light_response(time, light, lon, lat, ...)
}

#' @method print tf_light_response
#' @export
print.tf_light_response <- function(x, ...) {
  cat("\nTag light response (logistic in solar zenith)\n")
  cat("=============================================\n")
  cat(sprintf("dark level (floor) : %.1f\n", x$floor))
  cat(sprintf("clear-sky amplitude: %.1f\n", x$amp))
  cat(sprintf("half-amplitude at  : %.1f deg zenith\n", x$z50))
  cat(sprintf("transition width   : %.1f deg (90%% to 10%%)\n", x$width_deg))
  cat(sprintf("observations used  : %d\n", x$n_used))
  cat("---------------------------------------------\n")
  cat("A logger linear in irradiance gives a width near 11 degrees; a\n")
  cat("log-scaled channel gives several tens of degrees.\n")
  cat("=============================================\n")
  invisible(x)
}


#' Build a non-parametric light response for an engine
#'
#' Turns one or more fitted responses into a lookup table of expected light
#' against solar zenith, on the baselined scale the engines work in. Pass the
#' result as `calibration`.
#'
#' Use it when neither parametric form fits the channel. On 29 northern
#' elephant-seal deployments the measured clear-sky envelope holds a daytime
#' plateau, collapses across about 15 degrees of zenith, then settles on a floor
#' of roughly a quarter of the range. The two-parameter clamped line cannot
#' place a floor at all, and forcing it to zero is the single largest
#' misspecification in that data. A logistic can place one, but is still
#' committed to a symmetric transition. The table commits to nothing, at the
#' cost of `length(seq(from, to, by))` numbers instead of two or four.
#'
#' The table is built from the envelope as MEASURED, not the monotone-constrained
#' one used to fit `z50` and `scale`, because the constraint is wrong at night:
#' the dark reading drifts slightly upward with zenith and `cummin` pins it to
#' the running minimum.
#'
#' @param x A `tf_light_response`, or a list of them to pool across tags of the
#'   same model. Pooling is by median at each zenith, after each tag's envelope
#'   is put on its own baselined scale, so tags of different brightness combine
#'   without one dominating.
#' @param from,to,by The zenith grid, in degrees. The default spans full day to
#'   full night at one-degree steps.
#' @param min_tags When pooling, the number of tags that must contribute to a
#'   zenith before it is used. Bins below this are filled by interpolation from
#'   the ones that qualify.
#' @param flat_outside Optional `c(lo, hi)` zeniths beyond which the response is
#'   held constant at its value on the boundary. Pass `c(x$saturate_at,
#'   x$zero_at)` to reuse the window the clamped-linear response already
#'   implies, which introduces no new tuning constant.
#'
#'   This matters more than it looks. A response that varies with zenith
#'   everywhere lets every observation carry information about position, and on
#'   elephant-seal records 73 per cent of observations sit outside the
#'   transition, where the channel is flat, the residual is largest (sd 0.169 of
#'   the range in daylight against 0.094 at night) and any error in the curve
#'   converts straight into a position error. The clamped-linear response gets
#'   its accuracy from being FLAT there, not from being right: scored against
#'   Argos on 29 deployments it beats a logistic 246 km to 634 despite
#'   describing the measured curve five times worse. Flattening the table
#'   outside the transition keeps the measured shape where it carries signal and
#'   discards it where it carries noise.
#'
#' @param max_light The tag's own `max_light`. Supply it whenever more than one
#'   fit is pooled: the shape is then pooled in units of each contributing tag's
#'   own range and rescaled to this one, so a shared shape does not impose a
#'   shared gain. Build one table per tag this way, as
#'   [pool_light_responses()] gives every tag shared geometry and its own
#'   intensity scale. `NULL` pools in absolute units, which is only right for a
#'   single fit.
#'
#' @return A numeric vector `c(-1, from, by, y...)`, the leading `-1` marking
#'   the table form to the engines. Attribute `"n_tags"` records how many fits
#'   contributed.
#'
#' @seealso [fit_light_response()], [pool_light_responses()]
#' @export
light_response_table <- function(x, from = 30, to = 140, by = 1, min_tags = 1,
                                 flat_outside = NULL, max_light = NULL) {
  if (inherits(x, "tf_light_response")) x <- list(x)
  x <- Filter(Negate(is.null), x)
  if (!length(x)) stop("no fitted responses supplied")
  bad <- !vapply(x, inherits, logical(1), "tf_light_response")
  if (any(bad)) stop("`x` must contain tf_light_response objects")
  if (!is.finite(by) || by <= 0) stop("`by` must be positive")
  grid <- seq(from, to, by = by)
  if (length(grid) < 2) stop("the zenith grid needs at least two points")

  # Each tag's measured envelope, moved onto the baselined scale the engine sees
  # and interpolated onto the shared grid. Held constant beyond each tag's own
  # zenith range rather than extrapolated, since a fitted response says nothing
  # about zeniths it never observed.
  #
  # When `max_light` is given, pool in units of each tag's OWN range and rescale
  # at the end. Tags of the same model differ in gain: across 29 elephant-seal
  # deployments `max_light` runs 147 to 175. Pooling in absolute units hands
  # every tag the same curve, which for the dimmest tags tops out above their
  # range (so the day saturates and flattens) and for the brightest never
  # reaches it (so the engine cannot explain their brightest observations as
  # clear sky). This is the same separation `pool_light_responses()` makes
  # between shared geometry and per-tag intensity, and the same error that cost
  # 473 km against 363 when the tangent's slope was tied to `amp`.
  norm <- !is.null(max_light)
  if (norm && (!is.numeric(max_light) || length(max_light) != 1L ||
               !is.finite(max_light) || max_light <= 0))
    stop("`max_light` must be a single positive number")
  m <- vapply(x, function(r) {
    e <- r$envelope
    y <- if (!is.null(e$measured)) e$measured else e$light
    v <- pmax(y - r$baseline, 0)
    if (norm) v <- v / r$max_light
    stats::approx(e$zenith, v, grid, rule = 2)$y
  }, numeric(length(grid)))
  if (is.null(dim(m))) m <- matrix(m, ncol = 1)

  n_at <- rowSums(is.finite(m))
  yy <- apply(m, 1, stats::median, na.rm = TRUE)
  ok <- n_at >= min_tags & is.finite(yy)
  if (sum(ok) < 2) stop("too few zeniths met `min_tags`")
  if (any(!ok)) yy <- stats::approx(grid[ok], yy[ok], grid, rule = 2)$y

  if (!is.null(flat_outside)) {
    if (length(flat_outside) != 2L || !all(is.finite(flat_outside)) ||
        flat_outside[1] >= flat_outside[2])
      stop("`flat_outside` must be c(lo, hi) with lo < hi")
    lo <- which.min(abs(grid - flat_outside[1]))
    hi <- which.min(abs(grid - flat_outside[2]))
    yy[seq_len(lo)] <- yy[lo]
    yy[hi:length(yy)] <- yy[hi]
  }

  if (norm) yy <- yy * max_light

  out <- c(-1, from, by, pmax(yy, 0))
  attr(out, "n_tags") <- length(x)
  attr(out, "flat_outside") <- flat_outside
  attr(out, "max_light") <- max_light
  out
}


# R-side twin of the Rust `expected_light()`. The two must agree exactly, since
# the clock calibration profiles the same likelihood the engines evaluate.
# Length 2 = clamped linear (the original model); length 4 = logistic; a
# negative leading element = the lookup table `c(-1, z_min, dz, y...)`.
.tf_expected_light <- function(z, calibration, max_light) {
  if (calibration[1] < 0 && length(calibration) >= 5) {
    z_min <- calibration[2]
    dz <- calibration[3]
    if (abs(dz) < 1e-9) dz <- 1e-9
    y <- calibration[-(1:3)]
    n <- length(y)
    pos <- pmin(pmax((z - z_min) / dz, 0), n - 1)
    i <- pmin(floor(pos), n - 2)
    f <- pos - i
    return(pmin(pmax(y[i + 1] * (1 - f) + y[i + 2] * f, 0), max_light))
  }
  if (length(calibration) >= 4) {
    s <- calibration[4]
    if (abs(s) < 1e-9) s <- 1e-9
    pmin(pmax(calibration[1] +
                calibration[2] / (1 + exp((z - calibration[3]) / s)), 0), max_light)
  } else {
    pmin(pmax(calibration[1] - calibration[2] * z, 0), max_light)
  }
}
