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
#' @param time Observation times (POSIXct, or numeric Unix seconds).
#' @param light Observed light, on the tag's own scale. Do not baseline it: the
#'   fitted `floor` is the baseline.
#' @param lon,lat Known location for this stretch of record. Either scalars (a
#'   stationary calibration period) or vectors as long as `time` (a known or
#'   first-pass track).
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
#'   Also returned
#'   are `calibration_logistic` (`c(floor, amp, z50, scale)`), `width_deg` (the
#'   90 to 10 per cent zenith span), `saturate_at` and `zero_at`, and `envelope`,
#'   the binned envelope the fit was made to. `NULL` if the record is too short
#'   or too dark to support a fit.
#'
#' @seealso [refine_light_response()], [TwilightFreeGrid()]
#' @examples
#' # A simulated logger that is linear in irradiance: the fit recovers a narrow
#' # transition, as it should.
#' t <- seq(as.POSIXct("2021-06-01", tz = "UTC"), by = "5 min", length.out = 10 * 288)
#' z <- solar_zenith(as.numeric(t), rep(150, length(t)), rep(-45, length(t)))
#' lig <- pmin(pmax(558.5 - 5.818 * z, 0), 64)
#' fit <- fit_light_response(t, lig, 150, -45)
#' round(c(width = fit$width_deg, slope = fit$slope), 2)   # true slope 5.818
#' @importFrom stats quantile median optim
#' @export
fit_light_response <- function(time, light, lon, lat,
                               env_q = 0.95, bin = 1, night_zenith = 108,
                               min_bin = 3) {
  tt <- as.numeric(time)
  n <- length(tt)
  lon <- if (length(lon) == 1L) rep(lon, n) else lon
  lat <- if (length(lat) == 1L) rep(lat, n) else lat
  ok <- is.finite(tt) & is.finite(light) & is.finite(lon) & is.finite(lat)
  if (sum(ok) < 200) return(NULL)
  tt <- tt[ok]; li <- as.numeric(light)[ok]
  q05 <- as.numeric(stats::quantile(li, 0.05))
  rng <- as.numeric(stats::quantile(li, 0.95)) - q05

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
  imax <- which.max(env)
  keep <- keep[imax:length(keep)]; env <- cummin(env[imax:length(env)])

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
  lin_slope <- amp / (4 * scale)
  lin_zero  <- z50 + 2 * scale
  # `baseline` is what to subtract from the light before fitting a track, and it
  # is the low quantile, NOT the fitted dark level. Subtracting the dark level
  # clips every deep-dive observation to zero and destroys the gradient the
  # engine reads: scored against Argos it roughly doubles the error (832 km
  # against 421 km on the same tags). `floor` is still reported, as the measured
  # dark reading, because it is a useful diagnostic of the channel.
  out <- list(calibration = c(lin_slope * lin_zero, lin_slope),
              calibration_logistic = c(floor_l, amp, z50, scale),
              baseline = q05, max_light = rng, dark_level = floor_l,
              slope = lin_slope, zero_at = lin_zero,
              saturate_at = z50 - 2 * scale,
              width_deg = 4.394 * scale,
              z50 = z50, scale = scale, floor = floor_l, amp = amp,
              envelope = data.frame(zenith = keep, light = env),
              n_used = sum(ok))
  class(out) <- "tf_light_response"
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


# R-side twin of the Rust `expected_light()`. The two must agree exactly, since
# the clock calibration profiles the same likelihood the engines evaluate.
# Length 2 = clamped linear (the original model); length 4 = logistic.
.tf_expected_light <- function(z, calibration, max_light) {
  if (length(calibration) >= 4) {
    s <- calibration[4]
    if (abs(s) < 1e-9) s <- 1e-9
    calibration[1] + calibration[2] / (1 + exp((z - calibration[3]) / s))
  } else {
    pmin(pmax(calibration[1] - calibration[2] * z, 0), max_light)
  }
}
