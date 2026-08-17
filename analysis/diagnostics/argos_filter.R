# A speed filter for Argos fixes: delete points implying impossible movement.
#
# WHY 3 m/s. Not chosen freely -- the fvilches delivery arrived already filtered at
# exactly 3 m/s (19 of 19 tags cap there, none above), while the 2021 delivery is raw
# with 15-25% of steps above 5 m/s and a maximum of 159,722 m/s. Matching 3 m/s makes
# the reference CONSISTENT ACROSS DELIVERIES, which matters more than the exact value:
# a systematic difference in ground-truth quality between two thirds and one third of
# the dataset is a confound in every comparison that mixes them.
#
# For scale, a northern elephant seal transits at 0.7-1.0 m/s and sustains about
# 2 m/s, so 3 m/s is permissive rather than aggressive -- it deletes the impossible,
# not the merely fast.
#
# THE ALGORITHM, and why it is not just "drop the fast steps". A single bad fix makes
# TWO steps look fast: the one into it and the one out of it. Dropping every fast
# step would delete the good neighbour as often as the outlier. So a fix is only
# removed when it is inconsistent with BOTH of its neighbours, which is the signature
# of the fix itself being wrong rather than of a genuine fast run. Iterated to a
# fixed point, because removing one outlier can expose another.
#
# This is deliberately not a state-space model. It deletes points; it does not
# reconstruct a track, and it makes no claim to.

# Great-circle km; local copy so this file stands alone.
.gc <- function(l1, p1, l2, p2, R = 6371.0088) {
  d <- pi/180
  2*R*asin(pmin(1, sqrt(sin((p2-p1)*d/2)^2 +
                        cos(p1*d)*cos(p2*d)*sin((l2-l1)*d/2)^2)))
}

#' Remove Argos fixes implying impossible movement
#'
#' @param a data.table/data.frame with `time` (POSIXct), `lon`, `lat`.
#' @param vmax Maximum plausible speed, m/s (default 3).
#' @param max_iter Safety bound on the iteration (default 20).
#' @return The input with offending rows removed, plus attributes `n_removed`,
#'   `frac_removed` and `iterations`.
argos_speed_filter <- function(a, vmax = 3, max_drop_frac = 0.5) {
  a <- as.data.frame(a)
  a <- a[order(a$time), , drop = FALSE]
  a <- a[!duplicated(a$time), , drop = FALSE]   # duplicate times -> infinite speed
  n0 <- nrow(a)
  if (n0 < 3) return(structure(a, n_removed = 0L, frac_removed = 0, iterations = 0L))

  # Speeds of every step, m/s.
  spd <- function(x) {
    n <- nrow(x); if (n < 2) return(numeric(0))
    d <- .gc(x$lon[-n], x$lat[-n], x$lon[-1], x$lat[-1])
    d / pmax(as.numeric(diff(as.numeric(x$time))), 1) * 1000
  }

  # Remove ONE fix at a time, always the worst offender, until no step exceeds vmax.
  #
  # The earlier "flag any fix whose two adjacent steps are both too fast" rule was
  # not enough: two consecutive bad fixes displaced the same way have a SLOW step
  # between them, so neither is flagged, and endpoints were never tested at all. It
  # removed 9.3% of the 2021 fixes and still left a maximum of 23,485 m/s.
  #
  # For the fastest step, one of its two endpoints is wrong. Choose the one that is
  # also inconsistent with its OTHER neighbour -- a bad fix is fast in both
  # directions, a good fix next to a bad one is fast in only the one. At a track end
  # there is no other neighbour, so the endpoint itself is the suspect.
  it <- 0L
  repeat {
    n <- nrow(a); if (n < 3) break
    v <- spd(a)
    if (!length(v) || max(v) <= vmax) break
    if ((n0 - n) >= max_drop_frac * n0) {
      warning(sprintf("argos_speed_filter: stopped after removing %.0f%% of fixes",
                      100 * (n0 - n) / n0), call. = FALSE)
      break
    }
    i <- which.max(v)          # offending step joins fixes i and i+1
    left_other  <- if (i > 1)      v[i - 1] else Inf   # Inf: an endpoint has no alibi
    right_other <- if (i + 1 < n)  v[i + 1] else Inf
    drop <- if (left_other >= right_other) i else i + 1L
    a <- a[-drop, , drop = FALSE]
    it <- it + 1L
  }
  structure(a, n_removed = n0 - nrow(a),
            frac_removed = (n0 - nrow(a)) / n0, iterations = it)
}

#' Remove Argos fixes that form implausible out-and-back excursions
#'
#' A speed filter permits `vmax` TIMES THE GAP, which on a sparse track is a great
#' deal of rope: at 3 m/s a 13-hour gap allows 140 km, so a fix displaced 65 km passes
#' the speed test and leaves a triangular detour -- out and straight back. Measured
#' after speed filtering, 0.2 to 1.6 per cent of triples in every deployment are such
#' spikes, worst on the sparsest tracks (cor(median gap, rate) = +0.45).
#'
#' The statistic is the detour ratio
#'   \code{r = (d(i-1,i) + d(i,i+1)) / d(i-1,i+1)}
#' which is about 1 for travel along a path and large for a spike. It is scale-free,
#' so it catches what a speed threshold cannot.
#'
#' BOTH conditions are required. A large ratio alone is not evidence: an animal
#' genuinely milling on station produces tight turns with big ratios over small
#' distances, and Argos noise alone does the same. Requiring the excursion to exceed
#' `min_km` as well confines removal to displacements too large to be either.
#'
#' @param a data.frame with `time`, `lon`, `lat`.
#' @param max_detour Detour ratio above which a fix is a candidate (default 10).
#' @param min_km Excursion, km, that a candidate must also exceed (default 15, above
#'   the error of a poor Argos class).
#' @param max_drop_frac Safety bound (default 0.2).
#' @return Filtered data, with attributes `n_removed` and `frac_removed`.
argos_detour_filter <- function(a, max_detour = 10, min_km = 15, max_drop_frac = 0.2) {
  a <- as.data.frame(a); a <- a[order(a$time), , drop = FALSE]
  n0 <- nrow(a); if (n0 < 4) return(structure(a, n_removed = 0L, frac_removed = 0))
  repeat {
    n <- nrow(a); if (n < 4) break
    lo <- a$lon; la <- a$lat
    d_in  <- .gc(lo[1:(n-2)], la[1:(n-2)], lo[2:(n-1)], la[2:(n-1)])
    d_out <- .gc(lo[2:(n-1)], la[2:(n-1)], lo[3:n],     la[3:n])
    d_dir <- .gc(lo[1:(n-2)], la[1:(n-2)], lo[3:n],     la[3:n])
    r   <- (d_in + d_out) / pmax(d_dir, 0.5)
    exc <- pmin(d_in, d_out)
    bad <- r > max_detour & exc > min_km
    if (!any(bad) || (n0 - n) >= max_drop_frac * n0) break
    # remove the single worst, then recompute: neighbouring triples share fixes, so
    # removing several at once can delete a good fix adjacent to a bad one
    a <- a[-(which.max(ifelse(bad, r * exc, -Inf)) + 1L), , drop = FALSE]
  }
  structure(a, n_removed = n0 - nrow(a), frac_removed = (n0 - nrow(a)) / n0)
}

#' Speed filter then detour filter
#' @param ... passed to both stages via their own arguments
argos_clean <- function(a, vmax = 3, max_detour = 10, min_km = 15) {
  s <- argos_speed_filter(a, vmax = vmax)
  d <- argos_detour_filter(s, max_detour = max_detour, min_km = min_km)
  structure(d, n_speed = attr(s, "n_removed"), n_detour = attr(d, "n_removed"),
            frac_removed = 1 - nrow(d) / nrow(as.data.frame(a)))
}

#' Drop Argos fixes that no neighbour can corroborate
#'
#' A fix with a large gap on BOTH sides is unverifiable. Neither a speed nor a shape
#' criterion can reject it: 2023041 carries one at 2023-06-23 with gaps of 213 h
#' before and 352 h after, implying 2.36 and 2.01 m/s, both legal, because at 3 m/s
#' those gaps permit 3804 km of travel. Interpolating through it dragged the track
#' 1800 km east and back, which is the triangular detour visible by eye.
#'
#' This does not decide whether the animal was there. It says the record cannot tell
#' us, and a fix that cannot be corroborated should not be allowed to define a
#' reference position.
#'
#' @param a data.frame with `time`, `lon`, `lat`.
#' @param max_isolation_h Hours. A fix is dropped when the gaps on both sides exceed
#'   this (default 48).
argos_isolation_filter <- function(a, max_isolation_h = 48) {
  a <- as.data.frame(a)[order(as.data.frame(a)$time), , drop = FALSE]
  n <- nrow(a); if (n < 3) return(structure(a, n_removed = 0L))
  g <- as.numeric(diff(as.numeric(a$time))) / 3600
  before <- c(Inf, g); after <- c(g, Inf)      # ends are isolated on one side
  bad <- pmin(before, after) > max_isolation_h
  structure(a[!bad, , drop = FALSE], n_removed = sum(bad))
}

#' Speed, detour and isolation filters in sequence
argos_clean2 <- function(a, vmax = 3, max_detour = 10, min_km = 15,
                         max_isolation_h = 48) {
  n0 <- nrow(as.data.frame(a))
  s <- argos_speed_filter(a, vmax = vmax)
  d <- argos_detour_filter(s, max_detour = max_detour, min_km = min_km)
  i <- argos_isolation_filter(d, max_isolation_h = max_isolation_h)
  structure(i, n_speed = attr(s, "n_removed"), n_detour = attr(d, "n_removed"),
            n_isolated = attr(i, "n_removed"), frac_removed = 1 - nrow(i) / n0)
}
