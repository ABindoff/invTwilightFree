use extendr_api::prelude::*;

/// Calculate the solar zenith angle (in degrees)
///
/// @param unix_time Numeric vector of Unix timestamps (seconds since 1970-01-01 UTC)
/// @param lon Numeric vector of longitudes (decimal degrees)
/// @param lat Numeric vector of latitudes (decimal degrees)
/// @return Numeric vector of solar zenith angles (degrees; 0 = overhead, 90 = horizon)
/// @export
#[extendr]
fn solar_zenith(unix_time: &[f64], lon: &[f64], lat: &[f64]) -> Vec<f64> {
    let n = unix_time.len();
    let mut out = Vec::with_capacity(n);
    
    for i in 0..n {
        let t = unix_time[i];
        let l = lon[i];
        let phi = lat[i] * std::f64::consts::PI / 180.0;
        
        // Days since J2000.0 (2000-01-01 12:00 UTC = 946728000)
        let d = (t - 946728000.0) / 86400.0;
        
        // Mean anomaly of the sun
        let g = (357.529 + 0.98560028 * d).to_radians();
        
        // Mean longitude of the sun
        let q = 280.459 + 0.98564736 * d;
        
        // Ecliptic longitude of the sun
        let l_sun = (q + 1.915 * g.sin() + 0.020 * (2.0 * g).sin()).to_radians();
        
        // Obliquity of the ecliptic
        let e = (23.439 - 0.00000036 * d).to_radians();
        
        // Declination of the sun
        let delta = (e.sin() * l_sun.sin()).asin();
        
        // Right ascension of the sun
        let ra = f64::atan2(e.cos() * l_sun.sin(), l_sun.cos());
        
        // GMST (Greenwich Mean Sidereal Time)
        let gmst = (18.697374558 + 24.06570982441908 * d) * 15.0_f64.to_radians();
        
        // Local Mean Sidereal Time
        let lmst = gmst + l * std::f64::consts::PI / 180.0;
        
        // Hour angle
        let h = lmst - ra;
        
        // Zenith angle
        let cos_z = phi.sin() * delta.sin() + phi.cos() * delta.cos() * h.cos();
        out.push(cos_z.acos() * 180.0 / std::f64::consts::PI);
    }
    
    out
}

/// Normalising constant of the asymmetric-exponential spike over the tag's
/// support [0, max_light], given the expected (clamped) light `mu`. Without this
/// the spike is improper and its mass depends on `mu`, which (because `mu`
/// depends on latitude through the solar geometry) tilts the per-cell likelihood
/// by latitude and biases the latitude posterior. Confirmed by SBC; see
/// notes/topology/sbc_design.md.
#[inline]
fn spike_normaliser(mu: f64, lambda: f64, max_light: f64) -> f64 {
    let lo = 1.0 - (-lambda * mu).exp();                              // mass below mu
    let hi = 0.5 * (1.0 - (-2.0 * lambda * (max_light - mu)).exp());  // mass above mu
    (lo + hi).max(1e-12)
}

/// Copy a calibration slice into a fixed-size array so `Ctx` needs no lifetime
/// for it. Anything beyond the first four entries is ignored.
#[inline]
fn cal_array(c: &[f64]) -> [f64; 4] {
    let mut a = [0.0f64; 4];
    for i in 0..c.len().min(4) { a[i] = c[i]; }
    a
}

/// Clear-sky expected light as a function of solar zenith angle.
///
/// Two response models, selected by the length of `calibration`, so that
/// existing length-2 calls are bit-identical to before:
///
/// * length 2 — `c(intercept, slope)`: the original clamped-linear response,
///   `clamp(intercept - slope*z, 0, max_light)`. Correct for a logger whose
///   output is linear in irradiance, where light collapses within civil
///   twilight.
///
/// * length 4 — `c(floor, amp, z50, scale)`: a logistic response,
///   `floor + amp / (1 + exp((z - z50)/scale))`. `floor` is the sensor's reading
///   in darkness, `amp` the clear-sky amplitude above it, `z50` the zenith at
///   half amplitude, and `scale` sets the width of the transition (the 90% to
///   10% span is about 4.39*scale). This is the model for a log-scaled channel,
///   whose light declines smoothly over tens of degrees of zenith rather than
///   falling off a cliff, and it degenerates to the clamped line as `scale`
///   goes to zero. Archival tags differ in this response even within a batch,
///   so it is fitted per tag rather than assumed.
#[inline]
fn expected_light(z: f64, calibration: &[f64], max_light: f64) -> f64 {
    if calibration.len() >= 4 {
        let (floor, amp, z50, scale) = (calibration[0], calibration[1],
                                        calibration[2], calibration[3]);
        let s = if scale.abs() < 1e-9 { 1e-9 } else { scale };
        floor + amp / (1.0 + ((z - z50) / s).exp())
    } else {
        (calibration[0] - calibration[1] * z).max(0.0).min(max_light)
    }
}

/// Asymmetric-exponential spike density for the spike-and-slab light model,
/// normalised over [0, max_light]. obs <= expected: one-sided exponential
/// (shading); obs > expected: penalised with twice the decay rate (sensor
/// physics forbid observations brighter than the clear-sky maximum).
#[inline]
fn spike_density(obs: f64, expected: f64, lambda: f64, max_light: f64) -> f64 {
    let raw = if obs <= expected {
        lambda * (-lambda * (expected - obs)).exp()
    } else {
        lambda * (-lambda * 2.0 * (obs - expected)).exp()
    };
    raw / spike_normaliser(expected, lambda, max_light)
}

/// Calculate the log-likelihood of observed light given proposed tracks.
/// uses a spike-and-slab model for false light events.
///
/// @param obs_light Observed light values
/// @param expected_light Expected maximum light values (from solar zenith)
/// @param lambda Decay rate for the shading exponential distribution
/// @param max_light Maximum possible light value for the tag
/// @param prob_slab Probability of a false light event (the slab)
/// @export
#[extendr]
fn light_log_likelihood(
    obs_light: &[f64],
    expected_light: &[f64],
    lambda: f64,
    max_light: f64,
    prob_slab: f64,
) -> f64 {
    let mut log_lik = 0.0;
    
    // Density of the slab (uniform distribution over possible light values)
    let slab_density = 1.0 / max_light;
    
    for i in 0..obs_light.len() {
        let obs = obs_light[i];
        let exp = expected_light[i];
        
        let spike = spike_density(obs, exp, lambda, max_light);

        // Mixture model: (1 - pi) * True Light + pi * False Light
        let marginal_density = (1.0 - prob_slab) * spike + prob_slab * slab_density;
        
        // Accumulate log likelihood
        log_lik += marginal_density.ln();
    }
    
    log_lik
}

use rand::prelude::*;
use rand::rngs::StdRng;
use rand_distr::{Normal, Beta, Gamma};

#[derive(Clone, Copy, Debug)]
struct Particle {
    lat: f64,
    lon: f64,
    weight: f64,
    state: usize,
    prob_slab: f64,
}

/// Time-only solar ephemeris, shared across all particles at a given timestamp.
/// Returns (sin(delta), cos(delta), right ascension, GMST). The heavy trig here
/// depends only on the observation time, not on position, so it is computed once
/// per observation and reused for every particle.
#[inline]
fn solar_ephemeris(t: f64) -> (f64, f64, f64, f64) {
    let d = (t - 946728000.0) / 86400.0;
    let g = (357.529 + 0.98560028 * d).to_radians();
    let q = 280.459 + 0.98564736 * d;
    let l_sun = (q + 1.915 * g.sin() + 0.020 * (2.0 * g).sin()).to_radians();
    let e = (23.439 - 0.00000036 * d).to_radians();
    let delta = (e.sin() * l_sun.sin()).asin();
    let ra = f64::atan2(e.cos() * l_sun.sin(), l_sun.cos());
    let gmst = (18.697374558 + 24.06570982441908 * d) * 15.0_f64.to_radians();
    (delta.sin(), delta.cos(), ra, gmst)
}

/// Solar zenith (degrees) from a precomputed time-only ephemeris and a position.
/// Bit-identical to the inline computation in get_solar_zenith: the only change is
/// that the time-dependent terms are passed in rather than recomputed per particle.
#[inline]
fn zenith_from_ephem(sin_delta: f64, cos_delta: f64, ra: f64, gmst: f64, l: f64, phi_deg: f64) -> f64 {
    let phi = phi_deg.to_radians();
    let lmst = gmst + l.to_radians();
    let h = lmst - ra;
    let cos_z = phi.sin() * sin_delta + phi.cos() * cos_delta * h.cos();
    cos_z.acos().to_degrees()
}

#[allow(dead_code)] // kept as the reference implementation; hot loops use the split form
fn get_solar_zenith(t: f64, l: f64, phi_deg: f64) -> f64 {
    let (sin_delta, cos_delta, ra, gmst) = solar_ephemeris(t);
    zenith_from_ephem(sin_delta, cos_delta, ra, gmst, l, phi_deg)
}

fn interpolate_lon(lon1: f64, lon2: f64, f: f64) -> f64 {
    let mut dlon = lon2 - lon1;
    if dlon > 180.0 { dlon -= 360.0; }
    if dlon < -180.0 { dlon += 360.0; }
    let mut res = lon1 + f * dlon;
    if res > 180.0 { res -= 360.0; }
    if res < -180.0 { res += 360.0; }
    res
}

/// Run the Particle Filter Engine
/// 
/// Core SMC engine for TwilightFree track reconstruction.
/// 
/// @param unix_times Timestamps of observations
/// @param obs_light Observed light values
/// @param n_particles Number of particles to use
/// @param start_lat Starting Latitude
/// @param start_lon Starting Longitude
/// @param end_lat Ending Latitude (use NA if unknown)
/// @param end_lon Ending Longitude (use NA if unknown)
/// @param method Method to use: "guided", "ffbs", or "forward"
/// @param step_hours Hours between particle movement steps
/// @param diffusion Diffusion coefficient (kilometers per sqrt(day))
/// @param trans_prob Flattened row-major transition probability matrix for behavioral states
/// @param calibration Response parameters: `c(intercept, slope)` for the clamped-linear
///   model, or `c(floor, amp, z50, scale)` for the logistic model (see details in
///   `fit_light_response`). Length selects the model.
/// @param likelihood_params c(lambda, max_light, prob_slab)
/// @param mask_matrix Flattened spatial mask matrix (0 = impassable); empty for no mask
/// @param mask_extent c(xmin, xmax, ymin, ymax) extent of the mask raster
/// @param mask_nrow Number of rows in the mask raster
/// @param mask_ncol Number of columns in the mask raster
/// @param seed Integer seed for reproducibility; 0 means non-deterministic (uses entropy)
/// @param aux_logl_flat Flattened aux log-likelihood raster (k_steps * nrow * ncol, row-major); empty for no aux
/// @param aux_extent c(xmin, xmax, ymin, ymax) extent of the aux raster
/// @param aux_nrow Number of rows in the aux raster
/// @param aux_ncol Number of columns in the aux raster
/// @param flat_light_threshold Light-range threshold (max_obs - min_obs) below which a knot is
///   classified as informationally flat (no visible twilight). Diffusion is scaled by
///   `flat_light_scale` on such knots to prevent unconstrained wandering when the light
///   carries no positional signal. Default 10.0 (raw light units).
/// @param flat_light_scale Multiplicative scale applied to the diffusion coefficient on flat-light
///   knots (see `flat_light_threshold`). Values in (0, 1) constrain movement; 1.0 disables the
///   heuristic. Default 0.1 (10x tighter).
/// @name run_particle_filter
/// @export
#[extendr]
fn run_particle_filter(
    unix_times: &[f64],
    obs_light: &[f64],
    n_particles: i32,
    start_lat: f64,
    start_lon: f64,
    end_lat: f64,
    end_lon: f64,
    method: String,
    step_hours: f64,
    diffusion: Vec<f64>,
    trans_prob: Vec<f64>,
    calibration: Vec<f64>,
    likelihood_params: Vec<f64>,
    mask_matrix: Vec<f64>,
    mask_extent: Vec<f64>,
    mask_nrow: i32,
    mask_ncol: i32,
    seed: f64,
    aux_logl_flat: Vec<f64>,
    aux_extent: Vec<f64>,
    aux_nrow: i32,
    aux_ncol: i32,
    flat_light_threshold: f64,
    flat_light_scale: f64,
) -> List {
    let n = n_particles as usize;
    let num_obs = unix_times.len();
    let lambda = likelihood_params[0];
    let max_light = likelihood_params[1];
    let prob_slab = if likelihood_params.len() > 2 { likelihood_params[2] } else { 0.05 };
    let use_hyperprior = likelihood_params.len() >= 4;
    let prior_alpha = if use_hyperprior { likelihood_params[2] } else { 1.0 };
    let prior_beta = if use_hyperprior { likelihood_params[3] } else { 1.0 };
    let slab_density = 1.0 / max_light;
    let earth_radius = 6371.0; // km

    // Define Knots
    let mut k_steps = (((unix_times.last().unwrap() - unix_times[0]) / (step_hours * 3600.0)).ceil() as usize) + 1;
    let mut knot_times = vec![0.0; k_steps];
    let t_start = unix_times[0];
    let t_end = *unix_times.last().unwrap();
    let t_step = if k_steps > 1 { (t_end - t_start) / ((k_steps - 1) as f64) } else { 0.0 };
    
    for k in 0..k_steps {
        knot_times[k] = t_start + (k as f64) * t_step;
    }

    let mut rng: StdRng = if seed == 0.0 {
        StdRng::from_entropy()
    } else {
        StdRng::seed_from_u64(seed as u64)
    };
    
    // Forward History Storage (at knots)
    let mut hist_lat = vec![vec![0.0; n]; k_steps];
    let mut hist_lon = vec![vec![0.0; n]; k_steps];
    let mut hist_w = vec![vec![0.0; n]; k_steps];
    let mut hist_state = vec![vec![0; n]; k_steps];
    let mut hist_prob_slab = vec![vec![0.0; n]; k_steps];

    let mut particles = vec![
        Particle {
            lat: start_lat,
            lon: start_lon,
            weight: 1.0 / (n as f64),
            state: 0,
            prob_slab: 0.0,
        };
        n
    ];

    for i in 0..n {
        particles[i].prob_slab = if use_hyperprior {
            Beta::new(prior_alpha, prior_beta).unwrap().sample(&mut rng)
        } else {
            prob_slab
        };
        hist_lat[0][i] = start_lat;
        hist_lon[0][i] = start_lon;
        hist_w[0][i] = 1.0 / (n as f64);
        hist_state[0][i] = 0;
        hist_prob_slab[0][i] = particles[i].prob_slab;
    }

    let mut obs_idx = 0; // track which observation we are up to

    // Precompute the time-only solar ephemeris once per observation. Reused across
    // all particles in both the forward likelihood loop and the final diagnostics
    // pass; the per-particle code then only does the local hour-angle and cosine.
    let mut eph_sd = vec![0.0; num_obs];
    let mut eph_cd = vec![0.0; num_obs];
    let mut eph_ra = vec![0.0; num_obs];
    let mut eph_gmst = vec![0.0; num_obs];
    for j in 0..num_obs {
        let (sd, cd, ra, gmst) = solar_ephemeris(unix_times[j]);
        eph_sd[j] = sd; eph_cd[j] = cd; eph_ra[j] = ra; eph_gmst[j] = gmst;
    }

    // FORWARD FILTERING
    for k in 1..k_steps {
        let t_prev = knot_times[k-1];
        let t_curr = knot_times[k];
        let dt = (t_curr - t_prev) / 86400.0; // in days
        let time_remain = (t_end - t_prev) / 86400.0;
        let num_states = diffusion.len();
        
        let mut is_guided = false;

        if method == "guided" && !end_lat.is_nan() && !end_lon.is_nan() && time_remain > 0.0 {
            is_guided = true;
        }
        
        // Adaptive diffusion: if the light curve is flat (no twilight), constrain diffusion
        let mut min_obs = f64::INFINITY;
        let mut max_obs = f64::NEG_INFINITY;
        let mut obs_in_range = false;
        
        // Find observations in this knot segment
        let start_obs_idx = obs_idx;
        while obs_idx < num_obs && unix_times[obs_idx] <= t_curr {
            let o = obs_light[obs_idx];
            if o < min_obs { min_obs = o; }
            if o > max_obs { max_obs = o; }
            obs_in_range = true;
            obs_idx += 1;
        }
        
        let mut diff_scale = 1.0;
        if obs_in_range && (max_obs - min_obs) < flat_light_threshold {
            diff_scale = flat_light_scale;
        }
        let end_obs_idx = obs_idx;

        let mut log_weights = vec![0.0; n];
        let mut max_log_w = f64::NEG_INFINITY;

        for i in 0..n {
            let lat1 = particles[i].lat.to_radians();
            let lon1 = particles[i].lon.to_radians();
            
            let mut lat_mean = lat1;
            let mut lon_mean = lon1;
            
            // Sample new state
            let old_state = particles[i].state;
            let mut new_state = old_state;
            if num_states > 1 && trans_prob.len() >= num_states * num_states {
                let r: f64 = rng.gen();
                let mut cum_p = 0.0;
                let offset = old_state * num_states;
                for s in 0..num_states {
                    cum_p += trans_prob[offset + s];
                    if r <= cum_p {
                        new_state = s;
                        break;
                    }
                }
            }
            particles[i].state = new_state;

            if use_hyperprior {
                let jitter: f64 = rng.sample(Normal::new(0.0, 0.01).unwrap());
                particles[i].prob_slab = (particles[i].prob_slab + jitter).clamp(0.001, 0.999);
            }

            let sigma_u = diffusion[new_state] * dt.sqrt() * diff_scale;
            // Brownian-bridge variance: sigma_g^2 = sigma_u^2 * (time_remain - dt)/time_remain.
            // Degenerates to 0 when time_remain == dt (last step); fall back to unguided there.
            let sigma_g_sq = if is_guided {
                (diffusion[new_state] * diff_scale).powi(2)
                    * (dt * (time_remain - dt)).max(0.0) / time_remain
            } else {
                0.0
            };
            let do_guided = is_guided && sigma_g_sq > 1e-10;
            let sigma = if do_guided { sigma_g_sq.sqrt() } else { sigma_u };

            if do_guided {
                let end_lat_rad = end_lat.to_radians();
                let end_lon_rad = end_lon.to_radians();
                let a = ((end_lat_rad - lat1) / 2.0).sin().powi(2) +
                        lat1.cos() * end_lat_rad.cos() * ((end_lon_rad - lon1) / 2.0).sin().powi(2);
                let total_dist = 2.0 * earth_radius * a.sqrt().asin();
                let pull_frac = dt / time_remain;
                let pull_dist = total_dist * pull_frac;
                let pull_dr = pull_dist / earth_radius;

                let by = (end_lon_rad - lon1).sin() * end_lat_rad.cos();
                let bx = lat1.cos() * end_lat_rad.sin() - lat1.sin() * end_lat_rad.cos() * (end_lon_rad - lon1).cos();
                let bearing = by.atan2(bx);

                lat_mean = (lat1.sin() * pull_dr.cos() + lat1.cos() * pull_dr.sin() * bearing.cos()).asin();
                lon_mean = lon1 + (bearing.sin() * pull_dr.sin() * lat1.cos()).atan2(pull_dr.cos() - lat1.sin() * lat_mean.sin());
            }

            let theta: f64 = rng.gen::<f64>() * 2.0 * std::f64::consts::PI;
            let x: f64 = rng.sample(Normal::new(0.0, sigma).unwrap());
            let y: f64 = rng.sample(Normal::new(0.0, sigma).unwrap());
            let d = (x*x + y*y).sqrt();
            let dr = d / earth_radius;

            let lat2 = (lat_mean.sin() * dr.cos() + lat_mean.cos() * dr.sin() * theta.cos()).asin();
            let lon2 = lon_mean + (theta.sin() * dr.sin() * lat_mean.cos()).atan2(dr.cos() - lat_mean.sin() * lat2.sin());

            particles[i].lat = lat2.to_degrees();
            particles[i].lon = ((lon2.to_degrees() + 180.0) % 360.0 + 360.0) % 360.0 - 180.0;

            // Evaluate likelihood for all observations in segment
            let mut log_lik = 0.0;
            for j in start_obs_idx..end_obs_idx {
                let f = (unix_times[j] - t_prev) / (t_curr - t_prev);
                let p_lat = particles[i].lat * f + hist_lat[k-1][i] * (1.0 - f); // linear approx
                let p_lon = interpolate_lon(hist_lon[k-1][i], particles[i].lon, f);
                
                let zenith = zenith_from_ephem(eph_sd[j], eph_cd[j], eph_ra[j], eph_gmst[j], p_lon, p_lat);
                let expected = expected_light(zenith, &calibration, max_light);
                let obs = obs_light[j];
                let spike = spike_density(obs, expected, lambda, max_light);
                let current_prob_slab = particles[i].prob_slab;
                let den = (1.0 - current_prob_slab) * spike + current_prob_slab * slab_density;
                log_lik += den.ln();
            }

            if !aux_logl_flat.is_empty() && aux_extent.len() == 4 {
                let aux_ncols = aux_ncol as usize;
                let aux_nrows = aux_nrow as usize;
                let xmin = aux_extent[0]; let xmax = aux_extent[1];
                let ymax = aux_extent[3];
                let cell_w = (xmax - xmin) / aux_ncols as f64;
                let cell_h = (aux_extent[3] - aux_extent[2]) / aux_nrows as f64;
                let p_lon = particles[i].lon;
                let p_lat = particles[i].lat;
                let col_idx = ((p_lon - xmin) / cell_w).floor() as usize;
                let row_idx = ((ymax - p_lat) / cell_h).floor() as usize;
                let col_idx = col_idx.min(aux_ncols.saturating_sub(1));
                let row_idx = row_idx.min(aux_nrows.saturating_sub(1));
                let n_cells = aux_nrows * aux_ncols;
                log_lik += aux_logl_flat[k * n_cells + row_idx * aux_ncols + col_idx];
            }

            // IS correction for the guided proposal (Doob h-transform analog).
            // The guided kernel shifts mass toward the endpoint; without this correction
            // the weights contain only the likelihood and the estimate is biased.
            //
            // log(p_u(x_k|x_{k-1}) / q_g(x_k|x_{k-1},x_T))
            //   = log(σ_g²/σ_u²) + d_g²/(2σ_g²) − d_u²/(2σ_u²)
            //
            // d_g = d  (step from guided mean, km; constructed above)
            // d_u = great-circle distance from previous to new position (km)
            if do_guided {
                let cos_a = (lat2.sin() * lat1.sin()
                    + lat2.cos() * lat1.cos() * (lon2 - lon1).cos())
                    .clamp(-1.0, 1.0);
                let d_u = cos_a.acos() * earth_radius;
                let sigma_u_sq = sigma_u * sigma_u;
                let sigma_g_sq_v = sigma * sigma;
                log_lik += (sigma_g_sq_v / sigma_u_sq).ln()
                    + d * d / (2.0 * sigma_g_sq_v)
                    - d_u * d_u / (2.0 * sigma_u_sq);
            }

            log_weights[i] = log_lik;
            if log_lik > max_log_w {
                max_log_w = log_lik;
            }
        }

        let mut total_weight = 0.0;
        for i in 0..n {
            particles[i].weight *= (log_weights[i] - max_log_w).exp();
            total_weight += particles[i].weight;
        }

        if total_weight > 0.0 {
            for i in 0..n {
                particles[i].weight /= total_weight;
                hist_lat[k][i] = particles[i].lat;
                hist_lon[k][i] = particles[i].lon;
                hist_w[k][i] = particles[i].weight;
                hist_state[k][i] = particles[i].state;
                hist_prob_slab[k][i] = particles[i].prob_slab;
            }
        } else {
            for i in 0..n { 
                particles[i].weight = 1.0 / (n as f64);
                hist_lat[k][i] = particles[i].lat;
                hist_lon[k][i] = particles[i].lon; 
                hist_w[k][i] = 1.0 / (n as f64);
                hist_state[k][i] = particles[i].state;
                hist_prob_slab[k][i] = particles[i].prob_slab;
            }
        }

        let ess = 1.0 / particles.iter().map(|p| p.weight.powi(2)).sum::<f64>();
        if ess < (n as f64) / 2.0 {
            let mut new_particles = Vec::with_capacity(n);
            let mut sum = 0.0;
            let mut cum_weights = Vec::with_capacity(n);
            for p in &particles {
                sum += p.weight;
                cum_weights.push(sum);
            }
            for _ in 0..n {
                let r = rng.gen::<f64>();
                let idx = match cum_weights.binary_search_by(|w| w.partial_cmp(&r).unwrap()) {
                    Ok(i) => i, Err(i) => i,
                }.min(n-1);
                let mut p = particles[idx];
                p.weight = 1.0 / (n as f64);
                new_particles.push(p);
            }
            particles = new_particles;
            
            // Need to update the history at step k to match resampled particles
            for i in 0..n {
                hist_lat[k][i] = particles[i].lat;
                hist_lon[k][i] = particles[i].lon;
                hist_w[k][i] = particles[i].weight;
                hist_state[k][i] = particles[i].state;
                hist_prob_slab[k][i] = particles[i].prob_slab;
            }
        }
    }

    // SMOOTHING PASS (if ffbs)
    let mut smooth_lat = vec![vec![0.0; n]; k_steps];
    let mut smooth_lon = vec![vec![0.0; n]; k_steps];
    let mut smooth_state = vec![vec![0; n]; k_steps];
    let mut smooth_prob_slab = vec![vec![0.0; n]; k_steps];

    if method == "ffbs" {
        if !end_lat.is_nan() && !end_lon.is_nan() {
            for j in 0..n {
                smooth_lat[k_steps-1][j] = end_lat;
                smooth_lon[k_steps-1][j] = end_lon;
            }
        } else {
            let mut cum_w = Vec::with_capacity(n);
            let mut sum = 0.0;
            for i in 0..n {
                sum += hist_w[k_steps-1][i];
                cum_w.push(sum);
            }
            for j in 0..n {
                let r = rng.gen::<f64>();
                let idx = match cum_w.binary_search_by(|w| w.partial_cmp(&r).unwrap()) {
                    Ok(i) => i, Err(i) => i,
                }.min(n-1);
                smooth_lat[k_steps-1][j] = hist_lat[k_steps-1][idx];
                smooth_lon[k_steps-1][j] = hist_lon[k_steps-1][idx];
                smooth_state[k_steps-1][j] = hist_state[k_steps-1][idx];
                smooth_prob_slab[k_steps-1][j] = hist_prob_slab[k_steps-1][idx];
            }
        }

        for k in (0..k_steps-1).rev() {
            let dt = (knot_times[k+1] - knot_times[k]) / 86400.0;

            for j in 0..n {
                let next_lat = smooth_lat[k+1][j].to_radians();
                let next_lon = smooth_lon[k+1][j].to_radians();

                let mut back_w = Vec::with_capacity(n);
                let mut sum_w = 0.0;

                for i in 0..n {
                    let cur_lat = hist_lat[k][i].to_radians();
                    let cur_lon = hist_lon[k][i].to_radians();
                    let a = ((next_lat - cur_lat) / 2.0).sin().powi(2) +
                            cur_lat.cos() * next_lat.cos() * ((next_lon - cur_lon) / 2.0).sin().powi(2);
                    let dist = 2.0 * earth_radius * a.sqrt().asin();
                    
                    let st = smooth_state[k+1][j];
                    let sigma = diffusion[st] * dt.sqrt();
                    let var2 = 2.0 * sigma * sigma;
                    let mut trans_prob_val = (- (dist * dist) / var2).exp();
                    
                    if diffusion.len() > 1 && trans_prob.len() >= diffusion.len() * diffusion.len() {
                        let old_st = hist_state[k][i];
                        trans_prob_val *= trans_prob[old_st * diffusion.len() + st];
                    }
                    
                    let w = hist_w[k][i] * trans_prob_val;
                    sum_w += w;
                    back_w.push(sum_w);
                }

                if sum_w > 0.0 {
                    let r = rng.gen::<f64>() * sum_w;
                    let idx = match back_w.binary_search_by(|w| w.partial_cmp(&r).unwrap()) {
                        Ok(i) => i, Err(i) => i,
                    }.min(n-1);
                    smooth_lat[k][j] = hist_lat[k][idx];
                    smooth_lon[k][j] = hist_lon[k][idx];
                    smooth_state[k][j] = hist_state[k][idx];
                    smooth_prob_slab[k][j] = hist_prob_slab[k][idx];
                } else {
                    smooth_lat[k][j] = hist_lat[k][j];
                    smooth_lon[k][j] = hist_lon[k][j];
                    smooth_state[k][j] = hist_state[k][j];
                    smooth_prob_slab[k][j] = hist_prob_slab[k][j];
                }
            }
        }
    } else {
        smooth_lat = hist_lat;
        smooth_lon = hist_lon;
        smooth_state = hist_state;
        smooth_prob_slab = hist_prob_slab;
        if method == "guided" && !end_lat.is_nan() && !end_lon.is_nan() {
            for j in 0..n {
                smooth_lat[k_steps-1][j] = end_lat;
                smooth_lon[k_steps-1][j] = end_lon;
            }
        }
    }

    // Compute Knot Statistics
    let mut knot_lat = Vec::with_capacity(k_steps);
    let mut knot_lon = Vec::with_capacity(k_steps);
    let mut knot_lat_sd = Vec::with_capacity(k_steps);
    let mut knot_lon_sd = Vec::with_capacity(k_steps);
    let mut knot_prob_slab = Vec::with_capacity(k_steps);
    let mut knot_prob_state = vec![Vec::with_capacity(k_steps); diffusion.len()];

    for k in 0..k_steps {
        let mut m_lat = 0.0;
        let mut sum_x = 0.0;
        let mut sum_y = 0.0;
        let mut m_prob_slab = 0.0;
        let weight = if method == "ffbs" { 1.0 / (n as f64) } else { 0.0 }; // If guided/forward, we use hist_w

        for j in 0..n {
            let w = if method == "ffbs" { weight } else { hist_w[k][j] };
            m_lat += smooth_lat[k][j] * w;
            let lon_rad = smooth_lon[k][j].to_radians();
            sum_x += lon_rad.cos() * w;
            sum_y += lon_rad.sin() * w;
            m_prob_slab += smooth_prob_slab[k][j] * w;
        }
        let m_lon = sum_y.atan2(sum_x).to_degrees();
        
        knot_lat.push(m_lat);
        knot_lon.push(m_lon);
        knot_prob_slab.push(m_prob_slab);

        let mut s_lat = 0.0;
        let mut s_lon = 0.0;
        for j in 0..n {
            let w = if method == "ffbs" { weight } else { hist_w[k][j] };
            s_lat += (smooth_lat[k][j] - m_lat).powi(2) * w;
            let mut diff = smooth_lon[k][j] - m_lon;
            if diff > 180.0 { diff -= 360.0; }
            if diff < -180.0 { diff += 360.0; }
            s_lon += diff.powi(2) * w;
        }
        knot_lat_sd.push(s_lat.sqrt());
        knot_lon_sd.push(s_lon.sqrt());

        for s in 0..diffusion.len() {
            let mut p_state = 0.0;
            for j in 0..n {
                let w = if method == "ffbs" { weight } else { hist_w[k][j] };
                if smooth_state[k][j] == s {
                    p_state += w;
                }
            }
            knot_prob_state[s].push(p_state);
        }
    }

    // Final Pass: Compute Diagnostics at every observation point using the smoothed track
    let mut obs_zenith = Vec::with_capacity(num_obs);
    let mut obs_prob_false = Vec::with_capacity(num_obs);
    
    let mut k = 1;
    for j in 0..num_obs {
        while k < k_steps - 1 && unix_times[j] > knot_times[k] {
            k += 1;
        }
        let t_prev = knot_times[k-1];
        let t_curr = knot_times[k];
        let f = if t_curr > t_prev { (unix_times[j] - t_prev) / (t_curr - t_prev) } else { 0.0 };
        
        let mut p_false = 0.0;
        let mut mean_z = 0.0;
        
        for i in 0..n {
            let p_lat = smooth_lat[k][i] * f + smooth_lat[k-1][i] * (1.0 - f);
            let p_lon = interpolate_lon(smooth_lon[k-1][i], smooth_lon[k][i], f);
            let w = if method == "ffbs" { 1.0 / (n as f64) } else { hist_w[k][i] };
            
            let current_prob_slab = smooth_prob_slab[k][i] * f + smooth_prob_slab[k-1][i] * (1.0 - f);
            
            let z = zenith_from_ephem(eph_sd[j], eph_cd[j], eph_ra[j], eph_gmst[j], p_lon, p_lat);
            mean_z += z * w;
            
            let exp = expected_light(z, &calibration, max_light);
            let obs = obs_light[j];
            let spike = spike_density(obs, exp, lambda, max_light);
            let den = (1.0 - current_prob_slab) * spike + current_prob_slab * slab_density;
            let prob_f = (current_prob_slab * slab_density) / den;
            p_false += prob_f * w;
        }
        obs_zenith.push(mean_z);
        obs_prob_false.push(p_false);
    }

    let mut prob_state_list = List::new(diffusion.len());
    for s in 0..diffusion.len() {
        prob_state_list.set_elt(s, knot_prob_state[s].clone().into()).unwrap();
    }

    list!(
        knot_times = knot_times,
        lat = knot_lat, 
        lon = knot_lon,
        lat_sd = knot_lat_sd,
        lon_sd = knot_lon_sd,
        prob_state = prob_state_list,
        prob_slab = knot_prob_slab,
        obs_times = unix_times,
        obs_zenith = obs_zenith,
        prob_false = obs_prob_false
    )
}

/// Evaluate the spike-and-slab log-likelihood over a grid of locations
///
/// For each candidate location, returns the summed log-likelihood of the
/// observed light series under the continuous spike-and-slab model. Useful for
/// visualising the likelihood surface independently of the HMM smoother.
///
/// @param lon Longitudes of grid cells (degrees)
/// @param lat Latitudes of grid cells (degrees)
/// @param unix_times Observation timestamps (seconds since 1970-01-01)
/// @param obs_light Observed light values
/// @param calibration Response parameters: `c(intercept, slope)` for the clamped-linear
///   model, or `c(floor, amp, z50, scale)` for the logistic model (see details in
///   `fit_light_response`). Length selects the model.
/// @param likelihood_params c(lambda, max_light, prob_slab) or c(lambda, max_light, alpha, beta)
/// @return Numeric vector of log-likelihoods, one per grid cell
/// @name eval_logpk_grid
/// @export
#[extendr]
fn eval_logpk_grid(
    lon: &[f64],
    lat: &[f64],
    unix_times: &[f64],
    obs_light: &[f64],
    calibration: Vec<f64>,
    likelihood_params: Vec<f64>
) -> Vec<f64> {
    let n = lon.len();
    let num_obs = unix_times.len();
    let lambda = likelihood_params[0];
    let max_light = likelihood_params[1];
    let prob_slab = if likelihood_params.len() > 3 {
        let alpha = likelihood_params[2];
        let beta = likelihood_params[3];
        alpha / (alpha + beta)
    } else {
        likelihood_params[2]
    };
    let slab_density = 1.0 / max_light;

    let mut logl = vec![0.0; n];

    // Precompute time-only ephemeris once per observation, reused across grid cells.
    let eph: Vec<(f64, f64, f64, f64)> =
        (0..num_obs).map(|j| solar_ephemeris(unix_times[j])).collect();

    for i in 0..n {
        let mut sum_logl = 0.0;
        for j in 0..num_obs {
            let (sd, cd, ra, gmst) = eph[j];
            let zenith = zenith_from_ephem(sd, cd, ra, gmst, lon[i], lat[i]);
            let expected = expected_light(zenith, &calibration, max_light);
            let obs = obs_light[j];
            let spike = spike_density(obs, expected, lambda, max_light);
            let den = (1.0 - prob_slab) * spike + prob_slab * slab_density;
            sum_logl += den.ln();
        }
        logl[i] = sum_logl;
    }
    
    logl
}

#[extendr]
fn run_grid_hmm(
    lon: &[f64],
    lat: &[f64],
    knot_times: &[f64],
    obs_times: &[f64],
    obs_light: &[f64],
    fixed_idx: &[i32],
    fixed_lon: &[f64],
    fixed_lat: &[f64],
    diffusion: Vec<f64>,
    trans_prob: Vec<f64>,
    calibration: Vec<f64>,
    likelihood_params: Vec<f64>,
    aux_logl: Vec<f64>,
) -> List {
    let n = lon.len();
    let num_states = diffusion.len();
    let k_steps = knot_times.len();
    
    let lambda = likelihood_params[0];
    let max_light = likelihood_params[1];
    let prob_slab = if likelihood_params.len() > 3 {
        let alpha = likelihood_params[2];
        let beta = likelihood_params[3];
        alpha / (alpha + beta)
    } else {
        likelihood_params[2]
    };
    let slab_density = 1.0 / max_light;
    let earth_radius = 6371.0;

    let lon_rad: Vec<f64> = lon.iter().map(|x| x.to_radians()).collect();
    let lat_rad: Vec<f64> = lat.iter().map(|x| x.to_radians()).collect();

    let mut logpk = vec![vec![0.0; n]; k_steps];

    // Precompute time-only ephemeris once per observation, reused across grid cells.
    let eph: Vec<(f64, f64, f64, f64)> =
        (0..obs_times.len()).map(|j| solar_ephemeris(obs_times[j])).collect();

    // Auxiliary location terms (priors, SST, bathymetry, masks) precomputed in R
    // as a k_steps x n matrix, flattened row-major (index k*n + i). Seeded into
    // every cell at every knot, including knots with no light observations, so a
    // prior applies regardless of the light record. A -Inf entry hard-masks the
    // cell. This runs BEFORE the light loop so a hard-masked cell can skip the
    // per-observation light likelihood, which is the dominant cost.
    if !aux_logl.is_empty() {
        for k in 0..k_steps {
            for i in 0..n {
                logpk[k][i] = aux_logl[k * n + i];
            }
        }
    }

    for k in 0..k_steps {
        let t_curr = knot_times[k];
        let t_prev = if k == 0 { t_curr - (knot_times[1] - knot_times[0]) } else { knot_times[k-1] };

        let mut obs_in_k = Vec::new();
        for j in 0..obs_times.len() {
            if obs_times[j] > t_prev && obs_times[j] <= t_curr {
                obs_in_k.push(j);
            }
        }

        if obs_in_k.is_empty() { continue; }

        for i in 0..n {
            // Already ruled out by a hard auxiliary constraint: the light
            // likelihood cannot rescue a -inf, so do not compute it.
            if logpk[k][i] <= -1e29 { continue; }
            let mut sum_logl = 0.0;
            for &j in &obs_in_k {
                let (sd, cd, ra, gmst) = eph[j];
                let zenith = zenith_from_ephem(sd, cd, ra, gmst, lon[i], lat[i]);
                let expected = expected_light(zenith, &calibration, max_light);
                let obs = obs_light[j];
                let spike = spike_density(obs, expected, lambda, max_light);
                let den = (1.0 - prob_slab) * spike + prob_slab * slab_density;
                sum_logl += den.ln();
            }
            logpk[k][i] += sum_logl;
        }
    }

    for idx in 0..fixed_idx.len() {
        let k = fixed_idx[idx] as usize;
        let f_lon = fixed_lon[idx];
        let f_lat = fixed_lat[idx];
        
        let mut best_i = 0;
        let mut best_dist = f64::INFINITY;
        for i in 0..n {
            let a = ((lat_rad[i] - f_lat.to_radians()) / 2.0).sin().powi(2) +
                    f_lat.to_radians().cos() * lat_rad[i].cos() * ((lon_rad[i] - f_lon.to_radians()) / 2.0).sin().powi(2);
            let dist = 2.0 * earth_radius * a.sqrt().asin();
            if dist < best_dist {
                best_dist = dist;
                best_i = i;
            }
        }
        for i in 0..n {
            if i != best_i {
                logpk[k][i] = -1e30;
            }
        }
    }

    let mut alpha = vec![vec![-1e30; n * num_states]; k_steps];
    for i in 0..n {
        for s in 0..num_states {
            if logpk[0][i] > -1e29 {
                alpha[0][i * num_states + s] = logpk[0][i] - (n as f64 * num_states as f64).ln();
            }
        }
    }

    let mut dt_array = vec![0.0; k_steps];
    for k in 1..k_steps {
        dt_array[k] = (knot_times[k] - knot_times[k-1]) / 86400.0;
    }

    for k in 1..k_steps {
        let dt = dt_array[k];
        let mut max_sigma = 0.0;
        for s in 0..num_states {
            let sig = diffusion[s] * dt.sqrt();
            if sig > max_sigma { max_sigma = sig; }
        }
        let threshold_dist = max_sigma * 5.0; 
        
        for i in 0..n {
            if logpk[k][i] <= -1e29 { continue; } 
            
            for s in 0..num_states {
                let sigma = diffusion[s] * dt.sqrt();
                let var2 = 2.0 * sigma * sigma;
                let log_norm_const = - 2.0 * (sigma).ln(); 
                
                let mut max_val = -1e30;
                let mut sum_exp = 0.0;
                
                for j in 0..n {
                    let lat_diff = (lat[i] - lat[j]).abs();
                    if lat_diff > threshold_dist / 111.0 { continue; }
                    
                    let lon_diff = (lon[i] - lon[j]).abs();
                    // simple conservative wrap-around or local bound
                    let lon_diff_wrap = if lon_diff > 180.0 { 360.0 - lon_diff } else { lon_diff };
                    let max_cos = lat_rad[i].cos().min(lat_rad[j].cos()).max(0.1);
                    if lon_diff_wrap > threshold_dist / (111.0 * max_cos) { continue; }

                    let a = ((lat_rad[i] - lat_rad[j]) / 2.0).sin().powi(2) +
                            lat_rad[j].cos() * lat_rad[i].cos() * ((lon_rad[i] - lon_rad[j]) / 2.0).sin().powi(2);
                    let dist = 2.0 * earth_radius * a.sqrt().asin();
                    
                    if dist > threshold_dist { continue; }
                    
                    let log_spatial = - (dist * dist) / var2 + log_norm_const;
                    
                    for s_prev in 0..num_states {
                        let alpha_prev = alpha[k-1][j * num_states + s_prev];
                        if alpha_prev > -1e29 {
                            let log_t = if num_states > 1 { trans_prob[s_prev * num_states + s].ln() } else { 0.0 };
                            let val = alpha_prev + log_t + log_spatial;
                            
                            if val > max_val {
                                sum_exp = sum_exp * (max_val - val).exp() + 1.0;
                                max_val = val;
                            } else {
                                sum_exp += (val - max_val).exp();
                            }
                        }
                    }
                }
                
                if sum_exp > 0.0 {
                    alpha[k][i * num_states + s] = logpk[k][i] + max_val + sum_exp.ln();
                }
            }
        }
    }

    // Log marginal likelihood (model evidence): logsumexp of the final forward
    // column. alpha is the unnormalised log forward probability (the k=0 entries
    // already carry the uniform 1/(n*num_states) prior), so this is the log of
    // the total probability of the observations under the model and any priors
    // injected via aux_logl. Use it for model comparison, e.g. a hemisphere
    // Bayes factor: run twice under opposite priors and difference the logZ.
    let log_z = {
        let last = &alpha[k_steps - 1];
        let mut max_a = f64::NEG_INFINITY;
        for &a in last.iter() {
            if a > -1e29 && a > max_a { max_a = a; }
        }
        if max_a == f64::NEG_INFINITY {
            f64::NEG_INFINITY
        } else {
            let mut s = 0.0;
            for &a in last.iter() {
                if a > -1e29 { s += (a - max_a).exp(); }
            }
            max_a + s.ln()
        }
    };

    let mut beta = vec![vec![-1e30; n * num_states]; k_steps];
    for i in 0..n {
        for s in 0..num_states {
            if logpk[k_steps-1][i] > -1e29 {
                beta[k_steps-1][i * num_states + s] = 0.0;
            }
        }
    }

    for k in (0..k_steps-1).rev() {
        let dt = dt_array[k+1];
        let mut max_sigma = 0.0;
        for s in 0..num_states {
            let sig = diffusion[s] * dt.sqrt();
            if sig > max_sigma { max_sigma = sig; }
        }
        let threshold_dist = max_sigma * 5.0; 
        
        for i in 0..n {
            if logpk[k][i] <= -1e29 { continue; }
            
            for s in 0..num_states {
                let mut max_val = -1e30;
                let mut sum_exp = 0.0;
                
                for j in 0..n {
                    if logpk[k+1][j] <= -1e29 { continue; }
                    
                    let lat_diff = (lat[j] - lat[i]).abs();
                    if lat_diff > threshold_dist / 111.0 { continue; }
                    
                    let lon_diff = (lon[j] - lon[i]).abs();
                    let lon_diff_wrap = if lon_diff > 180.0 { 360.0 - lon_diff } else { lon_diff };
                    let max_cos = lat_rad[i].cos().min(lat_rad[j].cos()).max(0.1);
                    if lon_diff_wrap > threshold_dist / (111.0 * max_cos) { continue; }

                    let a = ((lat_rad[j] - lat_rad[i]) / 2.0).sin().powi(2) +
                            lat_rad[i].cos() * lat_rad[j].cos() * ((lon_rad[j] - lon_rad[i]) / 2.0).sin().powi(2);
                    let dist = 2.0 * earth_radius * a.sqrt().asin();
                    
                    if dist > threshold_dist { continue; }
                    
                    for s_next in 0..num_states {
                        let beta_next = beta[k+1][j * num_states + s_next];
                        if beta_next > -1e29 {
                            let sigma = diffusion[s_next] * dt.sqrt();
                            let var2 = 2.0 * sigma * sigma;
                            let log_norm_const = - 2.0 * (sigma).ln();
                            
                            let log_spatial = - (dist * dist) / var2 + log_norm_const;
                            let log_t = if num_states > 1 { trans_prob[s * num_states + s_next].ln() } else { 0.0 };
                            
                            let val = beta_next + log_t + log_spatial + logpk[k+1][j];
                            
                            if val > max_val {
                                sum_exp = sum_exp * (max_val - val).exp() + 1.0;
                                max_val = val;
                            } else {
                                sum_exp += (val - max_val).exp();
                            }
                        }
                    }
                }
                
                if sum_exp > 0.0 {
                    beta[k][i * num_states + s] = max_val + sum_exp.ln();
                }
            }
        }
    }

    let mut best_lat = vec![0.0; k_steps];
    let mut best_lon = vec![0.0; k_steps];
    let mut knot_prob_state = vec![Vec::with_capacity(k_steps); num_states];
    // Per-knot posterior over cells (marginalised over states), row-major k*n + i.
    // Exposed for uncertainty diagnostics and simulation-based calibration.
    let mut posterior = vec![0.0; k_steps * n];

    for k in 0..k_steps {
        let mut max_gamma = -1e30;
        let mut sum_gamma_exp = 0.0;
        let mut gamma = vec![0.0; n * num_states];
        
        for i in 0..n {
            for s in 0..num_states {
                let a_val = alpha[k][i * num_states + s];
                let b_val = beta[k][i * num_states + s];
                
                if a_val > -1e29 && b_val > -1e29 {
                    let g = a_val + b_val;
                    gamma[i * num_states + s] = g;
                    
                    if g > max_gamma {
                        sum_gamma_exp = sum_gamma_exp * (max_gamma - g).exp() + 1.0;
                        max_gamma = g;
                    } else {
                        sum_gamma_exp += (g - max_gamma).exp();
                    }
                } else {
                    gamma[i * num_states + s] = -1e30;
                }
            }
        }
        
        let mut prob_states = vec![0.0; num_states];
        let mut best_i = 0;
        let mut best_g = -1e30;

        if sum_gamma_exp > 0.0 {
            let log_denom = max_gamma + sum_gamma_exp.ln();
            for i in 0..n {
                for s in 0..num_states {
                    if gamma[i * num_states + s] > -1e29 {
                        let w = (gamma[i * num_states + s] - log_denom).exp();
                        prob_states[s] += w;
                        posterior[k * n + i] += w;        // marginal posterior of cell i at knot k

                        let g_state = gamma[i * num_states + s];
                        if g_state > best_g {
                            best_g = g_state;
                            best_i = i;
                        }
                    }
                }
            }
        } else {
            // Marginal collapsed (movement sigma too small for grid resolution).
            // Fall back to the light-only MAP cell so output is at least
            // in the right hemisphere rather than silently returning cell 0.
            let mut best_logl = -1e30;
            for i in 0..n {
                if logpk[k][i] > best_logl {
                    best_logl = logpk[k][i];
                    best_i = i;
                }
            }
        }

        best_lat[k] = lat[best_i];
        best_lon[k] = lon[best_i];
        for s in 0..num_states {
            knot_prob_state[s].push(prob_states[s]);
        }
    }

    let mut prob_state_list = List::new(num_states);
    for s in 0..num_states {
        prob_state_list.set_elt(s, knot_prob_state[s].clone().into()).unwrap();
    }

    list!(
        time = knot_times,
        lat = best_lat,
        lon = best_lon,
        prob_state = prob_state_list,
        log_z = log_z,
        posterior = posterior,
        n_cells = n as i32
    )
}

// =====================================================================
// Native single-track BLOCK + red-black-POLISH sampler kernel.
//
// Port of the validated R prototype (notes/da_block/da_hier.R). The coarse-HMM
// surrogate (per-knot Gaussian mean mu_k and 2x2 precision P_k) and the movement
// precision P_move are built in R and passed in; this runs the hot sweep loop
// natively. The delayed-acceptance correction's Gaussian normalisers cancel per
// knot between proposal and current, so only mu_k and P_k are needed (no
// determinants). Reuses solar_ephemeris/zenith_from_ephem/spike_density.
// =====================================================================

// 2x2 row-major helpers: m = [m00, m01, m10, m11]
#[inline]
fn m2_mulv(m: &[f64; 4], x: f64, y: f64) -> (f64, f64) { (m[0]*x + m[1]*y, m[2]*x + m[3]*y) }
#[inline]
fn m2_inv(m: &[f64; 4]) -> [f64; 4] {
    let idet = 1.0 / (m[0]*m[3] - m[1]*m[2]);
    [m[3]*idet, -m[1]*idet, -m[2]*idet, m[0]*idet]
}
#[inline]
fn m2_chol_lower(m: &[f64; 4]) -> [f64; 4] {  // L with L L^T = m; [l00, 0, l10, l11]
    let l00 = m[0].max(1e-12).sqrt();
    let l10 = m[2] / l00;
    let l11 = (m[3] - l10*l10).max(1e-12).sqrt();
    [l00, 0.0, l10, l11]
}
#[inline]
fn quad2(p: &[f64; 4], dx: f64, dy: f64) -> f64 { p[0]*dx*dx + 2.0*p[1]*dx*dy + p[3]*dy*dy }

// Great-circle (haversine) distance in km.
#[inline]
fn gcdist_km(lon1: f64, lat1: f64, lon2: f64, lat2: f64) -> f64 {
    let r = 6371.0_f64;
    let (p1, p2) = (lat1.to_radians(), lat2.to_radians());
    let (dlat, dlon) = ((lat2 - lat1).to_radians(), (lon2 - lon1).to_radians());
    let a = (dlat*0.5).sin().powi(2) + p1.cos()*p2.cos()*(dlon*0.5).sin().powi(2);
    2.0 * r * a.sqrt().asin()
}
// One edge's (great-circle^2 - planar^2) in km^2. planar^2 uses the tag's fixed
// reference metric (km_lon2, km_lat2) so it matches the linear proposal's density;
// the DA correction adds (this)/(-2 sig2) to upgrade the target prior to great-circle.
#[inline]
fn mv_edge(km_lon2: f64, km_lat2: f64, lon1: f64, lat1: f64, lon2: f64, lat2: f64) -> f64 {
    let gc = gcdist_km(lon1, lat1, lon2, lat2);
    let (dlon, dlat) = (lon2 - lon1, lat2 - lat1);
    gc*gc - (km_lon2*dlon*dlon + km_lat2*dlat*dlat)
}

// One innovation's contribution to (CRW - Brownian) in the log-prior numerator.
//
// The Brownian prior penalises |delta_t|^2; the correlated random walk penalises
// the INNOVATION |delta_t - rho*delta_{t-1}|^2. The difference is what has to be
// added to an acceptance ratio whose proposal already carries the Brownian term,
// which is how the spherical metric is handled too. Increments are in km via the
// tag's reference metric, so this composes with the spherical correction rather
// than duplicating it.
#[inline]
#[allow(clippy::too_many_arguments)]
fn crw_innov(km_lon2: f64, km_lat2: f64, rho: f64,
             lon0: f64, lat0: f64, lon1: f64, lat1: f64, lon2: f64, lat2: f64) -> f64 {
    let (d1x, d1y) = (lon1 - lon0, lat1 - lat0);
    let (d2x, d2y) = (lon2 - lon1, lat2 - lat1);
    let dot  = km_lon2*d2x*d1x + km_lat2*d2y*d1y;
    let prev = km_lon2*d1x*d1x + km_lat2*d1y*d1y;
    -2.0*rho*dot + rho*rho*prev
}

// Sum crw_innov over the innovations a move touches. Innovation t is built from
// knots t-2, t-1 and t, so moving knot m disturbs innovations m, m+1 and m+2.
// Defined for 2 ..= k-1.
fn crw_sum<F: Fn(usize) -> (f64, f64)>(ctx: &Ctx, pos: F, k: usize,
                                       lo: usize, hi: usize) -> f64 {
    let t0 = if lo < 2 { 2 } else { lo };
    let t1 = if hi > k - 1 { k - 1 } else { hi };
    let mut acc = 0.0;
    let mut t = t0;
    while t <= t1 {
        let a = pos(t - 2); let b = pos(t - 1); let c = pos(t);
        acc += crw_innov(ctx.km_lon2, ctx.km_lat2, ctx.rho.get(),
                         a.0, a.1, b.0, b.1, c.0, c.1);
        t += 1;
    }
    acc
}

// Dense upper Cholesky R (row-major n x n), R^T R = Q.
fn chol_upper(q: &[f64], n: usize) -> Vec<f64> {
    let mut r = vec![0.0f64; n*n];
    for j in 0..n {
        let mut d = q[j*n + j];
        for k in 0..j { d -= r[k*n + j]*r[k*n + j]; }
        let rjj = d.max(1e-12).sqrt();
        r[j*n + j] = rjj;
        for i in (j+1)..n {
            let mut s = q[j*n + i];
            for k in 0..j { s -= r[k*n + j]*r[k*n + i]; }
            r[j*n + i] = s / rjj;
        }
    }
    r
}
fn solve_lt(r: &[f64], b: &[f64], n: usize) -> Vec<f64> {  // R^T y = b (forward)
    let mut y = vec![0.0; n];
    for i in 0..n {
        let mut s = b[i];
        for k in 0..i { s -= r[k*n + i]*y[k]; }
        y[i] = s / r[i*n + i];
    }
    y
}
fn solve_ut(r: &[f64], y: &[f64], n: usize) -> Vec<f64> {  // R x = y (back)
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let mut s = y[i];
        for k in (i+1)..n { s -= r[i*n + k]*x[k]; }
        x[i] = s / r[i*n + i];
    }
    x
}

// Per-individual likelihood context. `aux` is an optional additive per-knot
// location log-likelihood (sensor terms: priors, SST, masks, ...), evaluated by
// nearest-cell lookup on a regular lon/lat grid shared by all of this tag's knots.
// koff maps a global knot index to this tag's local knot (for aux indexing).
struct Ctx<'a> {
    obs_start: &'a [i32], obs_len: &'a [i32],
    obs_light: &'a [f64], eph: &'a [(f64, f64, f64, f64)],
    // Response parameters, length 2 (clamped linear) or 4 (logistic); see
    // expected_light(). Carried as a fixed array so Ctx stays Copy-cheap.
    cal: [f64; 4], cal_len: usize,
    lambda: f64, max_light: f64, prob_slab: f64,
    koff: usize,
    aux: &'a [f64],       // this tag's K*ncell aux values (empty if no terms)
    aux_lon0: f64, aux_dlon: f64, aux_ncol: usize,
    aux_lat0: f64, aux_dlat: f64, aux_nrow: usize,
    spherical: bool, km_lon2: f64, km_lat2: f64,   // great-circle movement metric
    // Directional persistence of the correlated random walk: increments follow
    // an AR(1), delta_t = rho*delta_{t-1} + eps. rho = 0 is the memoryless
    // Brownian model and reproduces it exactly. A Cell because rho is resampled
    // every sweep while the context is borrowed immutably by the movers.
    rho: std::cell::Cell<f64>,
}
impl<'a> Ctx<'a> {
    #[inline]
    fn emit(&self, knot: usize, lon: f64, lat: f64) -> f64 {
        // Auxiliary terms FIRST. A hard constraint (a -inf from mask_rule() or a
        // hard floor_rule()) settles the emission on its own, so the light loop
        // over this knot's observations is skipped entirely. This is what makes
        // an infeasible location cheap rather than merely improbable: the light
        // likelihood is the dominant cost, and there is no point paying it for a
        // position the bathymetry or the coastline has already excluded.
        let mut ll = 0.0;
        if !self.aux.is_empty() {
            let li = (((lon - self.aux_lon0) / self.aux_dlon).round() as i64)
                .clamp(0, self.aux_ncol as i64 - 1) as usize;
            let lj = (((lat - self.aux_lat0) / self.aux_dlat).round() as i64)
                .clamp(0, self.aux_nrow as i64 - 1) as usize;
            let ncell = self.aux_ncol * self.aux_nrow;
            let local = knot - self.koff;
            ll = self.aux[local * ncell + lj * self.aux_ncol + li];
            if ll <= -1e29 { return f64::NEG_INFINITY; }
        }
        let s = self.obs_start[knot] as usize;
        let n = self.obs_len[knot] as usize;
        let slab = 1.0 / self.max_light;
        for j in s..(s + n) {
            let (sd, cd, ra, gmst) = self.eph[j];
            let z = zenith_from_ephem(sd, cd, ra, gmst, lon, lat);
            let expc = expected_light(z, &self.cal[..self.cal_len], self.max_light);
            let spike = spike_density(self.obs_light[j], expc, self.lambda, self.max_light);
            ll += ((1.0 - self.prob_slab)*spike + self.prob_slab*slab).ln();
        }
        ll
    }
}

// Track moves operate on a per-track local state (xlon/xlat/le indexed 0..k-1);
// koff maps a local knot t to its global index koff+t in the surrogate/obs arrays,
// so the same code serves one track or an individual inside the hierarchy.
// le[t] caches emit(koff+t, current x[t]); only proposals need a fresh emit.
fn single_site(ctx: &Ctx, xlon: &mut [f64], xlat: &mut [f64], le: &mut [f64], koff: usize,
               t: usize, k: usize, pm: &[f64; 4], rng: &mut StdRng, nrm: &Normal<f64>) {
    let (mx, my, prec) = if t == k - 1 {
        (xlon[t-1], xlat[t-1], *pm)
    } else {
        ((xlon[t-1]+xlon[t+1])*0.5, (xlat[t-1]+xlat[t+1])*0.5,
         [2.0*pm[0], 2.0*pm[1], 2.0*pm[2], 2.0*pm[3]])
    };
    let l = m2_chol_lower(&m2_inv(&prec));
    let (z0, z1): (f64, f64) = (rng.sample(nrm), rng.sample(nrm));
    let px = mx + l[0]*z0;
    let py = my + l[2]*z0 + l[3]*z1;
    let le_p = ctx.emit(koff + t, px, py);
    let mut corr = le_p - le[t];
    if ctx.spherical {
        let sig2 = ctx.km_lon2 / pm[0];
        let mut e = mv_edge(ctx.km_lon2, ctx.km_lat2, xlon[t-1], xlat[t-1], px, py)
                  - mv_edge(ctx.km_lon2, ctx.km_lat2, xlon[t-1], xlat[t-1], xlon[t], xlat[t]);
        if t < k - 1 {
            e += mv_edge(ctx.km_lon2, ctx.km_lat2, px, py, xlon[t+1], xlat[t+1])
               - mv_edge(ctx.km_lon2, ctx.km_lat2, xlon[t], xlat[t], xlon[t+1], xlat[t+1]);
        }
        corr += -e / (2.0 * sig2);
    }
    if ctx.rho.get() != 0.0 {
        let sig2 = ctx.km_lon2 / pm[0];
        let cur = |m: usize| (xlon[m], xlat[m]);
        let prp = |m: usize| if m == t { (px, py) } else { (xlon[m], xlat[m]) };
        let e = crw_sum(ctx, prp, k, t, t + 2) - crw_sum(ctx, cur, k, t, t + 2);
        corr += -e / (2.0 * sig2);
    }
    if rng.gen::<f64>().ln() < corr { xlon[t] = px; xlat[t] = py; le[t] = le_p; }
}

fn rb_polish(ctx: &Ctx, xlon: &mut [f64], xlat: &mut [f64], le: &mut [f64], koff: usize,
             k: usize, pm: &[f64; 4], rng: &mut StdRng, nrm: &Normal<f64>) {
    let l = m2_chol_lower(&m2_inv(&[2.0*pm[0], 2.0*pm[1], 2.0*pm[2], 2.0*pm[3]]));
    for parity in 0..2usize {
        for t in 1..(k - 1) {
            if t % 2 != parity { continue; }
            let mx = (xlon[t-1]+xlon[t+1])*0.5;
            let my = (xlat[t-1]+xlat[t+1])*0.5;
            let (z0, z1): (f64, f64) = (rng.sample(nrm), rng.sample(nrm));
            let px = mx + l[0]*z0;
            let py = my + l[2]*z0 + l[3]*z1;
            let le_p = ctx.emit(koff + t, px, py);
            let mut corr = le_p - le[t];
            if ctx.spherical {
                let sig2 = ctx.km_lon2 / pm[0];
                let e = mv_edge(ctx.km_lon2, ctx.km_lat2, xlon[t-1], xlat[t-1], px, py)
                      - mv_edge(ctx.km_lon2, ctx.km_lat2, xlon[t-1], xlat[t-1], xlon[t], xlat[t])
                      + mv_edge(ctx.km_lon2, ctx.km_lat2, px, py, xlon[t+1], xlat[t+1])
                      - mv_edge(ctx.km_lon2, ctx.km_lat2, xlon[t], xlat[t], xlon[t+1], xlat[t+1]);
                corr += -e / (2.0 * sig2);
            }
            if ctx.rho.get() != 0.0 {
                let sig2 = ctx.km_lon2 / pm[0];
                let cur = |m: usize| (xlon[m], xlat[m]);
                let prp = |m: usize| if m == t { (px, py) } else { (xlon[m], xlat[m]) };
                let e = crw_sum(ctx, prp, k, t, t + 2) - crw_sum(ctx, cur, k, t, t + 2);
                corr += -e / (2.0 * sig2);
            }
            if rng.gen::<f64>().ln() < corr { xlon[t] = px; xlat[t] = py; le[t] = le_p; }
        }
    }
}

#[allow(clippy::too_many_arguments)]
#[allow(clippy::too_many_arguments)]
fn block_update(ctx: &Ctx, xlon: &mut [f64], xlat: &mut [f64], le: &mut [f64], koff: usize,
                l: usize, r: usize, k: usize, surr_mu: &[f64], surr_p: &[f64], pm: &[f64; 4],
                rng: &mut StdRng, nrm: &Normal<f64>) -> bool {
    let b_len = r - l + 1;
    let d = 2 * b_len;
    let mut q = vec![0.0f64; d*d];
    let mut bb = vec![0.0f64; d];
    for a in 0..b_len {
        let gk = koff + l + a;
        let (i0, i1) = (2*a, 2*a + 1);
        let pk = [surr_p[4*gk], surr_p[4*gk+1], surr_p[4*gk+2], surr_p[4*gk+3]];
        q[i0*d + i0] += 2.0*pm[0] + pk[0]; q[i0*d + i1] += 2.0*pm[1] + pk[1];
        q[i1*d + i0] += 2.0*pm[2] + pk[2]; q[i1*d + i1] += 2.0*pm[3] + pk[3];
        if a > 0 {
            let (j0, j1) = (2*(a-1), 2*(a-1) + 1);
            q[i0*d + j0] -= pm[0]; q[i0*d + j1] -= pm[1];
            q[i1*d + j0] -= pm[2]; q[i1*d + j1] -= pm[3];
            q[j0*d + i0] -= pm[0]; q[j0*d + i1] -= pm[2];
            q[j1*d + i0] -= pm[1]; q[j1*d + i1] -= pm[3];
        }
        let (bx, by) = m2_mulv(&pk, surr_mu[2*gk], surr_mu[2*gk + 1]);
        bb[i0] += bx; bb[i1] += by;
    }
    let (elx, ely) = m2_mulv(pm, xlon[l-1], xlat[l-1]); bb[0] += elx; bb[1] += ely;
    let (erx, ery) = m2_mulv(pm, xlon[r+1], xlat[r+1]); bb[d-2] += erx; bb[d-1] += ery;
    let rmat = chol_upper(&q, d);
    let mu_blk = solve_ut(&rmat, &solve_lt(&rmat, &bb, d), d);
    let z: Vec<f64> = (0..d).map(|_| rng.sample(nrm)).collect();
    let xdraw = solve_ut(&rmat, &z, d);
    let prop: Vec<f64> = (0..d).map(|i| mu_blk[i] + xdraw[i]).collect();
    let mut corr = 0.0;
    let mut lep = vec![0.0f64; b_len];
    for a in 0..b_len {
        let gk = koff + l + a;
        let (px, py) = (prop[2*a], prop[2*a + 1]);
        let (cx, cy) = (xlon[l + a], xlat[l + a]);
        let pk = [surr_p[4*gk], surr_p[4*gk+1], surr_p[4*gk+2], surr_p[4*gk+3]];
        let (mkx, mky) = (surr_mu[2*gk], surr_mu[2*gk + 1]);
        lep[a] = ctx.emit(gk, px, py);
        corr += (lep[a] - (-0.5*quad2(&pk, px-mkx, py-mky))) - (le[l+a] - (-0.5*quad2(&pk, cx-mkx, cy-mky)));
    }
    let posp = |m: usize| -> (f64, f64) {
        if m >= l && m <= r { (prop[2*(m-l)], prop[2*(m-l)+1]) } else { (xlon[m], xlat[m]) }
    };
    if ctx.spherical {
        let sig2 = ctx.km_lon2 / pm[0];
        let mut e = 0.0;
        for m in (l-1)..=r {   // affected edges (m, m+1), m = l-1..r
            let (pa, pb) = (posp(m), posp(m + 1));
            e += mv_edge(ctx.km_lon2, ctx.km_lat2, pa.0, pa.1, pb.0, pb.1)
               - mv_edge(ctx.km_lon2, ctx.km_lat2, xlon[m], xlat[m], xlon[m+1], xlat[m+1]);
        }
        corr += -e / (2.0 * sig2);
    }
    if ctx.rho.get() != 0.0 {
        let sig2 = ctx.km_lon2 / pm[0];
        let cur = |m: usize| (xlon[m], xlat[m]);
        let e = crw_sum(ctx, &posp, k, l, r + 2) - crw_sum(ctx, cur, k, l, r + 2);
        corr += -e / (2.0 * sig2);
    }
    if rng.gen::<f64>().ln() < corr {
        for a in 0..b_len { xlon[l + a] = prop[2*a]; xlat[l + a] = prop[2*a + 1]; le[l + a] = lep[a]; }
        true
    } else { false }
}

// One full sweep of a single track: interior blocks (with single-site fallback on
// obs-free knots), single-site free endpoint, optional red-black polish.
// free_end: update the last knot (single-site) when retrieval is unknown; when a
// retrieval location is given the last knot is fixed and skipped.
#[allow(clippy::too_many_arguments)]
fn track_update(ctx: &Ctx, xlon: &mut [f64], xlat: &mut [f64], le: &mut [f64], koff: usize,
                k: usize, surr_mu: &[f64], surr_p: &[f64], pm: &[f64; 4], block_len: usize,
                polish: bool, free_end: bool, rng: &mut StdRng, nrm: &Normal<f64>) -> (u64, u64) {
    let phase = (rng.gen::<u64>() % block_len as u64) as usize;
    let mut blk_max = if block_len > phase { block_len - phase } else { 1 };
    let (mut acc, mut att) = (0u64, 0u64);
    let mut l = 1usize;
    while l <= k - 2 {
        if ctx.obs_len[koff + l] <= 0 {
            single_site(ctx, xlon, xlat, le, koff, l, k, pm, rng, nrm); l += 1; continue;
        }
        let mut r = l;
        while r < (l + blk_max - 1).min(k - 2) && ctx.obs_len[koff + r + 1] > 0 { r += 1; }
        blk_max = block_len;
        att += 1;
        if block_update(ctx, xlon, xlat, le, koff, l, r, k, surr_mu, surr_p, pm, rng, nrm) { acc += 1; }
        l = r + 1;
    }
    if free_end { single_site(ctx, xlon, xlat, le, koff, k - 1, k, pm, rng, nrm); }
    if polish { rb_polish(ctx, xlon, xlat, le, koff, k, pm, rng, nrm); }
    (acc, att)
}

/// Run the native single-track block + polish sampler.
///
/// @param knot_obs_start 0-based start index of each knot's contiguous obs block
/// @param knot_obs_len number of observations in each knot
/// @param obs_times observation timestamps (seconds since 1970)
/// @param obs_light observed light values (processed as by the R fit)
/// @param calibration Response parameters: `c(intercept, slope)` for the clamped-linear
///   model, or `c(floor, amp, z50, scale)` for the logistic model (see details in
///   `fit_light_response`). Length selects the model.
/// @param likelihood_params c(lambda, max_light, prob_slab)
/// @param surr_mu per-knot surrogate mean, length 2K (lon, lat interleaved)
/// @param surr_p per-knot surrogate 2x2 precision, length 4K (row-major)
/// @param p_move movement precision (2x2, row-major)
/// @param start_lon Fixed first-knot longitude
/// @param start_lat Fixed first-knot latitude
/// @param block_len Maximum block length (knots)
/// @param sweeps Total sweeps
/// @param burn Burn-in sweeps
/// @param thin Thinning interval
/// @param polish Whether to run the red-black single-site polish each sweep
/// @param seed RNG seed; 0 means entropy
/// @return List of per-knot posterior mean_lon, sd_lon, mean_lat, sd_lat, plus accept and n_kept
/// @name run_block_track
#[extendr]
fn run_block_track(
    knot_obs_start: &[i32], knot_obs_len: &[i32],
    obs_times: &[f64], obs_light: &[f64],
    calibration: Vec<f64>, likelihood_params: Vec<f64>,
    surr_mu: Vec<f64>, surr_p: Vec<f64>, p_move: Vec<f64>,
    start_lon: f64, start_lat: f64,
    block_len: i32, sweeps: i32, burn: i32, thin: i32, polish: bool, seed: f64,
) -> List {
    let k = knot_obs_start.len();
    let bl = block_len as usize;
    let pm = [p_move[0], p_move[1], p_move[2], p_move[3]];
    let eph: Vec<(f64, f64, f64, f64)> = obs_times.iter().map(|&t| solar_ephemeris(t)).collect();
    let ctx = Ctx {
        obs_start: knot_obs_start, obs_len: knot_obs_len, obs_light: &obs_light, eph: &eph,
        cal: cal_array(&calibration), cal_len: calibration.len().min(4),
        lambda: likelihood_params[0], max_light: likelihood_params[1], prob_slab: likelihood_params[2],
        koff: 0, aux: &[], aux_lon0: 0.0, aux_dlon: 1.0, aux_ncol: 0, aux_lat0: 0.0, aux_dlat: 1.0, aux_nrow: 0,
        spherical: false, km_lon2: 0.0, km_lat2: 0.0, rho: std::cell::Cell::new(0.0),
    };
    let informative: Vec<bool> = knot_obs_len.iter().map(|&n| n > 0).collect();

    let mut rng: StdRng = if seed == 0.0 { StdRng::from_entropy() } else { StdRng::seed_from_u64(seed as u64) };
    let nrm = Normal::new(0.0, 1.0).unwrap();

    // init at surrogate mean (start-fixed first knot)
    let mut xlon = vec![0.0; k];
    let mut xlat = vec![0.0; k];
    for i in 0..k {
        if informative[i] { xlon[i] = surr_mu[2*i]; xlat[i] = surr_mu[2*i + 1]; }
        else { xlon[i] = start_lon; xlat[i] = start_lat; }
    }
    xlon[0] = start_lon; xlat[0] = start_lat;
    let mut le: Vec<f64> = (0..k).map(|i| ctx.emit(i, xlon[i], xlat[i])).collect();

    let (mut mlon, mut m2lon) = (vec![0.0; k], vec![0.0; k]);
    let (mut mlat, mut m2lat) = (vec![0.0; k], vec![0.0; k]);
    let mut n_kept = 0usize;
    let (mut acc, mut att) = (0u64, 0u64);

    let _ = &informative;  // (informative == obs_len>0; track_update reads ctx.obs_len directly)
    for sweep in 0..sweeps {
        let (a, t) = track_update(&ctx, &mut xlon, &mut xlat, &mut le, 0, k,
                                  &surr_mu, &surr_p, &pm, bl, polish, true, &mut rng, &nrm);
        acc += a; att += t;
        if sweep >= burn && (sweep - burn) % thin == 0 {
            n_kept += 1;
            let inv = 1.0 / n_kept as f64;
            for i in 0..k {
                let dl = xlon[i] - mlon[i]; mlon[i] += dl*inv; m2lon[i] += dl*(xlon[i] - mlon[i]);
                let da = xlat[i] - mlat[i]; mlat[i] += da*inv; m2lat[i] += da*(xlat[i] - mlat[i]);
            }
        }
    }
    let den = if n_kept > 1 { (n_kept - 1) as f64 } else { 1.0 };
    let slon: Vec<f64> = m2lon.iter().map(|v| (v/den).sqrt()).collect();
    let slat: Vec<f64> = m2lat.iter().map(|v| (v/den).sqrt()).collect();
    list!(mean_lon = mlon, sd_lon = slon, mean_lat = mlat, sd_lat = slat,
          accept = if att > 0 { acc as f64 / att as f64 } else { 0.0 }, n_kept = n_kept as i32)
}

/// Native hierarchical partial-pooling sampler over a panel of tracks.
///
/// Runs the whole Gibbs sweep natively: each individual's track is updated with
/// the block+polish kernel conditional on its movement variance sig2_i, then
/// sig2_i and the population scale beta are drawn from their conjugate
/// full-conditionals. All arrays concatenate the individuals; knots_per_ind gives
/// each individual's knot count, and knot_obs_start indexes the GLOBAL obs arrays.
///
/// @param n_ind Number of individuals
/// @param knots_per_ind Knot count per individual (length n_ind)
/// @param knot_obs_start Per-knot global obs start index (length sum knots_per_ind)
/// @param knot_obs_len Per-knot obs count
/// @param obs_times Global concatenated observation timestamps
/// @param obs_light Global concatenated observed light
/// @param cal Per-individual response parameters, packed contiguously: length
///   2*n_ind for the clamped-linear model or 4*n_ind for the logistic one.
/// @param lp Per-individual c(lambda, max_light, prob_slab), length 3*n_ind
/// @param surr_mu Per-knot surrogate mean (2 per knot)
/// @param surr_p Per-knot surrogate 2x2 precision (4 per knot)
/// @param cinv Per-individual movement shape precision Cinv (2x2), length 4*n_ind
/// @param start_lon Per-individual fixed first-knot (deploy) longitude
/// @param start_lat Per-individual fixed first-knot (deploy) latitude
/// @param end_lon Per-individual fixed last-knot (retrieval) longitude; NaN = free
/// @param end_lat Per-individual fixed last-knot (retrieval) latitude; NaN = free
/// @param aux_flat Concatenated per-tag additive location log-likelihood (sensor
///   terms), each tag a K*ncell block on its own regular lon/lat grid; empty for none
/// @param aux_ncol,aux_nrow Per-tag aux grid dimensions (0 columns/rows = no terms)
/// @param aux_lon0,aux_dlon,aux_lat0,aux_dlat Per-tag aux grid origin and spacing
/// @param a_pop InvGamma shape for sig2_i
/// @param g0 Gamma shape hyperprior for beta
/// @param h0 Gamma rate hyperprior for beta
/// @param block_len Maximum block length
/// @param sweeps Total sweeps
/// @param burn Burn-in sweeps
/// @param thin Thinning interval
/// @param polish Whether to run the red-black polish each sweep
/// @param spherical Great-circle movement metric (exp(-gcdist^2/2 sig2)) if TRUE,
///   else the flat tangent-plane metric at each tag's reference latitude
/// @param crw Correlated random walk: model the increments as an AR(1),
///   delta_t = rho*delta_{t-1} + eps, with a per-tag persistence rho sampled
///   alongside the movement variance. FALSE gives the memoryless Brownian walk
///   and is bit-identical to the previous behaviour.
/// @param rho_prior_sd Standard deviation of the mean-zero normal prior on each
///   tag's persistence; rho is confined to (-0.99, 0.99)
/// @param seed RNG seed; 0 means entropy
/// @return List with beta and sig2 (kept draws), plus per-knot track posterior
///   mean_lon/sd_lon/mean_lat/sd_lat (concatenated across individuals, same order
///   as knots_per_ind), and rho (per-tag persistence draws, zero when crw is FALSE)
/// @name run_block_hier
#[extendr]
fn run_block_hier(
    n_ind: i32, knots_per_ind: &[i32],
    knot_obs_start: &[i32], knot_obs_len: &[i32],
    obs_times: &[f64], obs_light: &[f64],
    cal: Vec<f64>, lp: Vec<f64>,
    surr_mu: Vec<f64>, surr_p: Vec<f64>, cinv: Vec<f64>,
    start_lon: Vec<f64>, start_lat: Vec<f64>,
    end_lon: Vec<f64>, end_lat: Vec<f64>,   // retrieval; NaN => last knot free
    aux_flat: Vec<f64>, aux_ncol: &[i32], aux_nrow: &[i32],   // sensor terms (0 cells => none)
    aux_lon0: Vec<f64>, aux_dlon: Vec<f64>, aux_lat0: Vec<f64>, aux_dlat: Vec<f64>,
    a_pop: f64, g0: f64, h0: f64,
    block_len: i32, sweeps: i32, burn: i32, thin: i32, polish: bool, spherical: bool,
    crw: bool, rho_prior_sd: f64, seed: f64,
) -> List {
    let n = n_ind as usize;
    let bl = block_len as usize;
    let ki: Vec<usize> = knots_per_ind.iter().map(|&x| x as usize).collect();
    // knot offset per individual
    let mut koff = vec![0usize; n];
    for i in 1..n { koff[i] = koff[i-1] + ki[i-1]; }
    let eph: Vec<(f64, f64, f64, f64)> = obs_times.iter().map(|&t| solar_ephemeris(t)).collect();

    let mut rng: StdRng = if seed == 0.0 { StdRng::from_entropy() } else { StdRng::seed_from_u64(seed as u64) };
    let nrm = Normal::new(0.0, 1.0).unwrap();

    // per-individual state: track (local), emit cache, Cinv
    let mut xlon: Vec<Vec<f64>> = Vec::with_capacity(n);
    let mut xlat: Vec<Vec<f64>> = Vec::with_capacity(n);
    let mut le:   Vec<Vec<f64>> = Vec::with_capacity(n);
    let mut sig2 = vec![0.0f64; n];
    let cinv_i: Vec<[f64;4]> = (0..n).map(|i| [cinv[4*i], cinv[4*i+1], cinv[4*i+2], cinv[4*i+3]]).collect();

    // aux (sensor terms) offsets: aux_off[i] = sum_{j<i} K_j * ncell_j
    let ancol: Vec<usize> = aux_ncol.iter().map(|&x| x as usize).collect();
    let anrow: Vec<usize> = aux_nrow.iter().map(|&x| x as usize).collect();
    let mut aux_off = vec![0usize; n];
    for i in 1..n { aux_off[i] = aux_off[i-1] + ki[i-1] * ancol[i-1] * anrow[i-1]; }
    // Per-individual calibration is packed contiguously; the stride says which
    // response model it is (2 = clamped linear, 4 = logistic), so no extra
    // argument is needed and existing length-2*n callers are unaffected.
    let cal_stride = if n > 0 { cal.len() / n } else { 2 };
    // build a Ctx per individual (borrows global obs arrays + per-individual cal/lp/aux)
    let ctxs: Vec<Ctx> = (0..n).map(|i| {
        let ncell = ki[i] * ancol[i] * anrow[i];
        let aslice: &[f64] = if ncell > 0 { &aux_flat[aux_off[i]..aux_off[i] + ncell] } else { &[] };
        Ctx {
            obs_start: knot_obs_start, obs_len: knot_obs_len, obs_light: &obs_light, eph: &eph,
            cal: cal_array(&cal[cal_stride*i .. cal_stride*(i+1)]), cal_len: cal_stride,
            lambda: lp[3*i], max_light: lp[3*i+1], prob_slab: lp[3*i+2],
            koff: koff[i], aux: aslice,
            aux_lon0: aux_lon0[i], aux_dlon: aux_dlon[i], aux_ncol: ancol[i],
            aux_lat0: aux_lat0[i], aux_dlat: aux_dlat[i], aux_nrow: anrow[i],
            spherical, km_lon2: cinv_i[i][0], km_lat2: cinv_i[i][3], rho: std::cell::Cell::new(0.0),
        }
    }).collect();

    let free_end: Vec<bool> = (0..n).map(|i| !end_lon[i].is_finite()).collect();
    for i in 0..n {
        let k = ki[i];
        let mut xl = vec![0.0; k]; let mut xt = vec![0.0; k];
        for t in 0..k {
            let gk = koff[i] + t;
            if knot_obs_len[gk] > 0 { xl[t] = surr_mu[2*gk]; xt[t] = surr_mu[2*gk + 1]; }
            else { xl[t] = start_lon[i]; xt[t] = start_lat[i]; }
        }
        xl[0] = start_lon[i]; xt[0] = start_lat[i];
        if !free_end[i] { xl[k-1] = end_lon[i]; xt[k-1] = end_lat[i]; }
        let lei: Vec<f64> = (0..k).map(|t| ctxs[i].emit(koff[i] + t, xl[t], xt[t])).collect();
        // data-driven sig2 start from initial increments (great-circle^2 if spherical)
        let mut ss = 0.0;
        for t in 1..k {
            ss += if spherical { let d = gcdist_km(xl[t-1], xt[t-1], xl[t], xt[t]); d*d }
                  else { quad2(&cinv_i[i], xl[t]-xl[t-1], xt[t]-xt[t-1]) };
        }
        sig2[i] = (ss / (k as f64 - 1.0)).max(1e-6);
        xlon.push(xl); xlat.push(xt); le.push(lei);
    }
    let mut beta = { let s: f64 = sig2.iter().sum(); (s / n as f64) * (a_pop - 1.0) };
    // Per-tag directional persistence. Zero unless the correlated random walk is
    // requested, in which case it is sampled each sweep from its conjugate
    // normal full conditional.
    let mut rho = vec![0.0f64; n];

    let nkeep = ((sweeps - burn + thin - 1) / thin).max(0) as usize;
    let mut beta_draws: Vec<f64> = Vec::with_capacity(nkeep);
    let mut sig2_draws: Vec<f64> = Vec::with_capacity(nkeep * n);
    let mut rho_draws: Vec<f64> = Vec::with_capacity(nkeep * n);
    let total_knots: usize = ki.iter().sum();
    let (mut tmlon, mut tm2lon) = (vec![0.0; total_knots], vec![0.0; total_knots]);
    let (mut tmlat, mut tm2lat) = (vec![0.0; total_knots], vec![0.0; total_knots]);
    let mut rec = 0usize;

    for sweep in 0..sweeps {
        for i in 0..n {
            let k = ki[i];
            let pm = [cinv_i[i][0]/sig2[i], cinv_i[i][1]/sig2[i], cinv_i[i][2]/sig2[i], cinv_i[i][3]/sig2[i]];
            track_update(&ctxs[i], &mut xlon[i], &mut xlat[i], &mut le[i], koff[i], k,
                         &surr_mu, &surr_p, &pm, bl, polish, free_end[i], &mut rng, &nrm);
            // conjugate sig2_i | track ~ InvGamma(a_pop + (k-1), beta + 0.5 SS)
            // SS is sum of squared movement (great-circle when spherical, planar otherwise)
            // Under the correlated random walk the residual is the INNOVATION
            // delta_t - rho*delta_{t-1}, not the step itself; the first
            // increment has no predecessor and enters as-is.
            let mut ss = 0.0;
            for t in 1..k {
                let (dx, dy) = (xlon[i][t]-xlon[i][t-1], xlat[i][t]-xlat[i][t-1]);
                ss += if crw && t >= 2 {
                    let (px, py) = (xlon[i][t-1]-xlon[i][t-2], xlat[i][t-1]-xlat[i][t-2]);
                    quad2(&cinv_i[i], dx - rho[i]*px, dy - rho[i]*py)
                } else if spherical {
                    let d = gcdist_km(xlon[i][t-1], xlat[i][t-1], xlon[i][t], xlat[i][t]); d*d
                } else {
                    quad2(&cinv_i[i], dx, dy)
                };
            }
            let g = Gamma::new(a_pop + (k as f64 - 1.0), 1.0 / (beta + 0.5*ss)).unwrap();
            sig2[i] = 1.0 / rng.sample(g);

            // rho_i | track, sig2_i is the coefficient of a linear regression of
            // each increment on the one before it, so it has a conjugate normal
            // full conditional under a mean-zero normal prior. Confined to
            // (-0.99, 0.99), beyond which the walk is not stationary.
            if crw && k >= 3 {
                let (mut num, mut den) = (0.0, 0.0);
                for t in 2..k {
                    let (dx, dy) = (xlon[i][t]-xlon[i][t-1], xlat[i][t]-xlat[i][t-1]);
                    let (px, py) = (xlon[i][t-1]-xlon[i][t-2], xlat[i][t-1]-xlat[i][t-2]);
                    num += cinv_i[i][0]*dx*px + cinv_i[i][3]*dy*py;
                    den += cinv_i[i][0]*px*px + cinv_i[i][3]*py*py;
                }
                let prior_prec = if rho_prior_sd > 0.0 { 1.0/(rho_prior_sd*rho_prior_sd) } else { 0.0 };
                let prec = den / sig2[i] + prior_prec;
                if prec > 0.0 && den.is_finite() {
                    let mean = (num / sig2[i]) / prec;
                    let sd = (1.0 / prec).sqrt();
                    let mut r = mean + sd * rng.sample(&nrm);
                    if !r.is_finite() { r = 0.0; }
                    rho[i] = r.max(-0.99).min(0.99);
                    ctxs[i].rho.set(rho[i]);
                }
            }
        }
        // conjugate beta | . ~ Gamma(g0 + n*a_pop, h0 + sum 1/sig2)
        let inv_sum: f64 = sig2.iter().map(|s| 1.0/s).sum();
        let gb = Gamma::new(g0 + n as f64 * a_pop, 1.0 / (h0 + inv_sum)).unwrap();
        beta = rng.sample(gb);

        if sweep >= burn && (sweep - burn) % thin == 0 {
            beta_draws.push(beta);
            for i in 0..n { sig2_draws.push(sig2[i]); }
            for i in 0..n { rho_draws.push(rho[i]); }
            rec += 1; let inv = 1.0 / rec as f64;
            for i in 0..n { for t in 0..ki[i] {
                let gk = koff[i] + t;
                let dl = xlon[i][t] - tmlon[gk]; tmlon[gk] += dl*inv; tm2lon[gk] += dl*(xlon[i][t] - tmlon[gk]);
                let da = xlat[i][t] - tmlat[gk]; tmlat[gk] += da*inv; tm2lat[gk] += da*(xlat[i][t] - tmlat[gk]);
            } }
        }
    }
    let nk = beta_draws.len() as i32;
    let den = if rec > 1 { (rec - 1) as f64 } else { 1.0 };
    let sdlon: Vec<f64> = tm2lon.iter().map(|v| (v/den).sqrt()).collect();
    let sdlat: Vec<f64> = tm2lat.iter().map(|v| (v/den).sqrt()).collect();
    list!(beta = beta_draws, sig2 = sig2_draws, rho = rho_draws, n_kept = nk, n_ind = n_ind,
          mean_lon = tmlon, sd_lon = sdlon, mean_lat = tmlat, sd_lat = sdlat)
}

extendr_module! {
    mod invtwilightfree;
    fn solar_zenith;
    fn light_log_likelihood;
    fn run_particle_filter;
    fn eval_logpk_grid;
    fn run_grid_hmm;
    fn run_block_track;
    fn run_block_hier;
}
