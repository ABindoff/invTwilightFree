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
        
        // Probability under the spike (true light, subject to shading)
        let spike_density = if obs <= exp {
            // Shading reduces light. Exponential decay from expected max light.
            lambda * (-lambda * (exp - obs)).exp()
        } else {
            // If observed is brighter than expected, we penalize it.
            // A sharp decay prevents numerical underflow to 0 while enforcing physics.
            lambda * (-lambda * 10.0 * (obs - exp)).exp() 
        };
        
        // Mixture model: (1 - pi) * True Light + pi * False Light
        let marginal_density = (1.0 - prob_slab) * spike_density + prob_slab * slab_density;
        
        // Accumulate log likelihood
        log_lik += marginal_density.ln();
    }
    
    log_lik
}

use rand::prelude::*;
use rand::rngs::StdRng;
use rand_distr::{Normal, Beta};

#[derive(Clone, Copy, Debug)]
struct Particle {
    lat: f64,
    lon: f64,
    weight: f64,
    state: usize,
    prob_slab: f64,
}

fn get_solar_zenith(t: f64, l: f64, phi_deg: f64) -> f64 {
    let phi = phi_deg.to_radians();
    let d = (t - 946728000.0) / 86400.0;
    let g = (357.529 + 0.98560028 * d).to_radians();
    let q = 280.459 + 0.98564736 * d;
    let l_sun = (q + 1.915 * g.sin() + 0.020 * (2.0 * g).sin()).to_radians();
    let e = (23.439 - 0.00000036 * d).to_radians();
    let delta = (e.sin() * l_sun.sin()).asin();
    let ra = f64::atan2(e.cos() * l_sun.sin(), l_sun.cos());
    let gmst = (18.697374558 + 24.06570982441908 * d) * 15.0_f64.to_radians();
    let lmst = gmst + l.to_radians();
    let h = lmst - ra;
    let cos_z = phi.sin() * delta.sin() + phi.cos() * delta.cos() * h.cos();
    cos_z.acos().to_degrees()
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
/// @param calibration c(intercept, slope)
/// @param likelihood_params c(lambda, max_light, prob_slab)
/// @param mask_matrix Flattened spatial mask matrix (0 = impassable); empty for no mask
/// @param mask_extent c(xmin, xmax, ymin, ymax) extent of the mask raster
/// @param mask_nrow Number of rows in the mask raster
/// @param mask_ncol Number of columns in the mask raster
/// @param seed Integer seed for reproducibility; 0 means non-deterministic (uses entropy)
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
) -> List {
    let n = n_particles as usize;
    let num_obs = unix_times.len();
    let intercept = calibration[0];
    let slope = calibration[1];
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
        if obs_in_range && (max_obs - min_obs) < 10.0 {
            diff_scale = 0.1; // 10x tighter diffusion
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

            let mut sigma = diffusion[new_state] * dt.sqrt() * diff_scale;
            if is_guided {
                sigma = diffusion[new_state] * ((dt * (time_remain - dt)) / time_remain).max(0.0).sqrt() * diff_scale;
            }

            if is_guided {
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
                
                let zenith = get_solar_zenith(unix_times[j], p_lon, p_lat);
                let expected = (intercept - slope * zenith).max(0.0).min(max_light);
                let obs = obs_light[j];
                let spike_density = if obs <= expected {
                    lambda * (-lambda * (expected - obs)).exp()
                } else {
                    lambda * (-lambda * 2.0 * (obs - expected)).exp()
                };
                let current_prob_slab = particles[i].prob_slab;
                let den = (1.0 - current_prob_slab) * spike_density + current_prob_slab * slab_density;
                log_lik += den.ln();
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
            
            let z = get_solar_zenith(unix_times[j], p_lon, p_lat);
            mean_z += z * w;
            
            let exp = (intercept - slope * z).max(0.0).min(max_light);
            let obs = obs_light[j];
            let spike_density = if obs <= exp {
                lambda * (-lambda * (exp - obs)).exp()
            } else {
                lambda * (-lambda * 2.0 * (obs - exp)).exp()
            };
            let den = (1.0 - current_prob_slab) * spike_density + current_prob_slab * slab_density;
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
/// @param calibration c(intercept, slope)
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
    let intercept = calibration[0];
    let slope = calibration[1];
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
    
    for i in 0..n {
        let mut sum_logl = 0.0;
        for j in 0..num_obs {
            let zenith = get_solar_zenith(unix_times[j], lon[i], lat[i]);
            let expected = (intercept - slope * zenith).max(0.0).min(max_light);
            let obs = obs_light[j];
            let spike = if obs <= expected {
                lambda * (-lambda * (expected - obs)).exp()
            } else {
                lambda * (-lambda * 2.0 * (obs - expected)).exp()
            };
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
) -> List {
    let n = lon.len();
    let num_states = diffusion.len();
    let k_steps = knot_times.len();
    
    let intercept = calibration[0];
    let slope = calibration[1];
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
            let mut sum_logl = 0.0;
            for &j in &obs_in_k {
                let zenith = get_solar_zenith(obs_times[j], lon[i], lat[i]);
                let expected = (intercept - slope * zenith).max(0.0).min(max_light);
                let obs = obs_light[j];
                let spike = if obs <= expected {
                    lambda * (-lambda * (expected - obs)).exp()
                } else {
                    lambda * (-lambda * 2.0 * (obs - expected)).exp()
                };
                let den = (1.0 - prob_slab) * spike + prob_slab * slab_density;
                sum_logl += den.ln();
            }
            logpk[k][i] = sum_logl;
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
                        
                        let g_state = gamma[i * num_states + s];
                        if g_state > best_g {
                            best_g = g_state;
                            best_i = i;
                        }
                    }
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
        prob_state = prob_state_list
    )
}

extendr_module! {
    mod invtwilightfree;
    fn solar_zenith;
    fn light_log_likelihood;
    fn run_particle_filter;
    fn eval_logpk_grid;
    fn run_grid_hmm;
}
