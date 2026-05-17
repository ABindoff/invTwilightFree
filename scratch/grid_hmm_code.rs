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
    let prob_slab = likelihood_params[2];
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
                let log_norm_const = - (sigma).ln(); 
                
                let mut max_val = -1e30;
                let mut sum_exp = 0.0;
                
                for j in 0..n {
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
                    
                    let a = ((lat_rad[j] - lat_rad[i]) / 2.0).sin().powi(2) +
                            lat_rad[i].cos() * lat_rad[j].cos() * ((lon_rad[j] - lon_rad[i]) / 2.0).sin().powi(2);
                    let dist = 2.0 * earth_radius * a.sqrt().asin();
                    
                    if dist > threshold_dist { continue; }
                    
                    for s_next in 0..num_states {
                        let beta_next = beta[k+1][j * num_states + s_next];
                        if beta_next > -1e29 {
                            let sigma = diffusion[s_next] * dt.sqrt();
                            let var2 = 2.0 * sigma * sigma;
                            let log_norm_const = - (sigma).ln();
                            
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
                let mut sum_g_cell = -1e30;
                for s in 0..num_states {
                    if gamma[i * num_states + s] > -1e29 {
                        let w = (gamma[i * num_states + s] - log_denom).exp();
                        prob_states[s] += w;
                        
                        let g_state = gamma[i * num_states + s];
                        if g_state > sum_g_cell {
                            sum_g_cell = sum_g_cell * (g_state - sum_g_cell).exp() + g_state; // Approx
                        } else {
                            sum_g_cell = sum_g_cell * (g_state - sum_g_cell).exp() + sum_g_cell; 
                        }
                    }
                }
                if sum_g_cell > best_g {
                    best_g = sum_g_cell;
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
        prob_state = prob_state_list
    )
}
