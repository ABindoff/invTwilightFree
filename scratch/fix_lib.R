lines_lib <- readLines("src/rust/src/lib.rs")

# Find the start of the corrupted list! in run_particle_filter
idx_list <- which(grepl("obs_zenith = obs_zenith", lines_lib))
# Wait, obs_zenith was deleted. Let's find `obs_times = unix_times`
idx_list <- which(grepl("obs_times = unix_times,", lines_lib))

if (length(idx_list) > 0) {
  idx <- idx_list[1]
  # Keep up to here
  lines_lib <- lines_lib[1:idx]
  # Re-add obs_zenith and prob_false, then close run_particle_filter
  lines_lib <- c(lines_lib, 
                 "        obs_zenith = obs_zenith,",
                 "        prob_false = obs_prob_false",
                 "    )",
                 "}",
                 "")
}

# Re-add eval_logpk_grid
eval_logpk_grid_code <- c(
"/// @name eval_logpk_grid",
"/// @export",
"#[extendr]",
"fn eval_logpk_grid(",
"    lon: &[f64],",
"    lat: &[f64],",
"    unix_times: &[f64],",
"    obs_light: &[f64],",
"    calibration: Vec<f64>,",
"    likelihood_params: Vec<f64>",
") -> Vec<f64> {",
"    let n = lon.len();",
"    let num_obs = unix_times.len();",
"    let intercept = calibration[0];",
"    let slope = calibration[1];",
"    let lambda = likelihood_params[0];",
"    let max_light = likelihood_params[1];",
"    let prob_slab = likelihood_params[2];",
"    let slab_density = 1.0 / max_light;",
"",
"    let mut logl = vec![0.0; n];",
"    ",
"    for i in 0..n {",
"        let mut sum_logl = 0.0;",
"        for j in 0..num_obs {",
"            let zenith = get_solar_zenith(unix_times[j], lon[i], lat[i]);",
"            let expected = (intercept - slope * zenith).max(0.0).min(max_light);",
"            let obs = obs_light[j];",
"            let spike = if obs <= expected {",
"                lambda * (-lambda * (expected - obs)).exp()",
"            } else {",
"                lambda * (-lambda * 2.0 * (obs - expected)).exp()",
"            };",
"            let den = (1.0 - prob_slab) * spike + prob_slab * slab_density;",
"            sum_logl += den.ln();",
"        }",
"        logl[i] = sum_logl;",
"    }",
"    ",
"    logl",
"}",
""
)

grid_hmm_code <- readLines("scratch/grid_hmm_code.rs")

# Write everything out
new_lib <- c(lines_lib, eval_logpk_grid_code, grid_hmm_code)
writeLines(new_lib, "src/rust/src/lib.rs")

cat("Fixed lib.rs!\n")
