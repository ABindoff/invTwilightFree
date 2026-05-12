# invTwilightFree: High-Fidelity Solar Geolocation

`invTwilightFree` is the high-performance successor to the original `TwilightFree` package. It re-imagines solar geolocation by inverting the traditional workflow: instead of identifying discrete "twilights," it models the **continuous likelihood of the entire light curve**.

Powered by a high-speed **Rust engine** via `extendr`, it offers $50\times$ to $100\times$ speedups over traditional state-space models while providing superior robustness to sensor shading and non-solar light noise.

## Key Features

*   **Truly Twilight-Free:** No manual editing or thresholding of twilight events required. Processes raw light data directly.
*   **Rust-Powered SMC:** Implements Sequential Monte Carlo (Particle Filtering) in Rust for extreme performance.
*   **Equinox & Polar Support:** Handles the "blind spots" of equinoxes and the "midnight sun" of polar summers using continuous likelihood and Guided Brownian Bridges.
*   **Anomaly Robustness:** Automatically detects and down-weights artificial light at night (ALAN) and sensor shading using a built-in spike-and-slab likelihood model.

## Installation

`invTwilightFree` requires a working [Rust toolchain](https://rustup.rs/) and the `rextendr` package.

1.  **Install Rust:**
    Follow the instructions at [rustup.rs](https://rustup.rs/) for your OS.
2.  **Install from GitHub:**
    ```r
    # if (!requireNamespace("remotes")) install.packages("remotes")
    remotes::install_github("ABindoff/invTwilightFree")
    ```

## Quick Start

Reconstructing a track is as simple as providing your light data and start/end locations.

```r
library(invTwilightFree)

# 1. Load your data (Date/Time and Light intensity)
# track_data <- read.csv("your_data.csv")

# 2. Fit the model using Sequential Monte Carlo
fit <- TwilightFreeSMC(
  date_time = track_data$time,
  light = track_data$light,
  start_lat = -45.0, start_lon = 140.0,
  end_lat = -45.0, end_lon = 140.0,
  method = "guided",      # "guided" for speed, "ffbs" for maximum accuracy
  n_particles = 1000
)

# 3. Visualize the results
plot(fit, type = "track")       # Map of the estimated path
plot(fit, type = "diagnostics") # Check for anomalous light events
```

## Citation

If you use this method, please cite the original foundational work:

> Bindoff AD, Wotherspoon SJ, Guinet C, Hindell MA. Twilight-free geolocation from noisy light data. *Methods Ecol Evol.* 2017;00:1–9. https://doi.org/10.1111/2041-210X.12953

## License

This package is licensed under the GPL-3 License.
