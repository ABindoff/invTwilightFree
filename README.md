# invTwilightFree: High-Fidelity Solar Geolocation

`invTwilightFree` is the modern, high-performance successor to the original `TwilightFree` method. Building on the foundational approach of modeling the **continuous likelihood of the entire light curve** (avoiding the need for manual twilight identification), `invTwilightFree` introduces a re-engineered engine designed for the next generation of animal tracking.

Powered by a high-speed **Rust engine** via `extendr`, it provides the same "twilight-free" workflow with up to **~195× speedups** vs. SGAT/FLightR on 180-day deployments (41 s vs. ~132 min), 9–10× better accuracy under adverse shading conditions, and significantly enhanced robustness to sensor noise.

## Key Advancements

*   **Rust-Powered SMC:** A completely new Sequential Monte Carlo (Particle Filtering) engine implemented in Rust for extreme performance and scalability.
*   **Built-in Anomaly Detection:** Uses a sophisticated spike-and-slab likelihood model to automatically handle Artificial Light at Night (ALAN) and erratic sensor shading without manual preprocessing.
*   **Optimized Smoothers:** Includes both fast "Guided" Brownian Bridges for real-time processing and rigorous "FFBS" (Forward-Filtering Backward-Smoothing) for publication-quality tracks.
*   **Equinox & Polar Support:** Expertly handles the "blind spots" of equinoxes and the "midnight sun" of polar summers.

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
  start_time = as.POSIXct("2024-01-05"), # Optional: ignore pre-deployment data
  end_time = as.POSIXct("2024-06-15"),   # Optional: ignore post-retrieval data
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

> Bindoff AD, Wotherspoon SJ, Guinet C, Hindell MA. Twilight-free geolocation from noisy light data. *Methods Ecol Evol.* 2018;9(5):1190-1198. https://doi.org/10.1111/2041-210X.12953

## License

This package is licensed under the GPL-3 License.
