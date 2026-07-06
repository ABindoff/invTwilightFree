# build_vignettes_all.R
#
# Knits all vignettes and populates inst/doc/ for CRAN submission.
# Run from the project root, from RStudio / Positron (pandoc required).
#
# Usage:
#   source("build_vignettes_all.R")        # from RStudio console
#   ! Rscript build_vignettes_all.R        # from Claude Code terminal
#
# Benchmark vignettes need caches first — run scratch/run_benchmark2_prep.R
# once (requires SGAT, FLightR, TwGeos installed) before running this.

Sys.setenv(NOT_CRAN = "true")

if (!rmarkdown::pandoc_available()) {
  stop("Pandoc not found. Run this from RStudio or Positron, which bundle pandoc.")
}

stopifnot(file.exists("DESCRIPTION"))   # must be run from project root

timer <- function(expr) {
  t <- system.time(force(expr))
  cat(sprintf("  done in %.0f s\n", t["elapsed"]))
}

# ---------------------------------------------------------------------------
# 1. Reinstall the package so the latest source is used during rendering.
# ---------------------------------------------------------------------------
cat("Installing package from source...\n")
tryCatch(
  devtools::install(".", quiet = TRUE, dependencies = FALSE),
  error = function(e) {
    warning(
      "Could not reinstall package (DLL probably locked by another R session).\n",
      "Continuing with the currently installed version.\n",
      "If R or Rust source has changed, close all other R sessions first and re-run.",
      call. = FALSE
    )
  }
)

# ---------------------------------------------------------------------------
# 2. Official vignettes -> inst/doc/  (what CRAN sees)
#
#    Three files per vignette go into inst/doc/:
#      .html  — prebuilt output (avoids CRAN rebuilding from source)
#      .Rmd   — source copy (required by R)
#      .R     — knitr::purl() tangle (what R CMD check runs)
# ---------------------------------------------------------------------------
dir.create("inst/doc", recursive = TRUE, showWarnings = FALSE)

official <- c(
  "vignettes/sensor_fusion.Rmd",
  "vignettes/calibration_tuning.Rmd"
)

for (v in official) {
  cat("\nKnitting", basename(v), "...\n")
  timer({
    rmarkdown::render(v, output_dir = "inst/doc/", quiet = TRUE)
    knitr::purl(
      v,
      output        = file.path("inst/doc", sub("\\.Rmd$", ".R", basename(v))),
      documentation = 0L,
      quiet         = TRUE
    )
    file.copy(v, "inst/doc/", overwrite = TRUE)
  })
}

# ---------------------------------------------------------------------------
# 3. Benchmark vignettes -> vignettes/  (local inspection; excluded from build)
# ---------------------------------------------------------------------------
benchmarks <- c(
  "vignettes/benchmark_comparison.Rmd",
  "vignettes/benchmark_comparison2.Rmd",
  "vignettes/comprehensive_benchmark.Rmd"
)

has_data  <- file.exists("scratch/simulated_seal_light_scenarios.rds")
has_cache <- length(list.files("vignettes", pattern = "^cache_flightr_seal_cloudy")) > 0

if (!has_data) {
  cat("\nSKIPPING benchmark vignettes — scratch/simulated_seal_light_scenarios.rds not found.\n")
} else if (!has_cache) {
  cat("\nSKIPPING benchmark vignettes — FLightR/SGAT caches not found.\n")
  cat("Generate them first: source('scratch/run_benchmark2_prep.R')\n")
} else {
  for (v in benchmarks) {
    if (!file.exists(v)) next
    cat("\nKnitting", basename(v), "...\n")
    timer(
      rmarkdown::render(v, output_dir = "vignettes/", quiet = TRUE)
    )
  }
}

# ---------------------------------------------------------------------------
# 4. Quick sanity check (no vignette re-run).
# ---------------------------------------------------------------------------
cat("\nRunning devtools::check() (skipping vignettes)...\n")
devtools::check(vignettes = FALSE, document = FALSE)

cat("\n=== Done ===\n")
cat("  Official vignettes: inst/doc/\n")
cat("  Benchmark output:   vignettes/*.html  (if caches were present)\n")
