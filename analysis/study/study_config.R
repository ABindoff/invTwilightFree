# =============================================================================
# FROZEN STUDY CONFIGURATION
#
# Every constant the reporting run depends on, in one place. Nothing downstream
# may define one of these locally. If a value changes here, the run is a
# different study and STUDY_VERSION must be bumped -- the output files carry it,
# so results from two versions can never be pooled by accident.
#
# The decisions recorded here were settled on 2026-08-26 and are documented with
# their evidence in PREREGISTRATION.md. Do not change one without changing that.
# =============================================================================

STUDY_VERSION <- "v1"

# ---- geography --------------------------------------------------------------
# Ano Nuevo, in the 0-360 east convention this study works in throughout. Two
# animals cross the antimeridian, so a signed convention corrupts their steps.
COLONY <- c(lon = 237.67, lat = 37.11)

# Search domain, defined from the species' range and NOT from the observed
# tracks. Verified in preflight to contain every Argos fix with margin to spare.
DOMAIN <- c(xmin = 150, xmax = 250, ymin = 20, ymax = 70)
GRID_CELL_DEG <- 1

# ---- calibration (settled 2026-08-26) ---------------------------------------
# These tags were deployed to record depth. There was no field calibration: no
# period on a post at a known position, before or after. The pre-departure
# haul-out is used as a SUBSTITUTE for that period, and the manuscript says so.
#
# CAL_GEOM_MAX_DAYS is set a priori to match the standard field protocol (about a
# week at a known position) and is NOT tuned against the Argos truth. The
# departure truncation is what keeps the window honest: departure ranges 0.0 to
# 16.9 days across this panel (median 6.8), so at day 7 only 13 of 29 animals are
# still ashore and the 90th percentile is already 537 km from the colony.
CAL_GEOM_MAX_DAYS <- 7      # cap on the geometry window
CAL_GEOM_MIN_OBS  <- 200    # below this a tag contributes nothing and uses the pool
CAL_SCALE_DAYS    <- 15     # intensity-scale window, per tag, unchanged from the campaign
DEPARTURE_DIVE_M  <- 50     # a day with >50% of bins deeper than this is "at sea"

# ---- ground truth (R2, R3) --------------------------------------------------
# ONE scorer for the whole study: argos_at() from nes_common.R, which returns NA
# when the bracketing Argos gap exceeds ARGOS_MAX_GAP_H. The permissive
# approx(rule = 2) variant that several campaign scripts defined locally is
# forbidden here -- it extrapolates past both track ends and adds the
# interpolator's error to the model's, preferentially where Argos is worst.
ARGOS_MAX_GAP_H <- 24
ARGOS_KEEP_LC   <- c("3", "2", "1", "0", "A", "B")
ARGOS_VMAX_KMH  <- 10       # sustained NES travel is 3-4 km/h

# ONE endpoint rule: both ends pinned at the colony for every deployment. This is
# prior knowledge for a central-place forager, not borrowed truth. The Argos
# deploy and recover positions are recorded as diagnostics and are never fitted.
ENDPOINT_RULE <- "colony"
ENDPOINT_MAX_DEVIATION_KM <- 100   # preflight flags any deployment beyond this

# ---- engines ----------------------------------------------------------------
STEP_HOURS <- 12
DIFFUSION  <- 110           # km/sqrt(day)
SEED       <- 42

PROB_SLAB   <- 0.10
SHADE_RATIO <- 2
AREA_CORRECTION <- TRUE     # cell-area Jacobian; area_prior() must NOT also be passed

# hierarchical
HIER_SWEEPS <- 4000L
HIER_BURN   <- 1500L
HIER_THIN   <- 5L
HIER_COARSE_RES <- 2
HIER_SURROGATE_DIFFUSION <- 110

# FFBS subset: the tags with the densest Argos coverage, where a continuous-space
# engine can be compared against the discretised ones most precisely.
FFBS_N <- 3L
FFBS_PARTICLES <- 1000L

# ---- sensor fusion ----------------------------------------------------------
# SST is OFF. Section 3.6 of the analysis established that Mk9 temperature
# degrades latitude when used as an SST proxy on these tags.
USE_SST <- FALSE
BATHY_MARGIN_M <- 300       # subtract from the knot's deepest dive before constraining
BATHY_SOFT_SD  <- 300
SOFT_MASK_PENALTY <- -20

# ---- paths ------------------------------------------------------------------
CACHE_DIR <- "analysis/cache/nes"
OUT_DIR   <- file.path("analysis/output/study", STUDY_VERSION)
