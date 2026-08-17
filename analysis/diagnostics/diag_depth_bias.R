# Is the residual latitude bias driven by DEPTH ATTENUATION of the light?
#
# The bias is in the light, not the movement model or the sampler (eight
# mechanisms tested and excluded). The largest measured driver of light error in
# this dataset is depth, and the observation model currently treats attenuation
# as noise: the spike absorbs it. But attenuation is not noise. It is roughly
# Beer-Lambert, exponential in depth, and THE DEPTH IS MEASURED.
#
# The mechanism predicts a specific signature. Attenuation makes the model see
# less light than the true surface irradiance, so it infers a lower sun, which
# shortens the apparent day. In northern summer a shorter day reads as further
# SOUTH; in northern winter, further NORTH. Depth varies seasonally with
# foraging, so the bias should be seasonally varying and tag-specific -- which is
# the structure observed, including 2021032's 12-degree sign flip.
#
# THE TEST, within tag and within month so that between-animal differences and
# the common seasonal curve cannot produce it:
#
#   does month-to-month variation in dive depth predict month-to-month
#   variation in latitude bias, in the direction the physics requires?
#
# Prediction: deeper months should be biased SOUTH in summer (delta > 0) and
# NORTH in winter, i.e. the interaction depth x sign(declination) is what
# carries the signal, not depth alone. Testing depth alone would miss it and
# testing it pooled across seasons would cancel it.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))

declination <- function(t) {
  doy <- as.numeric(format(t, "%j"))
  g <- 2*pi/365 * (doy - 1)
  (0.006918 - 0.399912*cos(g) + 0.070257*sin(g) - 0.006758*cos(2*g) +
     0.000907*sin(2*g) - 0.002697*cos(3*g) + 0.00148*sin(3*g)) * 180/pi
}

k <- as.data.table(read.csv("analysis/output/nes/per_knot_errors.csv"))
k <- k[engine == "hier" & batch == "light" & is.finite(err_lat)]
k[, `:=`(id = as.character(id), time = as.POSIXct(time, tz = "UTC"))]
k[, `:=`(mon = format(time, "%Y-%m"), decl = declination(time))]

# depth per knot window, from the archive
arch <- lapply(readRDS("analysis/cache/nes/archives_v3_30min.rds"), `[[`, "main")
dep <- rbindlist(lapply(names(arch), function(id) {
  x <- as.data.table(arch[[id]])
  x[, .(id = id, mon = format(time, "%Y-%m"),
        med_depth = median(depth_max, na.rm = TRUE),
        frac_deep = mean(depth_max > 400, na.rm = TRUE),
        # how shallow does it get? the max light in a window comes from here
        med_shallow = median(depth_min, na.rm = TRUE),
        frac_surf = mean(depth_min < 10, na.rm = TRUE)), by = format(time, "%Y-%m")][, -1]
}))
dep <- dep[, lapply(.SD, mean, na.rm = TRUE), by = .(id, mon)]

M <- merge(k[, .(bias = mean(err_lat), decl = mean(decl), n = .N), by = .(id, mon)],
           dep, by = c("id", "mon"))
M <- M[n >= 20]
cat(sprintf("%d tag-months across %d tags\n\n", nrow(M), uniqueN(M$id)))

# within-tag centring: remove each animal's own mean, so between-animal
# differences cannot create the association
M[, `:=`(bias_c = bias - mean(bias), depth_c = med_depth - mean(med_depth),
         surf_c = frac_surf - mean(frac_surf)), by = id]
M[, season := ifelse(decl > 0, "summer (delta>0)", "winter (delta<0)")]

cat("=== within-tag: does a deeper month shift the bias, and which way? ===\n")
print(as.data.frame(M[, .(tag_months = .N,
  r_depth = round(cor(depth_c, bias_c), 3),
  r_surf  = round(cor(surf_c, bias_c), 3)), by = season]), row.names = FALSE)

cat("\nPrediction: r_depth NEGATIVE in summer (deeper -> apparent shorter day ->\n")
cat("biased south) and POSITIVE in winter. Opposite signs are the signature;\n")
cat("a common sign in both seasons would be something else.\n")

fit <- lm(bias_c ~ depth_c : I(sign(decl)) + depth_c, data = M)
cat("\n=== the interaction test ===\n")
print(round(summary(fit)$coefficients, 5))
cat(sprintf("\nR2 = %.3f on %d tag-months\n", summary(fit)$r.squared, nrow(M)))

cat("\n=== for scale: how much attenuation is there to model? ===\n")
allx <- rbindlist(lapply(names(arch), function(id) {
  x <- as.data.table(arch[[id]]); x[, .(id = id, depth_min, depth_max, light)] }))
cat(sprintf("median depth of the SHALLOWEST point in a 30-min window: %.0f m\n",
            median(allx$depth_min, na.rm = TRUE)))
cat(sprintf("fraction of windows reaching within 10 m of the surface: %.2f\n",
            mean(allx$depth_min < 10, na.rm = TRUE)))
cat("If most windows never reach the surface, the retained maximum light is\n")
cat("attenuated even after decimation, and the attenuation is measured.\n")
