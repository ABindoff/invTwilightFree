suppressMessages(library(data.table))
o <- fread("scratch/nes_calibration/all29_tags.csv")
cat("=== does borrowed geometry explain the seasonal degradation? ===\n")
print(as.data.frame(o[, .(n = .N, median_km = round(mean(median_km)),
  rmse_lat = round(mean(rmse_lat), 2), cover_lat = round(mean(cover_lat), 2)),
  by = .(season, own_geometry)][order(season, own_geometry)]), row.names = FALSE)

cat("\n=== overall, by whether the tag had its own haul-out ===\n")
print(as.data.frame(o[, .(n = .N, median_km = round(mean(median_km)),
  rmse_lat = round(mean(rmse_lat), 2), cover_lat = round(mean(cover_lat), 2)),
  by = own_geometry]), row.names = FALSE)
w <- wilcox.test(median_km ~ own_geometry, data = o)
cat(sprintf("Wilcoxon own vs borrowed geometry: p = %.3f\n", w$p.value))

cat("\n=== worst and best deployments ===\n")
setorder(o, -median_km)
print(as.data.frame(o[, .(id, season, family, own_geometry, median_km,
                          rmse_lat, cover_lat, n_scored)][c(1:5, (.N-4):.N)]),
      row.names = FALSE)

cat("\n=== is 2021 reproduced? (was 208-209 km in the earlier analysis) ===\n")
cat(sprintf("2021 here: %.0f km over %d deployments\n",
            mean(o[season == "2021"]$median_km), nrow(o[season == "2021"])))

cat("\n=== per-tag latitude scatter, the quantity that varies 4x ===\n")
k <- fread("scratch/nes_calibration/all29_knots.csv")
s <- k[is.finite(err_lat), .(scatter = round(sd(err_lat), 2),
                             offset = round(mean(err_lat), 2)), by = id]
cat(sprintf("scatter ranges %.2f to %.2f across %d deployments (was 1.47-6.16 at n=10)\n",
            min(s$scatter), max(s$scatter), nrow(s)))
cat(sprintf("offset ranges %.2f to %.2f\n", min(s$offset), max(s$offset)))
