suppressMessages(library(data.table))
o <- fread("scratch/nes_calibration/residence_results.csv")
o[, weight := sprintf("%.2f", weight)]
cat("=== per tag: light only against a gentle occupancy prior ===\n")
w <- merge(o[weight == "0.00", .(id, km0 = median_km, lat0 = rmse_lat, cov0 = cover_lat)],
           o[weight == "0.25", .(id, km25 = median_km, lat25 = rmse_lat, cov25 = cover_lat)],
           by = "id")
w[, `:=`(d_km = km25 - km0, d_lat = round(lat25 - lat0, 2), d_cov = round(cov25 - cov0, 2))]
print(as.data.frame(w[order(d_km)]), row.names = FALSE)
cat(sprintf("\nimproved on median error: %d of %d tags\n", sum(w$d_km < 0), nrow(w)))
cat(sprintf("improved on latitude RMSE: %d of %d\n", sum(w$d_lat < 0), nrow(w)))
cat(sprintf("paired Wilcoxon on median error: p = %.3g\n",
            wilcox.test(w$km0, w$km25, paired = TRUE)$p.value))
