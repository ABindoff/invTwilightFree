suppressMessages(library(data.table))
o <- fread("scratch/nes_calibration/anisotropy_results.csv")
cat(sprintf("%d fits; configs present: %s\n\n", nrow(o),
            paste(unique(o$config), collapse = " | ")))
cat("=== pooled (complete configs in bold territory; partial marked) ===\n")
p <- o[, .(tags = .N, median_km = round(mean(median_km)),
           rmse_lat = round(mean(rmse_lat), 2), rmse_lon = round(mean(rmse_lon), 2),
           cover_lat = round(mean(cover_lat), 2), cover_lon = round(mean(cover_lon), 2),
           pi_lat = round(mean(pi_lat), 3), pi_lon = round(mean(pi_lon), 3),
           log_z = round(mean(log_z))), by = config]
p[, complete := tags == max(tags)]
print(as.data.frame(p[order(ns <- config)]), row.names = FALSE)

# like-for-like on the tags that have both, so a partial config is not compared
# against a different set of animals
base <- "iso 110/110 (current)"
for (cf in setdiff(unique(o$config), base)) {
  m <- merge(o[config == base, .(id, km0 = median_km, cl0 = cover_lat, co0 = cover_lon)],
             o[config == cf, .(id, km = median_km, cl = cover_lat, co = cover_lon)],
             by = "id")
  if (!nrow(m)) next
  cat(sprintf("\n%s vs isotropic, on the %d tags with both:\n", cf, nrow(m)))
  cat(sprintf("  median km  %4.0f -> %4.0f\n", mean(m$km0), mean(m$km)))
  cat(sprintf("  cover_lat  %.2f -> %.2f   (target 0.95, was under-covering)\n",
              mean(m$cl0), mean(m$cl)))
  cat(sprintf("  cover_lon  %.2f -> %.2f   (target 0.95, was over-covering)\n",
              mean(m$co0), mean(m$co)))
}
