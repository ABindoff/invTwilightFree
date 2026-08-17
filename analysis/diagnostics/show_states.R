suppressMessages(library(data.table))
o <- fread("scratch/nes_calibration/states_results.csv")
cat(sprintf("%d fits on disk across %d configs\n", nrow(o), uniqueN(o$config)))
full <- o[, .N, by = config][N == uniqueN(o$id)]$config
cat("complete configs:", paste(full, collapse = " | "), "\n\n")
p <- o[config %in% full, .(tags = .N,
  median_km = round(mean(median_km)), worst_km = max(median_km),
  rmse_lat = round(mean(rmse_lat), 2),
  post_sd_lat = round(mean(post_sd_lat), 2),
  post_sd_lon = round(mean(post_sd_lon), 2),
  cover_lat = round(mean(cover_lat), 2), cover_lon = round(mean(cover_lon), 2),
  mean_log_z = round(mean(log_z), 1), mins = round(sum(secs)/60, 1)), by = config]
print(as.data.frame(p[order(median_km)]), row.names = FALSE)

cat("\n=== log_z: do the DATA prefer switching? (paired, against 1-state D=110) ===\n")
base <- o[config == "1-state D=110 (report)", .(id, lz0 = log_z, km0 = median_km,
                                                cov0 = cover_lat)]
for (cf in setdiff(full, "1-state D=110 (report)")) {
  m <- merge(base, o[config == cf, .(id, lz = log_z, km = median_km, cov = cover_lat)],
             by = "id")
  cat(sprintf("  %-26s d_log_z %+8.1f  better in %d/%d   d_km %+6.1f  d_cover %+.2f\n",
              cf, mean(m$lz - m$lz0), sum(m$lz > m$lz0), nrow(m),
              mean(m$km - m$km0), mean(m$cov - m$cov0)))
}
cat("\n=== per tag, latitude coverage ===\n")
w <- dcast(o[config %in% full], id ~ config, value.var = "cover_lat")
print(as.data.frame(w), row.names = FALSE)
