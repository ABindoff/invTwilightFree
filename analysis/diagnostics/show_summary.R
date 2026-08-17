suppressMessages(library(data.table))
o <- fread("scratch/nes_calibration/summary_stat_results.csv")
cat(sprintf("%d tags complete\n\n", uniqueN(o$id)))
cat("=== pooled over tags ===\n")
print(as.data.frame(o[, .(tags = .N, median_km = round(mean(median_km)),
  q90_km = round(mean(q90_km)), bias_lat = round(mean(bias_lat), 3),
  rmse_lat = round(mean(rmse_lat), 3), rmse_lon = round(mean(rmse_lon), 3)),
  by = summary][order(median_km)]), row.names = FALSE)

cat("\n=== paired against the mode (the current behaviour) ===\n")
base <- o[summary == "mode", .(id, km0 = median_km, la0 = rmse_lat, lo0 = rmse_lon)]
for (s in c("mean", "med", "trim")) {
  m <- merge(base, o[summary == s, .(id, km = median_km, la = rmse_lat, lo = rmse_lon)],
             by = "id")
  if (!nrow(m)) next
  p <- tryCatch(wilcox.test(m$km0, m$km, paired = TRUE)$p.value, error = function(e) NA)
  cat(sprintf("  %-5s d_median %+6.1f km   better %d/%d   d_rmse_lat %+.3f   d_rmse_lon %+.3f   p=%.3f\n",
              s, mean(m$km - m$km0), sum(m$km < m$km0), nrow(m),
              mean(m$la - m$la0), mean(m$lo - m$lo0), p))
}
cat("\nNegative favours the alternative over the mode.\n")
