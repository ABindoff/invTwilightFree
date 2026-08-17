SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
suppressMessages(library(data.table))
o <- readRDS(file.path(SCRATCH, "domain.rds"))
cat("=== pooled over tags ===\n")
print(as.data.frame(o[, .(tags = .N, median_km = round(mean(median_km)),
  worst_km = max(median_km), q90_km = round(mean(q90_km)),
  bias_lat = round(mean(bias_lat), 3), rmse_lat = round(mean(rmse_lat), 3),
  lat_min = min(lat_min), lat_max = max(lat_max),
  frac_at_lat_edge = round(mean(at_lat_edge), 4)),
  by = domain][order(domain)]), row.names = FALSE)

cat("\n=== are the fits literally identical across domains? ===\n")
w <- dcast(o, id ~ domain, value.var = "median_km")
print(as.data.frame(w), row.names = FALSE)
d <- dcast(o, id ~ domain, value.var = "rmse_lat")
cat("\nrmse_lat:\n"); print(as.data.frame(d), row.names = FALSE)
nm <- setdiff(names(w), "id")
cat(sprintf("\nmax absolute difference in median_km between domains: %.3f\n",
            max(abs(w[[nm[1]]] - w[[nm[2]]]), abs(w[[nm[1]]] - w[[nm[3]]]))))
cat(sprintf("max absolute difference in rmse_lat: %.4f\n",
            max(abs(d[[nm[1]]] - d[[nm[2]]]), abs(d[[nm[1]]] - d[[nm[3]]]))))

cat("\n=== how close does any fit get to the latitude bounds? ===\n")
print(as.data.frame(o[, .(id, domain, lat_min, lat_max, at_lat_edge)][order(id, domain)]),
      row.names = FALSE)
