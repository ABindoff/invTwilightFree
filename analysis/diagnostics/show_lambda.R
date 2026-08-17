suppressMessages(library(data.table))
o <- fread("scratch/nes_calibration/lambda_results.csv")
cat(sprintf("%d fits, multipliers present: %s\n", nrow(o),
            paste(sort(unique(o$mult)), collapse = ", ")))
full <- o[, .N, by = mult][N == uniqueN(o$id)]$mult
cat("complete multipliers:", paste(sort(full), collapse = ", "), "\n\n")

cat("=== pooled, complete multipliers only ===\n")
print(as.data.frame(o[mult %in% full, .(tags = .N,
  median_km = round(mean(median_km)), rmse_lat = round(mean(rmse_lat), 2),
  post_sd_lat = round(mean(post_sd_lat), 2),
  cover_lat = round(mean(cover_lat), 2), cover_lon = round(mean(cover_lon), 2),
  mean_log_z = round(mean(log_z), 1)), by = mult][order(mult)]), row.names = FALSE)

cat("\n=== log_z per tag, relative to the smallest multiplier so far ===\n")
w <- dcast(o, id ~ mult, value.var = "log_z")
mm <- setdiff(names(w), "id")
base <- w[[mm[1]]]
for (m in mm) w[[m]] <- round(w[[m]] - base, 1)
print(as.data.frame(w), row.names = FALSE)
cat(sprintf("\n(columns are multipliers; values are log_z gain over multiplier %s)\n", mm[1]))

cat("\n=== does the posterior width respond to lambda at all? ===\n")
print(as.data.frame(dcast(o[mult %in% full], id ~ mult, value.var = "post_sd_lat")),
      row.names = FALSE)
cat("\n=== and coverage ===\n")
print(as.data.frame(dcast(o[mult %in% full], id ~ mult, value.var = "cover_lat")),
      row.names = FALSE)
