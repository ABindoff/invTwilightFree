data_dir <- list.files("scratch/probGLS_repo", pattern="probGLS-master", full.names=TRUE)[1]
data_path <- file.path(data_dir, "data")

cat("Loading BBA_lux.RData...\n")
load(file.path(data_path, "BBA_lux.RData"))
cat("BBA_lux structure:\n")
str(BBA_lux)
cat("\nSaving BBA_lux to scratch/albatross_gls.rds...\n")
saveRDS(BBA_lux, "scratch/albatross_gls.rds")

cat("\nLoading BBA_deg.RData...\n")
load(file.path(data_path, "BBA_deg.RData"))
cat("BBA_deg structure:\n")
str(BBA_deg)
cat("\nSaving BBA_deg to scratch/albatross_gps.rds...\n")
saveRDS(BBA_deg, "scratch/albatross_gps.rds")

cat("\nSuccess!\n")
