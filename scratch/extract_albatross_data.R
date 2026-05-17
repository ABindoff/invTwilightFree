cat("Downloading probGLS repository archive...\n")

temp <- tempfile()
download.file("https://github.com/benjamin-merkel/probGLS/archive/refs/heads/master.zip", temp, mode="wb")

cat("Extracting archive...\n")
unzip(temp, exdir = "scratch/probGLS_repo")
unlink(temp)

# Find the walb data file. It might be .rda or .RData
data_dir <- list.files("scratch/probGLS_repo", pattern="probGLS-master", full.names=TRUE)[1]
data_path <- file.path(data_dir, "data")
walb_files <- list.files(data_path, pattern="walb", full.names=TRUE)

if (length(walb_files) > 0) {
  walb_file <- walb_files[1]
  cat(sprintf("Found data file: %s\n", walb_file))
  
  # Load the environment
  e <- new.env()
  load(walb_file, envir = e)
  
  if ("walb" %in% ls(e)) {
    walb <- e$walb
    
    cat("Saving GPS and GLS components...\n")
    saveRDS(walb$gls, "scratch/albatross_gls.rds")
    saveRDS(walb$gps, "scratch/albatross_gps.rds")
    
    cat("Success!\n")
    cat("GLS (Light) Observations:", nrow(walb$gls), "\n")
    cat("GPS Ground-truth Observations:", nrow(walb$gps), "\n")
  } else {
    cat("Error: 'walb' object not found in the loaded file.\n")
    print(ls(e))
  }
} else {
  cat("Error: Could not find walb data file in the repository.\n")
  print(list.files(data_path))
}
