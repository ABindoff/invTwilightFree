# Download Dual-Tagged Wandering Albatross Dataset from probGLS repository
cat("Downloading dual-tagged Wandering Albatross dataset...\n")

# URL to the raw .rda file in the probGLS GitHub repo
url <- "https://github.com/benjamin-merkel/probGLS/raw/master/data/walb.rda"
dest_file <- "scratch/walb.rda"

# Download the file
download.file(url, destfile = dest_file, mode = "wb")

# Load the dataset
load(dest_file)

# walb should now be in the environment. It is a list containing:
# walb$gls: The light-level geolocation data (twilight events)
# walb$gps: The high-resolution GPS track

cat("Download complete. Extracting GPS and GLS components...\n")

str(walb)

# Let's save them as clean CSVs or RDS files for easy access in our package
saveRDS(walb$gls, "scratch/albatross_gls.rds")
saveRDS(walb$gps, "scratch/albatross_gps.rds")

cat("\nSuccess! Saved to scratch/albatross_gls.rds and scratch/albatross_gps.rds\n")
cat("GPS points:", nrow(walb$gps), "\n")
cat("GLS twilight events:", nrow(walb$gls), "\n")
