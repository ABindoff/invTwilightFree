# The new delivery ships message-level RawArgos rather than processed Locations,
# so before writing the ingest: how many rows per pass, which columns are
# actually populated, what location classes appear, and whether the mirror
# solution (Latitude2/Longitude2) ever differs from the primary.
suppressMessages(library(data.table))
Sys.setlocale("LC_TIME", "C")
D <- "fvilches/extracted/Argos raw"
f <- list.files(D, pattern = "RawArgos[.]csv$", full.names = TRUE)
cat(sprintf("%d RawArgos files\n\n", length(f)))

a <- fread(f[grep("2023030", f)], showProgress = FALSE)
cat("=== column fill rates (2023030) ===\n")
fill <- sapply(a, function(x) mean(!is.na(x) & x != ""))
print(round(fill[fill > 0], 3))

cat("\n=== rows per pass ===\n")
a[, passkey := paste(PTT, PassDate, PassTime)]
rp <- a[, .N, by = passkey]
cat(sprintf("passes: %d   rows: %d   rows/pass: median %.0f, max %d\n",
            nrow(rp), nrow(a), median(rp$N), max(rp$N)))

cat("\n=== location classes ===\n")
print(table(a$Class, useNA = "ifany"))

cat("\n=== does the position vary within a pass? ===\n")
v <- a[, .(nlat = uniqueN(Latitude), nlon = uniqueN(Longitude)), by = passkey]
cat(sprintf("passes with >1 distinct latitude: %d of %d\n",
            sum(v$nlat > 1), nrow(v)))

cat("\n=== primary vs mirror solution ===\n")
u <- unique(a, by = "passkey")
u <- u[!is.na(Latitude) & !is.na(Latitude2)]
cat(sprintf("passes where the mirror differs from the primary: %d of %d\n",
            sum(abs(u$Latitude - u$Latitude2) > 1e-6 |
                abs(u$Longitude - u$Longitude2) > 1e-6), nrow(u)))

cat("\n=== error ellipse availability ===\n")
cat(sprintf("passes with an error radius: %.2f\n",
            mean(!is.na(u$`Error radius`))))
cat(sprintf("passes with a semi-major axis: %.2f\n",
            mean(!is.na(u$`Error Semi-major axis`))))

cat("\n=== first and last fixes (deployment / recovery) ===\n")
u[, time := as.POSIXct(paste(PassDate, PassTime), format = "%d-%b-%Y %H:%M:%S", tz = "UTC")]
setorder(u, time)
print(head(u[, .(time, Class, Latitude, Longitude)], 4))
print(tail(u[, .(time, Class, Latitude, Longitude)], 3))
cat(sprintf("\nspan: %s to %s (%.0f days)\n", min(u$time), max(u$time),
            as.numeric(difftime(max(u$time), min(u$time), units = "days"))))

cat("\n=== do all files parse the same way? ===\n")
chk <- rbindlist(lapply(f, function(p) {
  x <- fread(p, showProgress = FALSE, nrows = 5000)
  data.table(file = basename(p), cols = ncol(x),
             has_class = "Class" %in% names(x),
             has_latlon = all(c("Latitude", "Longitude") %in% names(x)))
}))
print(chk[, .(files = .N), by = .(cols, has_class, has_latlon)])
