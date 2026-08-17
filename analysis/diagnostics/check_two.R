# Two deployments survived the clip badly. Both need a decision, not a guess.
#   2023034 (PTT 214400): zero fixes inside its own archive span
#   2023039 (PTT 214397): first in-window fix 1315 km from the colony
# PTTs are reused across seasons, so the question for both is whether the file
# actually contains this deployment's fixes at all.
suppressMessages(library(data.table))
Sys.setlocale("LC_TIME", "C")
lon360 <- function(x) (x %% 360 + 360) %% 360
gc_km <- function(l1, p1, l2, p2, R = 6371.0088) { d <- pi/180
  2*R*asin(pmin(1, sqrt(sin((p2-p1)*d/2)^2 + cos(p1*d)*cos(p2*d)*sin((l2-l1)*d/2)^2))) }
COL <- c(237.67, 37.11)

arch <- readRDS("analysis/cache/nes/fvilches_archives_v1_30min.rds")
for (tg in c("2023034", "2023039")) {
  f <- list.files("fvilches/extracted/Argos raw", pattern = paste0("^", tg, "_"),
                  full.names = TRUE)
  a <- fread(f[1], showProgress = FALSE)
  a[, time := as.POSIXct(paste(PassDate, PassTime), format = "%d-%b-%Y %H:%M:%S", tz = "UTC")]
  a <- a[!is.na(time) & !is.na(Latitude)]
  d <- arch[[tg]]$main
  cat(sprintf("\n=== %s  (%s) ===\n", tg, basename(f[1])))
  cat(sprintf("archive span : %s .. %s\n", min(d$time), max(d$time)))
  cat(sprintf("argos  span  : %s .. %s  (%d fixes)\n", min(a$time), max(a$time), nrow(a)))
  ov <- a[time >= min(d$time) - 3*86400 & time <= max(d$time) + 3*86400]
  cat(sprintf("fixes inside the archive span: %d\n", nrow(ov)))
  if (nrow(ov)) {
    setorder(ov, time)
    ov[, km := gc_km(lon360(Longitude), Latitude, COL[1], COL[2])]
    cat("first 3 in-window fixes:\n")
    print(head(ov[, .(time, Class, lat = round(Latitude,3),
                      lon = round(lon360(Longitude),3), km = round(km))], 3))
    cat(sprintf("days from archive start to first fix: %.1f\n",
                as.numeric(difftime(min(ov$time), min(d$time), units = "days"))))
  }
  # what years does the file actually cover?
  cat("fixes by year in the file:\n")
  print(table(format(a$time, "%Y")))
}
