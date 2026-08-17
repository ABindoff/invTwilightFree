# Ingest of the 2022 and 2023 seasons (F. Vilches).
#
# The TDR archives are byte-for-byte the same format as the 2021 delivery, so
# the validated decimation is reused unchanged. The Argos side is NOT: this
# delivery ships message-level RawArgos rather than processed Locations, which
# needs collapsing to one fix per pass and gives no `Type == "User"` row for the
# deployment position. It also carries error ellipses on every pass, which the
# processed files did not, so truth can be weighted by its own uncertainty later.
#
# Provenance is recorded per deployment (season, source folder, tag serial, PTT)
# because the physical tags are redeployed across years and the analysis will
# need deployments nested within tags.
#
# WHAT THIS DOES NOT DO: it does not re-read the 2021 season, which is already
# cached from the earlier delivery in the identical format. The manifest records
# which delivery each deployment came from.
suppressMessages({ library(data.table) })
Sys.setlocale("LC_TIME", "C")
TDR <- "fvilches/extracted/TDR raw"
ARG <- "fvilches/extracted/Argos raw"
CACHE <- "analysis/cache/nes"
DECIMATE_MIN <- 30
# The 18A family names the temperature column "External Temperature" where the
# Mk9 219xxxx family uses "External Temp". Matching the name exactly silently
# drops four deployments that DO carry light, so match on what is needed and
# resolve the temperature column by either name.
LIGHT_NEED <- c("Time", "Depth", "Light Level")
TEMP_NAMES <- c("External Temp", "External Temperature")

lon360 <- function(x) (x %% 360 + 360) %% 360
dlon <- function(a, b) ((a - b + 180) %% 360) - 180
gc_km <- function(l1, p1, l2, p2, R = 6371.0088) { d <- pi/180
  2*R*asin(pmin(1, sqrt(sin((p2-p1)*d/2)^2 + cos(p1*d)*cos(p2*d)*sin((l2-l1)*d/2)^2))) }

# ---- the validated decimation, unchanged from the 2021 analysis -------------
# Maximum light over a window longer than the dive cycle, timestamped at the
# sample the value came from. Both details are load-bearing: subsampling samples
# the dive state rather than the light (r 0.43 against 0.81), and timestamping
# at the window start injects a pure longitude bias of (minutes/2)/4 degrees.
decimate_dive <- function(d, minutes) {
  d[, .b := floor(as.numeric(time) / (minutes * 60))]
  out <- d[, {
    j <- which.max(light)
    i <- which.min(depth)
    .(time = time[j],
      light = as.numeric(light[j]),
      depth_max = as.numeric(max(depth, na.rm = TRUE)),
      depth_min = as.numeric(min(depth, na.rm = TRUE)),
      temp_surf = if (length(i)) as.numeric(temp[i]) else NA_real_)
  }, by = .b]
  out[, .b := NULL]
  setorder(out, time)
  as.data.frame(out[is.finite(light)])
}

read_archive <- function(path) {
  hdr <- names(fread(path, nrows = 0L))
  if (!all(LIGHT_NEED %in% hdr)) return(NULL)          # genuinely no light channel
  tcol <- intersect(TEMP_NAMES, hdr)
  cols <- c("Time", "Depth", if (length(tcol)) tcol[1] else NULL, "Light Level")
  d <- fread(path, select = cols, showProgress = FALSE)
  setnames(d, cols, c("time", "depth", if (length(tcol)) "temp" else NULL, "light"))
  if (!length(tcol)) d[, temp := NA_real_]
  d[, time := as.POSIXct(time, format = "%H:%M:%S %d-%b-%Y", tz = "UTC")]
  d <- d[!is.na(time) & is.finite(light)]
  if (!nrow(d)) return(NULL)
  decimate_dive(d, DECIMATE_MIN)
}

# ---- RawArgos -> one fix per pass -------------------------------------------
read_rawargos <- function(path) {
  a <- fread(path, showProgress = FALSE)
  need <- c("PTT", "Class", "PassDate", "PassTime", "Latitude", "Longitude")
  if (!all(need %in% names(a))) return(NULL)
  a[, time := as.POSIXct(paste(PassDate, PassTime),
                         format = "%d-%b-%Y %H:%M:%S", tz = "UTC")]
  a <- a[!is.na(time) & !is.na(Latitude) & !is.na(Longitude)]
  if (!nrow(a)) return(NULL)
  # Rank the classes so that collapsing a multi-message pass keeps its best
  # solution rather than whichever row happened to come first.
  rank <- c("3" = 1, "2" = 2, "1" = 3, "0" = 4, "A" = 5, "B" = 6, "Z" = 9)
  a[, cr := rank[as.character(Class)]]
  a[is.na(cr), cr := 9L]
  setorder(a, time, cr)
  a <- unique(a, by = "time")
  er <- if ("Error radius" %in% names(a)) a[["Error radius"]] else NA_real_
  data.table(time = a$time, ptt = as.character(a$PTT),
             class = as.character(a$Class),
             lon = lon360(a$Longitude), lat = a$Latitude,
             err_radius_m = as.numeric(er))[class %in% c("3","2","1","0","A","B")]
}

# Forward speed filter, as in the 2021 analysis: sustained travel is 3-4 km/h,
# so a step implying more than 10 km/h is a bad solution rather than a fast seal.
speed_filter <- function(d, vmax = 10) {
  setorder(d, time)
  if (nrow(d) < 3) return(d)
  keep <- rep(TRUE, nrow(d)); last <- 1L
  for (i in 2:nrow(d)) {
    dt <- as.numeric(difftime(d$time[i], d$time[last], units = "hours"))
    if (dt <= 0) { keep[i] <- FALSE; next }
    v <- gc_km(d$lon[last], d$lat[last], d$lon[i], d$lat[i]) / dt
    if (v > vmax) keep[i] <- FALSE else last <- i
  }
  d[keep]
}

# ---- inventory --------------------------------------------------------------
tdr <- data.table(path = list.files(TDR, pattern = "Archive[.]csv$", full.names = TRUE))
tdr[, file := basename(path)]
tdr[, topp := sub("^([0-9]+)_.*$", "\\1", file)]
tdr[, serial := sub("^[0-9]+_(.+)-out-Archive[.]csv$", "\\1", file)]
tdr[, season := substr(topp, 1, 4)]
tdr <- tdr[season %in% c("2022", "2023")]          # 2021 is already cached

arg <- data.table(path = list.files(ARG, pattern = "RawArgos[.]csv$", full.names = TRUE))
arg[, file := basename(path)]
arg[, topp := sub("^([0-9]+)_.*$", "\\1", file)]
arg[, ptt := sub("^[0-9]+_(.+)-RawArgos[.]csv$", "\\1", file)]

cat(sprintf("to ingest: %d TDR archives (2022, 2023); %d Argos files available\n\n",
            nrow(tdr), nrow(arg)))

# ---- run --------------------------------------------------------------------
# Incremental: keep whatever is already cached so a fix to one tag family does
# not cost a reprocess of the whole delivery.
ap <- file.path(CACHE, "fvilches_archives_v1_30min.rds")
gp <- file.path(CACHE, "fvilches_argos_v1.rds")
archives <- if (file.exists(ap)) readRDS(ap) else list()
argos    <- if (file.exists(gp)) readRDS(gp) else list()
mp <- "scratch/nes_calibration/fvilches_manifest.csv"
man <- if (file.exists(mp)) split(as.data.frame(fread(mp)), fread(mp)$topp) else list()
for (i in seq_len(nrow(tdr))) {
  tg <- tdr$topp[i]
  if (!is.null(archives[[tg]])) { message(sprintf("[%d/%d] %s cached", i, nrow(tdr), tg)); next }
  message(sprintf("[%d/%d] %s (%s)", i, nrow(tdr), tg, tdr$serial[i]))
  a_path <- arg[topp == tg]$path
  if (!length(a_path)) { message("   no Argos, skipped"); next }
  fx <- read_rawargos(a_path[1])
  if (is.null(fx) || nrow(fx) < 50) { message("   too few fixes, skipped"); next }
  fx <- speed_filter(fx)

  d <- read_archive(tdr$path[i])
  if (is.null(d) || nrow(d) < 1000) { message("   no usable light, skipped"); next }

  # Deployment position: there is no Type == "User" row in RawArgos, so take the
  # median of the good-quality fixes in the first 48 h, when the animal is still
  # at the colony. Recovery: the same at the other end.
  t0 <- min(fx$time); t1 <- max(fx$time)
  dep <- fx[time <= t0 + 48*3600 & class %in% c("3","2","1")]
  if (!nrow(dep)) dep <- fx[time <= t0 + 48*3600]
  rec <- fx[time >= t1 - 72*3600]
  archives[[tg]] <- list(main = d)
  argos[[tg]] <- as.data.frame(fx)
  man[[tg]] <- data.frame(
    topp = tg, season = tdr$season[i], serial = tdr$serial[i],
    ptt = arg[topp == tg]$ptt[1], source = "fvilches 2026-08-09",
    n_light = nrow(d), n_argos = nrow(fx),
    first_light = format(min(d$time)), last_light = format(max(d$time)),
    days = round(as.numeric(difftime(t1, t0, units = "days")), 1),
    deploy_lon = round(median(dep$lon), 4), deploy_lat = round(median(dep$lat), 4),
    recover_lon = round(median(rec$lon), 4), recover_lat = round(median(rec$lat), 4),
    med_err_m = round(median(fx$err_radius_m, na.rm = TRUE)),
    row.names = NULL)
}
M <- rbindlist(man, fill = TRUE)
cat("\n=== ingested ===\n"); print(as.data.frame(M), row.names = FALSE)

saveRDS(archives, file.path(CACHE, "fvilches_archives_v1_30min.rds"))
saveRDS(argos,    file.path(CACHE, "fvilches_argos_v1.rds"))
fwrite(M, "scratch/nes_calibration/fvilches_manifest.csv")
cat(sprintf("\nwrote %d deployments to %s\n", length(archives), CACHE))

cat("\n=== sanity: is the deployment position the colony? ===\n")
COL <- c(237.67, 37.11)
M[, km_from_colony := round(gc_km(deploy_lon, deploy_lat, COL[1], COL[2]))]
print(as.data.frame(M[, .(topp, season, deploy_lon, deploy_lat, km_from_colony,
                          days, n_light, med_err_m)]), row.names = FALSE)
