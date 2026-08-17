# The deployment positions from the first ingest are wrong for six deployments,
# and the sanity check caught it: 2022042/043/046 place the animal at 357.2 E,
# 56.3 N -- Scotland, 8172 km from Ano Nuevo -- and 2022040/041 report trips of
# 613 days.
#
# CAUSE. The RawArgos files are named by TOPP id but keyed by PTT, and PTTs are
# REDEPLOYED across seasons just as the TDR tags are (214390 appears in both
# 2021023 and 2022040; 214397 in 2021031, 2022041 and 2023039). Each file
# therefore contains every fix that PTT ever produced, including other seasons'
# deployments and pre-deployment test transmissions from wherever the tag was
# programmed. Taking "the first 48 hours of the file" then lands on whatever
# came first in the tag's life, not on this deployment.
#
# FIX. The TDR archive IS deployment-specific -- one file per deployment -- so
# clip the Argos to the archive's own time span before deriving anything. The
# archives are already cached and correct, so only the fixes and the manifest
# need redoing, which is cheap.
suppressMessages(library(data.table))
Sys.setlocale("LC_TIME", "C")
ARG <- "fvilches/extracted/Argos raw"
CACHE <- "analysis/cache/nes"
COL <- c(237.67, 37.11)

lon360 <- function(x) (x %% 360 + 360) %% 360
gc_km <- function(l1, p1, l2, p2, R = 6371.0088) { d <- pi/180
  2*R*asin(pmin(1, sqrt(sin((p2-p1)*d/2)^2 + cos(p1*d)*cos(p2*d)*sin((l2-l1)*d/2)^2))) }

read_rawargos <- function(path) {
  a <- fread(path, showProgress = FALSE)
  need <- c("PTT", "Class", "PassDate", "PassTime", "Latitude", "Longitude")
  if (!all(need %in% names(a))) return(NULL)
  # The delivery mixes date formats: most files use "17-May-2023", at least one
  # uses "24-May-23" with a two-digit year. Assuming a single format silently
  # produced year 0023 and dropped a whole deployment. Rather than guess again,
  # try every plausible format and keep the one that parses the most rows to a
  # sane year -- and fail loudly if none does, instead of returning nothing.
  ts <- paste(a$PassDate, a$PassTime)
  fmts <- c("%d-%b-%Y %H:%M:%S", "%d-%b-%y %H:%M:%S",
            "%Y-%m-%d %H:%M:%S", "%y-%m-%d %H:%M:%S", "%d/%m/%Y %H:%M:%S")
  parsed <- lapply(fmts, function(f) as.POSIXct(ts, format = f, tz = "UTC"))
  score <- vapply(parsed, function(t) {
    y <- suppressWarnings(as.integer(format(t, "%Y")))
    mean(!is.na(y) & y >= 2000 & y <= 2100)
  }, 0)
  if (max(score) < 0.5)
    stop("no date format parses ", basename(path), ": first value '", ts[1], "'")
  a[, time := parsed[[which.max(score)]]]
  a <- a[!is.na(time) & !is.na(Latitude) & !is.na(Longitude)]
  if (!nrow(a)) return(NULL)
  rank <- c("3"=1,"2"=2,"1"=3,"0"=4,"A"=5,"B"=6,"Z"=9)
  a[, cr := rank[as.character(Class)]]; a[is.na(cr), cr := 9L]
  setorder(a, time, cr); a <- unique(a, by = "time")
  er <- if ("Error radius" %in% names(a)) a[["Error radius"]] else NA_real_
  data.table(time = a$time, ptt = as.character(a$PTT), class = as.character(a$Class),
             lon = lon360(a$Longitude), lat = a$Latitude,
             err_radius_m = as.numeric(er))[class %in% c("3","2","1","0","A","B")]
}
speed_filter <- function(d, vmax = 10) {
  setorder(d, time); if (nrow(d) < 3) return(d)
  keep <- rep(TRUE, nrow(d)); last <- 1L
  for (i in 2:nrow(d)) {
    dt <- as.numeric(difftime(d$time[i], d$time[last], units = "hours"))
    if (dt <= 0) { keep[i] <- FALSE; next }
    if (gc_km(d$lon[last], d$lat[last], d$lon[i], d$lat[i]) / dt > vmax) keep[i] <- FALSE
    else last <- i
  }
  d[keep]
}

archives <- readRDS(file.path(CACHE, "fvilches_archives_v1_30min.rds"))
arg <- data.table(path = list.files(ARG, pattern = "RawArgos[.]csv$", full.names = TRUE))
arg[, file := basename(path)]
arg[, topp := sub("^([0-9]+)_.*$", "\\1", file)]
arg[, ptt := sub("^[0-9]+_(.+)-RawArgos[.]csv$", "\\1", file)]

argos <- list(); man <- list()
for (tg in names(archives)) {
  d <- archives[[tg]]$main
  win <- range(d$time)
  p <- arg[topp == tg]$path
  if (!length(p)) next
  raw <- read_rawargos(p[1])
  if (is.null(raw)) next
  n_all <- nrow(raw)
  # clip to THIS deployment, using the archive's own span with a small margin
  fx <- raw[time >= win[1] - 3*86400 & time <= win[2] + 3*86400]
  if (nrow(fx) < 50) { message(tg, ": only ", nrow(fx), " fixes in window, skipped"); next }
  # Anchor the speed filter on a GOOD fix. A class-0 solution at the head of the
  # record (2023039 had one 1315 km out, duplicated three times) otherwise
  # becomes the reference every later fix is judged against.
  first_good <- which(fx$class %in% c("3", "2", "1"))[1]
  if (!is.na(first_good) && first_good > 1) fx <- fx[first_good:.N]
  fx <- speed_filter(fx)
  t0 <- min(fx$time); t1 <- max(fx$time)
  # Deployment position from good classes only, widening the window rather than
  # falling back to poor fixes: a wrong endpoint is worse than a late one.
  dep <- fx[time <= t0 + 48*3600 & class %in% c("3","2","1")]
  if (!nrow(dep)) dep <- fx[time <= t0 + 7*86400 & class %in% c("3","2","1","0")]
  if (!nrow(dep)) dep <- fx[time <= t0 + 48*3600]
  rec <- fx[time >= t1 - 72*3600]
  argos[[tg]] <- as.data.frame(fx)
  man[[tg]] <- data.frame(topp = tg, season = substr(tg, 1, 4),
    ptt = arg[topp == tg]$ptt[1], source = "fvilches 2026-08-09",
    n_argos_file = n_all, n_argos_used = nrow(fx),
    dropped_outside = n_all - nrow(fx),
    n_light = nrow(d),
    days = round(as.numeric(difftime(t1, t0, units = "days")), 1),
    deploy_lon = round(median(dep$lon), 4), deploy_lat = round(median(dep$lat), 4),
    recover_lon = round(median(rec$lon), 4), recover_lat = round(median(rec$lat), 4),
    med_err_m = round(median(fx$err_radius_m, na.rm = TRUE)), row.names = NULL)
}
M <- rbindlist(man)
M[, km_from_colony := round(gc_km(deploy_lon, deploy_lat, COL[1], COL[2]))]
M[, recover_km := round(gc_km(recover_lon, recover_lat, COL[1], COL[2]))]
setorder(M, season, topp)
cat("=== after clipping the Argos to each deployment's own archive span ===\n")
print(as.data.frame(M[, .(topp, season, ptt, n_argos_file, n_argos_used,
                          dropped_outside, days, km_from_colony, recover_km)]),
      row.names = FALSE)

bad <- M[km_from_colony > 100 | days > 320 | days < 60]
cat(sprintf("\ndeployments still implausible: %d\n", nrow(bad)))
if (nrow(bad)) print(as.data.frame(bad[, .(topp, days, km_from_colony)]), row.names = FALSE)

saveRDS(argos, file.path(CACHE, "fvilches_argos_v1.rds"))
fwrite(M, "scratch/nes_calibration/fvilches_manifest.csv")
cat(sprintf("\nrewrote Argos and manifest for %d deployments\n", length(argos)))
