# Shared scaffolding for the calibration diagnostics: the crosswalk, the Argos
# reader and the scoring, lifted from the analysis so the tests run against the
# same data the report does.
suppressMessages({ library(data.table); library(invTwilightFree); library(terra) })
Sys.setlocale("LC_TIME", "C")

DATA_DIR <- "data/nes_untracked"
lon360 <- function(x) (x %% 360 + 360) %% 360
dlon <- function(a, b) ((a - b + 180) %% 360) - 180
gc_km <- function(l1, p1, l2, p2, R = 6371.0088) { d <- pi/180
  2*R*asin(pmin(1, sqrt(sin((p2-p1)*d/2)^2 + cos(p1*d)*cos(p2*d)*sin((l2-l1)*d/2)^2))) }

read_xlsx_grid <- function(path, sheet = 1L) {
  ns <- c(d = "http://schemas.openxmlformats.org/spreadsheetml/2006/main")
  tmp <- file.path(tempdir(), "xlsx"); unlink(tmp, recursive = TRUE)
  utils::unzip(path, exdir = tmp)
  ss <- character(0)
  sfile <- file.path(tmp, "xl", "sharedStrings.xml")
  if (file.exists(sfile))
    ss <- vapply(xml2::xml_find_all(xml2::read_xml(sfile), "d:si", ns),
                 function(si) paste0(xml2::xml_text(xml2::xml_find_all(si, ".//d:t", ns)),
                                     collapse = ""), "")
  ws <- xml2::read_xml(file.path(tmp, "xl", "worksheets", sprintf("sheet%d.xml", sheet)))
  rows <- xml2::xml_find_all(ws, "//d:row", ns)
  cell_col <- function(ref) {
    L <- utf8ToInt(toupper(gsub("[0-9]", "", ref))) - 64L
    sum(L * 26^rev(seq_along(L) - 1))
  }
  parsed <- Filter(Negate(is.null), lapply(rows, function(r) {
    cs <- xml2::xml_find_all(r, "d:c", ns)
    if (!length(cs)) return(NULL)
    list(idx = vapply(xml2::xml_attr(cs, "r"), cell_col, 1),
         val = vapply(seq_along(cs), function(i) {
           cc <- cs[[i]]; ty <- xml2::xml_attr(cc, "t")
           v <- xml2::xml_find_first(cc, "d:v", ns)
           if (is.na(v)) {
             is_el <- xml2::xml_find_first(cc, "d:is", ns)
             if (is.na(is_el)) "" else xml2::xml_text(is_el)
           } else if (identical(ty, "s")) ss[as.integer(xml2::xml_text(v)) + 1L]
           else xml2::xml_text(v)
         }, ""))
  }))
  nc <- max(unlist(lapply(parsed, function(p) p$idx)))
  out <- matrix("", nrow = length(parsed), ncol = nc)
  for (i in seq_along(parsed)) out[i, parsed[[i]]$idx] <- parsed[[i]]$val
  out
}

nes_meta <- function() {
  xlsx <- list.files(DATA_DIR, pattern = "\\.xlsx$", full.names = TRUE)[1]
  g <- read_xlsx_grid(xlsx)
  tbl <- as.data.frame(g[-(1:2), , drop = FALSE], stringsAsFactors = FALSE)
  names(tbl) <- make.names(g[2, ], unique = TRUE)
  m <- data.frame(id = tbl$TOPP.ID, ptt = tbl$Sat.Tag.PTT,
                  deploy_date = as.Date(suppressWarnings(as.numeric(tbl$Date)),
                                        origin = "1899-12-30"),
                  stringsAsFactors = FALSE)
  m[!is.na(m$id) & nzchar(m$id), ]
}

nes_argos <- function() {
  f <- list.files(file.path(DATA_DIR, "argos"), pattern = "-Locations[.]csv$",
                  recursive = TRUE, full.names = TRUE)
  A <- rbindlist(lapply(f, fread, showProgress = FALSE,
                        colClasses = list(character = c("Ptt", "Quality"))), fill = TRUE)
  A[, time := as.POSIXct(Date, format = "%H:%M:%S %d-%b-%Y", tz = "UTC")]
  A[, `:=`(lon = lon360(Longitude), lat = Latitude)]
  dep <- A[Type == "User", .(deploy_lon = median(lon), deploy_lat = median(lat),
                             user_date = min(as.Date(time))), by = Ptt]
  fixes <- A[Type == "Argos" & Quality %in% c("3","2","1","0","A","B") &
               !is.na(time) & !is.na(lon)]
  list(deploy = dep, fixes = fixes)
}

# Argos position interpolated to arbitrary times; gaps over 24 h are left NA
# rather than bridged, so a scoring comparison never leans on invented truth.
argos_at <- function(a, tm) {
  setorder(a, time); at <- as.numeric(a$time); qt <- as.numeric(tm)
  j <- findInterval(qt, at); ok <- j >= 1 & j < length(at)
  lo <- la <- rep(NA_real_, length(qt)); jj <- j[ok]
  good <- (at[jj+1] - at[jj]) / 3600 <= 24
  w <- (qt[ok] - at[jj]) / (at[jj+1] - at[jj]); idx <- which(ok)[good]
  lo[idx] <- lon360(a$lon[jj[good]] + w[good]*dlon(a$lon[jj[good]+1], a$lon[jj[good]]))
  la[idx] <- a$lat[jj[good]] + w[good]*(a$lat[jj[good]+1] - a$lat[jj[good]])
  list(lon = lo, lat = la)
}

# Departure from the dive record alone: a hauled-out seal does not dive, so the
# first day after the animal is last ashore is when it left.
departure_time <- function(d, dive_m = 50) {
  dd <- as.data.table(d)[, .b := as.integer(floor(as.numeric(time)/86400))]
  s <- dd[, .(deep = mean(depth_max > dive_m, na.rm = TRUE)), by = .b][order(.b)]
  at_sea <- s$deep > 0.5
  if (!any(at_sea)) return(NA_real_)
  last_ashore <- if (any(!at_sea)) max(which(!at_sea)) else 0L
  i <- if (last_ashore >= length(at_sea)) which(at_sea)[1] else last_ashore + 1L
  s$.b[i] * 86400
}

# Response geometry from `rg`, light scale from `rs`.
graft_response <- function(rg, rs) {
  if (is.null(rg) || is.null(rs)) return(NULL)
  slope <- rs$max_light / (4 * rg$scale)
  rg$calibration <- c(slope * (rg$z50 + 2*rg$scale), slope)
  rg$baseline <- rs$baseline; rg$max_light <- rs$max_light
  rg$slope <- slope; rg$zero_at <- rg$z50 + 2*rg$scale
  rg
}
