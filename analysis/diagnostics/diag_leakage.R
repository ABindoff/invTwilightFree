# Audit every place the ANALYSIS could be using Argos to help itself, rather
# than only to score itself. Three candidates:
#
#   1. the grid extent          -- built from the Argos envelope (confirmed)
#   2. the recovery endpoint    -- median of the last 72 h of Argos fixes
#   3. the hierarchical mesh    -- deploy point +/- a hand-chosen pad
#
# 2 is only leakage if the recovery position is not simply the colony. These
# animals are recaptured at the colony to get the tag back, so a GLS-only study
# would know it; but that has to be checked, not assumed.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))

meta <- nes_meta(); ar <- nes_argos()
arch <- lapply(readRDS("analysis/cache/nes/archives_v3_30min.rds"), `[[`, "main")

rows <- list()
for (id in names(arch)) {
  m <- meta[match(id, meta$id), ]
  if (is.na(m$ptt) || !nzchar(m$ptt)) next
  df <- ar$deploy[Ptt == m$ptt]; a <- ar$fixes[Ptt == m$ptt]
  if (!nrow(df) || nrow(a) < 50) next
  tf <- a[time >= max(time) - 72*3600]
  rec <- c(median(tf$lon), median(tf$lat))
  rows[[length(rows)+1]] <- data.frame(
    id = id,
    deploy_lon = round(df$deploy_lon, 2), deploy_lat = round(df$deploy_lat, 2),
    recover_lon = round(rec[1], 2), recover_lat = round(rec[2], 2),
    recover_km_from_deploy = round(gc_km(rec[1], rec[2], df$deploy_lon, df$deploy_lat)),
    n_tail_fix = nrow(tf), row.names = NULL)
}
r <- do.call(rbind, rows)
cat("=== is the recovery endpoint just the colony? ===\n")
print(r, row.names = FALSE)
cat(sprintf("\nmedian distance of the recovery fix from the deployment fix: %.0f km\n",
            median(r$recover_km_from_deploy)))
cat("If that is small, pinning the end at the colony is legitimate prior\n")
cat("knowledge (the animal is recaptured there) rather than borrowed truth.\n")

cat("\n=== how much of the Argos envelope does the mesh pad cover? ===\n")
all_lon <- unlist(lapply(names(arch), function(id) {
  m <- meta[match(id, meta$id), ]
  if (is.na(m$ptt) || !nzchar(m$ptt)) return(NULL)
  ar$fixes[Ptt == m$ptt]$lon }))
all_lat <- unlist(lapply(names(arch), function(id) {
  m <- meta[match(id, meta$id), ]
  if (is.na(m$ptt) || !nzchar(m$ptt)) return(NULL)
  ar$fixes[Ptt == m$ptt]$lat }))
COL <- c(median(r$deploy_lon), median(r$deploy_lat))
cat(sprintf("colony %.2f E %.2f N\n", COL[1], COL[2]))
cat(sprintf("Argos envelope reaches %.1f deg west and %.1f deg north of it\n",
            COL[1] - min(all_lon), max(all_lat) - COL[2]))
cat(sprintf("mesh_pad_lon = 74, mesh_pad_lat = 26 -> covers to %.1f E, %.1f N\n",
            COL[1] - 74, COL[2] + 26))
