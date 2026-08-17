# Inventory of the new delivery: what is actually NEW, and what is a repeat of
# the ten deployments already analysed.
#
# Files are named <TOPPID>_<tagserial>-out-Archive.csv and
# <TOPPID>_<PTT>-RawArgos.csv, so both identifiers are recoverable from the
# names. Three things decide what the delivery is worth:
#   1. how many TOPP IDs are new (a reprocessed deployment adds no information)
#   2. how many of those have BOTH light and Argos
#   3. whether the same physical tag was redeployed across years -- which would
#      turn "pooling across tags" into a two-level hierarchy, deployments nested
#      within tags, and would let the response's stability be tested directly
suppressMessages(library(data.table))
NEW <- "fvilches/extracted"

tdr <- data.table(file = list.files(file.path(NEW, "TDR raw"), pattern = "Archive[.]csv$"))
tdr[, topp := sub("^([0-9]+)_.*$", "\\1", file)]
tdr[, serial := sub("^[0-9]+_(.+)-out-Archive[.]csv$", "\\1", file)]
tdr[, season := substr(topp, 1, 4)]

arg <- data.table(file = list.files(file.path(NEW, "Argos raw"), pattern = "RawArgos[.]csv$"))
arg[, topp := sub("^([0-9]+)_.*$", "\\1", file)]
arg[, ptt := sub("^[0-9]+_(.+)-RawArgos[.]csv$", "\\1", file)]

depth_only <- list.files(file.path(NEW, "TDR raw"), pattern = "compressed_tdr[.]txt$")

# the ten already analysed
have <- names(readRDS("analysis/cache/nes/archives_v3_30min.rds"))

cat("=== files delivered ===\n")
cat(sprintf("TDR archives (light + depth): %d\n", nrow(tdr)))
cat(sprintf("depth-only .txt (no light, cannot be geolocated): %d  [%s]\n",
            length(depth_only), paste(sub("_.*", "", depth_only), collapse = ", ")))
cat(sprintf("Argos raw: %d\n\n", nrow(arg)))

both <- merge(tdr[, .(topp, serial, season)], arg[, .(topp, ptt)], by = "topp")
both[, already := topp %in% have]
cat("=== double-tagged deployments, by season ===\n")
print(as.data.frame(both[, .(deployments = .N,
  already_analysed = sum(already), new = sum(!already)), by = season][order(season)]),
  row.names = FALSE)

cat(sprintf("\nALREADY ANALYSED: %d   NEW: %d   TOTAL AVAILABLE: %d\n",
            sum(both$already), sum(!both$already),
            length(have) + sum(!both$already)))
cat(sprintf("(the ten in the current analysis, plus the new ones)\n"))

cat("\n=== is the same physical tag redeployed across years? ===\n")
reuse <- both[, .(deployments = .N, seasons = paste(sort(season), collapse = "+"),
                  topps = paste(sort(topp), collapse = " ")), by = serial][deployments > 1]
if (nrow(reuse)) {
  print(as.data.frame(reuse[order(-deployments)]), row.names = FALSE)
  cat(sprintf("\n%d tag serials were deployed more than once, covering %d deployments.\n",
              nrow(reuse), sum(reuse$deployments)))
  cat("That makes deployments nested within physical tags: the response can be\n")
  cat("tested for stability across years on the SAME sensor, which is the\n")
  cat("assumption the pooling rests on.\n")
} else cat("no serial reuse\n")

cat("\n=== tag families present (pooling must be within model) ===\n")
both[, family := fifelse(grepl("^219", serial), "Mk9 219xxxx",
                  fifelse(grepl("^18A", serial), "18Axxxx",
                   fifelse(grepl("^19A", serial), "19Axxxx", "other")))]
print(as.data.frame(both[, .(deployments = .N, serials = uniqueN(serial)),
                         by = family][order(-deployments)]), row.names = FALSE)

cat("\n=== deployments with Argos but no light archive ===\n")
cat(paste(setdiff(arg$topp, tdr$topp), collapse = ", "), "\n")
cat("=== light archive but no Argos ===\n")
cat(paste(setdiff(tdr$topp, arg$topp), collapse = ", "), "\n")
