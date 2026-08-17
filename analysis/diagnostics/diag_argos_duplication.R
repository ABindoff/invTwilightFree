# Are Argos series being shared between deployments?
#
# The speed diagnostic showed several TOPP IDs with byte-identical fix counts, median
# speeds and maxima -- across DIFFERENT YEARS. If two deployments are being handed the
# same Argos series then at least one of them has the wrong ground truth, and every
# error, bias and coverage figure computed for it is meaningless.
#
# The likely mechanism is the PTT lookup: `ar$fixes[Ptt == m$ptt]` returns every fix
# ever transmitted by that satellite tag. If a PTT was redeployed on another animal in
# a later season, that query returns BOTH deployments' fixes, and the deployment is
# then time-windowed against a track that includes another animal entirely.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))
new_argos <- readRDS("analysis/cache/nes/fvilches_argos_v1.rds")
ar <- nes_argos(); meta <- nes_meta()

cat("=== does one PTT serve more than one TOPP deployment? ===\n")
m <- as.data.table(meta)[!is.na(ptt) & nzchar(ptt)]
dup <- m[, .N, by = ptt][N > 1]
cat(sprintf("  %d PTTs mapped to %d deployments; %d PTTs used more than once\n",
            uniqueN(m$ptt), nrow(m), nrow(dup)))
if (nrow(dup)) for (i in seq_len(nrow(dup)))
  cat(sprintf("    PTT %-8s -> %s\n", dup$ptt[i],
              paste(m[ptt == dup$ptt[i]]$id, collapse = ", ")))

cat("\n=== per PTT: does the fix series span more than one deployment season? ===\n")
for (p in unique(m$ptt)) {
  f <- ar$fixes[Ptt == p]
  if (nrow(f) < 50) next
  yrs <- sort(unique(format(f$time, "%Y")))
  gap <- if (nrow(f) > 2) max(as.numeric(diff(sort(f$time)), units = "days")) else 0
  if (length(yrs) > 1 || gap > 60)
    cat(sprintf("  PTT %-8s (%s): %d fixes, %s to %s, years %s, largest gap %.0f d\n",
                p, paste(m[ptt == p]$id, collapse = "/"), nrow(f),
                format(min(f$time), "%Y-%m-%d"), format(max(f$time), "%Y-%m-%d"),
                paste(yrs, collapse = "+"), gap))
}

cat("\n=== are any two deployments' assembled series literally identical? ===\n")
A <- list()
for (tg in names(new_argos)) A[[tg]] <- as.data.table(new_argos[[tg]])[order(time)]
for (i in seq_len(nrow(m))) {
  f <- ar$fixes[Ptt == m$ptt[i]][order(time)]
  if (nrow(f) >= 50) A[[m$id[i]]] <- f
}
key <- vapply(A, function(a) paste(nrow(a), round(sum(a$lat), 4), round(sum(lon360(a$lon)), 4),
                                   as.numeric(min(a$time)), as.numeric(max(a$time)), sep = "|"), "")
tb <- table(key); shared <- names(tb)[tb > 1]
cat(sprintf("  %d deployments, %d distinct Argos series\n", length(A), length(unique(key))))
if (length(shared)) {
  for (k in shared) {
    ids <- names(key)[key == k]
    a <- A[[ids[1]]]
    cat(sprintf("    IDENTICAL: %s  (%d fixes, %s to %s)\n", paste(ids, collapse = " = "),
                nrow(a), format(min(a$time), "%Y-%m-%d"), format(max(a$time), "%Y-%m-%d")))
  }
} else cat("  none identical\n")

cat("\n=== does each deployment's Argos actually overlap its own light record? ===\n")
old <- lapply(readRDS("analysis/cache/nes/archives_v3_30min.rds"), `[[`, "main")
new <- lapply(readRDS("analysis/cache/nes/fvilches_archives_v1_30min.rds"), `[[`, "main")
arch <- c(old, new)
bad <- 0
for (tg in names(A)) {
  if (is.null(arch[[tg]])) next
  lt <- range(arch[[tg]]$time); at <- range(A[[tg]]$time)
  ov <- as.numeric(difftime(min(lt[2], at[2]), max(lt[1], at[1]), units = "days"))
  span_l <- as.numeric(diff(lt), units = "days"); span_a <- as.numeric(diff(at), units = "days")
  if (ov < 0.5 * span_l || span_a > 1.8 * span_l) {
    bad <- bad + 1
    cat(sprintf("  %-9s light %.0f d, Argos %.0f d, overlap %.0f d %s\n", tg, span_l, span_a, ov,
                if (span_a > 1.8 * span_l) "<- Argos spans far longer than the light record" else ""))
  }
}
cat(sprintf("  %d of %d deployments look mismatched\n", bad, length(A)))
