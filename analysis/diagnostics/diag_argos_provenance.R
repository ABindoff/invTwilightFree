# Re-check Argos provenance using the EXACT assembly the production scripts use.
#
# My first duplication diagnostic was wrong: it built one lookup table for all 29
# deployments and let a PTT lookup into the 2021 delivery overwrite the fvilches
# tags' own Argos. That manufactured 8 groups of "identical" series and 9 negative
# overlaps. fit_rescore29.R does not do that -- it takes 2021 tags from
# `ar$fixes[Ptt == ...]` and fvilches tags from `new_argos[[topp]]`, which are
# different sources. Redo it that way and report what is actually true.
SCRATCH <- "analysis/diagnostics"   # repo-relative; run from the package root
source(file.path(SCRATCH, "nes_common.R"))

old <- lapply(readRDS("analysis/cache/nes/archives_v3_30min.rds"), `[[`, "main")
new <- lapply(readRDS("analysis/cache/nes/fvilches_archives_v1_30min.rds"), `[[`, "main")
new_argos <- readRDS("analysis/cache/nes/fvilches_argos_v1.rds")
meta <- nes_meta(); ar <- nes_argos()

REC <- list()
for (tg in names(old)) {                       # 2021 delivery: PTT lookup
  m <- meta[match(tg, meta$id), ]
  if (is.na(m$ptt) || !nzchar(m$ptt)) next
  a <- ar$fixes[Ptt == m$ptt]; if (nrow(a) < 50) next
  REC[[tg]] <- list(src = "2021/ptt", a = a[order(time)], l = as.data.table(old[[tg]]))
}
for (tg in names(new)) {                       # fvilches: keyed by topp
  a <- as.data.table(new_argos[[tg]]); if (is.null(a) || !nrow(a)) next
  REC[[tg]] <- list(src = "fvilches", a = a[order(time)], l = as.data.table(new[[tg]]))
}
cat(sprintf("%d deployments assembled the way production does\n\n", length(REC)))

cat("=== does each deployment's Argos overlap its OWN light record? ===\n")
bad <- character(0)
for (tg in names(REC)) {
  r <- REC[[tg]]
  lt <- range(r$l$time); at <- range(r$a$time)
  ov <- as.numeric(difftime(min(lt[2], at[2]), max(lt[1], at[1]), units = "days"))
  sl <- as.numeric(diff(lt), units = "days")
  if (ov < 0.5 * sl) {
    bad <- c(bad, tg)
    cat(sprintf("  %-9s [%s] light %s..%s | Argos %s..%s | overlap %.0f d\n", tg, r$src,
                format(lt[1], "%Y-%m-%d"), format(lt[2], "%Y-%m-%d"),
                format(at[1], "%Y-%m-%d"), format(at[2], "%Y-%m-%d"), ov))
  }
}
cat(sprintf("  %d of %d with poor overlap%s\n", length(bad), length(REC),
            if (!length(bad)) " -- all deployments are matched to their own Argos" else ""))

cat("\n=== are any two deployments sharing an identical series? ===\n")
key <- vapply(REC, function(r) paste(nrow(r$a), round(sum(r$a$lat), 4),
                                     as.numeric(min(r$a$time)), as.numeric(max(r$a$time)),
                                     sep = "|"), "")
tb <- table(key); sh <- names(tb)[tb > 1]
if (length(sh)) for (k in sh)
  cat(sprintf("  IDENTICAL: %s\n", paste(names(key)[key == k], collapse = " = ")))
cat(sprintf("  %d deployments, %d distinct series\n", length(REC), length(unique(key))))

cat("\n=== speed, by delivery ===\n")
sp <- function(a) {
  a <- unique(a[order(time)], by = "time"); n <- nrow(a); if (n < 3) return(NULL)
  d <- gc_km(lon360(a$lon[-n]), a$lat[-n], lon360(a$lon[-1]), a$lat[-1])
  dt <- as.numeric(difftime(a$time[-1], a$time[-n], units = "hours"))
  v <- d[dt > 0] / dt[dt > 0] * 1000 / 3600; v[is.finite(v)]
}
S <- rbindlist(lapply(names(REC), function(tg) {
  v <- sp(REC[[tg]]$a); if (is.null(v)) return(NULL)
  data.table(id = tg, src = REC[[tg]]$src, n = length(v), v_med = median(v),
             v_max = max(v), f2 = mean(v > 2), f5 = mean(v > 5))
}))
for (s in unique(S$src)) {
  x <- S[src == s]
  cat(sprintf("  %-9s n=%2d tags | median speed %.2f | max %.0f | steps >2 m/s %.1f%% | >5 m/s %.1f%%\n",
              s, nrow(x), median(x$v_med), max(x$v_max), 100*mean(x$f2), 100*mean(x$f5)))
}
cat("\n  per-tag maxima:\n")
for (i in seq_len(nrow(S)))
  cat(sprintf("    %-9s [%-8s] v_max %8.0f m/s | >5 m/s on %5.1f%% of steps\n",
              S$id[i], S$src[i], S$v_max[i], 100 * S$f5[i]))
