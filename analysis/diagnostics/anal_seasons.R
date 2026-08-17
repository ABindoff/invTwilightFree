# What seasons do these deployments actually cover, and does season explain the two
# lambda clusters?
#
# Northern elephant seals make TWO distinct annual migrations: a short post-breeding
# trip (roughly Feb-May) and a long post-moult trip (roughly May-Jan). They differ in
# duration, latitude range and -- critically for a light-based method -- in which part
# of the solar year they sample. A post-breeding trip sits either side of the March
# equinox; a post-moult trip spans the June solstice, the September equinox and often
# the December solstice.
#
# Per-logger optimal lambda splits into two clusters wanting OPPOSITE corrections
# (0.07-0.20 vs 1.8-18.4). If those line up with trip type, the "two branches" of the
# non-monotone bias curve are a seasonal phenomenon rather than a sensor one, and the
# per-tag problem becomes a per-TRIP problem, which is knowable without Argos.
#
# Aggregates and date ranges only. No per-observation rows printed.
suppressMessages(library(data.table))
# run from the package root (setwd removed for portability)

old <- lapply(readRDS("analysis/cache/nes/archives_v3_30min.rds"), `[[`, "main")
new <- lapply(readRDS("analysis/cache/nes/fvilches_archives_v1_30min.rds"), `[[`, "main")
arch <- c(old, new)

a1 <- data.table(file = list.files("data/nes_untracked", pattern = "^[0-9]+_.*\\.csv$"))
a1[, `:=`(topp = sub("^([0-9]+)_.*$", "\\1", file), serial = sub("^[0-9]+_(.+)\\.csv$", "\\1", file))]
b1 <- data.table(file = list.files("fvilches/extracted/TDR raw", pattern = "Archive[.]csv$"))
b1[, `:=`(topp = sub("^([0-9]+)_.*$", "\\1", file),
          serial = sub("^[0-9]+_(.+)-out-Archive[.]csv$", "\\1", file))]
SER <- unique(rbind(a1[, .(topp, serial)], b1[, .(topp, serial)]), by = c("topp", "serial"))

S <- rbindlist(lapply(names(arch), function(tg) {
  tm <- arch[[tg]]$time
  tm <- tm[is.finite(as.numeric(tm))]
  if (!length(tm)) return(NULL)
  data.table(id = tg, start = min(tm), end = max(tm),
             days = as.numeric(difftime(max(tm), min(tm), units = "days")))
}))
S <- merge(S, SER, by.x = "id", by.y = "topp", all.x = TRUE)
S[, `:=`(yr = format(start, "%Y"), m0 = as.integer(format(start, "%m")),
         m1 = as.integer(format(end, "%m")))]
# post-breeding trips are short and start Feb-Apr; post-moult are long and start May-Aug
S[, trip := fifelse(days < 120, "post-breeding (short)", "post-moult (long)")]

cat("=== temporal coverage ===\n")
cat(sprintf("deployments %d | overall span %s to %s\n", nrow(S),
            format(min(S$start), "%Y-%m-%d"), format(max(S$end), "%Y-%m-%d")))
cat(sprintf("deployment years: %s\n",
            paste(sprintf("%s (n=%d)", names(table(S$yr)), as.integer(table(S$yr))), collapse = ", ")))
cat(sprintf("duration days: median %.0f | range %.0f-%.0f\n",
            median(S$days), min(S$days), max(S$days)))
cat(sprintf("start months: %s\n",
            paste(sprintf("%s:%d", month.abb[as.integer(names(table(S$m0)))],
                          as.integer(table(S$m0))), collapse = " ")))
cat(sprintf("end months:   %s\n\n",
            paste(sprintf("%s:%d", month.abb[as.integer(names(table(S$m1)))],
                          as.integer(table(S$m1))), collapse = " ")))

cat("=== trip type ===\n")
print(as.data.frame(S[, .(n = .N, med_days = round(median(days)),
                          starts = paste(sort(unique(month.abb[m0])), collapse = "/"),
                          ends = paste(sort(unique(month.abb[m1])), collapse = "/")),
                      by = trip]), row.names = FALSE)

# how many calendar months does each deployment cover, and do they see a solstice?
S[, months_covered := mapply(function(a, b)
    length(unique(format(seq(a, b, by = "week"), "%Y-%m"))), start, end)]
cat(sprintf("\ncalendar months covered per deployment: median %.0f, range %d-%d\n",
            median(S$months_covered), min(S$months_covered), max(S$months_covered)))

# ---- does trip type explain the lambda clusters? ---------------------------
f <- "scratch/nes_calibration/bias_lambda_wide.csv"
if (file.exists(f)) {
  W <- fread(f)[band == "all" & !is.na(serial)]
  xing <- function(lam, y) {
    x <- log2(lam); o <- order(x); x <- x[o]; y <- y[o]
    k <- which(y[-1] * y[-length(y)] < 0); if (!length(k)) return(NA_real_)
    i <- k[1]; 2^(x[i] - y[i] * (x[i+1] - x[i]) / (y[i+1] - y[i]))
  }
  # per DEPLOYMENT this time, so it can be joined to that deployment's season
  pd <- W[, .(bias = mean(bias, na.rm = TRUE)), by = .(id, lam_mult)]
  lam <- pd[, .(lam_opt = xing(lam_mult, bias)), by = id]
  lam[, id := as.character(id)]; S[, id := as.character(id)]
  M <- merge(S, lam, by = "id")
  M <- M[is.finite(lam_opt)]
  cat(sprintf("\n=== lambda cluster vs trip type (%d deployments with a crossing) ===\n", nrow(M)))
  M[, cluster := fifelse(lam_opt < 1, "loose (<1x)", "tight (>1x)")]
  print(as.data.frame(dcast(M[, .N, by = .(trip, cluster)], trip ~ cluster,
                            value.var = "N", fill = 0)), row.names = FALSE)
  cat(sprintf("\n  median lam_opt: post-breeding %.2f | post-moult %.2f\n",
              median(M[trip == "post-breeding (short)"]$lam_opt),
              median(M[trip == "post-moult (long)"]$lam_opt)))
  if (uniqueN(M$trip) > 1 && nrow(M) >= 6) {
    p <- suppressWarnings(wilcox.test(log2(lam_opt) ~ trip, data = M)$p.value)
    cat(sprintf("  Wilcoxon on log2(lam_opt) by trip type: p = %.4f\n", p))
  }
  ct <- suppressWarnings(cor.test(M$days, log2(M$lam_opt)))
  cat(sprintf("  cor(deployment length, log2 lam_opt) = %+.3f, p = %.3f\n",
              ct$estimate, ct$p.value))
  ct2 <- suppressWarnings(cor.test(M$m0, log2(M$lam_opt)))
  cat(sprintf("  cor(start month,      log2 lam_opt) = %+.3f, p = %.3f\n",
              ct2$estimate, ct2$p.value))
}
