# ARM A vs ARM B, complete at 29 tags, on BOTH yardsticks.
#
# The pre-registered tests are per-knot accuracy, bias, coverage and the seasonal
# and latitude slopes. Those are reported first because they are what was written
# down in advance.
#
# The second block asks a question that was not pre-registered and should be
# labelled as such: does the depth gate transfer to HABITAT OCCUPANCY? The gate
# collapses the seasonal emission drift by 92% and does not transfer to
# kilometres, but kilometres is the wrong yardstick for what these tags are
# deployed to answer. Population occupancy overlap and, more to the point, the
# month x band overlap, are the metrics the ecology rests on. A change that does
# nothing for per-fix accuracy could still matter there, and the reverse.
suppressMessages({ library(data.table) })

D <- fread("scratch/nes_calibration/final29_tags.csv", colClasses = list(character = "id"))
K <- fread("scratch/nes_calibration/final29_knots.csv", colClasses = list(character = "id"))
A <- D[arm == "A_shipped"]; B <- D[arm == "B_gate"]
w <- merge(A[, .(id, kA = median_km, bA = bias_lat, rA = rmse_lat, sA = lat_sd, cA = cover_lat, nA = n)],
           B[, .(id, kB = median_km, bB = bias_lat, rB = rmse_lat, sB = lat_sd, cB = cover_lat, nB = n)],
           by = "id")
cat("=== PRE-REGISTERED: paired A_shipped vs B_gate, n =", nrow(w), "tags ===\n")
cat(sprintf("  %-22s A %8.0f   B %8.0f   delta %+8.0f   B better %2d/%2d\n",
            "median error (km)", median(w$kA), median(w$kB), median(w$kB) - median(w$kA),
            sum(w$kB < w$kA), nrow(w)))
f <- function(lab, a, b, better) cat(sprintf(
  "  %-22s A %8.3f   B %8.3f   delta %+8.3f   B better %2d/%2d\n",
  lab, mean(a), mean(b), mean(b) - mean(a), better, length(a)))
f("bias_lat (deg)", w$bA, w$bB, sum(abs(w$bB) < abs(w$bA)))
f("|bias_lat| (deg)", abs(w$bA), abs(w$bB), sum(abs(w$bB) < abs(w$bA)))
f("rmse_lat (deg)", w$rA, w$rB, sum(w$rB < w$rA))
f("posterior lat_sd (deg)", w$sA, w$sB, sum(w$sB < w$sA))
cat(sprintf("  %-22s A %8.3f   B %8.3f   delta %+8.3f   B better %2d/%2d\n",
            "coverage (lat)", mean(w$cA), mean(w$cB), mean(w$cB) - mean(w$cA),
            sum(w$cB > w$cA), nrow(w)))

cat("\n=== PRE-REGISTERED: per-knot structure (tags present in both arms) ===\n")
KK <- K[id %in% w$id & is.finite(err_lat)]
for (a in c("A_shipped", "B_gate")) {
  s <- KK[arm == a]
  sw <- s[, .(b = median(err_lat)), by = .(db = round(decl / 5) * 5)]
  cat(sprintf("  %-10s n=%6d | mean bias %+.3f | seasonal slope %+.5f (swing %.2f) | latitude slope %+.5f\n",
              a, nrow(s), mean(s$err_lat), coef(lm(err_lat ~ decl, data = s))[2],
              diff(range(sw$b)), coef(lm(err_lat ~ true_lat, data = s))[2]))
}

cat("\n=== NOT PRE-REGISTERED: habitat occupancy, the ecological yardstick ===\n")
KK[, est_lat := true_lat + err_lat]
KK[, month := as.integer(format(as.POSIXct(time, tz = "UTC"), "%m"))]
brk <- seq(10, 80, by = 5)
ov <- function(a, b) { pa <- prop.table(table(cut(a, brk))); pb <- prop.table(table(cut(b, brk)))
  sum(pmin(pa, pb)) }
ovmb <- function(d) { pa <- prop.table(table(cut(d$true_lat, brk), d$month))
  pb <- prop.table(table(cut(d$est_lat, brk), d$month)); sum(pmin(pa, pb)) }
for (a in c("A_shipped", "B_gate")) {
  s <- KK[arm == a]
  per <- s[, .(o = ov(true_lat, est_lat)), by = id]
  cat(sprintf("  %-10s pooled occupancy overlap %.3f | per-tag median %.3f | month x band %.3f | mean lat shift %+.2f deg\n",
              a, ov(s$true_lat, s$est_lat), median(per$o), ovmb(s),
              mean(s$est_lat) - mean(s$true_lat)))
}
cat("\n  overlap: 1.0 = occupancy distribution reproduced exactly.\n")
cat("  month x band is the 'when and where' metric; it is the one the ecology needs.\n")

cat("\n=== seasonal shift by month, both arms (n >= 500 months only) ===\n")
ms <- KK[, .(n = .N, shift = mean(err_lat)), by = .(arm, month)][n >= 500]
print(as.data.frame(dcast(ms, month ~ arm, value.var = "shift")[
  , lapply(.SD, function(x) if (is.numeric(x)) round(x, 2) else x)]), row.names = FALSE)
