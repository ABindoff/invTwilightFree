# =============================================================================
# REPORTING -- turns the run's output into the tables the manuscript needs.
#
#   Rscript analysis/study/report_study.R          # the study
#   STUDY_DIR=analysis/output/study/v1_smoke Rscript analysis/study/report_study.R
#
# Every rule here is fixed by PREREGISTRATION.md section 3 and none of them may
# be chosen after seeing the numbers:
#
#   * every within-panel contrast is PAIRED by deployment and tested with a
#     Wilcoxon signed-rank test; marginal summaries are reported beside the test,
#     never in place of it
#   * when the mean and the median of per-tag values disagree in DIRECTION, both
#     are printed and the disagreement is flagged
#   * the primary significance statement clusters at the LOGGER (15 of them for
#     29 deployments); deployment-level figures are printed beside it and are the
#     optimistic ones
#   * `frac_scored` travels with every error
#
# Writes the tables to CSV as well as printing them, so the manuscript quotes a
# file rather than a screen.
# =============================================================================

suppressMessages(library(data.table))
source("analysis/study/study_config.R")

STUDY_DIR <- Sys.getenv("STUDY_DIR", OUT_DIR)
if (!file.exists(file.path(STUDY_DIR, "study_tags.csv")))
  stop("no study output in ", STUDY_DIR, " -- run run_study.R first")

TAGS <- fread(file.path(STUDY_DIR, "study_tags.csv"), colClasses = list(character = c("id", "serial")))
TIME <- fread(file.path(STUDY_DIR, "study_timing.csv"), colClasses = list(character = "id"))
PROV <- readRDS(file.path(STUDY_DIR, "study_provenance.rds"))
REP  <- file.path(STUDY_DIR, "report")
dir.create(REP, recursive = TRUE, showWarnings = FALSE)

rule <- function(title) cat(sprintf("\n%s\n%s\n", title, strrep("=", nchar(title))))
put  <- function(dt, name) { fwrite(dt, file.path(REP, paste0(name, ".csv"))); print(as.data.frame(dt), row.names = FALSE) }

# =============================================================================
rule("0. PROVENANCE")
cat(sprintf("study %s | started %s | finished %s (%.2f h)\n",
            PROV$study_version, format(PROV$started), format(PROV$finished),
            as.numeric(difftime(PROV$finished, PROV$started, units = "hours"))))
cat(sprintf("git %s | %s\n", substr(PROV$git[1], 1, 8), PROV$session$R.version$version.string))
cat(sprintf("machine probe: %.2f s at start, %.2f s at end (drift %+.0f%%) on %d cores\n",
            PROV$machine$start$median, PROV$machine$end$median,
            100 * (PROV$machine$end$median / PROV$machine$start$median - 1),
            PROV$machine$start$cores))
cat("The engine is single-threaded; the probe is the run's own certificate that\n")
cat("the machine stayed quiescent. A large drift invalidates the timing table only.\n")

# =============================================================================
rule("1. THE PANEL")
P <- as.data.table(PROV$panel)
put(P[, .(deployments = .N, loggers = uniqueN(serial),
          days_median = round(median(days)), days_min = round(min(days)),
          days_max = round(max(days)),
          tag_days = round(sum(days))), by = season][order(season)], "panel_by_season")
cat(sprintf("\ntotal: %d deployments, %d distinct loggers, %.0f tag-days\n",
            nrow(P), uniqueN(P$serial), sum(P$days)))
reuse <- P[, .N, by = serial][N > 1]
cat(sprintf("loggers deployed more than once: %d (covering %d deployments)\n",
            nrow(reuse), sum(reuse$N)))
cat("Deployments are NOT independent at the logger level, and every claim in this\n")
cat("study is about the sensor. Clustered figures below are the defensible ones.\n")

# =============================================================================
rule("2. PRIMARY ENDPOINT -- grid_light, distribution across deployments")
A <- TAGS[arm == "grid_light" & calibrated == TRUE]
put(A[, .(n = .N,
          median_km = round(median(median_km)), IQR_km = round(IQR(median_km)),
          min_km = round(min(median_km)), max_km = round(max(median_km)),
          bias_lat = round(median(bias_lat), 2),
          rmse_lat = round(median(rmse_lat), 2),
          cover_lat = round(median(cover_lat), 2),
          frac_scored = round(median(frac_scored), 3))], "primary_endpoint")
cat("\nReport the distribution, not the headline alone. Both endpoints are pinned at\n")
cat("one colony, which suppresses the sloppy latitude direction about twentyfold,\n")
cat("so these errors are a floor rather than a typical expectation.\n")

cat("\nscored fraction, which travels with the error:\n")
put(TAGS[arm == "grid_light", .(min = round(min(frac_scored), 3),
                                med = round(median(frac_scored), 3),
                                max = round(max(frac_scored), 3),
                                n_under_half = sum(frac_scored < 0.5)), by = season][order(season)],
    "scored_fraction")

# =============================================================================
rule("3. PAIRED CONTRASTS")
ser <- P[, .(id, serial, season)]

paired <- function(arm_a, arm_b, what = "median_km") {
  x <- merge(TAGS[arm == arm_a & calibrated == TRUE, c("id", what), with = FALSE],
             TAGS[arm == arm_b & calibrated == TRUE, c("id", what), with = FALSE],
             by = "id", suffixes = c("_a", "_b"))
  if (!nrow(x)) return(NULL)
  a <- x[[paste0(what, "_a")]]; b <- x[[paste0(what, "_b")]]
  # deployment level
  w1 <- suppressWarnings(wilcox.test(a, b, paired = TRUE)$p.value)
  # logger level: average within logger first, so a thrice-used logger counts once
  y <- merge(x, ser[, .(id, serial)], by = "id")
  z <- y[, .(a = mean(get(paste0(what, "_a"))), b = mean(get(paste0(what, "_b")))), by = serial]
  w2 <- suppressWarnings(wilcox.test(z$a, z$b, paired = TRUE)$p.value)
  d_mean <- mean(a) - mean(b); d_med <- median(a) - median(b)
  data.table(contrast = paste(arm_a, "-", arm_b), metric = what,
             n_deploy = nrow(x), n_logger = nrow(z),
             mean_a = round(mean(a), 3), mean_b = round(mean(b), 3),
             median_a = round(median(a), 3), median_b = round(median(b), 3),
             paired_median_diff = round(median(a - b), 3),
             better_a = sum(a < b), p_deployment = signif(w1, 3), p_logger = signif(w2, 3),
             summaries_disagree = sign(d_mean) != sign(d_med))
}

CONTRASTS <- rbindlist(Filter(Negate(is.null), c(
  lapply(c("median_km", "rmse_lat", "bias_lat", "cover_lat"),
         function(m) paired("grid_light", "grid_cal_untrunc", m)),
  lapply(c("median_km", "rmse_lat", "cover_lat"),
         function(m) paired("grid_cal_untrunc", "grid_cal_pertag", m)),
  lapply(c("median_km", "rmse_lat", "cover_lat"),
         function(m) paired("grid_fusion", "grid_light", m)),
  lapply(c("median_km", "cover_lat"),
         function(m) paired("hier_light", "grid_light", m)),
  lapply(c("median_km", "cover_lat"),
         function(m) paired("ffbs_light", "grid_light", m)),
  lapply(c("median_km", "cover_lat"),
         function(m) paired("ffbs_fusion", "ffbs_light", m)))))
put(CONTRASTS, "paired_contrasts")
if (any(CONTRASTS$summaries_disagree))
  cat("\n*** mean and median disagree in DIRECTION on the rows flagged above.\n",
      "    Report both, as the campaign's 242-vs-276 km reversal required.\n")

rule("3b. THE CALIBRATION DECOMPOSITION")
cat("grid_light -> grid_cal_untrunc isolates the departure truncation;\n")
cat("grid_cal_untrunc -> grid_cal_pertag isolates pooling within tag family.\n\n")
put(TAGS[arm %in% c("grid_light", "grid_cal_untrunc", "grid_cal_pertag"),
         .(n_calibrated = sum(calibrated), n_failed = sum(!calibrated),
           median_km = round(median(median_km, na.rm = TRUE)),
           mean_km = round(mean(median_km, na.rm = TRUE)),
           bias_lat = round(median(bias_lat, na.rm = TRUE), 2),
           cover_lat = round(median(cover_lat, na.rm = TRUE), 2)), by = arm],
    "calibration_decomposition")
cat("\nA control arm that fails to calibrate a deployment records that failure.\n")
cat("n_failed is a result, not a bug: it counts the animals whose own haul-out\n")
cat("cannot calibrate them.\n")

# =============================================================================
rule("4. SEASON -- reported as a test, not as a gradient")
S <- TAGS[arm == "grid_light" & calibrated == TRUE]
put(S[, .(n = .N, median_km = round(median(median_km)),
          IQR_km = round(IQR(median_km)),
          bias_lat = round(median(bias_lat), 2),
          cover_lat = round(median(cover_lat), 2)), by = season][order(season)],
    "season_summary")
season_p <- function(v) {
  if (uniqueN(S$season) < 2) return(NA_real_)
  kruskal.test(S[[v]] ~ factor(S$season))$p.value
}
cat("\n")
if (uniqueN(S$season) < 2) {
  cat("  only one season present -- the season test does not apply (smoke run?)\n")
} else {
  for (v in c("median_km", "rmse_lat", "bias_lat", "cover_lat")) {
    r2 <- summary(lm(S[[v]] ~ factor(S$season)))$r.squared
    cat(sprintf("  %-10s Kruskal-Wallis p = %.3f | R2 = %.3f\n", v, season_p(v), r2))
  }
}
cat("\nWithin-season IQR against the between-season shift is the number that matters.\n")
cat("Three groups of about ten cannot resolve a difference smaller than that spread.\n")

# =============================================================================
rule("5. COVERAGE -- nominal 0.95")
CV <- TAGS[calibrated == TRUE, .(n = .N,
      median = round(median(cover_lat), 3), min = round(min(cover_lat), 3),
      max = round(max(cover_lat), 3),
      n_under_0.5 = sum(cover_lat < 0.5), n_over_0.9 = sum(cover_lat > 0.9),
      lat_sd = round(mean(lat_sd), 3), rmse_lat = round(mean(rmse_lat), 3),
      inflation = round(mean(rmse_lat) / mean(lat_sd), 2)), by = arm]
put(CV, "coverage")
cat("\n`inflation` is the factor by which the reported latitude sd would have to be\n")
cat("multiplied to match the realised error. SBC passes on both arms, so a\n")
cat("real-data shortfall is misspecification against real light residuals rather\n")
cat("than an implementation error.\n")

# =============================================================================
rule("6. WHAT EXPLAINS THE DEPLOYMENT-LEVEL SCATTER?")
vc <- function(y, g) {
  g <- factor(g); k <- nlevels(g); n <- length(y)
  if (k < 2 || n <= k) return(c(between = NA, within = NA, icc = NA))
  ni <- table(g); ybar <- tapply(y, g, mean)
  msb <- sum(ni * (ybar - mean(y))^2) / (k - 1)
  msw <- sum((y - ybar[g])^2) / (n - k)
  n0  <- (n - sum(ni^2) / n) / (k - 1)
  s2b <- max(0, (msb - msw) / n0)
  c(between = s2b, within = msw, icc = s2b / (s2b + msw))
}
G <- merge(TAGS[arm == "grid_light" & calibrated == TRUE], ser[, .(id, serial)], by = "id")
VC <- rbindlist(lapply(c("bias_lat", "rmse_lat", "cover_lat", "median_km"), function(v) {
  z <- vc(G[[v]], G$serial)
  data.table(metric = v, between_logger = round(z[[1]], 4),
             within_logger = round(z[[2]], 4), ICC = round(z[[3]], 3))
}))
put(VC, "variance_components")
cat(sprintf("\n%d loggers, %d deployments. A one-way estimate truncated at zero is weak\n",
            uniqueN(G$serial), nrow(G)))
cat("evidence of absence, not evidence of no effect. Report the ICC with that caveat.\n")

# =============================================================================
rule("7. TIMING -- work-normalised, never a bare wall clock")
put(TIME[, .(n_fits = .N, engine = engine[1],
             total_min = round(sum(seconds) / 60, 1),
             sec_per_fit_median = round(median(seconds)),
             sec_per_knot = round(median(sec_per_knot), 3),
             sec_per_tag_day = round(median(sec_per_tag_day), 4)), by = arm],
    "timing")
gl <- TIME[arm == "grid_light", sum(seconds)]
gf <- TIME[arm == "grid_fusion", sum(seconds)]
if (length(gl) && length(gf) && gl > 0)
  cat(sprintf("\nHARD constraint versus light only: %.0f s -> %.0f s (%+.1f%% wall clock)\n",
              gl, gf, 100 * (gf / gl - 1)))
cat("The claim is narrow and testable: a hard auxiliary constraint speeds the grid\n")
cat("HMM roughly in proportion to the cells it excludes, because the light\n")
cat("likelihood is never evaluated there. A soft constraint carrying the same\n")
cat("information buys none of it.\n")

# =============================================================================
rule("8. PRE-REGISTERED PREDICTIONS -- scored")
verdict <- function(label, passed, detail) cat(sprintf("  [%s] %s\n        %s\n",
                                                       if (isTRUE(passed)) "PASS" else if (is.na(passed)) " -- " else "FAIL",
                                                       label, detail))
p1 <- CONTRASTS[contrast == "grid_light - grid_cal_untrunc" & metric == "median_km"]
verdict("1. truncation improves accuracy, paired p < 0.05",
        nrow(p1) && p1$paired_median_diff < 0 && p1$p_deployment < 0.05,
        if (nrow(p1)) sprintf("paired median diff %+.1f km, better on %d/%d, p = %.4f",
                              p1$paired_median_diff, p1$better_a, p1$n_deploy, p1$p_deployment) else "not run")
p2 <- CONTRASTS[contrast == "grid_cal_untrunc - grid_cal_pertag" & metric == "median_km"]
verdict("2. pooling beats per-tag calibration",
        nrow(p2) && p2$paired_median_diff < 0,
        if (nrow(p2)) sprintf("paired median diff %+.1f km, p = %.4f",
                              p2$paired_median_diff, p2$p_deployment) else "not run")
nf <- TAGS[arm == "grid_cal_pertag", sum(!calibrated)]
verdict("3. per-tag calibration fails outright on some deployments", nf > 0,
        sprintf("%d of %d deployments got no response at all", nf, TAGS[arm == "grid_cal_pertag", .N]))
p4a <- if (length(gl) && length(gf) && gl > 0) gf < gl else NA
p4b <- CONTRASTS[contrast == "grid_fusion - grid_light" & metric == "median_km"]
verdict("4. hard fusion is FASTER and not detectably more accurate",
        isTRUE(p4a) && nrow(p4b) && p4b$p_deployment > 0.05,
        sprintf("wall clock %+.1f%%; accuracy p = %s", 100 * (gf / gl - 1),
                if (nrow(p4b)) format(p4b$p_deployment) else "NA"))
cv <- CV[arm == "grid_light"]
verdict("5. coverage well below 0.95 and bimodal across tags",
        nrow(cv) && cv$median < 0.85,
        sprintf("median %.2f, range %.2f-%.2f, %d tags below 0.5",
                cv$median, cv$min, cv$max, cv$n_under_0.5))
neg <- A[, sum(bias_lat < 0)]
verdict("6. latitude bias negative on most deployments, scatter unexplained",
        neg >= 0.75 * nrow(A),
        sprintf("%d of %d negative, median %+.2f deg, range %.2f deg; logger ICC %s",
                neg, nrow(A), median(A$bias_lat), diff(range(A$bias_lat)),
                format(VC[metric == "bias_lat"]$ICC)))
kp <- season_p("median_km")
verdict("7. season is not significant", if (is.na(kp)) NA else kp > 0.05,
        if (is.na(kp)) "only one season present -- not testable"
        else sprintf("Kruskal-Wallis p = %.3f on median error", kp))

cat(sprintf("\n\ntables written to %s\n", REP))
