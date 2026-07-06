# =============================================================================
# ecdf_bands.R  -- reusable simultaneous ECDF confidence bands for SBC ranks
#
# Vendored from the sbc-validation skill for self-contained use in this package.
# Method: Sailynoja, Burkner & Vehtari (2022), "Graphical test for discrete
# uniformity ...", Statistics and Computing 32(2):32. A single pointwise tail
# level gamma is chosen by Monte Carlo so that, under the discrete-uniform null,
# the ENTIRE rank ECDF stays inside the band with probability `conf`. No external
# packages are required for the band itself; ggplot2 is used only for the
# optional plot and is loaded lazily.
#
# Conventions:
#   * ranks are integers in {0, ..., L}, where L = number of (thinned) posterior
#     draws used to rank each truth. rank = #{draws < truth}.
#   * a parameter PASSES if its rank-ECDF minus the uniform diagonal stays within
#     the simultaneous band across the whole [0,1] range.
# =============================================================================

sbc_sim_band <- function(N, L, conf = 0.95, M = 4000L, grid = NULL) {
  if (N < 1) return(NULL)
  if (is.null(grid)) grid <- (1:L) / (L + 1)
  thr <- grid * (L + 1) - 1
  null_counts <- function() {
    R <- matrix(sample(0:L, N * M, replace = TRUE), nrow = M)
    vapply(thr, function(t) rowSums(R <= t), numeric(M))   # M x G counts
  }
  cal <- null_counts(); ev <- null_counts()
  cover_at <- function(gamma) {
    lo <- apply(cal, 2, quantile, probs = gamma / 2,     type = 1)
    hi <- apply(cal, 2, quantile, probs = 1 - gamma / 2, type = 1)
    inside <- rowMeans((ev >= matrix(lo, nrow(ev), length(lo), byrow = TRUE)) &
                       (ev <= matrix(hi, nrow(ev), length(hi), byrow = TRUE))) == 1
    list(cov = mean(inside), lo = lo, hi = hi)
  }
  g_lo <- 1e-4; g_hi <- 0.5    # coverage is monotone decreasing in gamma
  for (it in 1:40) {
    g <- (g_lo + g_hi) / 2
    if (cover_at(g)$cov >= conf) g_lo <- g else g_hi <- g
  }
  b <- cover_at(g_lo)
  data.frame(p = grid, lo = b$lo / N, hi = b$hi / N, gamma = g_lo)
}

sbc_ecdf_diff <- function(ranks, L, band) {
  ranks <- ranks[!is.na(ranks)]
  thr <- band$p * (L + 1) - 1
  obs <- vapply(thr, function(t) mean(ranks <= t), numeric(1))
  data.frame(p = band$p, diff = obs - band$p, lo = band$lo - band$p, hi = band$hi - band$p)
}

sbc_ecdf_panel <- function(ranks, L, title = "", conf = 0.95, plot = TRUE) {
  N <- sum(!is.na(ranks))
  band <- sbc_sim_band(N, L, conf = conf)
  if (is.null(band)) return(list(plot = NULL, pass = NA, gamma = NA, N = N))
  dd <- sbc_ecdf_diff(ranks, L, band)
  pass <- all(dd$diff >= dd$lo & dd$diff <= dd$hi)
  p <- NULL
  if (plot && requireNamespace("ggplot2", quietly = TRUE)) {
    p <- ggplot2::ggplot(dd, ggplot2::aes(p)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), fill = "grey80") +
      ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
      ggplot2::geom_step(ggplot2::aes(y = diff), linewidth = 0.6) +
      ggplot2::labs(x = "Fractional rank", y = "ECDF - uniform",
                    title = tryCatch(parse(text = title), error = function(e) title)) +
      ggplot2::theme_minimal(base_size = 11)
  }
  list(plot = p, pass = pass, gamma = band$gamma[1], N = N)
}

sbc_ecdf_report <- function(df, specs, L, conf = 0.95, fig = NULL, tag = "SBC") {
  res <- lapply(names(specs), function(col) {
    pn <- sbc_ecdf_panel(df[[col]], L, specs[[col]], conf = conf, plot = !is.null(fig))
    data.frame(param = specs[[col]], column = col, N = pn$N,
               pass = pn$pass, gamma = pn$gamma, stringsAsFactors = FALSE,
               row.names = NULL)
  })
  out <- do.call(rbind, res)
  cat(sprintf("\n[%s] simultaneous ECDF bands (%.0f%%):\n", tag, 100 * conf))
  for (i in seq_len(nrow(out)))
    cat(sprintf("  %-22s %s (N=%d, gamma=%.4f)\n", out$param[i],
                ifelse(is.na(out$pass[i]), "n/a", ifelse(out$pass[i], "PASS", "FAIL")),
                out$N[i], out$gamma[i]))
  if (!is.null(fig) && requireNamespace("ggplot2", quietly = TRUE) &&
      requireNamespace("patchwork", quietly = TRUE)) {
    panels <- lapply(names(specs), function(col)
      sbc_ecdf_panel(df[[col]], L, specs[[col]], conf = conf, plot = TRUE)$plot)
    panels <- Filter(Negate(is.null), panels)
    if (length(panels)) {
      combined <- Reduce(`+`, panels) + patchwork::plot_layout(nrow = 1)
      ggplot2::ggsave(fig, combined, width = 3.7 * length(panels), height = 3.6, dpi = 200)
      cat(sprintf("[%s] figure -> %s\n", tag, fig))
    }
  }
  invisible(out)
}
