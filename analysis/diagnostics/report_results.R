r <- read.csv("analysis/output/nes/per_tag_results.csv")
k <- read.csv("analysis/output/nes/per_knot_errors.csv")

agg <- do.call(rbind, lapply(split(k, list(k$engine, k$batch), drop = TRUE), function(d) {
  d <- d[is.finite(d$err_km), ]
  data.frame(engine = d$engine[1], batch = d$batch[1], n = nrow(d),
             median_km = round(median(d$err_km)),
             q90_km = round(unname(quantile(d$err_km, .9))),
             bias_lon = round(mean(d$err_lon), 2),
             bias_lat = round(mean(d$err_lat), 2),
             rmse_lat = round(sqrt(mean(d$err_lat^2)), 2), row.names = NULL)
}))
cat("=== batches, pooled over knots ===\n")
print(agg[order(agg$median_km), ], row.names = FALSE)

cat("\n=== coverage, mean over tags ===\n")
cv <- aggregate(cbind(cover95_lon, cover95_lat) ~ engine + batch, r,
                function(x) round(mean(x), 2))
print(cv, row.names = FALSE)

cat("\n=== per tag, hier/light ===\n")
h <- r[r$engine == "hier" & r$batch == "light",
       c("id", "median_km", "bias_lat", "rmse_lat_deg", "cover95_lat")]
h[] <- lapply(h, function(x) if (is.numeric(x)) round(x, 2) else x)
print(h[order(h$median_km), ], row.names = FALSE)

cat("\n=== calibration recipes (grid HMM) ===\n")
cal <- r[r$batch %in% c("naive", "pooled"), ]
if (nrow(cal)) {
  w <- merge(cal[cal$batch == "naive",  c("id", "median_km", "rmse_lat_deg")],
             cal[cal$batch == "pooled", c("id", "median_km", "rmse_lat_deg")],
             by = "id", suffixes = c("_naive", "_pooled"))
  w[] <- lapply(w, function(x) if (is.numeric(x)) round(x) else x)
  print(w, row.names = FALSE)
  ck <- k[k$batch %in% c("naive", "pooled"), ]
  print(do.call(rbind, lapply(split(ck, ck$batch), function(d) {
    d <- d[is.finite(d$err_km), ]
    data.frame(recipe = d$batch[1], median_km = round(median(d$err_km)),
               q90_km = round(unname(quantile(d$err_km, .9))),
               bias_lat = round(mean(d$err_lat), 2),
               rmse_lat = round(sqrt(mean(d$err_lat^2)), 2), row.names = NULL)
  })), row.names = FALSE)
}
