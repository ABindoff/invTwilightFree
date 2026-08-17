cl <- readRDS("analysis/cache/nes/clocks_v2.rds")
cat("class:", class(cl), " length:", length(cl), "\n")
out <- do.call(rbind, lapply(names(cl), function(id) {
  x <- cl[[id]]
  if (is.null(x)) return(data.frame(id = id, offset_min = NA_real_, note = "NULL"))
  o <- if (!is.null(x$offset)) x$offset else if (!is.null(x$offset_sec)) x$offset_sec else NA
  data.frame(id = id, offset_min = round(as.numeric(o) / 60, 1),
             note = paste(names(x), collapse = ","), row.names = NULL)
}))
print(out, row.names = FALSE)
cat("\nfields of the first entry:\n")
str(cl[[1]], max.level = 1)
