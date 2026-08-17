a <- readRDS("analysis/cache/nes/fvilches_archives_v1_30min.rds")
g <- readRDS("analysis/cache/nes/fvilches_argos_v1.rds")
o <- readRDS("analysis/cache/nes/archives_v3_30min.rds")
cat(sprintf("new archives: %d   new argos: %d   original 2021: %d\n",
            length(a), length(g), length(o)))
cat(sprintf("combined: %d   overlap (should be 0): %d\n",
            length(a) + length(o), length(intersect(names(a), names(o)))))
cat("\nnew ids:", paste(sort(names(a)), collapse = " "), "\n")
cat("\nboth archive and argos:", length(intersect(names(a), names(g))), "\n")
