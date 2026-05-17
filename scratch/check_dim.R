c_sgat <- readRDS('vignettes/cache_sgat_seal_ideal.rds')
obj <- c_sgat$fit
cat('Dim:', dim(obj$x[[1]]), '\n')
