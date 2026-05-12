
library(SGAT)
cached <- readRDS('c:/Users/bindoffa/antigravity projects/invTwilightFree/vignettes/cache_sgat.rds')
fit_sgat <- cached$fit

cat("Structure of fit_sgat:\n")
str(fit_sgat, max.level = 1)

if(!is.null(fit_sgat$x)) {
  cat("Structure of fit_sgat$x:\n")
  str(fit_sgat$x, max.level = 2)
  
  sgat_chain <- fit_sgat$x[[1]]
  cat("Class of sgat_chain:", class(sgat_chain), "\n")
  cat("Dimensions of sgat_chain:", paste(dim(sgat_chain), collapse="x"), "\n")
  
  if(is.array(sgat_chain)) {
    sgat_lon <- apply(sgat_chain[,1,], 1, mean)
    cat("First 5 Lon:", head(sgat_lon), "\n")
  }
}
