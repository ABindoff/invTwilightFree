#' Find max light observations over a coarse time grid then interpolate at finer grid
#' @param d data frame with Light, Date, and optionally Depth columns
#' @param depth1 depth filter applied to the coarse grid pass; FALSE to skip
#' @param depth2 depth filter applied to the fine grid pass; FALSE to skip
#' @param period1 time period for the coarse grid (e.g. "30 minutes")
#' @param period2 time period for the fine grid (e.g. "4 minutes")
#' @export
#' @return see help(max_light_delta) for relevant details
#' @importFrom stats splinefun
interpolate_max_light <- function(d, depth1 = FALSE, depth2 = FALSE, period1 = "30 minutes", period2 = "4 minutes"){
  if(depth1){
    b <- max_light(d[d$Depth < depth1,], period1)
  } else {
    b <- max_light(d, period1)
  }

  fn <- splinefun(b$Date, b$Light, method = "periodic")  # build spline function on max_light

  if(depth2){
    b <- max_light(d[d$Depth < depth2,],  period2)
  } else {
    b <- max_light(d, period2)
  }

  b$Light <- fn(b$Date)      # interpolate on finer scale
  b$Light[1] <- b$Light0[1]  # correct the interpolated first observation
  b$Light[b$Light0 > b$Light] = b$Light0[b$Light0 > b$Light]  # correct under-interpolated observations
  b$Light[b$Light > quantile(b$Light0, 0.99)] = quantile(b$Light0, 0.99)  # correct over-interpolated observations

  b$delta <- c(0, b$Light[2:nrow(b)] - b$Light[1:(nrow(b)-1)])# calculate delta-Light for latent class models
  return(b)
}
