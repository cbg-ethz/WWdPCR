#' loglik_binom_n
#'
#' @param y empirical success frequency
#' @param p probability of success
#'
#' @return log-likelihood kernel
#' @export
loglik_binom_n <- function(y,p){
  y*log(p/(1-p))+log(1-p)
}
attributes(loglik_binom_n) <- list("overdispersion"=function(ydata, weights, fitted_values){
  (ydata*weights - fitted_values*weights)^2 / (weights*fitted_values*(1-fitted_values))
})
