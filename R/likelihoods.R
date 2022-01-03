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




#' loglik_trinom_prof
#'
#' @param y data matrix of counts
#' @param r expected ratio of concentrations
#'
#' @return log-likelihood kernel
#' @export
loglik_trinom_prof <- function(y, r){
  if(is.null(dim(y))){
    x_0 <- y[1]
    x_1 <- y[2]
    x_2 <- y[3]
  }else{
    x_0 <- y[,1]
    x_1 <- y[,2]
    x_2 <- y[,3]
  }
  n <-  x_0 + x_1 + x_2
  x_0*log(x_0/n) +
    ifelse(x_1 == 0, yes = 0, no = x_1*log(-x_0/n + exp((1-r)*log(x_0/n)))) +
    ifelse(x_2 == 0, yes = 0, no = x_2*log(1 - exp((1-r)*log(x_0/n))))
}
attributes(loglik_trinom_prof) <- list("overdispersion"=function(ydata, weights, r_fitted_values, l_fitted_values=NULL){
  n <- apply(ydata, 1, sum)
  if(is.null(l_fitted_values)){
    l_fitted_values <- -log(ydata[,1]/n)
  }
  p1 <- exp(-(1-r_fitted_values)*l_fitted_values) - exp(-l_fitted_values)
  p2 <- 1 - exp(-(1-r_fitted_values)*l_fitted_values)
  var_ratio_1 <- (ydata[,2] - p1*n)^2 / (n*p1*(1-p1))
  var_ratio_2 <- (ydata[,3] - p2*n)^2 / (n*p2*(1-p2))
  return(list(var_ratio_1, var_ratio_2))
})


