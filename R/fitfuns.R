#' pl2model
#'
#' @param x time
#' @param a rate parameter
#' @param b bias parameter
#'
#' @return a two parametric logistic curve
#' @details returns a logistic curve with bottom asymptote at 0, such that: \eqn{r(x;a,b) = \frac{\exp(ax+b)}{1+exp(a*x+b)}}
#' @export
pl2model <- function(x,a,b){
  (exp(a*x+b)/(1+exp(a*x+b)))
}

#' pl3model
#'
#' @param x time
#' @param a rate parameter
#' @param b bias parameter
#' @param c background parameter
#'
#' @return a three parametric logistic curve
#' @details returns a logistic curve with bottom asymptote not at 0 but at parameter c, such that: \eqn{r(x;a,b,c) = c + (1-c) \frac{\exp(ax+b)}{1+exp(a*x+b)}}
#' @export
pl3model <- function(x,a,b,c){
  c + (1-c)*(exp(a*x+b)/(1+exp(a*x+b)))
}


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


fit_pl <- function(ydata, xdata, weights, fitfun=pl3model, likelihood=loglik_binom_n,
                   starting_values=c("a"=0.5, "b"=-10, "c"=0.2),
                   lower=c(-Inf, -Inf, 0), upper=c(Inf, Inf, 1)){

  opt1 <- optim(starting_values,
                function(x){
                  fit <- do.call(fitfun, c(list(xdata), as.list(x)))
                  sum(likelihood(ydata, fit)*weights)
                },
                control=list(fnscale=-1),
                method="L-BFGS-B",
                lower=lower, upper=upper,
                hessian=T)
  out <- list(
    estimates=opt1$par,
    se=sqrt(diag(solve(-1*opt1$hessian))),
    fitted=do.call(fitfun, c(list(xdata), as.list(opt1$par))),
    ydata=ydata,
    xdata=xdata,
    weights=weights,
    opt=opt1
  )
  class(out) = "pl_fit"
  return(out)
}


