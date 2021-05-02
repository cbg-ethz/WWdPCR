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


coef.pl_fit <- function(x){
  x$estimates
}

se.pl_fit <- function(x){
  x$se
}

fitted.pl_fit <- function(x){
  x$fitted
}

# confint.pl_fit <- function(x){
#   x$
# }
