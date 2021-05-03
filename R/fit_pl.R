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
    opt=opt1,
    infl_factor=1,
    fitfun=fitfun
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

confint.pl_fit <- function(x, level=0.95){
  pars_est_df <- cbind(x$estimates,
                       x$estimates - qnorm(1-(1-level)/2) * x$se,
                       x$estimates + qnorm(1-(1-level)/2) * x$se)
  colnames(pars_est_df) <- c("estimate", "lower", "upper")
  pars_est_df
}

predict.pl_fit <- function(x, newdata = NULL, scale=c("linear", "logit"), se=TRUE){
  if(is.null(newdata)){
    newdata <- x$xdata
  }
  pred <- do.call(x$fitfun, c(list(newdata), as.list(x$estimates)))
  if(scale=="logit"){
    pred=logit(pred)
  }
  out = list("pred"=pred)

  if(se==TRUE){
    # choose appropriate jacobian for delta method
    grad <- switch(scale,
                   "linear"=attributes(x$fitfun)$grad,
                   "logit"=attributes(x$fitfun)$logit_grad
                   )
    # compute standard errors of the fit
    stand.errors <- apply(as.matrix(newdata), 1, function(t){
      grad_eval <- do.call(grad, c(list(t), as.list(x$estimates)))
      t(grad_eval) %*% (solve(-1*x$opt$hessian)) %*% grad_eval
    })
    out$se <- stand.errors * x$infl_factor
  }
  out$bla <- 1
  return(out)
}

