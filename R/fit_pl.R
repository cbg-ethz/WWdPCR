#' fit_pl
#'
#' @param ydata response variable
#' @param xdata independent variable
#' @param weights weights vector (optional).
#' @param fitfun fitting function
#' @param likelihood likelihood function
#' @param quasi compute overdispersion and use quasilikelihood model to build confidence intervals
#' @param starting_values (named) vector of starting values for optimization
#' @param lower vector of lower bounds for parameters
#' @param upper vector of upper bounds for parameters
#' @param ... additional argument(s)
#'
#' @return pl_fit object
#' @export

fit_pl <- function(ydata, xdata, weights=1, fitfun=pl3model, likelihood=loglik_binom_n,
                   quasi=F,
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

  fitted_values <- do.call(fitfun, c(list(xdata), as.list(opt1$par)))
  if(quasi==T){
    overdispersion <- attributes(likelihood)$overdispersion(ydata, weights, fitted_values)
    infl_factor <- 1/(length(unlist(overdispersion))-1) * sum(unlist(overdispersion))
  }else{
    overdispersion <- NULL
    infl_factor <- 1
  }
  out <- list(
    estimates=opt1$par,
    se=sqrt(diag(solve(-1*opt1$hessian)))*infl_factor,
    fitted=fitted_values,
    ydata=ydata,
    xdata=xdata,
    weights=weights,
    opt=opt1,
    infl_factor=infl_factor,
    overdispersion=overdispersion,
    likelihood=opt1$value,
    fitfun=fitfun
  )
  class(out) <- "pl_fit"
  return(out)
}

#' @method coef pl_fit
#' @rdname fit_pl
#' @param object an object of class \code{pl_fit}
#' @export
coef.pl_fit <- function(object, ...){
  object$estimates
}

#' @method se pl_fit
#' @rdname fit_pl
#' @param x an object of class \code{pl_fit}
#' @export
se.pl_fit <- function(x){
  x$se
}

#' @method fitted pl_fit
#' @rdname fit_pl
#' @param object an object of class \code{pl_fit}
#' @export
fitted.pl_fit <- function(object, ...){
  object$fitted
}


#' @method confint pl_fit
#' @rdname fit_pl
#' @param object an object of class \code{pl_fit}
#' @param parm a specification of which parameters are to be given confidence
#' intervals, either a vector of numbers or a vector of names. If missing,
#' all parameters are considered.
#' @param level confidence level 1-alpha
#'
#' @export
confint.pl_fit <- function(object, parm, level=0.95, ...){
  pars_est_df <- cbind(object$estimates,
                       object$estimates - qnorm(1-(1-level)/2) * object$se,
                       object$estimates + qnorm(1-(1-level)/2) * object$se)
  colnames(pars_est_df) <- c("estimate", "lower", "upper")
  # if(!is.null(parm)){
  #   pars_est_df[parm,]
  if(!missing(parm)){
    pars_est_df[parm,]
  }else{
    pars_est_df
  }
}

#' @method predict pl_fit
#' @rdname fit_pl
#' @param object an object of class \code{pl_fit}
#' @param newdata an optional vector of observed independent variables with which to predict
#' @param scale scale on which to predict, either linear or logit scale
#' @param se whether to return the standard errors of the predictions or not
#'
#' @export
predict.pl_fit <- function(object, newdata = NULL, scale=c("linear", "logit"), se=TRUE, ...){
  if(is.null(newdata)){
    newdata <- object$xdata
  }
  pred <- do.call(object$fitfun, c(list(newdata), as.list(object$estimates)))
  if(scale=="logit"){
    pred=logit(pred)
  }
  out = list("pred"=pred)

  if(se==TRUE){
    # choose appropriate jacobian for delta method
    grad <- switch(scale,
                   "linear"=attributes(object$fitfun)$grad,
                   "logit"=attributes(object$fitfun)$logit_grad
                   )
    # compute standard errors of the fit
    stand.errors <- apply(as.matrix(newdata), 1, function(t){
      grad_eval <- do.call(grad, c(list(t), as.list(object$estimates)))
      sqrt(t(grad_eval) %*% (solve(-1*object$opt$hessian) * object$infl_factor) %*% grad_eval)
    })
    out$se <- stand.errors
  }
  return(out)
}

#' logist_confint
#'
#' @param x an object of class \code{pl_fit}
#' @param ... additional argument(s) for methods
#'
#' @return list of confidence intervals
#' @export
logist_confint <- function(x, ...) {
  UseMethod("logist_confint")
}

#' @method logist_confint pl_fit
#' @rdname fit_pl
#' @param x an object of class \code{pl_fit}
#' @param newdata an optional vector of observed independent variables with which to predict
#' @param level confidence level 1-alpha
#' @export
logist_confint.pl_fit <- function(x, newdata = NULL, level = 0.95, ...){
  predicted <-  predict(x, scale="logit", se=TRUE)
  return(
    list("lower" = logit_inv(predicted$pred - qnorm(1-(1-level)/2)*predicted$se),
         "upper" = logit_inv(predicted$pred + qnorm(1-(1-level)/2)*predicted$se))
         )

}

#' logist_predint
#'
#' @param x an object of class \code{pl_fit}
#' @param ... additional argument(s) for methods
#'
#' @return list of prediction intervals
#' @export
logist_predint <- function(x, ...) {
  UseMethod("logist_predint", ...)
}


#' @method logist_predint pl_fit
#' @rdname fit_pl
#' @param x an object of class \code{pl_fit}
#' @param newdata an optional vector of observed independent variables with which to predict
#' @param level confidence level 1-alpha for prediction
#' @export
logist_predint.pl_fit <- function(x, newdata = NULL, level = 0.95, likelihood=loglik_binom_n, ...){

  weighted.var <- function(x, w){
    if(length(w) == 1){w=rep(w, length(x))}
    sum((x - weighted.mean(x, w))^2 * w) / sum(w)
  }

  predicted <-  predict(x, scale="logit", se=TRUE)

  if(is.null(dim(x$ydata))){ #hack to select binom other trinom
    se <- sqrt(weighted.var(logit((x$ydata * x$weights +1) / (x$weights+2)) - logit(x$fitted), x$weights) + predicted$se^2)
  }else{
    se <- sqrt(weighted.var(logit(with(x$ydata+1, 1-log((x0+x1)/(x0+x1+x2))/log(x0/(x0+x1+x2)))) - logit(c(x$fitted)), rowSums(x$ydata)) + predicted$se^2)
  }

  return(
    list("lower" = logit_inv(predicted$pred - qnorm(1-(1-level)/2)*se),
         "upper" = logit_inv(predicted$pred + qnorm(1-(1-level)/2)*se))
  )

}
