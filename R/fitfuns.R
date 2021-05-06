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
attributes(pl2model) <- list(
  "grad"=function(x,a,b){
    c(x*exp(a*x+b)/(1+exp(a*x+b))^2,
      exp(a*x+b)/(1+exp(a*x+b))^2)
    },
  "logit_grad"=function(x,a,b){
    c(x, 1)
    })


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
attributes(pl3model) <- list(
  "grad"=function(x,a,b,c){
  c((1-c)*x*(exp(a*x+b)/(1+exp(a*x+b))),
    (1-c)*(exp(a*x+b)/(1+exp(a*x+b))),
    1-(exp(a*x+b)/(1+exp(a*x+b)))
  )},
  "logit_grad"=function(x,a,b,c){
    c(x*exp(a*x + b)/(c + exp(a*x + b)),
      exp(a*x + b)/(c + exp(a*x + b)),
      -(exp(a*x + b) + 1)/((c - 1)*(c + exp(a*x + b)))
    )})




# functions to compute forward and backward reparametrizations
logit <- function(r){log(r/(1-r))}
logit_inv <- function(psi){1/(1+exp(-psi))}





