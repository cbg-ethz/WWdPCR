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
