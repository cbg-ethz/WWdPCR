r_hat <- function(x){
  if(is.null(dim(x))){
    x_0 <- x[1]
    x_1 <- x[2]
    x_2 <- x[3]
  }else{
    x_0 <- x[,1]
    x_1 <- x[,2]
    x_2 <- x[,3]
  }
  n <- x_0 + x_1 + x_2
  l1 <- -log((x_0+x_1)/n)
  l2 <- -log(x_0/n)
  1-(l1/l2)
}


lower_e <- function(x, level=0.95){
  if(is.nan(r_hat(x))){
    return(NaN)
  }else if(r_hat(x) == 0){
    return(0)
  }else{
    uniroot(function(y){2*(loglik_trinom_prof(x,y) - loglik_trinom_prof(x,r_hat(x))) + qchisq(level, 1)},
            interval=c(0, r_hat(x)))$root
  }
}

upper_e <- function(x, level=0.95){
  if (is.nan(r_hat(x))){
    return(NaN)
  }else if (r_hat(x) == 1){
    return(1)
  }else{
    uniroot(function(y){2*(loglik_trinom_prof(x,y) - loglik_trinom_prof(x,r_hat(x))) + qchisq(level, 1)},
            interval=c(r_hat(x), 1), f.lower = qchisq(level, 1))$root
  }
}

lrt_confint <- function(x, level=0.95){
  tmp_res <- apply(x, 1, function(i){
    list(lower=lower_e(i, level = level), upper=upper_e(i, level = level))
  })
  Reduce(function(...) rbind.data.frame(...), tmp_res)
}

