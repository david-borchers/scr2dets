# function to compute edr of listening post given maximum call detection distance
edr = function(distances, max_dist, true_sigma = FALSE){
  # find mle of detection function parameter
  if(true_sigma == FALSE){
    mle = optimize(neg_log_lk, interval = c(0, 50000), dists = distances, 
                   n = length(distances), max_dist = max_dist)
  } else {
    mle = list(minimum = true_sigma)
  }
 
  
  
  # fitted half normal detection function
  p = function(x){
    return(exp(-x^2 / (2*mle$minimum^2)))
  }
  
  # compute edr
  edr = integrate(p, 0, max_dist)$value
  attr(edr, "par") = mle$minimum
  return(edr)
}
