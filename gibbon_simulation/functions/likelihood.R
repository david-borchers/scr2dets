# negative log-likelihood 
neg_log_lk = function(sigma, dists, n, max_dist){
  out = 0
  f = function(x){return(exp(-x^2 / (2 * sigma^2)) * 2 * x / max_dist^2)}
  ep = integrate(f, 0, max_dist)
  for(i in 1:n){
    out = out + log(f(dists[i]) /  ep$value)
  }
  return(-out)
}
