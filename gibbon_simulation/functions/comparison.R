require(secr)
source("./functions/likelihood.R")
source("./functions/edr.R")
source("./functions/circle.R")
source("./functions/generate_data.R")
source("./functions/effective_area.R")

comparison = function(x, density, sigma, spacing, max_dist, multi = TRUE, ratio = FALSE){
  tri_output = NA
  # avoid NaN output in extreme data cases
  while(is.nan(tri_output) | is.na(tri_output)){
    # generate data
    traps = generate_traps(spacing = spacing)
    capt_hist = capture_history(traps, density, sigma, spacing)
    
    # secr estimate
    fitted_model = secr.fit(capt_hist, fixed = list(g0 = 1), detectfn = , buffer = max_dist, trace = FALSE)
    
    # triangulation estimate
    distances = trial_distances(capt_hist, traps, density, sigma, spacing, multi = multi)
    area = effective_area(distances, traps, max_dist, union = FALSE, edr = TRUE)
    tri_output = 10000 * length(distances[1, ]) / area
  }
  if(ratio == TRUE){
    return(c(tri_output, exp(fitted_model$fit$estimate)[1]) / density)
  } else {
  return(c(tri_output, exp(fitted_model$fit$estimate)[1]))
  }
}
