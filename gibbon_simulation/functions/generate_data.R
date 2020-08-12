require(secr)
require(spatstat)

# wrapper function for secr make.circle
generate_traps = function(spacing){
  xrange = yrange = spacing * 30
  return(make.circle(n = 3, spacing = spacing, detector = "proximity",  originxy = c(xrange / 2, yrange / 2)))
}

capture_history = function(traps, density, sigma, spacing, spatial = FALSE, noccasions = 1){
  xrange = yrange = spacing * 30 
  # simulate the population with some reasonable guess for density, plot(pop) to verify
  
  pop = sim.popn(D = density, expand.grid(x = c(0, xrange), y = c(0, xrange)), 0)
  window = owin(xrange = c(0, xrange), yrange = c(0, xrange))
  
  if(spatial != FALSE){
    if(spatial == 'hardcore'){
      hardcore = list( beta = 1e-5, hc = 400)
      sim_location = rmh(model = list(cif = spatial, par = hardcore, w = window),
                        start = list(n.start = nrow(pop)),
                        control = list(p = 1), verbose = FALSE)
    } else if(spatial == 'straush'){
      strausshard = list( beta = 1e-5, gamma=0.5, r = 700, hc = 400)
      sim_location = rmh(model = list(cif = spatial, par = strausshard, w = window),
                        start = list(n.start = nrow(pop)),
                        control = list(p = 1), verbose = FALSE)
    }
    pop$x <- sim_location$x
    pop$y <- sim_location$y
  }
  
  out = sim.capthist(traps, popn = pop, noccasions = noccasions, detectpar = list(g0 = 1, sigma = sigma), 
                     savepopn = TRUE, renumber = FALSE, p.available = 1)
  while(is.null(dim(out[, 1, ]))){
    out = sim.capthist(traps, popn = pop, noccasions = noccasions, detectpar = list(g0 = 1, sigma = sigma), 
                       savepopn = TRUE, renumber = FALSE, p.available = 1)
  }
  
  attr(out, 'loc') = pop
  return(out)
}

# function to return distances from detected animals to traps
trial_distances = function(capt_hist, traps, density, sigma, spacing, multi){
  # calculate the number of times eacapt_hist individual was detected
  detect.ns = apply(capt_hist[, 1, ], 1, sum)
    
  # get the locations of all detected individuals
  detect.IDs = as.integer(row.names(capt_hist[,1,]))
  detect.XYs = attr(capt_hist,"popn")[detect.IDs,]
    
  # get distances from all traps to all detected individuals
  detect.dists = edist(traps, detect.XYs)
  colnames(detect.dists) = detect.IDs

  # get distances from all traps to those individuals detected more than once
  if(multi == TRUE){
    return(detect.dists[, detect.ns > 1])
  } else if(multi == FALSE){
    return(detect.dists)
  }
}

# function to repeatedly generate population until sufficient distance data generated 
# by omitting cases with only 1 gibbon group detected
generate_distances = function(traps, density, sigma, spacing, multi = multi, spatial = FALSE){
  distances = c()
  while(length(distances) < 3 | is.null(dim(distances))){
    capt_hist = capture_history(traps, density, sigma, spacing, spatial = spatial)
    distances = trial_distances(capt_hist, traps, density, sigma, spacing, multi = multi)
  }
  return(distances)
}
