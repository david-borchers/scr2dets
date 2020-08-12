
options(warn=-1)
setwd("~/Dropbox/academic/phd/gibbons/gibbon_simulation")
source("./functions/secr_sim.R")

# load simulation parameters
source("./simulations/sim_parameters.R")

# run simulation and print summary statistics
result = sapply(rep(0, secr_n_sims), secr_sim, density = true_density, 
                sigma = sigma, spacing = spacing, max_dist = max_dist * 2, spatial = FALSE, noccasions = 5)
print("Results for secr analysis")
print(list("mean" = mean(result), "sd" = sd(result), "true_density" = true_density))
