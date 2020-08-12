rm(list = ls())
options(warn=-1)
setwd("~/Dropbox/academic/phd/gibbons/gibbon_simulation")
source("./functions/comparison.R")

# load simulation parameters
source("./sim_parameters.R")

# run simulation and print summary statistics
result = sapply(rep(0, comp_n_sims), comparison, density = true_density, 
                sigma = sigma, spacing = spacing, max_dist = max_dist, 
                multi = TRUE, ratio = FALSE)

print("Comparing methods on identical simulated data sets")
print("Triangulation:")
print(list("mean" = mean(result[1, ]), "sd" = sd(result[1, ]), "true_density" = true_density))

print("SECR:")
print(list("mean" = mean(result[2, ]), "sd" = sd(result[2, ]), "true_density" = true_density))
