rm(list = ls())
options(warn=-1)
setwd("~/Dropbox/academic/phd/gibbons/gibbon_simulation")

# load simulation parameters
source("./simulations/sim_parameters.R")

# run simulation and print summary statistics
source("./functions/tri_sim.R")

result = sapply(rep(0, tri_n_sims), tri_sim, density = true_density, sigma = sigma, 
                spacing = spacing, multi = TRUE, max_dist = max_dist, union = FALSE, 
                edr_val = TRUE)
print("Results for original triangulation method")
print(list("mean" = mean(result), "sd" = sd(result), "true_density" = true_density))

# run simulation computing edr including calls heard by one listener
result = sapply(rep(0, tri_n_sims), tri_sim, density = true_density, sigma = sigma, 
                spacing = spacing, multi = FALSE, max_dist = max_dist, union = FALSE, 
                edr_val = TRUE)
print("Results for modified edr triangulation method")
print(list("mean" = mean(result), "sd" = sd(result), "true_density" = true_density))

# run simulation computing effective area with union rather than intersection
result = sapply(rep(0, tri_n_sims), tri_sim, density = true_density, sigma = sigma, 
                spacing = spacing, multi = FALSE, max_dist = max_dist, union = TRUE,
                edr_val = TRUE)
print("Results for union triangulation method")
print(list("mean" = mean(result), "sd" = sd(result), "true_density" = true_density))

# run original simulation with fixed edr = 1000 (1km)
result = sapply(rep(0, tri_n_sims), tri_sim, density = true_density, sigma = sigma, 
                spacing = spacing, multi = TRUE, max_dist = max_dist, union = FALSE, 
                edr_val = 1000)
print("Results for original triangulation method with fixed edr")
print(list("mean" = mean(result), "sd" = sd(result), "true_density" = true_density))

# run original simulation with fixed edr = 1000 (1km)
result = sapply(rep(0, tri_n_sims), tri_sim, density = true_density, sigma = sigma, 
                spacing = spacing, multi = TRUE, max_dist = max_dist, union = FALSE, 
                edr_val = TRUE, true_sigma = sigma)
print("Results for with true sigma for EDR computation:")
print(list("mean" = mean(result), "sd" = sd(result), "true_density" = true_density))

# run original simulation with fixed edr = 1000 (1km)
result = sapply(rep(0, tri_n_sims), tri_sim, density = true_density, sigma = sigma, 
                spacing = spacing, multi = TRUE, max_dist = max_dist, union = TRUE, 
                edr_val = TRUE, true_sigma = sigma)
print("Results for with true sigma for EDR computation and union:")
print(list("mean" = mean(result), "sd" = sd(result), "true_density" = true_density))

# original simulation with repelling data (hardcore)
result = sapply(rep(0, tri_n_sims), tri_sim, density = true_density, sigma = sigma, 
                spacing = spacing, multi = TRUE, max_dist = max_dist, union = FALSE, 
                edr_val = TRUE, true_sigma = sigma, spatial = 'hardcore')
print("Results for with true sigma for hardcore process:")
print(list("mean" = mean(result), "sd" = sd(result), "true_density" = true_density))

# original simulation with repelling data (strauss hardcore)
result = sapply(rep(0, tri_n_sims), tri_sim, density = true_density, sigma = sigma, 
                spacing = spacing, multi = TRUE, max_dist = max_dist, union = FALSE, 
                edr_val = TRUE, true_sigma = sigma, spatial = 'straush')
print("Results for with true sigma for strauss hardcore process:")
print(list("mean" = mean(result), "sd" = sd(result), "true_density" = true_density))