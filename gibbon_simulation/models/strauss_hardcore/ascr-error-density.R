rm(list=setdiff(ls(), c("edout", "esout", "edout_union")))
options(warn=-1)
setwd("~/Dropbox/academic/phd/gibbons/gibbon_simulation")

# run simulation and print summary statistics
source("./functions/tri_sim.R")

# specify parameters
n_sims = 100
sigma = 500
spacing = 400
true_density = 0.005
num_instances = 100
output = matrix(0, nrow = num_instances, ncol = n_sims)

for(i in 1:num_instances){
  output[i, ] = sapply(output[i, ], ascr_sim, density = true_density, sigma = sigma, 
                       spacing = spacing) / true_density
  true_density = true_density + 0.0005
  print(i)
}

# plot mean error ratios
plot(y = rowMeans(output), x = seq(0.001, 0.0505, 0.0005), xlab = "Density / Hectare", ylab = "Estimate / True Density", 
     xlim = c(0, 0.05), ylim = c(0, 6))

lines(seq(0.001, 0.0505, 0.0005), rowMeans(output)-apply(output, 1, sd), col = "BLUE")
lines(seq(0.001, 0.0505, 0.0005), rowMeans(output)+apply(output, 1, sd), col = "BLUE")

ascr_edout = output