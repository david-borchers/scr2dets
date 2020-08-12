rm(list = ls())
options(warn=-1)
setwd("~/Dropbox/academic/phd/gibbons/gibbon_simulation")

# run simulation and print summary statistics
source("./functions/tri_sim.R")

# specify parameters
tri_n_sims = 100
sigma = 400
spacing = 400
max_dist = sigma * 3
true_density = 0.001
num_instances = 49

output = matrix(0, nrow = num_instances, ncol = tri_n_sims)
for(i in 1:num_instances){
  output[i, ] = sapply(output[i, ], tri_sim, density = true_density, sigma = sigma, 
                  spacing = spacing, multi = TRUE, max_dist = max_dist, union = FALSE, 
                  edr_val = TRUE, spatial = 'hardcore') / true_density
  true_density = true_density + 0.001
}

# plot mean error ratios
plot(y = rowMeans(output), x = seq(0.001, 0.049, 0.001), xlab = "Density / Hectare", ylab = "Estimate / True Density", 
     xlim = c(0, 0.05), ylim = c(1, 10), main = "Effect of density on error for Triangulation method")

# fit reciprocal quadratic model and plot
qr.model = lm(rowMeans(output) ~ I(seq(0.001, 0.049, 0.001)^(-1)))
qr.preds = predict(qr.model, list(density = seq(0.001, 0.049, 0.001)))
lines(seq(0.001, 0.049, 0.001), qr.preds, lwd=2, col = "red")
