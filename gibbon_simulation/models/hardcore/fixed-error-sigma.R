rm(list = ls())
options(warn=-1)
setwd("~/Dropbox/academic/phd/gibbons/gibbon_simulation")

# run simulation and print summary statistics
source("./functions/tri_sim.R")

# specify parameters
tri_n_sims = 100
sigma = 300
spacing = 400
max_dist = sigma * 3
true_density = 0.04
num_instances = 200

output = matrix(0, nrow = num_instances, ncol = tri_n_sims)
for(i in 1:num_instances){
  output[i, ] = sapply(output[i, ], tri_sim, density = true_density, sigma = sigma, 
                       spacing = spacing, multi = TRUE, max_dist = max_dist, union = FALSE, 
                       edr = 1000) / true_density
  sigma = sigma + 10
  max_dist = sigma * 3
  print(i)
}

par(mfrow = c(2, 1))
# plot mean error ratios
plot(y = rowMeans(output), x = seq(400, 2390, 10), xlab = "Sigma", ylab = "Estimate / True Density", 
     xlim = c(390, 2500), ylim = c(-1,10), main = "Effect of sigma on bias with fixed edr of 1000m")

# fit polynomial regression model and plot
qr.model = lm(rowMeans(output) ~ seq(400, 2390, 10) + I(seq(400, 2390, 10)^2) + I(seq(400, 2390, 10)^3))
qr.preds = predict(qr.model, list(density = seq(400, 2390, 10)))
lines(seq(400, 2390, 10), qr.preds, lwd=2, col = "red")

plot(y = rowMeans(output), x = seq(400, 2390, 10), xlab = "Sigma", ylab = "Estimate / True Density", 
     xlim = c(500, 1000), ylim = c(0,2), main = "Effect of sigma on bias with fixed edr of 1000m")
abline(a = 1, b = 0, col = "BLUE")
