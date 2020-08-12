rm(list = ls())
options(warn=-1)
setwd("~/Dropbox/academic/phd/gibbons/gibbon_simulation")

# run simulation and print summary statistics
source("./functions/tri_sim.R")

# specify parameters
tri_n_sims = 100
spacing = 400
true_density = 0.04
num_instances = 200
sigma = 300
sigma_step = 10
xsigma = seq(sigma, sigma + (num_instances - 1) * sigma_step, sigma_step)
max_dist = sigma * 3

output = matrix(0, nrow = num_instances, ncol = tri_n_sims)
for(i in 1:num_instances){
  output[i, ] = sapply(output[i, ], tri_sim, density = true_density, sigma = sigma, 
                       spacing = spacing, multi = TRUE, max_dist = max_dist, union = FALSE, 
                       edr = TRUE, true_sigma = sigma) / true_density
  sigma = sigma + sigma_step
  max_dist = sigma * 3
  print(i)
}

par(mfrow = c(1, 1))
# plot mean error ratios

plot(y = rowMeans(output), x = xsigma, xlab = "Sigma", ylab = "Estimate / True Density", 
     xlim = c(300, 300 + num_instances * sigma_step), ylim = c(0.8, 2.5), main = "Effect of sigma on bias with
     true_sigma used in computation of EDR (density 0.04)")

# fit polynomial regression model and plot
qr.model = lm(rowMeans(output) ~ xsigma + I(xsigma^2) + I(xsigma^3) + I(xsigma^4)+ I(xsigma^5))
qr.preds = predict(qr.model, list(density = xsigma))
lines(xsigma, qr.preds, lwd=2, col = "red")
abline(a = 1, b = 0, col = "BLUE")
