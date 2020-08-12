rm(list=setdiff(ls(), c("edout", "esout", "edout_union")))
options(warn=-1)
setwd("~/Dropbox/academic/phd/gibbons/gibbon_simulation")

# run simulation and print summary statistics
source("./functions/tri_sim.R")

# specify parameters
tri_n_sims = 100
sigma = 500
spacing = 400
max_dist = sigma * 3
true_density = 0.005
num_instances = 100

output = matrix(0, nrow = num_instances, ncol = tri_n_sims)
for(i in 1:num_instances){
  output[i, ] = sapply(output[i, ], tri_sim, density = true_density, sigma = sigma, 
                  spacing = spacing, multi = TRUE, max_dist = max_dist, union = FALSE, 
                  edr_val = TRUE, spatial = 'straush') / true_density
  true_density = true_density + 0.0005
  print(i)
}

# plot mean error ratios
plot(y = rowMeans(output), x = seq(0.001, 0.0505, 0.0005), xlab = "Density / Hectare", ylab = "Estimate / True Density", 
     xlim = c(0, 0.05), ylim = c(0, 6))

# fit reciprocal quadratic model and plot
# qr.model = lm(rowMeans(output) ~ I(seq(0.001, 0.0505, 0.0005)^(-1)))
# qr.preds = predict(qr.model, list(density = seq(0.001, 0.0505, 0.0005)))
# lines(seq(0.001, 0.0505, 0.0005), qr.preds, lwd=2, col = "red")
lines(seq(0.001, 0.0505, 0.0005), rowMeans(output)-apply(output, 1, sd), col = "BLUE")
lines(seq(0.001, 0.0505, 0.0005), rowMeans(output)+apply(output, 1, sd), col = "BLUE")

edout = output