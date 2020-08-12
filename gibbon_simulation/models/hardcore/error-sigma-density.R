rm(list = ls())
options(warn=-1)
setwd("~/Dropbox/academic/phd/gibbons/gibbon_simulation")

# run simulation and print summary statistics
source("./functions/tri_sim.R")

# custom parameters
tri_n_sims = 100
sigma = 400
spacing = 400
max_dist = sigma * 3
true_density = 0.01
num_instances = 100
num_dens = 3

output = array(0, dim = c(num_instances, tri_n_sims, num_dens))
               
for(j in 1:num_dens){
  sigma = 400
  for(i in 1:num_instances){
    output[i, , j] = sapply(output[i, , j], tri_sim, density = true_density, sigma = sigma, 
                         spacing = spacing, multi = TRUE, max_dist = max_dist, union = FALSE, 
                         edr_val = TRUE) / true_density
    sigma = sigma + 10
    max_dist = sigma * 3
    print(i)
  }
  true_density = true_density + 0.01
}


# plot mean error ratios
par(mfrow = c(2, 2))
plot(y = rowMeans(output[,,1]), x = seq(400, 1390, 10), xlab = "Sigma", ylab = "Estimate / True Density", 
     xlim = c(390, 1400), ylim = c(0.8, 6.5), main = "Density = 0.01")
plot(y = rowMeans(output[,,2]), x = seq(400, 1390, 10), xlab = "Sigma", ylab = "Estimate / True Density", 
       xlim = c(390, 1400), ylim = c(0.8, 4), main = "Density = 0.02", col = "RED")
plot(y = rowMeans(output[,,3]), x = seq(400, 1390, 10), xlab = "Sigma", ylab = "Estimate / True Density", 
       xlim = c(390, 1400), ylim = c(0.8, 4), main = "Density = 0.03", col = "BLUE")
plot(y = rowMeans(output[,,1]), x = seq(400, 1390, 10), xlab = "Sigma", ylab = "Estimate / True Density", 
     xlim = c(390, 1400), ylim = c(0.8, 6.5), main = "Superimposed")
points(y = rowMeans(output[,,2]), x = seq(400, 1390, 10), xlab = "Sigma", ylab = "Estimate / True Density", 
     xlim = c(390, 1400), ylim = c(0.8, 4), col = "Red")
points(y = rowMeans(output[,,3]), x = seq(400, 1390, 10), xlab = "Sigma", ylab = "Estimate / True Density", 
     xlim = c(390, 1400), ylim = c(0.8, 4), col = "Blue")

# fit model and compare fit across densities through R-squared statistic
qr.model = lm(rowMeans(output[,, 1]) ~ seq(400, 1390, 10) + I(seq(400, 1390, 10)^2) + I(seq(400, 1390, 10)^3)
              + I(seq(400, 1390, 10)^4))
qr.preds = predict(qr.model, list(density = seq(400, 1390, 10)))
rs1 = 1 - (sum((output[,,1] - qr.preds)^2)/sum((output[,,1] - mean(output[,, 1]))^2))
rs2 = 1 - (sum((output[,,2] - qr.preds)^2)/sum((output[,,2] - mean(output[,, 2]))^2))
rs3 = 1 - (sum((output[,,3] - qr.preds)^2)/sum((output[,,3] - mean(output[,, 3]))^2))
print(c(rs1, rs2, rs3))
