rm(list = ls())
options(warn=-1)
setwd("~/Dropbox/academic/phd/gibbons/gibbon_simulation")

require(secr)
source("./functions/likelihood.R")
source("./functions/edr.R")
source("./functions/circle.R")
source("./functions/generate_data.R")
source("./functions/effective_area.R")

n_sims = 500
true_density = 0.04
spacing = 400


sigma = 300
max_dist = sigma * 3
est_sigma = rep(NA, n_sims)
true_sigma = rep(NA, n_sims)
for(i in 1:n_sims){
  traps = generate_traps(spacing = spacing)
  distances = generate_distances(traps, true_density, sigma, spacing, multi = TRUE)
  e = edr(distances, max_dist)
  est_sigma[i] = attr(e, "par")
  true_sigma[i] = sigma
  sigma = sigma + 4
  max_dist = sigma * 3
}

par(mfrow = c(1, 1))
plot(true_sigma, est_sigma, xlim = c(300, 2500), ylim = c(100, 4300), main = "True sigma vs estimated sigma used in EDR")
abline(a = 0, b = 1)

# fit polynomial regression line
qr.model = lm(est_sigma ~ true_sigma + I(true_sigma^2))
qr.preds = predict(qr.model, list(density = true_sigma))
lines(true_sigma, qr.preds, lwd=2, col = "red")
