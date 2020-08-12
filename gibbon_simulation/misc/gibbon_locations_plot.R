source('~/Dropbox/academic/phd/gibbons/gibbon_simulation/functions/generate_data.R')

# parameters 
sigma = 600
spacing = 400
max_dist = 3 * sigma
density = 0.03

# generate traps
traps = generate_traps(spacing)

# generate locations
loc_poiss = attr(capture_history(traps, density, sigma, spacing, spatial = FALSE), 'loc')
loc_hardc = attr(capture_history(traps, density, sigma, spacing, spatial = 'hardcore'), 'loc')
loc_strau = attr(capture_history(traps, density, sigma, spacing, spatial = 'straush'), 'loc')

# plots
par(mfrow = c(1, 1))


plot(loc_poiss, main = "Poisson Point Process: Density = 0.03 / Hectare", col = "BLUE", pch = 4)
plot(traps, add = TRUE)

plot(loc_hardc, main = "Hardcore Point Process: Density = 0.03 / Hectare", col = "BLUE", pch = 4)
plot(traps, add = TRUE)

plot(loc_strau, main = "Hardcore Strauss Process: Density = 0.03 / Hectare", col = "BLUE", pch = 4)
plot(traps, add = TRUE)
