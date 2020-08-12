# simulation to estimate gibbon density via trangulation method
setwd("~/Dropbox/academic/phd/gibbons/gibbon_simulation/")
require(secr)
require(ascr)
require(spatstat)
require(CircStats)
source("./functions/mod.sim.capt.R")
source("./functions/generate_data.R")

ascr_sim = function(x, true_density, sigma, spacing){
    # generate traps and mask/survey area
    t = as.matrix(generate_traps(spacing = spacing))
    m = create.mask(t, buffer = 5*sigma, nx = 64, ny = 64)
    
    # generate population
    popn = sim.popn(D = true_density, buffer = 0, core = data.frame(x = m[, 1], y = m[, 2]))
    strausshard = list( beta = 1e-5, gamma=0.5, r = 700, hc = 400)
    window = owin(xrange = range(m[, 1]), yrange = range(m[, 2]))
    sim_location = rmh(model = list(cif = 'straush', par = strausshard, w = window),
                       start = list(n.start = nrow(popn)),
                       control = list(p = 1), verbose = FALSE)
    popn$x <- sim_location$x
    popn$y <- sim_location$y
    popn = as.matrix(popn)
    
    # simulate capture history --- choose which kappa?
    output = NA
    ch = try(repel.sim.capt(traps = t, mask = m, popn = popn,
                          infotypes = c("bearing"), detfn = "hn",
                          pars = list(D = true_density, g0 = 1, sigma = sigma, kappa = 50)), silent = TRUE)
    test = try(coef(fit.ascr(ch, traps = t, mask = m, fix = list(g0 = 1)), "D"), silent = TRUE)
    if(is.numeric(test)){
      output = test
    } else {
      output = 0
    }
    print(x)
    output
}

Ds <- sapply(1:100, ascr_sim, true_density = 0.02, sigma = 500, spacing = 400)
boxplot(Ds)

