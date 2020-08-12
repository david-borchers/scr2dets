library(ascr)

# set parameters
spacing = 400
sigma = 500
max_dist = sigma * 3
nocc = 5

######################## script #############################

nsims = 100
spacings = seq(sigma, 3*sigma, 53)
out = array(dim = c(nsims, 1, length(spacings)))
for(i in 1:nsims){
  xrange = yrange = spacing * 30
  traps = make.circle(n = 3, spacing = spacings[i], detector = "proximity",  originxy = c(xrange / 2, yrange / 2))
    
  pop = sim.popn(D = 0.02, expand.grid(x = c(0, xrange), y = c(0, xrange)), 0)
    
  #window = owin(xrange = c(0, xrange), yrange = c(0, xrange))
  #strausshard = list( beta = 1e-5, gamma=0.5, r = 700, hc = 400)
  #sim_location = rmh(model = list(cif = 'straush', par = strausshard, w = window),
  #                   start = list(n.start = nrow(pop)),
  #                  control = list(p = 1), verbose = FALSE)
  #pop$x <- sim_location$x
  #pop$y <- sim_location$y
    
  ch = sim.capthist(traps, popn = pop, noccasions = nocc, detectpar = list(g0 = 1, sigma = sigma), 
                        savepopn = TRUE, renumber = FALSE, p.available = 1)
  msk = create.mask(traps, max_dist)
  ch = sim.capt(fit = NULL, traps = traps, mask = msk,
                infotypes = character(0), detfn = "hn", pars = NULL, ss.opts = NULL,
                cue.rates = NULL, survey.length = NULL, freq.dist = "edf",
                sound.speed = 330, test.detfn = FALSE, first.only = FALSE,
                keep.locs = FALSE, keep.ids = FALSE, ...)
      
 
  fitted_model = fit.ascr(ch, traps, msk)
  est = exp(fitted_model$fit$estimate)[1]
  out[i] = est
}

    print(i)
  }
}

out = out[1:100,,]
plot(apply(out, 3, mean))
