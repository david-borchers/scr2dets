
require(secr)
source("./functions/generate_data.R")

# simulation to estimate gibbon density via SECR
secr_sim = function(x, density, sigma, spacing, max_dist, spatial = FALSE, noccasions = 1){
  traps = generate_traps(spacing = spacing)
  capt_hist = capture_history(traps, density, sigma, spacing, spatial = spatial, noccasions = noccasions)
  fitted_model = secr.fit(capt_hist, fixed = list(g0 = 1), detectfn = 0, buffer = max_dist, trace = FALSE)
  return(exp(fitted_model$fit$estimate)[1])
}
