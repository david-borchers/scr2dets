# simulation to estimate gibbon density via trangulation method
require(secr)
source("./functions/likelihood.R")
source("./functions/edr.R")
source("./functions/circle.R")
source("./functions/generate_data.R")
source("./functions/effective_area.R")
source("./functions/refine.R")

tri_sim = function(x, density, sigma, spacing, multi, max_dist, union, true_sigma = FALSE, spatial = FALSE, edr_val = TRUE){
  output = NA
  # avoid NaN output in extreme data cases
  while(is.nan(output) | is.na(output)){
    traps = generate_traps(spacing = spacing)
    distances = generate_distances(traps, true_density, sigma, spacing, multi, spatial = spatial)
    area = effective_area(distances, traps, max_dist, union = union, true_sigma = true_sigma, edr_val = edr_val)

    valid_detections = refine(distances, attr(area, 'edr'))
    if(is.vector(valid_detections)){
      output = 10000 / area
    } else{
      output = 10000 * length(valid_detections[1, ]) / area
    }
  }
  return(output)
}
