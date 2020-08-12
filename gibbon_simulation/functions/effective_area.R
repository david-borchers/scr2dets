effective_area = function(distances, traps, max_dist, union, true_sigma = FALSE, edr_val = TRUE){
  # trap locations
  x1 = traps[1, 1]
  y1 = traps[1, 2]
  x2 = traps[2, 1]
  y2 = traps[2, 2]
  x3 = traps[3, 1]
  y3 = traps[3, 2]
  
  # effective distance radius
  if(edr_val == TRUE){
    r1 = edr(distances[1, ], max_dist, true_sigma)
    r2 = edr(distances[2, ], max_dist, true_sigma)
    r3 = edr(distances[3, ], max_dist, true_sigma)
  } else {
    r1 = r2 = r3 = edr_val
  }
  if(union == FALSE){
    ea = intersect_area(x1, y1, r1, x2, y2, r2) +
      intersect_area(x2, y2, r2, x3, y3, r3) + 
      intersect_area(x1, y1, r1, x3, y3, r3) -
      2 * common_area(x1, y1, r1, x2, y2, r2, x3, y3, r3)
    
  } else {
    ea = pi * (r1^2 + r2^2 + r3^2) -
      intersect_area(x1, y1, r1, x2, y2, r2) -
      intersect_area(x2, y2, r2, x3, y3, r3) -
      intersect_area(x1, y1, r1, x3, y3, r3) +
      common_area(x1, y1, r1, x2, y2, r2, x3, y3, r3)
  }
  attr(ea, 'edr') = list('t1' = r1, 't2' = r2, 't3' = r3)
  return(ea)
}
