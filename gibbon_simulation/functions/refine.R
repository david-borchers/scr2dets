refine = function(distances, radials){
  edr1 = radials$t1
  edr2 = radials$t2
  edr3 = radials$t3
  remove_cols = c()
  for(i in 1:ncol(distances)){
    count = 0
    if(distances[1, i] > edr1){
      count = count + 1
    }
    if(distances[2, i] > edr2){
      count = count + 1
    }
    if(distances[3, i] > edr3){
      count = count + 1
    }
    if(count == 3){
      remove_cols = c(remove_cols, i)
    }
  }
  if(length(remove_cols) > 0){
    return(distances[, -remove_cols])
  } else{
    return(distances)
  }
}
