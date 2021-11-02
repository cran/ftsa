semimetric.L2 <-
function(DATA1, DATA2)
{
  ###############################################################
  # Computes  L2 distance between curves.
  ################################################################
  #    "DATA1" matrix contains a first set of curves stored row by row
  #    "DATA2" matrix contains a second set of curves stored row by row
  # Returns a "semimetric" matrix containing the semimetrics
  # between all pairs (curve1, curve2) with curve1 in DATA1 and curve2
  # in DATA2 
  ###############################################################
  SEMIMETRIC = matrix(0, nrow(DATA1), nrow(DATA2))
  for(ii in 1:nrow(DATA2)){
    DIFF = t(DATA1) - DATA2[ii, ]
    SEMIMETRIC[, ii] = apply(DIFF^2, 2, sum) / (ncol(DATA2) - 1)
  }
  return(sqrt(SEMIMETRIC))
}
