nice_nrow_ncol <-
function(n, increasing=TRUE, maxsteps=20, maxratio=3){
  # Proposes ideal number of columns and rows for the total number of figures.
  #  | 1 | 2 |...|ncol|
  #  | 2 |   |   |    |
  #  |...|   |   |    |
  #  |nrow|  |   |    |
  #
  #  Total number of cells should be higher than the number of figures.

  sq = sqrt(n)
  sqlow = floor(sq)
  squpp = ceiling(sq)
  upper = c()
  lower = c()
  remainder = c()
  x = squpp
  y = squpp
  i = 0
  while((i < maxsteps) & (y >= 1) & (x/y <= maxratio)){
    while(y*x >= n){
      y = y-1
    }
    lower = c(lower, y+1)
    upper = c(upper, x)
    remainder = c(remainder, (y+1)*x-n)
    x = x+1
    i = i+1
  }
  
  # find the first pair [x, y] minimizing the remainder
  i = which.min(remainder)
  if(increasing){
    nrow = lower[i]
    ncol = upper[i]
  }else{
    nrow = upper[i]
    ncol = lower[i]
  }
  
  return(c(nrow, ncol))
}
