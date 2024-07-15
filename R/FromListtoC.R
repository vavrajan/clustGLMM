FromListtoC <-
function(chain,
         iterations,
         p,
         settings,
         yspecd1,
         yspecd2,
         family)
{
  v <- c(settings[p, ])
  Ys <- c()
  if(v$ynums){Ys <- c(Ys, names(family)[is.element(family, c("num", "gaussian"))])}
  if(v$ypois){Ys <- c(Ys, names(family)[is.element(family, c("poi", "poisson"))])}
  if(v$ybins){Ys <- c(Ys, names(family)[is.element(family, c("bin", "bernoulli"))])}
  if(v$yords){Ys <- c(Ys, names(family)[is.element(family, c("ord", "cumulative"))])}
  if(v$ycats){Ys <- c(Ys, names(family)[is.element(family, c("cat", "categorical"))])}
  
  if(missing(iterations)){iterations <- 1:v$iter}
  
  sumd1 <- ifelse(missing(yspecd1), ifelse(v$d1==0,1,v$d1)*length(Ys), sum(yspecd1))
  RET <- c()
  
  
  if(v$isspec){
    if(v$isy){
      if(length(Ys) > 0){
        if(v$D==0){
          for(i in iterations){
            for(k in 1:v$G){
              for(y in 1:length(Ys)){
                RET <- c(RET, as.data.frame(chain[[p]][[k]][i, Ys[y]]))
              }
            }
          }
        } # end of D==0
        
        if(v$D==1){
          for(i in iterations){
            for(k in 1:v$G){
              for(y in 1:length(Ys)){
                RET <- c(RET, chain[[p]][[k]][[Ys[y]]][i, ])
              }
            }
          }
        } # end of D==1
        
        if(v$D==2){
          for(i in iterations){
            for(k in 1:v$G){
              for(y in 1:length(Ys)){
                d1 = ifelse(missing(yspecd1), v$d1, yspecd1[Ys[y]])
                d2 = ifelse(missing(yspecd2), v$d2, yspecd2[Ys[y]])
                if(v$sym){
                  # matrix is symmetrical and only upper-right triangle is stored
                  # d1 = d2
                  if(v$diag){
                    for(col in 1:d1){
                      for(row in 1:col){
                        RET <- c(RET, chain[[p]][[k]][[Ys[y]]][i, row, col])
                      }
                    }
                  }else{
                    for(col in 2:d1){
                      for(row in 1:(col-1)){
                        RET <- c(RET, chain[[p]][[k]][[Ys[y]]][i, row, col])
                      }
                    }
                  } # end of else of diag
                }else{
                  # matrix is general rectangular --> all is stored
                  if(v$BYROW){
                    RET <- c(RET, t(chain[[p]][[k]][[Ys[y]]][i,,]))
                  }else{
                    RET <- c(RET, chain[[p]][[k]][[Ys[y]]][i,,])
                  }
                } # end of else of sym
              }
            }
          }
        } # end of D==2
      }
    }else{
      # is class-specific, however not y-specific

      if(v$D==0){
        RET <- c(t(chain[[p]][iterations, ]))
      } # end of D==0
      
      if(v$D==1){
        for(i in iterations){
          for(k in 1:v$G){
            RET <- c(RET, chain[[p]][[k]][i,])
          }
        }
      } # end of D==1
      
      if(v$D==2){
        for(i in iterations){
          for(k in 1:v$G){
            if(v$sym){
              # matrix is symmetrical and only upper-right triangle is stored
              # d1 = d2
              if(v$diag){
                for(col in 1:v$d1){
                  for(row in 1:col){
                    RET <- c(RET, chain[[p]][[k]][i, row, col])
                  }
                }
              }else{
                for(col in 2:v$d1){
                  for(row in 1:(col-1)){
                    RET <- c(RET, chain[[p]][[k]][i, row, col])
                  }
                }
              } # end of else of diag
            }else{
              # matrix is general rectangular --> all is stored
              if(v$BYROW){
                RET <- c(RET, t(chain[[p]][[k]][i,,]))
              }else{
                RET <- c(RET, chain[[p]][[k]][i,,])
              }
            } # end of else of sym
          }
        }
      } # end of D==2
      
    } # end of else of isy
  }else{
    # not class-specific
    if(v$isy){
      if(length(Ys) > 0){
        # not class-specific, but still y-specific
        if(v$D==0){
          RET <- c(t(chain[[p]]))
        } # end of D==0
      
        if(v$D==1){
          for(i in iterations){
            for(y in 1:length(Ys)){
              RET <- c(RET, chain[[p]][[Ys[y]]][i, ])
            }
          }
        } # end of D==1
        
        if(v$D==2){
          for(i in iterations){
            for(y in 1:length(Ys)){
              d1 = ifelse(missing(yspecd1), v$d1, yspecd1[Ys[y]])
              d2 = ifelse(missing(yspecd2), v$d2, yspecd2[Ys[y]])
              if(v$sym){
                # matrix is symmetrical and only upper-right triangle is stored
                # d1 = d2
                if(v$diag){
                  for(col in 1:d1){
                    for(row in 1:col){
                      RET <- c(RET, chain[[p]][[Ys[y]]][i, row, col])
                    }
                  }
                }else{
                  for(col in 2:d1){
                    for(row in 1:(col-1)){
                      RET <- c(RET, chain[[p]][[Ys[y]]][i, row, col])
                    }
                  }
                } # end of else of diag
              }else{
                # matrix is general rectangular --> all is stored
                if(v$BYROW){
                  RET <- c(RET, t(chain[[p]][[Ys[y]]][i,,]))
                }else{
                  RET <- c(RET, chain[[p]][[Ys[y]]][i,,])
                }
              } # end of else of sym
            }
          }
        } # end of D==2
      }
    }else{
      # is NEITHER class-specific, nor y-specific
      
      if(v$D==0){
        RET <- chain[[p]][iterations]
      } # end of D==0
      
      if(v$D==1){
        RET <- c(t(chain[[p]][iterations,]))
      } # end of D==1
      
      if(v$D==2){
        for(i in iterations){
          if(v$sym){
            # matrix is symmetrical and only upper-right triangle is stored
            # d1 = d2
            if(v$diag){
              for(col in 1:v$d1){
                for(row in 1:col){
                  RET <- c(RET, chain[[p]][i, row, col])
                }
              }
            }else{
              for(col in 2:v$d1){
                for(row in 1:(col-1)){
                  RET <- c(RET, chain[[p]][i, row, col])
                }
              }
            } # end of else of diag
          }else{
            # matrix is general rectangular --> all is stored
            if(v$BYROW){
              RET <- c(RET, t(chain[[p]][i,,]))
            }else{
              RET <- c(RET, chain[[p]][i,,])
            }
          } # end of else of sym
        }
      } # end of D==2
    } # end of else of isy
  } # end of else of isspec
  
  return(RET)

}
