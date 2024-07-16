from_C_to_list <-
function(values,
         p,
         settings,
         yspecd1,
         yspecd2,
         family)        # parameter to be listed
{
  # data frame settings consists of needed parameters:
  # save, isspec, G, isy, ynums, yords, ybins, iter, 
  # d1, d2, BYROW, sym, diag, diagval, D
  v <- c(settings[p, ])
  Ys <- c()
  if(v$ynums){Ys <- c(Ys, names(family)[is.element(family, c("num", "gaussian"))])}
  if(v$ypois){Ys <- c(Ys, names(family)[is.element(family, c("poi", "poisson"))])}
  if(v$ybins){Ys <- c(Ys, names(family)[is.element(family, c("bin", "bernoulli"))])}
  if(v$yords){Ys <- c(Ys, names(family)[is.element(family, c("ord", "cumulative"))])}
  if(v$ycats){Ys <- c(Ys, names(family)[is.element(family, c("cat", "categorical"))])}
  
  sumd1 <- ifelse(missing(yspecd1), v$d1*length(Ys), sum(yspecd1))
  #sumd1 <- v$d1*length(Ys)
  
  if(v$isspec){
    RET <- list()
    if(v$isy){
      
      if(v$D==0){
        for(k in 1:v$G){
          ## option with [[y]]
          #RET[[k]] <- list()
          #for(y in 1:length(Ys)){
          #RET[[k]][[Ys[y]]] <- values[seq(from = (k-1)*length(Ys)+y, 
          #                                to = G*length(Ys)*iter, 
          #                                by = G*length(Ys))]
          # }
          
          ## better --> no need for [[y]] --> can be stored in matrix with columns names Ys
          RET[[k]] <- matrix(0, nrow = v$iter, ncol = length(Ys))
          colnames(RET[[k]]) <- Ys
          for(y in 1:length(Ys)){
            RET[[k]][, Ys[y]] = values[seq(from = (k-1)*length(Ys)+y, 
                                           to = v$G*length(Ys)*v$iter, 
                                           by = v$G*length(Ys))]
          }
          
        }
      } # end of D==0
      
      if(v$D==1){
        for(k in 1:v$G){
          RET[[k]] <- list()
          dimy = 0
          for(y in 1:length(Ys)){
            d1 = ifelse(missing(yspecd1), v$d1, yspecd1[Ys[y]])
            #d1 = v$d1
            RET[[k]][[Ys[y]]] <- matrix(0, nrow = v$iter, ncol = d1)
            for(i in 1:v$iter){
              RET[[k]][[Ys[y]]][i, ] <- values[(i-1)*v$G*sumd1 + (k-1)*sumd1 + dimy + 1:d1]
            }
            dimy = dimy + d1
          }
        }
      } # end of D==1
      
      if(v$D==2){
        for(k in 1:v$G){
          RET[[k]] <- list()
          if(v$sym){
            if(v$diag){
              dimtot = ifelse(missing(yspecd1), length(Ys)*v$d1*(v$d1+1)/2,sum(yspecd1*(yspecd1+1)/2))
            }else{
              dimtot = ifelse(missing(yspecd1), length(Ys)*(v$d1-1)*v$d1/2,sum((yspecd1-1)*yspecd1/2))
            }
          }else{
            dimtot = ifelse(missing(yspecd1)|missing(yspecd2), length(Ys)*v$d1*v$d2, sum(yspecd1*yspecd2))
          }
          dimy = 0
          for(y in 1:length(Ys)){
            d1 = ifelse(missing(yspecd1), v$d1, yspecd1[Ys[y]])
            d2 = ifelse(missing(yspecd2), v$d2, yspecd2[Ys[y]])
            RET[[k]][[Ys[y]]] <- array(0, dim = c(v$iter, d1, d2))
            for(i in 1:v$iter){
              if(v$sym){
                # matrix is symmetrical and only upper-right triangle is stored
                # d1 = d2
                if(v$diag){
                  dd = d1*(d1+1)/2 # number of elements 
                  for(row in 1:d1){
                    for(col in row:d1){
                      RET[[k]][[Ys[y]]][i, row, col] <- RET[[k]][[Ys[y]]][i, col, row] <-
                        values[(i-1)*v$G*dimtot + (k-1)*dimtot + dimy + row + (col-1)*col/2]
                    }
                  }
                }else{
                  dd = (d1-1)*d1/2 # number of elements
                  # non-diagonal elements
                  for(row in 1:(d1-1)){
                    for(col in (row+1):d1){
                      RET[[k]][[Ys[y]]][i, row, col] <- RET[[k]][[Ys[y]]][i, col, row] <-
                        values[(i-1)*v$G*dimtot + (k-1)*dimtot + dimy + row + (col-2)*(col-1)/2]
                    }
                  }
                  # diagonal elements
                  for(row in 1:d1){
                    RET[[k]][[Ys[y]]][i, row, row] <- v$diagval
                  }
                } # end of else of diag
              }else{
                # matrix is general rectangular --> all is stored
                dd = d1*d2
                RET[[k]][[Ys[y]]][i,,] <- matrix(values[(i-1)*v$G*dimtot + (k-1)*dimtot + dimy + 1:dd],
                                                 nrow = d1,
                                                 ncol = d2,
                                                 byrow = v$BYROW)
              } # end of else of sym
            }
            dimy = dimy + dd
          }
        }
      } # end of D==2
      
    }else{
      # is class-specific, however not y-specific
      
      if(v$D==0){
        RET <- matrix(0, nrow = v$iter, ncol = v$G)
        for(k in 1:v$G){
          # no need to do [[k]] --> can be stored in matrix
          # column corresponds to k
          RET[, k] = values[seq(from = k,
                                to = v$G*v$iter,
                                by = v$G)]
          
        }
      } # end of D==0
      
      if(v$D==1){
        for(k in 1:v$G){
          RET[[k]] <- matrix(0, nrow = v$iter, ncol = v$d1)
          for(i in 1:v$iter){
            RET[[k]][i, ] <- values[(i-1)*v$G*v$d1 + (k-1)*v$d1 + 1:v$d1]          
          }
        }
      } # end of D==1
      
      if(v$D==2){
        for(k in 1:v$G){
          RET[[k]] <- array(0, dim = c(v$iter, v$d1, v$d2))
          for(i in 1:v$iter){
            if(v$sym){
              # matrix is symmetrical and only upper-right triangle is stored
              # d1 = d2
              if(v$diag){
                dd = v$d1*(v$d1+1)/2 # number of elements 
                for(row in 1:v$d1){
                  for(col in row:v$d1){
                    RET[[k]][i, row, col] <- RET[[k]][i, col, row] <-
                      values[(i-1)*v$G*dd + (k-1)*dd + row + (col-1)*col/2]
                  }
                }
              }else{
                dd = (v$d1-1)*v$d1/2 # number of elements
                # non-diagonal elements
                for(row in 1:(v$d1-1)){
                  for(col in (row+1):v$d1){
                    RET[[k]][i, row, col] <- RET[[k]][i, col, row] <-
                      values[(i-1)*v$G*dd + (k-1)*dd + row + (col-2)*(col-1)/2]
                  }
                }
                # diagonal elements
                for(row in 1:v$d1){
                  RET[[k]][i, row, row] <- v$diagval
                }
              } # end of else of diag
            }else{
              # matrix is general rectangular --> all is stored
              RET[[k]][i,,] <- matrix(values[(i-1)*v$G*v$d1*v$d2 + (k-1)*v$d1*v$d2 + 1:(v$d1*v$d2)],
                                      nrow = v$d1,
                                      ncol = v$d2,
                                      byrow = v$BYROW)
            } # end of else of sym
          }
        }
      } # end of D==2
      
    } # end of else of isy
  }else{
    # not class-specific
    if(v$isy){
      # not class-specific, but still y-specific
      if(v$D==0){
        ## option with [[y]]
        #RET <- list()
        #for(y in 1:length(Ys)){
        #RET[[Ys[y]]] <- values[seq(from = y, 
        #                                to = length(Ys)*iter, 
        #                                by = length(Ys))]
        # }
        
        ## better --> no need for [[y]] --> can be stored in matrix with columns names Ys
        RET <- matrix(0, nrow = v$iter, ncol = length(Ys))
        colnames(RET) <- Ys
        for(y in 1:length(Ys)){
          RET[, Ys[y]] = values[seq(from = y,
                                    to = length(Ys)*v$iter,
                                    by = length(Ys))]
        }
        
      } # end of D==0
      
      if(v$D==1){
        RET <- list()
        dimtot = ifelse(missing(yspecd1), length(Ys)*v$d1, sum(yspecd1))
        dimy = 0
        for(y in 1:length(Ys)){
          d1 = ifelse(missing(yspecd1), v$d1, yspecd1[Ys[y]])
          RET[[Ys[y]]] <- matrix(0, nrow = v$iter, ncol = d1)
          for(i in 1:v$iter){
            RET[[Ys[y]]][i, ] <- values[(i-1)*dimtot + dimy + 1:d1]
          }
          dimy = dimy + d1
        }
      } # end of D==1
      
      if(v$D==2){
        RET <- list()
        if(v$sym){
          if(v$diag){
            dimtot = ifelse(missing(yspecd1), length(Ys)*v$d1*(v$d1+1)/2, sum(yspecd1*(yspecd1+1)/2))
          }else{
            dimtot = ifelse(missing(yspecd1), length(Ys)*(v$d1-1)*v$d1/2, sum((yspecd1-1)*yspecd1/2))
          }
        }else{
          dimtot = ifelse(missing(yspecd1)|missing(yspecd2), length(Ys)*v$d1*v$d2, sum(yspecd1*yspecd2))
        }
        dimy = 0
        for(y in 1:length(Ys)){
          d1 = ifelse(missing(yspecd1), v$d1, yspecd1[Ys[y]])
          d2 = ifelse(missing(yspecd2), v$d2, yspecd2[Ys[y]])
          RET[[Ys[y]]] <- array(0, dim = c(v$iter, d1, d2))
          for(i in 1:v$iter){
            if(v$sym){
              # matrix is symmetrical and only upper-right triangle is stored
              # d1 = d2
              if(v$diag){
                dd = d1*(d1+1)/2 # number of elements 
                for(row in 1:d1){
                  for(col in row:d1){
                    RET[[Ys[y]]][i, row, col] <- RET[[Ys[y]]][i, col, row] <-
                      values[(i-1)*dimtot + dimy + row + (col-1)*col/2]
                  }
                }
              }else{
                dd = (d1-1)*d1/2 # number of elements
                # non-diagonal elements
                for(row in 1:(d1-1)){
                  for(col in (row+1):d1){
                    RET[[Ys[y]]][i, row, col] <- RET[[Ys[y]]][i, col, row] <-
                      values[(i-1)*dimtot  + dimy + row + (col-2)*(col-1)/2]
                  }
                }
                # diagonal elements
                for(row in 1:d1){
                  RET[[Ys[y]]][i, row, row] <- v$diagval
                }
              } # end of else of diag
            }else{
              # matrix is general rectangular --> all is stored
              dd = d1 * d2
              RET[[Ys[y]]][i,,] <- matrix(values[(i-1)*dimtot + dimy + 1:dd],
                                          nrow = d1,
                                          ncol = d2,
                                          byrow = v$BYROW)
            } # end of else of sym
          }
          dimy = dimy + dd
        }
      } # end of D==2
      
    }else{
      # is NEITHER class-specific, nor y-specific
      
      if(v$D==0){
        RET <- values
      } # end of D==0
      
      if(v$D==1){
        RET <-  matrix(values, nrow = v$iter, ncol = v$d1, byrow = v$BYROW)
      } # end of D==1
      
      if(v$D==2){
        RET <- array(0, dim = c(v$iter, v$d1, v$d2))
        for(i in 1:v$iter){
          if(v$sym){
            # matrix is symmetrical and only upper-right triangle is stored
            # d1 = d2
            if(v$diag){
              dd = v$d1*(v$d1+1)/2 # number of elements 
              for(row in 1:v$d1){
                for(col in row:v$d1){
                  RET[i, row, col] <- RET[i, col, row] <-
                    values[(i-1)*dd + row + (col-1)*col/2]
                }
              }
            }else{
              dd = (v$d1-1)*v$d1/2 # number of elements
              # non-diagonal elements
              for(row in 1:(v$d1-1)){
                for(col in (row+1):v$d1){
                  RET[i, row, col] <- RET[i, col, row] <-
                    values[(i-1)*dd + row + (col-2)*(col-1)/2]
                }
              }
              # diagonal elements
              for(row in 1:v$d1){
                RET[i, row, row] <- v$diagval
              }
            } # end of else of diag
          }else{
            # matrix is general rectangular --> all is stored
            RET[i,,] <- matrix(values[(i-1)*v$d1*v$d2 + 1:(v$d1*v$d2)],
                               nrow = v$d1,
                               ncol = v$d2,
                               byrow = v$BYROW)
          } # end of else of sym
        }
      } # end of D==2
    } # end of else of isy
  } # end of else of isspec
  
  return(RET)
}
