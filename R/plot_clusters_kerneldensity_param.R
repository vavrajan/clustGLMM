plot_clusters_kerneldensity_param <-
  function(mcmc, what, burnin = 0, thin = 1, ChainsToPlot){
    opar <- par(no.readonly = TRUE)
    on.exit(par(opar))
    v <- c(mcmc$settings[what, ])
    Ys <- c()
    if(v$ynums){Ys <- c(Ys, mcmc$Nums)}
    if(v$ypois){Ys <- c(Ys, mcmc$Pois)}
    if(v$ybins){Ys <- c(Ys, mcmc$Bins)}
    if(v$yords){Ys <- c(Ys, mcmc$Ords)}
    if(v$ycats){Ys <- c(Ys, mcmc$Cats)}
    
    if(missing(ChainsToPlot)){ChainsToPlot = 1:mcmc$chains}
    nonempty <- list()
    for(ch in 1:mcmc$chains){
      TAB <- table(mcmc$last[[ch]]$U)
      nonempty[[ch]] <- sort(as.numeric(names(TAB)[which(TAB>0)]))
    }
    
    if(v$dimswithG > 0){  
      if(v$sym){
        if(v$diag){
          if(v$isspec){
            # symmetric matrix with diagonal is group-specific
            for(ch in ChainsToPlot){
              par(mfrow = c(v$d1, v$d2), mar = c(4,4,1,1))
              for(i in 1:v$d1){
                for(j in 1:v$d2){
                  plot_clusters(mcmc = mcmc, what = what, 
                                dimspec = c(i,j), whichg = nonempty[[ch]], 
                                burnin = burnin, thin = thin, ChainsToPlot = ch,
                                setparmfrow = FALSE, doKern = TRUE, doECDF = FALSE, 
                                labcex = ifelse(v$d1 <= 3, 0.6, 0.5))
                }
              }
            }
          }else{
            # symmetric matrix with diagonal is not group-specific
            warning(paste0("Parameter ", what, " is not group-specific."))
          }
        }else{
          if(v$isspec){
            # symmetric matrix withOUT diagonal is group-specific
            for(ch in ChainsToPlot){
              par(mfrow = c(v$d1-1, v$d2-1), mar = c(4,4,1,1))
              for(i in 1:(v$d1-1)){
                for(j in 2:v$d2){
                  plot_clusters(mcmc = mcmc, what = what, 
                                dimspec = c(i,j), whichg = nonempty[[ch]], 
                                burnin = burnin, thin = thin, ChainsToPlot = ch,
                                setparmfrow = FALSE, doKern = TRUE, doECDF = FALSE, 
                                labcex = ifelse(v$d1 <= 3, 0.6, 0.5))
                }
              }
            }
          }else{
            # symmetric matrix withOUT diagonal is not group-specific
            warning(paste0("Parameter ", what, " is not group-specific."))
          }
        }
      }else{
        # parameter is not a symmetric matrix
        if(v$isspec){
          # group-specific parameter
          if(v$D == 0){
            # zero-dimensional parameter
            if(v$isy){
              # y-specific 
              par(mfrow = c(length(ChainsToPlot), length(Ys)), mar = c(4,4,1,1))
              for(ch in ChainsToPlot){
                for(y in Ys){
                  plot_clusters(mcmc = mcmc, what = what, 
                                yspec = y, whichg = nonempty[[ch]],
                                burnin = burnin, thin = thin, ChainsToPlot = ch,
                                setparmfrow = FALSE, doKern = TRUE, doECDF = FALSE, 
                                labcex = 0.6)
                }
              }
            }else{
              # not y-specific, but still g-specific
              par(mfrow = nice_nrow_ncol(length(ChainsToPlot)), mar = c(4,4,1,1))
              for(ch in ChainsToPlot){
                plot_clusters(mcmc = mcmc, what = what,
                              whichg = nonempty[[ch]],
                              burnin = burnin, thin = thin, ChainsToPlot = ch,
                              setparmfrow = FALSE, doKern = TRUE, doECDF = FALSE, 
                              labcex = 0.6)
              }
            }
          }
          
          if(v$D == 1){
            # one-dimensional parameter
            if(v$ydepd1){
              # the size of the dimension depends on outcomes
              for(ch in ChainsToPlot){
                for(y in Ys){
                  par(mfrow = nice_nrow_ncol(mcmc$yspecd1[[what]][[y]]), 
                      mar = c(4,4,1,1))
                  for(i in 1:mcmc$yspecd1[[what]][[y]]){
                    plot_clusters(mcmc = mcmc, what = what, 
                                  yspec = y, dimspec = i, whichg = nonempty[[ch]],
                                  burnin = burnin, thin = thin, ChainsToPlot = ch,
                                  setparmfrow = FALSE, doKern = TRUE, doECDF = FALSE, 
                                  labcex = 0.6)
                  }
                }
              }
            }else{
              # the size of the dimension depends on outcomes
              par(mfrow = c(length(ChainsToPlot), v$d1), mar = c(4,4,1,1))
              for(ch in ChainsToPlot){
                for(i in 1:v$d1){
                  plot_clusters(mcmc = mcmc, what = what, 
                                dimspec = i, whichg = nonempty[[ch]],
                                burnin = burnin, thin = thin, ChainsToPlot = ch,
                                setparmfrow = FALSE, doKern = TRUE, doECDF = FALSE, 
                                labcex = 0.6)
                }
              }
            }
          }
          
          if(v$D == 2){
            # two-dimensional parameter
            if((v$ydepd1) | (v$ydepd2)){
              # some size of the dimension depends on outcomes
              for(ch in ChainsToPlot){
                for(y in Ys){
                  if(v$ydepd1){
                    d1 = mcmc$yspecd1[[what]][[y]]
                  }else{
                    d1 = v$d1
                  }
                  if(v$ydepd2){
                    d2 = mcmc$yspecd2[[what]][[y]]
                  }else{
                    d2 = v$d2
                  }
                  par(mfrow = c(d1,d2), mar = c(4,4,1,1))
                  for(i in 1:d1){
                    for(j in 1:d2){
                      plot_clusters(mcmc = mcmc, what = what, 
                                    yspec = y, dimspec = c(i,j), whichg = nonempty[[ch]],
                                    burnin = burnin, thin = thin, ChainsToPlot = ch,
                                    setparmfrow = FALSE, doKern = TRUE, doECDF = FALSE, 
                                    labcex = 0.6)
                    }
                  }
                }
              }
            }else{
              # none dimension is y-specific
              for(ch in ChainsToPlot){
                par(mfrow = c(v$d1, v$d2), mar = c(4,4,1,1))
                for(i in 1:v$d1){
                  for(j in 1:v$d2){
                    plot_clusters(mcmc = mcmc, what = what, 
                                  dimspec = c(i,j), whichg = nonempty[[ch]],
                                  burnin = burnin, thin = thin, ChainsToPlot = ch,
                                  setparmfrow = FALSE, doKern = TRUE, doECDF = FALSE, 
                                  labcex = 0.6)
                  }
                }
              }
            }
          }
          
        }else{
          # parameter is not group-specific 
          warning(paste0("Parameter ", what, " is not group-specific."))
        }
      }
    }else{
      warning("There is nothing to be traceplotted.")
    }
  }
