plot_kerneldensity_param <-
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
          for(g in 1:v$G){
            par(mfrow = c(v$d1, v$d2), mar = c(4,4,1,1))
            for(i in 1:v$d1){
              for(j in 1:v$d2){
                if(i <= j){
                  plot_kerneldensity(mcmc = mcmc, what = what, 
                                     gspec = g, dimspec = c(i,j), 
                                     burnin = burnin, thin = thin, ChainsToPlot = ChainsToPlot,
                                     labcex = ifelse(v$d1 <= 3, 0.6, 0.5))
                }else{
                  pomdata <- FindAllData(mcmc = mcmc, what = what, 
                                         gspec = g, dimspec = c(j,i),
                                         burnin = burnin, thin = thin, ChainsToPlot = ChainsToPlot)
                  plot(x=c(0,1), y=c(0,1), type = "n", xlab = "", ylab = "",
                       xaxt = "n", yaxt = "n", bty = "n")
                  text(0.5,0.5, labels = format(mean(pomdata$y), digits = 2, nsmall = 2))
                }
              }
            }
          }
        }else{
          # symmetric matrix with diagonal is not group-specific
          par(mfrow = c(v$d1, v$d2), mar = c(4,4,1,1))
          for(i in 1:v$d1){
            for(j in 1:v$d2){
              if(i <= j){
                plot_kerneldensity(mcmc = mcmc, what = what, 
                                   dimspec = c(i,j),
                                   burnin = burnin, thin = thin, ChainsToPlot = ChainsToPlot,
                                   labcex = ifelse(v$d1 <= 3, 0.6, 0.5))
              }else{
                pomdata <- FindAllData(mcmc = mcmc, what = what, 
                                       dimspec = c(j,i),
                                       burnin = burnin, thin = thin, ChainsToPlot = ChainsToPlot)
                plot(x=c(0,1), y=c(0,1), type = "n", xlab = "", ylab = "",
                     xaxt = "n", yaxt = "n", bty = "n")
                text(0.5,0.5, labels = format(mean(pomdata$y), digits = 2, nsmall = 2))
              }
            }
          }
        }
      }else{
        if(v$isspec){
          # symmetric matrix withOUT diagonal is group-specific
          for(g in 1:v$G){
            par(mfrow = c(v$d1-1, v$d2-1), mar = c(4,4,1,1))
            for(i in 1:(v$d1-1)){
              for(j in 2:v$d2){
                if(i < j){
                  plot_kerneldensity(mcmc = mcmc, what = what, 
                                     gspec = g, dimspec = c(i,j), 
                                     burnin = burnin, thin = thin, ChainsToPlot = ChainsToPlot,
                                     labcex = ifelse(v$d1 <= 3, 0.6, 0.5))
                }else{
                  pomdata <- FindAllData(mcmc = mcmc, what = what, 
                                         gspec = g, dimspec = c(j,i),
                                         burnin = burnin, thin = thin, ChainsToPlot = ChainsToPlot)
                  plot(x=c(0,1), y=c(0,1), type = "n", xlab = "", ylab = "",
                       xaxt = "n", yaxt = "n", bty = "n")
                  text(0.5,0.5, labels = format(mean(pomdata$y), digits = 2, nsmall = 2))
                }
              }
            }
          }
        }else{
          # symmetric matrix withOUT diagonal is not group-specific
          par(mfrow = c(v$d1, v$d2), mar = c(4,4,1,1))
          for(i in 1:v$d1){
            for(j in 1:v$d2){
              if(i < j){
                plot_kerneldensity(mcmc = mcmc, what = what, 
                                   dimspec = c(i,j), 
                                   burnin = burnin, thin = thin, ChainsToPlot = ChainsToPlot,
                                   labcex = ifelse(v$d1 <= 3, 0.6, 0.5))
              }else{
                plot(x=c(0,1), y=c(0,1), type = "n", xlab = "", ylab = "",
                     xaxt = "n", yaxt = "n", bty = "n")
                if(i == j){
                  text(0.5,0.5, labels = "1")
                }else{
                  pomdata <- FindAllData(mcmc = mcmc, what = what, 
                                         dimspec = c(j,i),
                                         burnin = burnin, thin = thin, ChainsToPlot = ChainsToPlot)
                  text(0.5,0.5, labels = format(mean(pomdata$y), digits = 2, nsmall = 2))
                }
              }
            }
          }
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
            for(ch in ChainsToPlot){
              par(mfrow = c(length(nonempty[[ch]]), length(Ys)), mar = c(4,4,1,1))
              for(g in nonempty[[ch]]){
                for(y in Ys){
                  plot_kerneldensity(mcmc = mcmc, what = what,
                                     gspec = g, yspec = y,
                                     burnin = burnin, thin = thin, ChainsToPlot = ch, 
                                     labcex = 0.6)
                }
              }
            }
          }else{
            # not y-specific, but still g-specific
            for(ch in ChainsToPlot){
              par(mfrow = nice_nrow_ncol(length(nonempty[[ch]])), mar = c(4,4,1,1))
              for(g in nonempty[[ch]]){
                plot_kerneldensity(mcmc = mcmc, what = what,
                                   gspec = g,
                                   burnin = burnin, thin = thin, ChainsToPlot = ch, 
                                   labcex = 0.6)
              }
            }
          }
        }
        
        if(v$D == 1){
          # one-dimensional parameter
          if(v$ydepd1){
            # the size of the dimension depends on outcomes
            for(ch in ChainsToPlot){
              for(y in Ys){
                par(mfrow = c(length(nonempty[[ch]]),mcmc$yspecd1[[what]][[y]]), 
                    mar = c(4,4,1,1))
                for(g in nonempty[[ch]]){
                  for(i in 1:mcmc$yspecd1[[what]][[y]]){
                    plot_kerneldensity(mcmc = mcmc, what = what,
                                       gspec = g, yspec = y, dimspec = i,
                                       burnin = burnin, thin = thin, ChainsToPlot = ch, 
                                       labcex = 0.6)
                  }
                }
              }
            }
          }else{
            # the size of the dimension depends on outcomes
            for(ch in ChainsToPlot){
              par(mfrow = c(length(nonempty[[ch]]), v$d1), mar = c(4,4,1,1))
              for(g in nonempty[[ch]]){
                for(i in 1:v$d1){
                  plot_kerneldensity(mcmc = mcmc, what = what,
                                     gspec = g, dimspec = i,
                                     burnin = burnin, thin = thin, ChainsToPlot = ch, 
                                     labcex = 0.6)
                }
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
                for(g in nonempty[[ch]]){
                  par(mfrow = c(d1,d2), mar = c(4,4,1,1))
                  for(i in 1:d1){
                    for(j in 1:d2){
                      plot_kerneldensity(mcmc = mcmc, what = what,
                                         gspec = g, yspec = y, dimspec = c(i,j),
                                         burnin = burnin, thin = thin, ChainsToPlot = ch, 
                                         labcex = 0.6)
                    }
                  }
                }
              }
            }
          }else{
            # none dimension is y-specific
            for(ch in ChainsToPlot){
              for(g in nonempty[[ch]]){
                par(mfrow = c(v$d1, v$d2), mar = c(4,4,1,1))
                for(i in 1:v$d1){
                  for(j in 1:v$d2){
                    plot_kerneldensity(mcmc = mcmc, what = what,
                                       gspec = g, dimspec = c(i,j),
                                       burnin = burnin, thin = thin, ChainsToPlot = ch, 
                                       labcex = 0.6)
                  }
                }
              }
            }
          }
        }
        
      }else{
        # not group-specific parameter
        if(v$D == 0){
          # zero-dimensional parameter
          if(v$isy){
            # y-specific
            par(mfrow = nice_nrow_ncol(length(Ys)), mar = c(4,4,1,1))
            for(y in Ys){
              plot_kerneldensity(mcmc = mcmc, what = what,
                                 yspec = y,
                                 burnin = burnin, thin = thin, ChainsToPlot = ChainsToPlot, 
                                 labcex = 0.6)
            }
          }else{
            # not y-specific, neither g-specific
            par(mfrow = c(1,1), mar = c(4,4,1,1))
            plot_kerneldensity(mcmc = mcmc, what = what,
                               burnin = burnin, thin = thin, ChainsToPlot = ChainsToPlot,
                               labcex = 0.6)
          }
        }
        
        if(v$D == 1){
          # one-dimensional parameter
          if(v$ydepd1){
            # the size of the dimension depends on outcomes
            for(y in Ys){
              par(mfrow = nice_nrow_ncol(mcmc$yspecd1[[what]][[y]]), 
                  mar = c(4,4,1,1))
              for(i in 1:mcmc$yspecd1[[what]][[y]]){
                plot_kerneldensity(mcmc = mcmc, what = what,
                                   yspec = y, dimspec = i,
                                   burnin = burnin, thin = thin, ChainsToPlot = ChainsToPlot, 
                                   labcex = 0.6)
              }
            }
          }else{
            # the size of the dimension depends on outcomes
            par(mfrow = nice_nrow_ncol(v$d1), mar = c(4,4,1,1))
            for(i in 1:v$d1){
              plot_kerneldensity(mcmc = mcmc, what = what,
                                 dimspec = i,
                                 burnin = burnin, thin = thin, ChainsToPlot = ChainsToPlot, 
                                 labcex = 0.6)
            }
          }
        }
        
        if(v$D == 2){
          # two-dimensional parameter
          if((v$ydepd1) | (v$ydepd2)){
            # some size of the dimension depends on outcomes
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
                  plot_kerneldensity(mcmc = mcmc, what = what,
                                     yspec = y, dimspec = c(i,j),
                                     burnin = burnin, thin = thin, ChainsToPlot = ChainsToPlot, 
                                     labcex = 0.6)
                }
              }
            }
          }else{
            # none dimension is y-specific
            par(mfrow = c(v$d1, v$d2), mar = c(4,4,1,1))
            for(i in 1:v$d1){
              for(j in 1:v$d2){
                plot_kerneldensity(mcmc = mcmc, what = what,
                                   dimspec = c(i,j),
                                   burnin = burnin, thin = thin, ChainsToPlot = ChainsToPlot, 
                                   labcex = 0.6)
              }
            }
          }
        }
      }
    }
  }else{
    warning("There is nothing to be plotted.")
  }
}

