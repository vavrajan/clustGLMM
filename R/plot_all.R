plot_all <-
function(mcmc, scalar_samples, 
         what, gspec, dimspec, yspec,
         burnin = 0, thin = 1, lag.max = 30, iterations,
         chains, trueval, setparmfrow = TRUE,
         COL, move.width = 0.3, whatLAB){
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))
  if(setparmfrow){par(mfrow = c(2,2))}
  if(missing(chains)){chains = mcmc$chains}
  #if (missing(gspec)){gspec <- 1}
  #if (missing(yspec)){yspec <- ifelse(is.null(mcmc$Nums), mcmc$Ords[1], mcmc$Nums[1])}
  if(!missing(yspec)){
    if (!is.element(yspec, names(mcmc$family))){
      stop("Unexisting response variable specified.")
    }
  }
  if(missing(scalar_samples)){
    scalar_samples <- get_scalar_samples(mcmc = mcmc, what = what, 
                                         gspec = gspec, dimspec = dimspec, yspec = yspec,
                                         burnin = burnin, thin = thin, 
                                         iterations = iterations, chains = chains)
  }
  
  if(missing(whatLAB)){
    whatLAB <- what
    
    if(!missing(yspec)){
      whatLAB <- paste(whatLAB, yspec, sep = "_")
    }
    
    if(!missing(gspec)){
      whatLAB <- paste0(whatLAB, "(", gspec, ")")
    }
    
    if(missing(dimspec)){
      # no dimension supplied = no need for supplying -> 1-dim parameter
    }else{
      if(length(dimspec) == 1){
        whatLAB <- paste0(whatLAB,"[",dimspec[1],"]")
      }
      if(length(dimspec) == 2){
        whatLAB <- paste0(whatLAB,"[",dimspec[1],",",dimspec[2],"]")
      }
    }
  }
  
  plot_traceplots(mcmc=mcmc, scalar_samples = scalar_samples, 
                  what=what, gspec=gspec, dimspec=dimspec, yspec=yspec,
                  chains=chains, COL=COL, labcex = 0.8, whatLAB=whatLAB)
  if(!missing(trueval)){
    abline(h = trueval, col = "seagreen", lty = 2)
  }
  plot_ECDF(mcmc=mcmc, scalar_samples = scalar_samples,
            what=what, gspec=gspec, dimspec=dimspec, yspec=yspec,
            chains=chains, COL=COL, labcex = 0.8, whatLAB=whatLAB)
  if(!missing(trueval)){
    abline(v = trueval, col = "seagreen", lty = 2)
  }
  plot_kerneldensity(mcmc=mcmc, scalar_samples = scalar_samples,
                     what=what, gspec=gspec, dimspec=dimspec, yspec=yspec,
                     chains=chains, COL=COL, labcex = 0.8, whatLAB=whatLAB)
  if(!missing(trueval)){
    abline(v = trueval, col = "seagreen", lty = 2)
  }
  plot_ACF(mcmc=mcmc, scalar_samples = scalar_samples,
           what=what, gspec=gspec, dimspec=dimspec, yspec=yspec,
           chains=chains, COL=COL, labcex = 0.8, whatLAB=whatLAB,
           lag.max = lag.max, move.width = move.width)
}
