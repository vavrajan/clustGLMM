plot_diagnostics <-
function(mcmc, scalar_samples, 
         what, gspec, dimspec, yspec,
         burnin = 0, thin = 1, lag.max = 30, iterations,
         chains, trueval, setparmfrow = TRUE,
         COL, move.width = 0.3, whatLAB){
  if(setparmfrow){
      opar <- par(mfrow = c(2, 2))
      on.exit(par(opar))
  }
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
  
  if(missing(whatLAB)){whatLAB <- scalar_samples$whatLAB}
  if(what == "b"){
    awhat <- sub(",\\d", "", whatLAB)
  }else{
    awhat <- sub("\\[[0-9]+(,[0-9]+)?\\]", "", whatLAB)
    awhat <- sub("^c_ord", "a_ord", awhat)
    awhat <- sub("^pi_ord", "a_ord", awhat)
  }
  accept <- numeric(length(mcmc$chains))
  for(ich in 1:length(mcmc$chains)){
    if(is.element(awhat, colnames(mcmc$accept[[mcmc$chains[ich]]]))){
      accept[ich] <- 100*mean(mcmc$accept[[mcmc$chains[ich]]][,awhat], na.rm = TRUE)
    }else{
      accept[ich] <- NA
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
  
  if(missing(COL)){
    COL <- diverge_hcl(mcmc$nchains, c = 80, l = 70)
  }
  
  if(any(!is.na(accept))){
    legend("topright", legend = paste0(mcmc$chains, ": ", format(accept, digits = 4), "%"), 
           col = COL[mcmc$chains], 
           lty = 1, title = "Chain: Acceptance", cex = 0.8,
           y.intersp = 0.6, x.intersp = 0.5)
  }else{
    legend("topright", legend = mcmc$chains, 
           col = COL[mcmc$chains], 
           lty = 1, title = "Chain", cex = 0.8,
           y.intersp = 0.6, x.intersp = 0.5)
  }
  
}
