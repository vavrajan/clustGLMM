plot_ng_trace_chain_split <- function(mcmc, 
                                      burnin = 0, thin = 1, iterations, chains,
                                      setparmfrow = TRUE){
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))
  if(missing(chains)){chains = mcmc$chains}
  if(setparmfrow){par(mfrow = nice_nrow_ncol(length(chains)))}
  
  COL <- diverge_hcl(mcmc$nchains, c = 80, l = 70)
  ngdraws <- iters <- list()
  for(ch in chains){
    if(missing(iterations)){
      iters[[ch]] <- mcmc$draws[[ch]]$m
      iters[[ch]] <- iters[[ch]][iters[[ch]] > burnin]
      iters[[ch]] <- iters[[ch]][seq(1, length(iters[[ch]]), by=thin)]
    }else{
      iters[[ch]] <- iterations
    }
    if(mcmc$howsave == "matrix"){
      ngdraws[[ch]] <- mcmc$draws[[ch]][iterations,paste0("ng(", 1:mcmc$G[ch], ")")]
    }
    if(mcmc$howsave == "list"){
      ngdraws[[ch]] <- mcmc$draws[[ch]]$ng[iterations,]
    }
  }
  
  MIN <- min(sapply(chains, 
                    function(ch){min(ngdraws[[ch]], na.rm = TRUE)}))
  MIN <- max(sapply(chains, 
                    function(ch){max(ngdraws[[ch]], na.rm = TRUE)}))
  
  par(mfrow = nice_nrow_ncol(length(chains)), mar = c(4,4,0.5,0.5))
  for(ch in chains){
    plot(0, 0, type = "n", 
         xlim = range(iters[[ch]]), ylim = c(MIN, MAX),
         xlab = "Iteration", ylab = "Cluster occupancy number ng")
    for(g in 1:mcmc$G[ch]){
      lines(ngdraws[[ch]][,g] ~ iters[[ch]],
            lty = 1, lwd = 1, col = COL[ch])
    }
  }
}