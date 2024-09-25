from_matrix_to_coda <- function(mcmc, burnin=0, thin=1){
  codamcmc <- list()
  # There is slight possibility that it will produce some NA values
  # We also omit n allocation indicators U
  for(ch in 1:mcmc$nchains){
    codamcmc[[ch]] <- mcmc$draws[[ch]][mcmc$draws[[ch]]$m > burnin, ]
    codamcmc[[ch]] <- codamcmc[[ch]][seq(1, dim(codamcmc[[ch]])[1], by = thin),-1,drop=FALSE]
    if(dim(mcmc$draws[[ch]])[1] == mcmc$iter){
      # We work with original chain of length mcmc$iter
      codamcmc[[ch]] <- mcmc(codamcmc[[ch]], 
                             thin = thin, start = burnin+1, end = mcmc$iter)
    }else{
      # We work with a subset of a chain (restriction to Gplus, clustering probabilities, ...)
      # Start numbering from 1 and omit any information about thinning 
      codamcmc[[ch]] <- mcmc(codamcmc[[ch]], 
                             thin = 1, start = 1, end = dim(mcmc$draws[[ch]])[1])
    }
  }
  codamcmc <- as.mcmc.list(codamcmc)
  return(codamcmc)
}