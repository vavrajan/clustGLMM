from_matrix_to_coda <- function(mcmc, burnin=0, thin=1){
  codamcmc <- list()
  # There is slight possibility that it will produce some NA values
  # We also omit n allocation indicators U
  for(ch in 1:mcmc$nchains){
    codamcmc[[ch]] <- mcmc$draws[mcmc$draws$chain == ch, 2:dim(mcmc$draws)[2]]
    codamcmc[[ch]] <- mcmc(codamcmc[[ch]][(burnin+1):mcmc$iter,], 
                           thin = thin, start = burnin+1, end = mcmc$iter)
  }
  codamcmc <- as.mcmc.list(codamcmc)
  return(codamcmc)
}