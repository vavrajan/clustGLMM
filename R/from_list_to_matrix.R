from_list_to_matrix <- function(mcmc){
  params <- names(mcmc$save[mcmc$save])
  draws = mcmc$draws
  cmcmc <- list()
  retmcmc <- mcmc
  nparams <- sum(mcmc$settings[mcmc$settings$save, "dimswithG"])
  
  # insert the samples 
  for(ch in mcmc$chains){
    iters = draws[[ch]]$m
    retmcmc$draws[[ch]] <- matrix(0, nrow = length(iters), ncol = 1+nparams)
    retmcmc$draws[[ch]][,1] <- iters
    colnames(retmcmc$draws[[ch]]) <- c("m", unlist(mcmc$param_names))
    retmcmc$draws[[ch]] <- as.data.frame(retmcmc$draws[[ch]])
  
    for(p in params){
      # transfer from list to C
      cmcmc[[p]] <- from_list_to_C(draws = draws[[ch]],
                                   p = p,
                                   settings = mcmc$settings,
                                   yspecd1 = mcmc$yspecd1[[p]],
                                   yspecd2 = mcmc$yspecd2[[p]],
                                   d2spec = mcmc$d2spec[[p]],
                                   family = mcmc$family)
      # transfer from C to matrix
      aux = from_C_to_matrix(values = cmcmc[[p]],
                             p = p,
                             settings = mcmc$settings,
                             yspecd1 = mcmc$yspecd1[[p]],
                             yspecd2 = mcmc$yspecd2[[p]],
                             d2spec = mcmc$d2spec[[p]],
                             family = mcmc$family)
      retmcmc$draws[[ch]][, colnames(aux)] <- aux
    }
    # chain ch completed
  }
  
  retmcmc$howsave <- "matrix"
  class(retmcmc) <- "clustglmm"
  
  return(retmcmc)
}
