from_list_to_matrix <- function(mcmc){
  params <- names(mcmc$save[mcmc$save])
  draws = mcmc$draws
  cmcmc <- list()
  retmcmc <- mcmc
  nparams <- sum(mcmc$settings[mcmc$settings$save, "dimswithG"])
  
  # insert the samples 
  for(ch in 1:mcmc$nchains){
    iters = draws[[ch]]$m
    retmcmc$draws[[ch]] <- matrix(0, nrow = length(iters), ncol = 1+nparams)
    retmcmc$draws[[ch]][,1] <- iters
    colnames(retmcmc$draws[[ch]]) <- c("m", unlist(mcmc$param_names))
    retmcmc$draws[[ch]] <- as.data.frame(retmcmc$draws[[ch]])
  
    for(p in params){
      #if(settings[p,"ydepd1"] & settings[p,"ydepd2"]){
        # transfer from list to C
        cmcmc[[p]] <- from_list_to_C(draws = draws[[ch]],
                                     p = p,
                                     settings = mcmc$settings,
                                     yspecd1 = mcmc$yspecd1[[p]],
                                     yspecd2 = mcmc$yspecd2[[p]],
                                     family = mcmc$family)
        aux = from_C_to_matrix(values = cmcmc[[p]],
                               p = p,
                               settings = mcmc$settings,
                               yspecd1 = mcmc$yspecd1[[p]],
                               yspecd2 = mcmc$yspecd2[[p]],
                               family = mcmc$family)
      # }
      # if(settings[p,"ydepd1"] & !settings[p,"ydepd2"]){
      #   cmcmc[[p]] <- from_list_to_C(chain = draws[[ch]],
      #                                p = p,
      #                                settings = settings,
      #                                yspecd1 = yspecd1[[p]],
      #                                family = family)
      #   aux = from_C_to_matrix(values = cmcmc[[p]],
      #                          p = p,
      #                          settings = settings,
      #                          yspecd1 = yspecd1[[p]],
      #                          family = family)
      # }
      # if(!settings[p,"ydepd1"] & settings[p,"ydepd2"]){
      #   cmcmc[[p]] <- from_list_to_C(chain = draws[[ch]],
      #                                p = p,
      #                                settings = settings,
      #                                yspecd2 = yspecd2[[p]],
      #                                family = family)
      #   aux = from_C_to_matrix(values = cmcmc[[p]],
      #                          p = p,
      #                          settings = settings,
      #                          yspecd2 = yspecd2[[p]],
      #                          family = family)
      # }
      # if(!settings[p,"ydepd1"] & !settings[p,"ydepd2"]){
      #   cmcmc[[p]] <- from_list_to_C(chain = draws[[ch]],
      #                                p = p,
      #                                settings = settings,
      #                                family = family)
      #   aux = from_C_to_matrix(values = cmcmc[[p]],
      #                          p = p,
      #                          settings = settings,
      #                          family = family)
      # }
      retmcmc$draws[[ch]][, colnames(aux)] <- aux
    }
    # chain ch completed
  }
  
  retmcmc$howsave <- "matrix"
  class(retmcmc) <- "clustglmm"
  
  return(retmcmc)
}
