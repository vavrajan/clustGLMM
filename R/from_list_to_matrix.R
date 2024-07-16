from_list_to_matrix <- function(mcmc){
  params <- names(mcmc$save[mcmc$save])
  settings <- mcmc$settings
  yspecd1 <- mcmc$yspecd1
  yspecd2 <- mcmc$yspecd2
  family = mcmc$family
  draws = mcmc$draws
  cmcmc <- list()
  retmcmc <- mcmc
  nparams <- sum(mcmc$settings[mcmc$settings$save, "dimswithG"])
  retmcmc$draws <- matrix(0, nrow = mcmc$nchains*mcmc$iter, ncol = 1+nparams)
  retmcmc$draws[,1] <- rep(1:mcmc$nchains, each=mcmc$iter)
  colnames(retmcmc$draws) <- c("chain", unlist(mcmc$param_names))
  retmcmc$draws <- as.data.frame(retmcmc$draws)
  
  # insert the samples 
  for(ch in 1:mcmc$nchains){
    iters = ((ch-1)*mcmc$iter+1):(ch*mcmc$iter)
    
    for(p in params){
      if(settings[p,"ydepd1"] & settings[p,"ydepd2"]){
        # transfer from list to C
        cmcmc[[p]] <- from_list_to_C(chain = draws[[ch]],
                                     p = p,
                                     settings = settings,
                                     yspecd1 = yspecd1[[p]],
                                     yspecd2 = yspecd2[[p]],
                                     family = family)
        aux = from_C_to_matrix(values = cmcmc[[p]],
                               p = p,
                               settings = settings,
                               yspecd1 = yspecd1[[p]],
                               yspecd2 = yspecd2[[p]],
                               family = family)
      }
      if(settings[p,"ydepd1"] & !settings[p,"ydepd2"]){
        cmcmc[[p]] <- from_list_to_C(chain = draws[[ch]],
                                     p = p,
                                     settings = settings,
                                     yspecd1 = yspecd1[[p]],
                                     family = family)
        aux = from_C_to_matrix(values = cmcmc[[p]],
                               p = p,
                               settings = settings,
                               yspecd1 = yspecd1[[p]],
                               family = family)
      }
      if(!settings[p,"ydepd1"] & settings[p,"ydepd2"]){
        cmcmc[[p]] <- from_list_to_C(chain = draws[[ch]],
                                     p = p,
                                     settings = settings,
                                     yspecd2 = yspecd2[[p]],
                                     family = family)
        aux = from_C_to_matrix(values = cmcmc[[p]],
                               p = p,
                               settings = settings,
                               yspecd2 = yspecd2[[p]],
                               family = family)
      }
      if(!settings[p,"ydepd1"] & !settings[p,"ydepd2"]){
        cmcmc[[p]] <- from_list_to_C(chain = draws[[ch]],
                                     p = p,
                                     settings = settings,
                                     family = family)
        aux = from_C_to_matrix(values = cmcmc[[p]],
                               p = p,
                               settings = settings,
                               family = family)
      }
      retmcmc$draws[iters, colnames(aux)] <- aux
    }
    # chain ch completed
  }
  
  retmcmc$howsave <- "matrix"
  class(retmcmc) <- "clusglmm"
  
  return(retmcmc)
}
