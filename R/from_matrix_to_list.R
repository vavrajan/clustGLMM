from_matrix_to_list <- function(mcmc){
  params <- names(mcmc$save[mcmc$save])
  settings <- mcmc$settings
  yspecd1 <- mcmc$yspecd1
  yspecd2 <- mcmc$yspecd2
  family <- mcmc$family
  draws <- mcmc$draws
  cmcmc <- list()
  iterations <- 1:mcmc$iter
  retmcmc <- mcmc
  retmcmc$draws <- list()

  for(ch in 1:mcmc$nchains){
    cmcmc <- list()
    chain <- list()
    for(p in params){
      ### First from matrix to C
      cmcmc[[p]] <- from_matrix_to_C(draws, chain = ch, iterations, p, settings)
      
      ### From C to list
      if(settings[p, "save"] & settings[p,"dims"] > 0){
        if(settings[p,"ydepd1"] & settings[p,"ydepd2"]){
          chain[[p]] <- from_C_to_list(values = cmcmc[[p]],
                                       p = p,
                                       settings = settings,
                                       yspecd1 = yspecd1[[p]],
                                       yspecd2 = yspecd2[[p]],
                                       family = family)
        }
        if(settings[p,"ydepd1"] & !settings[p,"ydepd2"]){
          chain[[p]] <- from_C_to_list(values = cmcmc[[p]],
                                       p = p,
                                       settings = settings,
                                       yspecd1 = yspecd1[[p]],
                                       family = family)
        }
        if(!settings[p,"ydepd1"] & settings[p,"ydepd2"]){
          chain[[p]] <- from_C_to_list(values = cmcmc[[p]],
                                       p = p,
                                       settings = settings,
                                       yspecd2 = yspecd2[[p]],
                                       family = family)
        }
        if(!settings[p,"ydepd1"] & !settings[p,"ydepd2"]){
          chain[[p]] <- from_C_to_list(values = cmcmc[[p]],
                                       p = p,
                                       settings = settings,
                                       family = family)
        }
      }
    }
    retmcmc$draws[[ch]] <- chain
  }
  
  retmcmc$howsave <- "list"
  class(retmcmc) <- "clusglmm"
  
  return(retmcmc)

}
