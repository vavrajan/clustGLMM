FromMatrixtoList <- function(mcmc){
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

  for(ch in 1:mcmc$chains){
    cmcmc <- list()
    chain <- list()
    for(p in params){
      ### First from Matrix to C
      cmcmc[[p]] <- FromMatrixtoC(draws, chain = ch, iterations, p, settings)
      
      ### From C to List
      if(settings[p, "save"] & settings[p,"dims"] > 0){
        if(settings[p,"ydepd1"] & settings[p,"ydepd2"]){
          chain[[p]] <- FromCtoList(values = cmcmc[[p]],
                                    p = p,
                                    settings = settings,
                                    yspecd1 = yspecd1[[p]],
                                    yspecd2 = yspecd2[[p]],
                                    family = family)
        }
        if(settings[p,"ydepd1"] & !settings[p,"ydepd2"]){
          chain[[p]] <- FromCtoList(values = cmcmc[[p]],
                                    p = p,
                                    settings = settings,
                                    yspecd1 = yspecd1[[p]],
                                    family = family)
        }
        if(!settings[p,"ydepd1"] & settings[p,"ydepd2"]){
          chain[[p]] <- FromCtoList(values = cmcmc[[p]],
                                    p = p,
                                    settings = settings,
                                    yspecd2 = yspecd2[[p]],
                                    family = family)
        }
        if(!settings[p,"ydepd1"] & !settings[p,"ydepd2"]){
          chain[[p]] <- FromCtoList(values = cmcmc[[p]],
                                    p = p,
                                    settings = settings,
                                    family = family)
        }
      }
    }
    retmcmc$draws[[ch]] <- chain
  }
  
  retmcmc$howsave <- "list"
  class(retmcmc) <- "mbcglmms"
  
  return(retmcmc)

}
