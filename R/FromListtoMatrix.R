FromListtoMatrix <- function(mcmc){
  params <- names(mcmc$save[mcmc$save])
  settings <- mcmc$settings
  yspecd1 <- mcmc$yspecd1
  yspecd2 <- mcmc$yspecd2
  family = mcmc$family
  draws = mcmc$draws
  cmcmc <- list()
  retmcmc <- mcmc
  retmcmc$draws <- data.frame()
  
  for(ch in 1:mcmc$chains){
    AllData <- matrix(ch, nrow = mcmc$iter, ncol = 1)
    colnames(AllData) <- "chain"
    for(p in params){
      if(settings[p,"ydepd1"] & settings[p,"ydepd2"]){
        # transfer from list to C
        cmcmc[[p]] <- FromListtoC(chain = draws[[ch]],
                                  p = p,
                                  settings = settings,
                                  yspecd1 = yspecd1[[p]],
                                  yspecd2 = yspecd2[[p]],
                                  family = family)
        AllData <- cbind(AllData, FromCtoMatrix(values = cmcmc[[p]],
                                                p = p,
                                                settings = settings,
                                                yspecd1 = yspecd1[[p]],
                                                yspecd2 = yspecd2[[p]],
                                                family = family))
      }
      if(settings[p,"ydepd1"] & !settings[p,"ydepd2"]){
        cmcmc[[p]] <- FromListtoC(chain = draws[[ch]],
                                  p = p,
                                  settings = settings,
                                  yspecd1 = yspecd1[[p]],
                                  family = family)
        AllData <- cbind(AllData, FromCtoMatrix(values = cmcmc[[p]],
                                                p = p,
                                                settings = settings,
                                                yspecd1 = yspecd1[[p]],
                                                family = family))
      }
      if(!settings[p,"ydepd1"] & settings[p,"ydepd2"]){
        cmcmc[[p]] <- FromListtoC(chain = draws[[ch]],
                                  p = p,
                                  settings = settings,
                                  yspecd2 = yspecd2[[p]],
                                  family = family)
        AllData <- cbind(AllData, FromCtoMatrix(values = cmcmc[[p]],
                                                p = p,
                                                settings = settings,
                                                yspecd2 = yspecd2[[p]],
                                                family = family))
      }
      if(!settings[p,"ydepd1"] & !settings[p,"ydepd2"]){
        cmcmc[[p]] <- FromListtoC(chain = draws[[ch]],
                                  p = p,
                                  settings = settings,
                                  family = family)
        AllData <- cbind(AllData, FromCtoMatrix(values = cmcmc[[p]],
                                                p = p,
                                                settings = settings,
                                                family = family))
      }
    }
    # chain ch completed
    retmcmc$draws <- rbind(retmcmc$draws, AllData)
  }
  
  retmcmc$howsave <- "matrix"
  class(retmcmc) <- "mbcglmms"
  
  return(retmcmc)
}
