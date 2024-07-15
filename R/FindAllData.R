FindAllData <-
function(mcmc, what, gspec, dimspec, yspec, burnin = 0, thin = 1, ChainsToPlot,
         restrictedtoGplus = FALSE){
  if(thin < 1){
    stop("Thinning parameter has to be integer and thin >= 1.")
  }
  AllData <- data.frame()
  
  ### All is saved in structured list
  if(restrictedtoGplus){
    # we should work with data where spurious clusters were deleted and permutations were applied
    all <- NULL
    if(what == "w"){all <- mcmc$wGplus}
    if(what == "pUig"){all <- mcmc$pUigGplus}
    if(what == "U"){all <- mcmc$UGplus}
    if(is.element(what, rownames(mcmc$settings)[as.logical(mcmc$settings$isspec)])){all <- mcmc$clusterspec}
    
    if(is.null(all)){
      stop("Selected what parameter does not belong to the group of parameters that required restriction to non-spurious clusters.")
    }
    
    colname <- what
    
    if(!missing(yspec)){
      colname <- paste(colname, yspec, sep = "_")
    }
    
    if(!missing(gspec)){
      colname <- paste0(colname, "(", gspec, ")")
    }
    
    
    if(missing(dimspec)){
      # no dimension supplied = no need for supplying -> 1-dim parameter
    }else{
      if(length(dimspec) == 1){
        colname <- paste0(colname,"[",dimspec[1],"]")
      }
      if(length(dimspec) == 2){
        colname <- paste0(colname,"[",dimspec[1],",",dimspec[2],"]")
      }
    }
    
    if(!is.element(colname, colnames(all[[ChainsToPlot[1]]]))){
      stop(paste("Wrong specification, column ", colname, " does not exist in sampled data freed of spurious clusters!"))
    }
    
    AllData <- list()
    
    for(ch in ChainsToPlot){
      auxdata <- data.frame(y=unlist(all[[ch]][, colname]), 
                            chain = ch,
                            inds = mcmc$inds[[ch]])
      # burnin removal
      auxdata <- auxdata[auxdata$inds > burnin,]
      # thinning application
      auxdata <- auxdata[seq(1, dim(auxdata)[1], by = thin),]
      
      AllData[[ch]] <- auxdata
    }
    
    
  }else{
    states = seq(burnin+1, mcmc$iter, by = thin)
    # we work with the original data containing even spurious clusters  
    if(mcmc$howsave == "list"){
      
      v <- mcmc$settings[what,]
        
      for(ch in ChainsToPlot){
        if(v$isspec){
          if(v$isy){
            # group-specific and y-specific
            if(v$D == 0){
              AllData <- rbind(AllData, data.frame(y=mcmc$draws[[ch]][[what]][[gspec]][states, yspec], chain = ch))
            }
            if(v$D == 1){
              AllData <- rbind(AllData, data.frame(y=mcmc$draws[[ch]][[what]][[gspec]][[yspec]][states,dimspec[1]], chain = ch))
            }
            if(v$D == 2){
              AllData <- rbind(AllData, data.frame(y=mcmc$draws[[ch]][[what]][[gspec]][[yspec]][states,dimspec[1],dimspec[2]], chain = ch))
            }
          }else{
            # group-specific but not y-specific
            if(v$D == 0){
              AllData <- rbind(AllData, data.frame(y=mcmc$draws[[ch]][[what]][[gspec]][states], chain = ch))
            }
            if(v$D == 1){
              AllData <- rbind(AllData, data.frame(y=mcmc$draws[[ch]][[what]][[gspec]][states,dimspec[1]], chain = ch))
            }
            if(v$D == 2){
              AllData <- rbind(AllData, data.frame(y=mcmc$draws[[ch]][[what]][[gspec]][states,dimspec[1],dimspec[2]], chain = ch))
            }
          } # end of else of if is g-specific
        }else{
          if(v$isy){
            # not group-specific but y-specific
            if(v$D == 0){
              AllData <- rbind(AllData, data.frame(y=mcmc$draws[[ch]][[what]][states, yspec], chain = ch))
            }
            if(v$D == 1){
              AllData <- rbind(AllData, data.frame(y=mcmc$draws[[ch]][[what]][[yspec]][states,dimspec[1]], chain = ch))
            }
            if(v$D == 2){
              AllData <- rbind(AllData, data.frame(y=mcmc$draws[[ch]][[what]][[yspec]][states,dimspec[1],dimspec[2]], chain = ch))
            }
          }else{
            # neither group-specific nor y-specific
            if(v$D == 0){
              AllData <- rbind(AllData, data.frame(y=mcmc$draws[[ch]][[what]][states], chain = ch))
            }
            if(v$D == 1){
              AllData <- rbind(AllData, data.frame(y=mcmc$draws[[ch]][[what]][states,dimspec[1]], chain = ch))
            }
            if(v$D == 2){
              AllData <- rbind(AllData, data.frame(y=mcmc$draws[[ch]][[what]][states,dimspec[1],dimspec[2]], chain = ch))
            }
          }
        } 
      } # end of for ch in ChainsToPlot
    }# end of if mcmc$howsave == list
    
    ### All is saved in large matrix
    if(mcmc$howsave == "matrix"){
      #mcmc, burnin, iter, what, gspec, dimspec, yspec, thin, ChainsToPlot
      colname <- what
      
      if(!missing(yspec)){
        colname <- paste(colname, yspec, sep = "_")
      }
      
      if(!missing(gspec)){
        colname <- paste0(colname, "(", gspec, ")")
      }
      
      
      if(missing(dimspec)){
        # no dimension supplied = no need for supplying -> scalar parameter
      }else{
        if(length(dimspec) == 1){
          colname <- paste0(colname,"[",dimspec[1],"]")
        }
        if(length(dimspec) == 2){
          colname <- paste0(colname,"[",dimspec[1],",",dimspec[2],"]")
        }
      }
      
      if(!is.element(colname, colnames(mcmc$draws))){
        stop(paste("Wrong specification, column ", colname, " does not exist in mcmc$draws!"))
      }
      
      AllData <- data.frame()
      
      for(ch in ChainsToPlot){
        AllData <- rbind(AllData, data.frame(y=mcmc$draws[states + (ch-1)*mcmc$iter, colname], chain = ch))
      }
      
    } # end of if howsave="matrix"
  }
  return(AllData)
}
