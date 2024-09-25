get_scalar_samples <-
function(mcmc, what, gspec, dimspec, yspec, 
         burnin = 0, thin = 1, 
         iterations, chains,
         restrictedtoGplus = FALSE){
  if(thin < 1){
    stop("Thinning parameter has to be integer and thin >= 1.")
  }
  if(missing(chains)){chains = mcmc$chains}
  scalar_samples <- data.frame()
  misiter <- missing(iterations)
  
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
    
    if(!is.element(colname, colnames(all[[chains[1]]]))){
      stop(paste("Wrong specification, column ", colname, " does not exist in sampled data freed of spurious clusters!"))
    }
    
    scalar_samples <- list()
    
    for(ch in chains){
      auxdata <- data.frame(y=unlist(all[[ch]][, colname]), 
                            chain = ch,
                            m = mcmc$inds[[ch]])
      # burnin removal
      auxdata <- auxdata[auxdata$m > burnin,]
      # thinning application
      auxdata <- auxdata[seq(1, dim(auxdata)[1], by = thin),]
      
      scalar_samples[[ch]] <- auxdata
    }
    
    
  }else{
    # we work with the original data containing even spurious clusters
    scalar_samples <- list()
    for(ch in chains){
      if(misiter){
        # iterations = seq(burnin+1, mcmc$iter, by = thin)
        iterations <- mcmc$draws[[ch]]$m
        iterations <- iterations[iterations > burnin]
        iterations <- iterations[seq(1, length(iterations), by=thin)]
      }
      
      ind_iterations <- is.element(mcmc$draws[[ch]]$m, iterations)    
      v <- mcmc$settings[what,]
      
      ### Draws are save in a structured list
      out <- data.frame(chain = ch, m = iterations)
      if(mcmc$howsave == "list"){
        if(v$isspec){
          if(v$isy){
            # group-specific and y-specific
            if(v$D == 0){
              out$y = mcmc$draws[[ch]][[what]][[gspec]][ind_iterations, yspec]
            }
            if(v$D == 1){
              out$y = mcmc$draws[[ch]][[what]][[gspec]][[yspec]][ind_iterations,dimspec[1]]
            }
            if(v$D == 2){
              out$y = mcmc$draws[[ch]][[what]][[gspec]][[yspec]][ind_iterations,dimspec[1],dimspec[2]]
            }
          }else{
            # group-specific but not y-specific
            if(v$D == 0){
              out$y = mcmc$draws[[ch]][[what]][[gspec]][ind_iterations]
            }
            if(v$D == 1){
              out$y = mcmc$draws[[ch]][[what]][[gspec]][ind_iterations,dimspec[1]]
            }
            if(v$D == 2){
              out$y = mcmc$draws[[ch]][[what]][[gspec]][ind_iterations,dimspec[1],dimspec[2]]
            }
          } # end of else of if is g-specific
        }else{
          if(v$isy){
            # not group-specific but y-specific
            if(v$D == 0){
              out$y = mcmc$draws[[ch]][[what]][ind_iterations, yspec]
            }
            if(v$D == 1){
              out$y = mcmc$draws[[ch]][[what]][[yspec]][ind_iterations,dimspec[1]]
            }
            if(v$D == 2){
              out$y = mcmc$draws[[ch]][[what]][[yspec]][ind_iterations,dimspec[1],dimspec[2]]
            }
          }else{
            # neither group-specific nor y-specific
            if(v$D == 0){
              out$y = mcmc$draws[[ch]][[what]][ind_iterations]
            }
            if(v$D == 1){
              out$y = mcmc$draws[[ch]][[what]][ind_iterations,dimspec[1]]
            }
            if(v$D == 2){
              out$y = mcmc$draws[[ch]][[what]][ind_iterations,dimspec[1],dimspec[2]]
            }
          }
        } 
      } # end of if mcmc$howsave == list
    
      ### Draws are saved in large matrix
      if(mcmc$howsave == "matrix"){
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
        
        if(!is.element(colname, colnames(mcmc$draws[[ch]]))){
          stop(paste("Wrong specification, column ", colname, " does not exist in mcmc$draws[[",ch,"]]!"))
        }
        
        out$y = mcmc$draws[[ch]][ind_iterations, colname]
        
      } # end of if howsave="matrix"
      scalar_samples[[ch]] <- out
    }# end of for ch in chains
  }
  return(scalar_samples)
}
