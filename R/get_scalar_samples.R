get_scalar_samples <-
function(mcmc, what, gspec, dimspec, yspec, 
         burnin = 0, thin = 1, 
         iterations, chains){
  if(thin < 1){
    stop("Thinning parameter has to be integer and thin >= 1.")
  }
  if(missing(chains)){chains = mcmc$chains}

  scalar_samples <- list()
  for(ch in chains){
    if(missing(iterations)){
      # iterations = seq(burnin+1, mcmc$iter, by = thin)
      # iterations <- mcmc$draws[[ch]]$m
      if(missing(gspec)){
        iterations <- mcmc$draws[[ch]]$m
      }else{
        # to avoid NA values created by post_processing()
        iterations <- mcmc$iterations[[ch]]
      }
      iterations <- iterations[iterations > burnin]
      iterations <- iterations[seq(1, length(iterations), by=thin)]
    }
    
    ind_iterations <- is.element(mcmc$draws[[ch]]$m, iterations)    
    v <- mcmc$settings[[ch]][what,]
    
    ### Draws are save in a structured list
    out <- data.frame(chain = ch, m = iterations)
    if(mcmc$howsave == "list"){
      if(v$gspec){
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

  return(scalar_samples)
}
