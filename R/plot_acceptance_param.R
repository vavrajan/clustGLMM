plot_acceptance_param <- 
function(mcmc, what = "e0", 
         burnin = 0, thin = 1, iterations, chains, 
         setparmfrow = TRUE, allclusters = FALSE,
         adjust.ylim = 0, which.b = 1:min(mcmc$n, 20)){
  if (setparmfrow) {        
    opar <- par(no.readonly = TRUE)
    on.exit(par(opar))
  }
  
  if(missing(chains)){chains = mcmc$chains}
  if(!is.element(what, c("e0", "b", 
                         "beta_poi_fix", "beta_poi", 
                         "beta_bin_fix", "beta_bin", 
                         "beta_ord_fix", "beta_ord", 
                         "a_ord", "c_ord", "pi_ord",
                         "beta_cat_fix", "beta_cat"))){
    stop(paste0("Parameter ", what, " does not exist or is not sampled via Metropolis-Hastings proposal algorithm.
                Try one of 'e0', 'b', 'beta_poi_fix', 'beta_poi', 'beta_bin_fix', 'beta_bin', 'beta_ord_fix', 'beta_ord', 'a_ord', 'c_ord', 'pi_ord', 'beta_cat_fix', 'beta_cat'."))
  }
  
  if(what == "e0"){
    if(setparmfrow){par(mfrow = c(1, 1))}
    plot_acceptance(mcmc, "e0", 
                    burnin = burnin, thin = thin, 
                    iterations = iterations, chains = chains,
                    adjust.ylim = adjust.ylim,
                    labcex = 1.0)
  }
  
  if(what == "b"){
    if(setparmfrow){par(mfrow = nice_nrow_ncol(length(which.b)))}
    for(i in which.b){
      plot_acceptance(mcmc, paste0("b[", i, "]"), 
                      burnin = burnin, thin = thin, 
                      iterations = iterations, chains = chains,
                      adjust.ylim = adjust.ylim,
                      labcex = ifelse(length(which.b) > 1, 0.6, 1.0))
    }
  }
  
  fixs <- c("beta_poi_fix", "beta_bin_fix", "beta_ord_fix", "beta_cat_fix")
  grps <- c("beta_poi", "beta_bin", "beta_ord", "beta_cat")
  if(mcmc$varying["c_ord"]){
    grps <- c(grps, "a_ord", "c_ord", "pi_ord")
  }else{
    fixs <- c(fixs, "a_ord", "c_ord", "pi_ord")
  }
  
  awhat <- sub("^c_ord", "a_ord", what)
  awhat <- sub("^pi_ord", "a_ord", awhat)
  
  if(is.element(what, fixs)){
    wnames <- colnames(mcmc$accept[[chains[1]]])[grepl(paste0("^", awhat), colnames(mcmc$accept[[chains[1]]]))]
    if(length(wnames) > 0){
      if(setparmfrow){par(mfrow = nice_nrow_ncol(length(wnames)))}
      for(w in wnames){
        aw <- sub("^a_ord", what, w)
        plot_acceptance(mcmc, aw, 
                        burnin = burnin, thin = thin, 
                        iterations = iterations, chains = chains,
                        adjust.ylim = adjust.ylim,
                        labcex = ifelse(length(wnames) > 1, 0.6, 1.0))
      }
    }else{
      warning(paste0("There are no acceptances in mcmc$accept[[", chains[1], "]] for parameter ", what, "."))
    }
  }
  
  if(is.element(what, grps)){
    
    if(allclusters){
      nonempty <- list()
      for(ch in chains){
        nonempty[[ch]] <- 1:mcmc$G[ch] # plot all clusters, not only non-empty
      }
    }else{
      nonempty <- mcmc$clusters
    }
    sameclusters <- TRUE
    for(ch in chains[-1]){
      sameclusters <- sameclusters & setequal(nonempty[[chains[1]]], 
                                              nonempty[[ch]])
    }
    
    if(sameclusters){
      for_chains <- chains[1]
      plot_chains <- list()
      for(ch in for_chains){
        plot_chains[[ch]] <- chains
      }
    }else{
      for_chains <- chains
      plot_chains <- list()
      for(ch in for_chains){
        plot_chains[[ch]] <- ch
      }
    }
    
    
    for(ch in for_chains){
      whatnames <- colnames(mcmc$accept[[ch]])[grepl(paste0("^", awhat), colnames(mcmc$accept[[ch]]))]
      whatfixnames <- colnames(mcmc$accept[[ch]])[grepl(paste0("^", what, "_fix"), colnames(mcmc$accept[[ch]]))]
      wnames <- setdiff(whatnames, whatfixnames)
      wnames <- wnames[grepl(paste0("\\(", nonempty[[ch]][1], "\\)"), wnames)]
      if(length(wnames) > 0){
        if(setparmfrow){par(mfrow = nice_nrow_ncol(length(wnames)*length(nonempty[[ch]])))}
        for(w in wnames){
          for(g in nonempty[[ch]]){
            wg <- sub("\\(\\d+\\)", paste0("(", g, ")"), w)
            awg <- sub("^a_ord", what, wg)
            plot_acceptance(mcmc, awg, 
                            burnin = burnin, thin = thin, 
                            iterations = iterations, chains = plot_chains[[ch]],
                            adjust.ylim = adjust.ylim,
                            labcex = ifelse(length(wnames) > 1, 0.6, 1.0))
          }
        }
      }else{
        warning(paste0("There are no acceptances in mcmc$accept[[", ch, "]] for parameter ", what, "."))
      }
    }
    
  }
  
  
}