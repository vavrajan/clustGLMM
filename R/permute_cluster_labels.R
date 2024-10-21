permute_cluster_labels <-
  function(mcmc, perm){
    # perm contains list of permutations - each for one chain
    # perm[[ch]][1] contains which old cluster should be the new first one (for chain ch)
    # ... and so on
    
    for(ch in mcmc$chains){
      if((length(perm[[ch]]) != mcmc$G[ch]) | (length(unique(perm[[ch]])) < mcmc$G[ch])){
        stop(paste0("Vector given for chain=",ch," is not a permutation on 1:",mcmc$G[ch]))
      }
      
      
      ## permute everything in the original data  
      if(mcmc$howsave == "matrix"){
        # data are save in matrix format
        aux <- mcmc$draws[[ch]]
        
        # permute U
        if(mcmc$save["U"]){
          Unames <- grep("^U\\[", colnames(aux))
          aux[,Unames] <- apply(aux[,Unames], c(1,2), function(g){which(perm[[ch]] == g)})
        }
        
        # cluster-specific parameters
        cnames <- colnames(aux)[grep("\\(1\\)", colnames(aux))]
        all_cnames <- c()
        for(g in 1:mcmc$G[ch]){
          all_cnames <- c(all_cnames, gsub("\\(1\\)", paste0("\\(",g,"\\)"), cnames))
        }
        auxcs <- aux[,all_cnames]
        for(g in 1:mcmc$G[ch]){
          pnames <- gsub("\\(1\\)", paste0("\\(",perm[[ch]][g],"\\)"), cnames)
          gnames <- gsub("\\(1\\)", paste0("\\(",g,"\\)"), cnames)
          aux[,gnames] <- auxcs[,pnames]
        } 
        mcmc$draws[[ch]] <- aux
      }
      
      if(mcmc$howsave == "list"){
        # data are saved in lists
        # permute U
        if(mcmc$save["U"]){
          mcmc$draws[[ch]]$U <- apply(mcmc$draws[[ch]]$U, c(1,2), function(g){which(perm[[ch]] == g)})
        }
        
        # cluster-specific parameters
        for(p in rownames(mcmc$settings[[ch]])[as.logical(mcmc$settings[[ch]]$gspec)]){
          aux <- list()
          if(mcmc$settings[[ch]][p,"ydepd1"]){
            for(g in 1:mcmc$G[ch]){
              aux[[g]] <- list()
              for(y in names(mcmc$draws[[ch]][[p]][[g]])){
                aux[[g]][[y]] <- mcmc$draws[[ch]][[p]][[perm[[ch]][g]]][[y]]
              } 
            }
          }else{
            for(g in 1:mcmc$G[ch]){
              aux[[g]] <- mcmc$draws[[ch]][[p]][[perm[[ch]][g]]]
            } 
          }
          mcmc$draws[[ch]][[p]] <- aux
        }
      }
      
      ## Permute mcmc$last, the last values
      # data are saved in lists
      # permute U
      mcmc$last[[ch]]$U <- sapply(mcmc$last[[ch]]$U, function(g){which(perm[[ch]] == g)})
      
      # cluster-specific parameters
      params <- rownames(mcmc$settings[[ch]])[as.logical(mcmc$settings[[ch]]$gspec)]
      for(p in intersect(params, names(mcmc$last[[ch]]))){
        aux <- list()
        if(mcmc$settings[[ch]][p,"ydepd1"]){
          for(y in names(mcmc$last[[ch]][[p]])){
            aux[[y]] <- list()
            for(g in 1:mcmc$G[ch]){
              aux[[y]][[g]] <- mcmc$last[[ch]][[p]][[y]][[perm[[ch]][g]]]
            } 
          }
        }else{
          for(g in 1:mcmc$G[ch]){
            aux[[g]] <- mcmc$last[[ch]][[p]][[perm[[ch]][g]]]
          } 
        }
        mcmc$last[[ch]][[p]] <- aux
      }
    
      ## Permute other outputs saved in mcmc
      # clustering, clusters
      mcmc$clustering[, ch] <- sapply(mcmc$clustering[, ch], function(g){which(perm[[ch]] == g)})
      mcmc$clusters[[ch]] <- sort(sapply(mcmc$clusters[[ch]], function(g){which(perm[[ch]] == g)}))
    } # end for ch in mcmc$chains
    
    # sameclusters
    sameclusters <- TRUE
    for(ch in mcmc$chains[-1]){
      sameclusters <- (sameclusters & setequal(mcmc$clusters[[mcmc$chains[1]]], mcmc$clusters[[ch]]))
    }
    mcmc$sameclusters <- sameclusters
    
    return(mcmc)
  }
