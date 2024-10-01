PermuteClusterSpec <-
function(mcmc, perm){
  # perm contains list of permutations - each for one chain
  # perm[[ch]][1] contains which old cluster should be the new first one (for chain ch)
  # ... and so on
  
  ## permute only the restricted parameters according to given list of permutations
  for(ch in 1:mcmc$nchains){
    # permute w
    if((!is.null(mcmc$wGplus)) & (length(mcmc$wGplus) > 0)){
      auxwGplus <- matrix(-1, nrow = dim(mcmc$wGplus[[ch]])[1], ncol = 0)
      for(g in 1:mcmc$modeGplus[[ch]]){
        auxwGplus <- cbind(auxwGplus, mcmc$wGplus[[ch]][,perm[[ch]][g]])
      }
      colnames(auxwGplus) <- colnames(mcmc$wGplus[[ch]])
      mcmc$wGplus[[ch]] <- auxwGplus
    }

    # permute ng
    if((!is.null(mcmc$ngGplus)) & (length(mcmc$ngGplus) > 0)){
      auxngGplus <- matrix(-1, nrow = dim(mcmc$ngGplus[[ch]])[1], ncol = 0)
      for(g in 1:mcmc$modeGplus[[ch]]){
        auxngGplus <- cbind(auxngGplus, mcmc$ngGplus[[ch]][,perm[[ch]][g]])
      }
      colnames(auxngGplus) <- colnames(mcmc$ngGplus[[ch]])
      mcmc$ngGplus[[ch]] <- auxngGplus
    }

    # permute U
    if((!is.null(mcmc$UGplus)) & (length(mcmc$UGplus) > 0)){
      mcmc$UGplus[[ch]] <- apply(mcmc$UGplus[[ch]], c(1,2), function(g){which(perm[[ch]] == g)})
    }
    
    # permute pUig
    if((!is.null(mcmc$pUigGplus)) & (length(mcmc$pUigGplus) > 0)){
      auxpU <- matrix(-1, nrow = dim(mcmc$pUigGplus[[ch]])[1], ncol = 0)
      for(i in 1:mcmc$settings["U","dims"]){
        auxpUig <- matrix(-1, nrow = dim(mcmc$pUigGplus[[ch]])[1], ncol = 0)
        auxpUi <- mcmc$pUigGplus[[ch]][,paste0("pUig[",i,",",1:mcmc$modeGplus[[ch]], "]")]
        for(g in 1:mcmc$modeGplus[[ch]]){
          auxpUig <- cbind(auxpUig, auxpUi[,perm[[ch]][g]])
        }
        colnames(auxpUig) <- colnames(auxpUi)
        auxpU <- cbind(auxpU, auxpUig)
      }
      mcmc$pUigGplus[[ch]] <- auxpU
    }
    
    # permute clusterspec
    if((!is.null(mcmc$clusterspec)) & (length(mcmc$clusterspec) > 0)){
      aux <- mcmc$clusterspec[[ch]]
      cnames <- colnames(aux)[grep("\\(1\\)", colnames(aux))]
      for(g in 1:mcmc$modeGplus[[ch]]){
        pnames <- gsub("\\(1\\)", paste0("\\(",perm[[ch]][g],"\\)"), cnames)
        gnames <- gsub("\\(1\\)", paste0("\\(",g,"\\)"), cnames)
        aux[,gnames] <- mcmc$clusterspec[[ch]][,pnames]
      } 
      mcmc$clusterspec[[ch]] <- aux
    }
    
  } # end of for ch in 1:nchains
  
  
  ## permute everything in the original data  
  if(mcmc$howsave == "matrix"){
    # data are save in matrix mcmc$draws
    draws <- list()
    for(ch in 1:mcmc$nchains){
      aux <- mcmc$draws[[ch]]
      
      # permute U
      if(mcmc$save["U"]){
        Unames <- grep("^U\\[", colnames(aux))
        aux[,Unames] <- apply(aux[,Unames], c(1,2), function(g){which(perm[[ch]] == g)})
      }
      
      # # permute w
      # if(mcmc$save["w"]){
      #   wnames <- grep("^w\\[", colnames(aux))
      #   auxw <- matrix(-1, nrow = dim(aux)[1], ncol = 0)
      #   for(g in 1:mcmc$G){
      #     auxw <- cbind(auxw, aux[,paste0("w[", perm[[ch]][g], "]")])
      #   }
      #   colnames(auxw) <- wnames
      #   aux[,wnames] <- auxw
      # }
      # 
      # # permute ng
      # if(mcmc$save["ng"]){
      #   ngnames <- grep("^ng\\[", colnames(aux))
      #   auxng <- matrix(-1, nrow = dim(aux)[1], ncol = 0)
      #   for(g in 1:mcmc$G){
      #     auxng <- cbind(auxng, aux[,paste0("ng[", perm[[ch]][g], "]")])
      #   }
      #   colnames(auxng) <- ngnames
      #   aux[,ngnames] <- auxng
      # }
      # 
      # # permute pUig
      # if(mcmc$save["pUig"]){
      #   pUignames <- grep("^pUig\\[", colnames(aux))
      #   auxpU <- matrix(-1, nrow = dim(aux)[1], ncol = 0)
      #   for(i in 1:mcmc$settings["U","dims"]){
      #     auxpUig <- matrix(-1, nrow = dim(aux)[1], ncol = 0)
      #     auxpUi <- aux[,paste0("pUig[",i,",",1:mcmc$G, "]")]
      #     for(g in 1:mcmc$G){
      #       auxpUig <- cbind(auxpUig, auxpUi[,perm[[ch]][g]])
      #     }
      #     colnames(auxpUig) <- colnames(auxpUi)
      #     auxpU <- cbind(auxpU, auxpUig)
      #   }
      #   aux[,pUignames] <- auxpU[,pUignames]
      # }
      
      # cluster-specific parameters
      cnames <- colnames(aux)[grep("\\(1\\)", colnames(aux))]
      all_cnames <- c()
      for(g in 1:mcmc$G){
        all_cnames <- c(all_cnames, gsub("\\(1\\)", paste0("\\(",g,"\\)"), cnames))
      }
      auxcs <- aux[,all_cnames]
      for(g in 1:mcmc$G){
        pnames <- gsub("\\(1\\)", paste0("\\(",perm[[ch]][g],"\\)"), cnames)
        gnames <- gsub("\\(1\\)", paste0("\\(",g,"\\)"), cnames)
        aux[,gnames] <- auxcs[,pnames]
      } 
      draws[[ch]] <- aux
    }
    mcmc$draws <- draws
  }
  
  if(mcmc$howsave == "list"){
    # data are saved in lists
    for(ch in 1:mcmc$nchains){
      # permute U
      if(mcmc$save["U"]){
        mcmc$draws[[ch]]$U <- apply(mcmc$draws[[ch]]$U, c(1,2), function(g){which(perm[[ch]] == g)})
      }
      
      # # permute w
      # if(mcmc$save["w"]){
      #   auxw <- matrix(-1, nrow = dim(mcmc$draws[[ch]]$w)[1], ncol = 0)
      #   for(g in 1:mcmc$G){
      #     auxw <- cbind(auxw, mcmc$draws[[ch]]$w[,perm[[ch]][g]])
      #   }
      #   colnames(auxw) <- colnames(mcmc$draws[[ch]]$w)
      #   mcmc$draws[[ch]]$w <- auxw
      # }
      # 
      # # permute ng
      # if(mcmc$save["ng"]){
      #   auxng <- matrix(-1, nrow = dim(mcmc$draws[[ch]]$ng)[1], ncol = 0)
      #   for(g in 1:mcmc$G){
      #     auxng <- cbind(auxng, mcmc$draws[[ch]]$ng[,perm[[ch]][g]])
      #   }
      #   colnames(auxng) <- colnames(mcmc$draws[[ch]]$ng)
      #   mcmc$draws[[ch]]$ng <- auxng
      # }
      # 
      # # permute pUig
      # if(mcmc$save["pUig"]){
      #   auxpU <- matrix(-1, nrow = dim(mcmc$draws[[ch]]$pUig)[1], ncol = 0)
      #   for(i in 1:mcmc$settings["U","dims"]){
      #     auxpUig <- matrix(-1, nrow = dim(mcmc$draws[[ch]]$pUig)[1], ncol = 0)
      #     auxpUi <- mcmc$draws[[ch]]$pUig[,paste0("pUig[",i,",",1:mcmc$G, "]")]
      #     for(g in 1:mcmc$G){
      #       auxpUig <- cbind(auxpUig, auxpUi[,perm[[ch]][g]])
      #     }
      #     colnames(auxpUig) <- colnames(auxpUi)
      #     auxpU <- cbind(auxpU, auxpUig)
      #   }
      #   mcmc$draws[[ch]]$pUig <- auxpU
      # }
      
      # cluster-specific parameters
      for(p in rownames(mcmc$settings)[as.logical(mcmc$settings$gspec)]){
        aux <- list()
        if(mcmc$settings[p,"ydepd1"]){
          for(g in 1:mcmc$G){
            aux[[g]] <- list()
            for(y in names(mcmc$draws[[ch]][[p]][[g]])){
              aux[[g]][[y]] <- mcmc$draws[[ch]][[p]][[perm[[ch]][g]]][[y]]
            } 
          }
        }else{
          for(g in 1:mcmc$G){
            aux[[g]] <- mcmc$draws[[ch]][[p]][[perm[[ch]][g]]]
          } 
        }
        mcmc$draws[[ch]][[p]] <- aux
      }
    }
  }
  
  ## Permute mcmc$last, the last values
  # data are saved in lists
  for(ch in 1:mcmc$nchains){
    # permute U
    mcmc$last[[ch]]$U <- sapply(mcmc$last[[ch]]$U, function(g){which(perm[[ch]] == g)})
    
    # # permute w
    # mcmc$last[[ch]]$w <- mcmc$last[[ch]]$w[perm[[ch]]]
    # 
    # # permute ng
    # mcmc$last[[ch]]$ng <- mcmc$last[[ch]]$ng[perm[[ch]]]
    # 
    # # permute pUig
    # auxpUig <- matrix(-1, nrow = dim(mcmc$last[[ch]]$pUig)[1], ncol = 0)
    # for(g in 1:mcmc$G){
    #   auxpUig <- cbind(auxpUig, mcmc$last[[ch]]$pUig[,perm[[ch]][g]])
    # }
    # colnames(auxpUig) <- colnames(mcmc$last[[ch]]$pUig)
    # mcmc$last[[ch]]$pUig <- auxpUig
    
    # cluster-specific parameters
    params <- rownames(mcmc$settings)[as.logical(mcmc$settings$gspec)]
    for(p in intersect(params, names(mcmc$last[[ch]]))){
      aux <- list()
      if(mcmc$settings[p,"ydepd1"]){
        for(y in names(mcmc$last[[ch]][[p]])){
          aux[[y]] <- list()
          for(g in 1:mcmc$G){
            aux[[y]][[g]] <- mcmc$last[[ch]][[p]][[y]][[perm[[ch]][g]]]
          } 
        }
      }else{
        for(g in 1:mcmc$G){
          aux[[g]] <- mcmc$last[[ch]][[p]][[perm[[ch]][g]]]
        } 
      }
      mcmc$last[[ch]][[p]] <- aux
    }
  }
  
  ## Permute other outputs saved in mcmc
  # clustering, clusters
  for(ch in 1:mcmc$nchains){
    mcmc$clustering[, ch] <- sapply(mcmc$clustering[, ch], function(g){which(perm[[ch]] == g)})
    mcmc$clusters[[ch]] <- sort(sapply(mcmc$clusters[[ch]], function(g){which(perm[[ch]] == g)}))
  }
  # sameclusters
  sameclusters <- TRUE
  for(ch in 2:mcmc$nchains){
    if(ch <= mcmc$nchains){
      sameclusters <- (sameclusters & setequal(mcmc$clusters[[1]], mcmc$clusters[[ch]]))
    }
  }
  mcmc$sameclusters <- sameclusters
    
  return(mcmc)
}
