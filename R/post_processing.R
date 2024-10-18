post_processing <-
function(mcmc, cols){
  
  ### 0) The type of saving results has to be "matrix" for the code below
  backtolist <- FALSE
  if(mcmc$howsave=="list"){
    # Then transfer to howsave="matrix" -- create mcmc$draws matrix
    backtolist <- TRUE
    mcmc <- from_list_to_matrix(mcmc)
  }
  
  ### 1) Find the posterior draws with Gplus = modeGplus
  inds <- ms <- list()
  iters <- numeric(mcmc$nchains)
  for(ch in mcmc$chains){
    submcmc <- mcmc$draws[[ch]][, c("m", "Gplus")]
    ms[[ch]] <- submcmc[submcmc$Gplus == mcmc$modeGplus[ch], "m"]
    inds[[ch]] <- which(submcmc$Gplus == mcmc$modeGplus[ch])
    iters[ch] <- length(inds[[ch]])
  }
  message("modeGplus")
  print(mcmc$modeGplus)
  message("Frequency of the modeGplus value")
  print(iters)
  
  ### 2) Columns of mcmc$draws that are cluster-specific parameters
  allgcols <- colnames(mcmc$draws[[ch]])[grep("\\(\\d+\\)", colnames(mcmc$draws[[ch]]))]
  # but not clustering related variables (weights, sizes, probabilities)
  allgcols <- allgcols[!(grepl("^w", allgcols) | grepl("^ng", allgcols) | grepl("^pUig", allgcols) | grepl("^pUig_int", allgcols))]
  
  if(missing(cols)){
    cols <- allgcols
  }else{
    # use only selected cols for k-means clustering
  }
  message("columns used for k-means")
  print(cols)
  
  ### 3) Creating submatrices of sampled cluster-specific parameters
  clusterspec <- whichG <- list()
  for(ch in mcmc$chains){
    suball <- mcmc$draws[[ch]][inds[[ch]], 
                               c("m", "Gplus", paste0("ng(", 1:mcmc$G, ")"), cols)]
    clusterspec[[ch]] <- matrix(NA, nrow = iters[ch], ncol = 1 + length(cols)/mcmc$G*mcmc$modeGplus[ch])
    whichG[[ch]] <- matrix(NA, nrow = iters[ch], ncol = mcmc$modeGplus[ch])
    for(m in 1:iters[ch]){
      whichG[[ch]][m,] <- which(suball[m, paste0("ng(", 1:mcmc$G, ")")] > 0)
      #print(whichG)
      expr <- paste(paste0("\\(", whichG[[ch]][m,], "\\)"), collapse = "|")
      subcols <- cols[grep(expr, cols)]
      clusterspec[[ch]][m,] <- unlist(suball[m,c("m",subcols)])
    }
    expr <- paste(paste0("\\(", 1:mcmc$modeGplus[ch], "\\)"), collapse = "|")
    subcols <- cols[grep(expr, cols)]
    cnames <- c("m", subcols)
    colnames(clusterspec[[ch]]) <- cnames
  }
  message("Subset of states of cluster-specific parameters with Gplus=modeGplus 
created for each chain separately")
  
  
  ### 4) Prepare for k-means clustering
  forkmeans <- list()
  for(ch in mcmc$chains){
    forkmeans[[ch]] <- matrix(NA, nrow = 0, ncol = length(cols)/mcmc$G)
    for(g in 1:mcmc$modeGplus[ch]){
      forkmeans[[ch]] <- rbind(forkmeans[[ch]], clusterspec[[ch]][,grep(paste0("\\(", g, "\\)"), 
                                                                        colnames(clusterspec[[ch]]))])
    }
    cnames <- colnames(clusterspec[[ch]])[grep("\\(1\\)", 
                                               colnames(clusterspec[[ch]]))]
    colnames(forkmeans[[ch]]) <- gsub("\\(1\\)", "", cnames)
  }
  message("Group-specific parameters stacked below each other, prepared for k-means")
  
  ### 5) Perform K-means clustering + checking the permutations
  Mrhos <- numeric(mcmc$nchains)
  indmat <- rowuniques <- list()
  for(ch in mcmc$chains){
    # standardization
    forkmeans[[ch]] <- matrix(as.numeric(forkmeans[[ch]]), ncol = dim(forkmeans[[ch]])[2])
    colnames(forkmeans[[ch]]) <- paste0("V", 1:dim(forkmeans[[ch]])[2])
    #class(forkmeans[[ch]])
    #dim(forkmeans[[ch]])
    #summary(forkmeans[[ch]])
    forkmeans[[ch]] <- scale(as.matrix(forkmeans[[ch]]))
    km <- kmeans(forkmeans[[ch]], centers = modeGplus[ch])
    # create index matrix
    indmat[[ch]] <- matrix(km$cluster, ncol = modeGplus[ch], byrow = F)
    #head(indmat[[ch]])
    # number of not permutations = number of rows where is something repeated
    rowuniques[[ch]] <- apply(indmat[[ch]], 1, function(x){length(unique(x))})
    Mrhos[ch] <- sum(rowuniques[[ch]] < modeGplus[ch])
  }
  # mcmc$Mrhos <- Mrhos
  message("K-means performed")
  message("Percentage of non-permutations within K-means")
  print(paste0(Mrhos / iters * 100, " %"))
  
  ### 6) Update the clusterspec matrix according to suggested permutations
  for(ch in mcmc$chains){
    #head(clusterspec[[ch]])
    cnames <- colnames(clusterspec[[ch]])[grep("\\(1\\)", 
                                               colnames(clusterspec[[ch]]))]
    for(m in 1:iters[ch]){
      if(rowuniques[[ch]][m] == modeGplus[ch]){
        # it is a permutation --> permute accordingly
        auxrow <- unlist(clusterspec[[ch]][m,])
        for(g in 1:modeGplus[ch]){
          gnames <- gsub("\\(1\\)", paste0("\\(",g,"\\)"), cnames)
          pnames <- gsub("\\(1\\)", paste0("\\(",indmat[[ch]][m,g],"\\)"), cnames)
          clusterspec[[ch]][m,pnames] <- auxrow[gnames]
        }
      }else{
        # not a permutation --> replace with NA
        clusterspec[[ch]][m,] <- NA
      }
    }
  }
  message("Each row permuted according to the found permutation")
  
  ### 7) New draws
  draws <- list()
  for(ch in mcmc$chains){
    # Parameters that are not group-specific and do not help with clustering
    nongcols <- setdiff(colnames(mcmc$draws[[ch]]), allgcols)
    nongcols <- nongcols[!(grepl("^w", nongcols) | 
                             grepl("^ng",nongcols) | 
                             grepl("^pUig", nongcols)| 
                             grepl("^pUig_int", nongcols)|
                             grepl("^U", nongcols))]
    draws[[ch]] <- mcmc$draws[[ch]][inds[[ch]], nongcols]
    
    # Add permuted clusterspec - only cols
    draws[[ch]] <- cbind(draws[[ch]], clusterspec[[ch]])
    
    # Add other group-specific columns than those given
    othergcols <- setdiff(allgcols, cols)
    if(length(othergcols) > 0){
      suball <- mcmc$draws[[ch]][inds[[ch]], othergcols]
      cnames <- othergcols[grep("\\(1\\)", othergcols)]
      oGplus <- matrix(NA, nrow = iters[ch], ncol = modeGplus[ch]*length(cnames))
      allcnames <- c()
      for(g in 1:mcmc$modeGplus[ch]){
        allcnames <- c(allcnames, gsub("\\(1\\)", paste0("\\(",g,"\\)"), cnames))
      }
      colnames(oGplus) <- allcnames
      
      for(m in 1:iters[ch]){
        if(rowuniques[[ch]][m] == modeGplus[ch]){
          expr <- paste(paste0("\\(", whichG[[ch]][m,], "\\)"), collapse = "|")
          subcols <- othergcols[grep(expr, othergcols)]
          auxrow <- unlist(suball[m,subcols])
          for(g in 1:modeGplus[ch]){
            gnames <- gsub("\\(1\\)", paste0("\\(",g,"\\)"), cnames)
            pnames <- gsub("\\(1\\)", paste0("\\(",indmat[[ch]][m,g],"\\)"), cnames)
            oGplus[m,pnames] <- auxrow[gnames]
          }
        }else{
          # remains NA
        }
      }
      draws[[ch]] <- cbind(draws[[ch]], oGplus)
    }
    
    # Add w and ng
    if(is.element("w", rownames(mcmc$settings[mcmc$settings$save,]))){
      # if "ng" would not be saved, then it would already complain above
      suball <- mcmc$draws[[ch]][inds[[ch]], c(paste0("w(", 1:mcmc$G, ")"),
                                               paste0("ng(", 1:mcmc$G, ")"))]
      for(p in c("w", "ng")){
        pGplus <- matrix(NA, nrow = iters[ch], ncol = modeGplus[ch])
        colnames(pGplus) <- paste0(p,"(", 1:modeGplus[ch], ")")
        for(m in 1:iters[ch]){
          if(rowuniques[[ch]][m] == modeGplus[ch]){
            pcols <- paste0(p,"(", whichG[[ch]][m,], ")")
            auxrow <- unlist(suball[m,pcols])
            for(g in 1:modeGplus[ch]){
              gname <- paste0(p,"(",g,")")
              pname <- paste0(p,"(",indmat[[ch]][m,g],")")
              pGplus[m,pname] <- auxrow[gname]
            }
          }else{
            # remains NA
          }
        }
        draws[[ch]] <- cbind(draws[[ch]], pGplus)
      }
    } # end of adding w and ng
    
    # Add U
    if(is.element("U", rownames(mcmc$settings[mcmc$settings$save,]))){
      suball <- mcmc$draws[[ch]][inds[[ch]], paste0("U[", 1:mcmc$n, "]")]
      UGplus <- matrix(NA, nrow = iters[ch], ncol = mcmc$n)
      colnames(UGplus) <- paste0("U[", 1:mcmc$n, "]")
      for(m in 1:iters[ch]){
        if(rowuniques[[ch]][m] == modeGplus[ch]){
          for(i in 1:mcmc$n){
            u <- suball[m,paste0("U[",i,"]")] 
            whu <- which(whichG[[ch]][m,] == u) # order within non-zero clusters
            UGplus[m,i] <- indmat[[ch]][m,whu] # applying permutation
          }
        }else{
          # remain NA
        }
      }
      draws[[ch]] <- cbind(draws[[ch]], UGplus)
    } # end of adding U
    
    # Add pUig, pUig_int
    for(p in c("pUig", "pUig_int")){
      if(is.element(p, rownames(mcmc$settings[mcmc$settings$save,]))){
        pGplus <- matrix(NA, nrow = iters[ch], ncol = 0)
        for(i in 1:mcmc$n){
          suball <- mcmc$draws[[ch]][inds[[ch]], 
                                     c(paste0(p,"(",1:mcmc$G,")[",i,"]"),
                                       paste0("ng(",1:mcmc$G,")"))]
          pU <- matrix(NA, nrow = iters[ch], ncol = modeGplus[ch])
          colnames(pU) <- paste0(p,"(",1:modeGplus[ch],")[",i,"]")
          for(m in 1:iters[ch]){
            if(rowuniques[[ch]][m] == modeGplus[ch]){
              pcols <- paste0(p,"(",whichG[[ch]][m,],")[",i,"]")
              auxrow <- suball[m,pcols]
              for(g in 1:modeGplus[ch]){
                gname <- paste0(p,"(",g,")[",i,"]")
                pname <- paste0(p,"(",indmat[[ch]][m,g],")[",i,"]")
                pU[m,pname] <- auxrow[gname]
              }
            }else{
              # remain NA
            }
          }
          pGplus <- cbind(pGplus, pU)
        }
      }
    } # end adding pUig, pUig_int
  } # end of creating draws
  
  
  fun <- function(U, ch){
    tab = table(U)
    if(clusteringU_modeGplus){
      tab <- tab[as.character(clusters[[ch]])]
    }
    c(as.numeric(names(which.max(tab))), max(tab, na.rm=TRUE)/iter)
  }
  
  ### 8) make the output of the same class as input
  # + make necessary changes
  clusters <- list()
  clustering <- matrix(NA, nrow = as.numeric(n), ncol = mcmc$nchains)
  certainty <- matrix(NA, nrow = as.numeric(n), ncol = mcmc$nchains)
  for(ch in mcmc$chains){
    clusters[[ch]] <- 1:modeGplus[ch]
    # U samples
    if(mcmc$settings["U","save"]){
      Us <- draws[[ch]][,paste0("U[",1:mcmc$n,"]")]
      ccU <- apply(Us, 2, fun, ch=ch)
      clustering[,ch] <- ccU[1, ]
      certainty[,ch] <- ccU[2, ]
    }else{
      # we no longer have allocation indicators saved
      # remains as NA
      # TODO think about how it could be done? ... 
      #      problem is each m could be permuted differently
    }
  }
  rownames(clustering) <- rownames(certainty) <- unique_ids

  sameclusters <- TRUE
  for(ch in mcmc$chains[-1]){
    sameclusters <- (sameclusters & setequal(clusters[[mcmc$chains[1]]], clusters[[ch]]))
  }
  mcmc$sameclusters <- sameclusters
  
  ### Final output
  res <- mcmc
  res$clusters <- clusters
  res$clustering <- clustering
  res$certainty <- certainty
  res$sameclusters <- sameclusters
  res$draws <- draws
  # inits, last?
  # settings?
  # iter?
  # G?
  # different chain ... different G possible?
  # res$post_processed <- TRUE
  
  # call
  res$call <- c("###------------------------------------------###\n")
  res$call <- paste0(res$call, "### Post-processed MCMC samples of clustGLMM ###\n")
  res$call <- paste0(res$call, "###------------------------------------------###\n")
  # TODO what should be added to call?
  
  
  if(backtolist){
    res <- from_matrix_to_list(res)
  }
  class(res) <- "clustglmm"
  
  return(res)
}
