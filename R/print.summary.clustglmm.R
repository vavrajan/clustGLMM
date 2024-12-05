print.summary.clustglmm <- function(x, chains=c(1L), 
                                    which = c("mcmc", "inv_param", "clust_param"), ...){
  # cat("Call:\n")
  # cat(x$call)
  # cat("\n## Clustering based on sampled indicators:\n")
  # for(ch in chains){
  #   cat("-----------------------------------------------------------------------\n")
  #   cat(paste0("Chain ", ch, "\n"))
  #   cat("-----------------------------------------------------------------------\n")
  #   cat(paste0("The most-frequent number of non-empty clusters: ", x$modeGplus[ch], "\n"))
  #   cat(paste0("Empty clusters: ", paste0(setdiff(1:x$G[ch], x$clusters[[ch]]), collapse = ", "), "\n"))
  #   cat(paste0("Non-empty clusters: ", paste0(x$clusters[[ch]], collapse = ", "), "\n"))
  #   print(summary(factor(x$clustering[,ch])))
  #   cat(paste0("Average certainty of clustering into a cluster: \n"))
  #   print(tapply(x$certainty[,ch], x$clustering[,ch], mean))
  # }
  # 
  
  cat("\n## Summary of the sampled parameters:\n")
  for(ch in chains){
    ss <- as.matrix(x$coda_summary[[ch]]$statistics)
    sq <- as.matrix(x$coda_summary[[ch]]$quantiles)
    s <- cbind(sq[,c("2.5%", "25%", "50%")], ss[,"Mean",drop=F], sq[,c("75%", "97.5%")])
    cat("-----------------------------------------------------------------------\n")
    cat(paste0("Chain ", ch, "\n"))
    cat("-----------------------------------------------------------------------\n")
    if(is.element("mcmc", which)){
      cat(paste0("Iterations = ", x$coda_summary[[ch]]$start, ":", x$coda_summary[[ch]]$end, "\n"))
      cat(paste0("Thinning interval = ", x$coda_summary[[ch]]$thin, "\n"))
      cat(paste0("Number of chains = ", x$coda_summary[[ch]]$nchain, "\n"))
      cat(paste0("Sample size per chain = ", ceiling((x$coda_summary[[ch]]$end - x$coda_summary[[ch]]$start+1) / x$coda_summary[[ch]]$nchain), "\n"))
    }
    
    if(is.element("clustering", which)){
      cat(paste0("\nSummary of the clustering based on ",
                 ifelse(x$output == "mcmc", 
                        "sampled allocation indicators U",
                        "approximated clustering probabilities pUig_int"),
                 "\n"))
      cat("-----------------------------------------------------------------------\n")
      cat(paste0("The most-frequent number of non-empty clusters: ", x$modeGplus[ch], "\n"))
      cat(paste0("Empty clusters: ", paste0(setdiff(1:x$G[ch], x$clusters[[ch]]), collapse = ", "), "\n"))
      cat(paste0("Non-empty clusters: ", paste0(x$clusters[[ch]], collapse = ", "), "\n"))
      print(summary(factor(x$clustering[,ch])))
      cat(paste0("Average certainty of clustering into a cluster: \n"))
      print(tapply(x$certainty[,ch], x$clustering[,ch], mean))
    }
    
    if(length(intersect(which, 
                        c("inv_param", "clust_param", "latent_param", "group_param"))
              ) > 0){
      cat("\nSummary of the key parameter estimates\n")
      cat("-----------------------------------------------------------------------\n")

      if(is.element("inv_param", which)){
        if(length(x$fnames[[ch]]) > 0){
          cat("\nCluster-invariant parameters:\n")
          print(s[x$fnames[[ch]],,drop=FALSE])
        }
      }
      
      if(is.element("clust_param", which)){
        cint <- intersect(x$cnames[[ch]], rownames(s))
        if(length(cint) > 0){
          cat("\nClustering-related parameters:\n")
          print(s[cint,,drop=FALSE])
        }
      }
      
      if(is.element("latent_param", which)){
        iint <- intersect(x$inames[[ch]], rownames(s))
        if(length(iint) > 0){
          cat("\nLatent variables:\n")
          print(s[iint,,drop=FALSE])
        }
      }
      
      if(is.element("group_param", which)){
        if(length(x$gnames[[ch]]) > 0){
          for(g in x$clusters[[ch]]){
            gnames <- gsub("\\(1\\)", paste0("\\(",g,"\\)"), x$gnames[[ch]])
            cat(paste0("\nParameters specific for cluster g=", g, ":\n"))
            print(s[gnames,,drop=FALSE])
          }
        }
      }
    }
    
    cat(paste0("\n"))
  }
}