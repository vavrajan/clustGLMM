print.summary.clustglmm <- function(x, which = c("mcmc", "inv_param", "clust_param"), 
                                    chains = c(1L), 
                                    clusters,
                                    pattern = "", ...){
  # what clusters to print
  if(missing(clusters)){
    clusters <- x$clusters
  }else{
    if(!is.list(clusters)){
      lclusters <- list()
      for(ch in chains){
        lclusters[[ch]] <- clusters
      }
      clusters <- lclusters
    }
  }
  
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
        fnames <- intersect(x$fnames[[ch]], rownames(s))
        fnames <- fnames[grep(pattern, fnames)]
        if(length(fnames) > 0){
          cat("\nCluster-invariant parameters:\n")
          print(s[fnames,,drop=FALSE])
        }
      }
      
      if(is.element("clust_param", which)){
        cnames <- intersect(x$cnames[[ch]], rownames(s))
        cnames <- cnames[grep(pattern, cnames)]
        if(length(cnames) > 0){
          cat("\nClustering-related parameters:\n")
          print(s[cnames,,drop=FALSE])
        }
      }
      
      if(is.element("latent_param", which)){
        lnames <- intersect(x$lnames[[ch]], rownames(s))
        lnames <- lnames[grep(pattern, lnames)]
        if(length(lnames) > 0){
          cat("\nLatent variables:\n")
          print(s[lnames,,drop=FALSE])
        }
      }
      
      if(is.element("group_param", which)){
        gnames <- intersect(x$gnames[[ch]], rownames(s))
        gnames <- gnames[grep(pattern, gnames)]
        if(length(gnames) > 0){
          for(g in clusters[[ch]]){
            ggnames <- gsub("\\(1\\)", paste0("\\(",g,"\\)"), gnames)
            cat(paste0("\nParameters specific for cluster g=", g, ":\n"))
            print(s[ggnames,,drop=FALSE])
          }
        }
      }
    }
    
    cat(paste0("\n"))
  }
}