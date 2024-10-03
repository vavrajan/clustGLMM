print.summary.clustglmm <- function(x, chains=c(1L), ...){
  # cat("Call:\n")
  # cat(x$call)
  # cat("\n## Clustering based on sampled indicators:\n")
  # for(ch in chains){
  #   cat("-----------------------------------------------------------------------\n")
  #   cat(paste0("Chain ", ch, "\n"))
  #   cat("-----------------------------------------------------------------------\n")
  #   cat(paste0("The most-frequent number of non-empty clusters: ", x$modeGplus[ch], "\n"))
  #   cat(paste0("Empty clusters: ", paste0(setdiff(1:x$G, x$clusters[[ch]]), collapse = ", "), "\n"))
  #   cat(paste0("Non-empty clusters: ", paste0(x$clusters[[ch]], collapse = ", "), "\n"))
  #   print(summary(factor(x$clustering[,ch])))
  #   cat(paste0("Average certainty of clustering into a cluster: \n"))
  #   print(tapply(x$certainty[,ch], x$clustering[,ch], mean))
  # }
  # 
  
  cat("\n## Summary of the sampled parameters:\n")
  for(ch in chains){
    ss <- x$coda_summary[[ch]]$statistics
    sq <- x$coda_summary[[ch]]$quantiles
    s <- cbind(sq[,c("2.5%", "25%", "50%")], ss[,"Mean",drop=F], sq[,c("75%", "97.5%")])
    cat("-----------------------------------------------------------------------\n")
    cat(paste0("Chain ", ch, "\n"))
    cat("-----------------------------------------------------------------------\n")
    cat(paste0("Iterations = ", x$coda_summary[[ch]]$start, ":", x$coda_summary[[ch]]$end, "\n"))
    cat(paste0("Thinning interval = ", x$coda_summary[[ch]]$thin, "\n"))
    cat(paste0("Number of chains = ", x$coda_summary[[ch]]$nchain, "\n"))
    cat(paste0("Sample size per chain = ", ceiling((x$coda_summary[[ch]]$end - x$coda_summary[[ch]]$start+1) / x$coda_summary[[ch]]$nchain), "\n"))
    
    cat(paste0("\nSummary of the clustering based on ",
               ifelse(x$output == "mcmc", 
                      "sampled allocation indicators U",
                      "approximated clustering probabilities pUig_int"),
               "\n"))
    cat("-----------------------------------------------------------------------\n")
    cat(paste0("The most-frequent number of non-empty clusters: ", x$modeGplus[ch], "\n"))
    cat(paste0("Empty clusters: ", paste0(setdiff(1:x$G, x$clusters[[ch]]), collapse = ", "), "\n"))
    cat(paste0("Non-empty clusters: ", paste0(x$clusters[[ch]], collapse = ", "), "\n"))
    print(summary(factor(x$clustering[,ch])))
    cat(paste0("Average certainty of clustering into a cluster: \n"))
    print(tapply(x$certainty[,ch], x$clustering[,ch], mean))
    
    cat("\nSummary of the key parameter estimates\n")
    cat("-----------------------------------------------------------------------\n")
    
    if(length(x$fnames) > 0){
      cat("\nCluster-invariant parameters:\n")
      print(s[x$fnames,,drop=FALSE])
    }
    
    cint <- intersect(x$cnames, rownames(s))
    if(length(cint) > 0){
      cat("\nClustering-related parameters:\n")
      print(s[cint,,drop=FALSE])
    }
    
    if(length(x$gnames) > 0){
      for(g in x$clusters[[ch]]){
        gnames <- gsub("\\(1\\)", paste0("\\(",g,"\\)"), x$gnames)
        cat(paste0("\nParameters specific for cluster g=", g, ":\n"))
        print(s[gnames,,drop=FALSE])
      }
    }
    
    cat(paste0("\n"))
  }
}