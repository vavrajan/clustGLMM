print.summary.clustglmm <- function(x, chains=c(1), type="statistics", ...){
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
  wnames <- paste0("w[",1:x$G,"]")
  
  cat("\n## Summary of the sampled parameters:\n")
  for(ch in chains){
    s <- x$coda_summary[[ch]][[type]]
    cat("-----------------------------------------------------------------------\n")
    cat(paste0("Chain ", ch, "\n"))
    cat("-----------------------------------------------------------------------\n")
    cat(paste0("Iterations = ", x$coda_summary[[ch]]$start, ":", x$coda_summary[[ch]]$end, "\n"))
    cat(paste0("Thinning interval = ", x$coda_summary[[ch]]$thin, "\n"))
    cat(paste0("Number of chains = ", x$coda_summary[[ch]]$nchain, "\n"))
    cat(paste0("Sample size per chain = ", ceiling((x$coda_summary[[ch]]$end - x$coda_summary[[ch]]$start+1) / x$coda_summary[[ch]]$nchain), "\n"))
    cat("\nCluster-invariant parameters:\n")
    print(s[x$fnames,,drop=FALSE])
    if(all(is.element(wnames, rownames(s)))){
      cat("\nCluster weights:\n")
      print(s[wnames,,drop=FALSE])
    }
    for(g in x$clusters[[ch]]){
      gnames <- gsub("\\(1\\)", paste0("\\(",g,"\\)"), x$gnames)
      cat(paste0("\nParameters specific for cluster g=", g, ":\n"))
      print(s[gnames,,drop=FALSE])
    }
  }
}