print.clusglmm <- function(x, ...){
  cat("Call:\n\n")
  cat(x$call)
  cat("\nClustering based on sampled indicators:\n\n")
  for(ch in x$chains){
    cat(paste0("\nChain ", ch, "\n"))
    cat(paste0("The most-frequent number of non-empty clusters: ", x$modeGplus[ch], "\n"))
    cat(paste0("Empty clusters: ", paste0(setdiff(1:x$G, x$clusters[[ch]]), collapse = ", "), "\n"))
    cat(paste0("Non-empty clusters: ", paste0(x$clusters[[ch]], collapse = ", "), "\n"))
    print(summary(factor(x$clustering[,ch])))
    cat(paste0("Average certainty of clustering into a cluster: \n"))
    print(tapply(mcmc$certainty[,ch], mcmc$clustering[,ch], mean))
  }
}
