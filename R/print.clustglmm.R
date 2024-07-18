print.clustglmm <- function(x, ...){
  cat("Call:\n\n")
  cat(x$call)
  
  if(x$G > 1){
    cat("\n## Clustering based on sampled indicators:\n")
    for(ch in x$chains){
      cat("-----------------------------------------------------------------------\n")
      cat(paste0("Chain ", ch, "\n"))
      cat("-----------------------------------------------------------------------\n")
      cat(paste0("The most-frequent number of non-empty clusters: ", x$modeGplus[ch], "\n"))
      cat(paste0("Empty clusters: ", paste0(setdiff(1:x$G, x$clusters[[ch]]), collapse = ", "), "\n"))
      cat(paste0("Non-empty clusters: ", paste0(x$clusters[[ch]], collapse = ", "), "\n"))
      print(summary(factor(x$clustering[,ch])))
      cat(paste0("Average certainty of clustering into a cluster: \n"))
      print(tapply(x$certainty[,ch], x$clustering[,ch], mean))
    }
  }
}
