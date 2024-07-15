summary.clusglmm <- function(object, ...){
  # computes necessary outputs for summary
  
  res <- list(call=object$call,
              clustering=object$clustering)
  
  class(res) <- "summary.clusglmm"
  return(res)
}