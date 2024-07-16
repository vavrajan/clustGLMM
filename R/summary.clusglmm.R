summary.clusglmm <- function(object, ...){
  coda = as_coda(object, ...)
  
  # computes necessary outputs for summary
  coda_summary <- list()
  for(ch in 1:object$nchains){
    coda_summary[[ch]] <- summary(coda[[ch]])
  }
  
  fset = object$settings[object$settings$save & !object$settings$isspec & object$settings$dims>0, ]
  fset = fset[setdiff(rownames(fset), c("loglik", "ng", "w", "U", "pUig", "Gplus", "e0")), ]
  fnames <- c()
  for(p in rownames(fset)){
    fnames <- c(fnames, object$param_names[[p]])
  }
  gnames <- grep("\\(1\\)", rownames(coda_summary[[1]]$statistics))
  
  res <- list(call=object$call,
              G=object$G,
              modeGplus=object$modeGplus,
              clustering=object$clustering,
              certainty=object$certainty,
              clusters=object$clusters,
              coda_summary=coda_summary,
              fnames=fnames,
              gnames=gnames)
  
  class(res) <- "summary.clusglmm"
  return(res)
}