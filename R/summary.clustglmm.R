summary.clustglmm <- function(object, ...){
  coda = as_coda(object, ...)
  
  # computes necessary outputs for summary
  coda_summary <- list()
  for(ch in 1:object$nchains){
    coda_summary[[ch]] <- summary(coda[[ch]])
  }
  
  output <- ifelse(length(intersect(rownames(object$settings), c("deviance", "dev_i", "pUig_int"))) > 0,
                   "probs", "mcmc")

  if(output == "probs"){
    fnames <- c("deviance")
    gnames <- cnames <- c()
  }else{
    fset = object$settings[object$settings$save & !object$settings$gspec & object$settings$dims>0, ]
    fset = fset[setdiff(rownames(fset), c("loglik", "U", "Gplus", "e0")), ]
    fnames <- c()
    for(p in rownames(fset)){
      fnames <- c(fnames, object$param_names[[p]])
    }
    
    gnames <- rownames(coda_summary[[1]]$statistics)[
      grep("\\(1\\)", rownames(coda_summary[[1]]$statistics))]
    gnames <- gnames[!(grepl("^w", gnames) | grepl("^ng", gnames) | grepl("^pUig", gnames))]
    
    cnames <- c("Gplus", "e0", paste0("w(",1:x$G,")"))
  }
  
  res <- list(call=object$call,
              G=object$G,
              modeGplus=object$modeGplus,
              clustering=object$clustering,
              certainty=object$certainty,
              clusters=object$clusters,
              coda_summary=coda_summary,
              output=output,
              fnames=fnames,
              gnames=gnames,
              cnames=cnames)
  
  class(res) <- "summary.clustglmm"
  return(res)
}