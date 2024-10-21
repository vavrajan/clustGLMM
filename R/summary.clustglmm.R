summary.clustglmm <- function(object, ...){
  coda = as_coda(object, ...)
  
  # computes necessary outputs for summary
  coda_summary <- fnames <- gnames <- cnames <- list()
  for(ch in object$chains){
    coda_summary[[ch]] <- summary(coda[[ch]])
  
    settings <- object$settings[[ch]]
    output <- ifelse(length(intersect(rownames(settings), c("deviance", "dev_i", "pUig_int"))) > 0,
                     "probs", "mcmc")
  
    if(output == "probs"){
      fnames[[ch]] <- c("deviance")
      gnames[[ch]] <- cnames[[ch]] <- c()
    }else{
      fset = settings[settings$save & !settings$gspec & settings$dims>0, ]
      fset = fset[setdiff(rownames(fset), c("loglik", "U", "Gplus", "e0")), ]
      fnames[[ch]] <- c()
      for(p in rownames(fset)){
        fnames[[ch]] <- c(fnames[[ch]], object$param_names[[p]])
      }
      
      gnames[[ch]] <- rownames(coda_summary[[1]]$statistics)[
        grep("\\(1\\)", rownames(coda_summary[[1]]$statistics))]
      gnames[[ch]] <- gnames[[ch]][!(grepl("^w", gnames[[ch]]) | 
                                       grepl("^ng", gnames[[ch]]) | 
                                       grepl("^pUig", gnames[[ch]]))]
      
      cnames[[ch]] <- c("Gplus", "e0", paste0("w(",1:object$G,")"))
    }
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