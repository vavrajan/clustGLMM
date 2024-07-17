predict.clusglmm <- 
function(object, newdata, type="f", posterior="E", chain=1, gs=object$clusters[[chain]]){
  coda_type <- switch(posterior,
                      "E" = "statistics",
                      "median" = "quantiles")
  coda_name <- switch(posterior,
                      "E" = "Mean",
                      "median" = "50%")
  object_as_coda <- as_coda(object)
  coda_summary <- summary(object_as_coda[[ch]])[[coda_type]]
  
  res <- list()
  if(type == "f"){
    for(y in names(object$family)){
      fX <- model.matrix(object$formula[[y]]$fixed, newdata)
      fX <- fX[, object$lfixnames[[y]]]
      
      frows <- grep(paste0("beta_",object$family[[y]],"_fix_", y, "\\["), rownames(coda_summary))
      beta_fix <- coda_summary[frows, coda_name]
      res[[y]] <- fX %*% beta_fix
    }
  }
  
  if(type == "g"){
    for(y in names(object$family)){
      gX <- model.matrix(object$formula[[y]]$group, newdata)
      gX <- gX[, object$lgrpnames[[y]]]
      
      res[[y]] <- matrix(0, nrow=dim(newdata)[1], ncol=length(gs))
      
      for(ig in 1:length(gs)){
        g <- gs[ig]
        grows <- grep(paste0("beta_",object$family[[y]],"_", y, "\\(", g, "\\)"), rownames(coda_summary))
        beta_grp <- coda_summary[grows, coda_name]
        res[[y]][, ig] <- gX %*% beta_grp
      }
      colnames(res[[y]]) <- paste0("(", gs, ")")
    }
  }
  
  if(type == "fg"){
    for(y in names(object$family)){
      fX <- model.matrix(object$formula[[y]]$fixed, newdata)
      gX <- model.matrix(object$formula[[y]]$group, newdata)
      fX <- fX[, object$lfixnames[[y]]]
      gX <- gX[, object$lgrpnames[[y]]]
      
      frows <- grep(paste0("beta_",object$family[[y]],"_fix_", y, "\\["), rownames(coda_summary))
      beta_fix <- coda_summary[frows, coda_name]
      res_fix <- fX %*% beta_fix
      res[[y]] <- matrix(0, nrow=dim(newdata)[1], ncol=length(gs))
      
      for(ig in 1:length(gs)){
        g <- gs[ig]
        grows <- grep(paste0("beta_",object$family[[y]],"_", y, "\\(", g, "\\)"), rownames(coda_summary))
        beta_grp <- coda_summary[grows, coda_name]
        res[[y]][, ig] <- res_fix + gX %*% beta_grp
      }
      colnames(res[[y]]) <- paste0("(", gs, ")")
    }
  }

  if(type == "as_parametric_function"){
    for(y in names(object$family)){
      fX <- model.matrix(object$formula[[y]]$fixed, newdata)
      gX <- model.matrix(object$formula[[y]]$group, newdata)
      fX <- fX[, object$lfixnames[[y]]]
      gX <- gX[, object$lgrpnames[[y]]]
      
      beta_fix_name <- paste0("beta_",object$family[[y]],"_fix")
      beta_grp_name <- paste0("beta_",object$family[[y]])
          
      fnames <- rownames(coda_summary)[grep(paste0(beta_fix_name, "_", y, "\\["), rownames(coda_summary))]
      if(object$howsave == "list"){
        beta_fix <- object$draws[[chain]][[beta_fix_name]][[y]]
      }
      if(object$howsave == "matrix"){
        beta_fix <- object$draws[object$draws$chain == chain, fnames]
      }
      
      res_fix <- fX %*% t(beta_fix)
      
      res[[y]] <- matrix(0, nrow=dim(newdata)[1], ncol=length(gs))
      colnames(res[[y]]) <- paste0("(",gs,")")
      
      for(g in gs){
        gnames <- rownames(coda_summary)[grep(paste0(beta_grp_name, "_", y, "\\(", g, "\\)"), rownames(coda_summary))]
        if(object$howsave == "list"){
          beta_grp <- object$draws[[chain]][[beta_grp_name]][[g]][[y]]
        }
        if(object$howsave == "matrix"){
          beta_grp <- object$draws[object$draws$chain == chain, gnames]
        }
        
        samples <- res_fix + gX %*% t(beta_grp)
        res[[y]][, paste0("(",g,")")] <- switch(coda_name,
                                                "Mean" = apply(samples, 1, mean, na.rm=TRUE),
                                                "50%" = apply(samples, 1, median, na.rm=TRUE))
      }
    }
  }
}