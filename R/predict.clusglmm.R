predict.clusglmm <- 
function(object, newdata, 
         type="lp-fg", how="samples", fun_posterior=mean, 
         level=0.95, interval="ET",
         chain=1, burnin=0, thin=1, gs=object$clusters[[chain]]){
  # we rather work with "matrix" format here
  if(object$howsave == "list"){
    object <- from_list_to_matrix(object)
  }
  
  # allow type lp without specification
  if(type=="lp"){
    type <- "lp-fg"
  }
  
  alpha_half <- (1-level)/2
  params <- colnames(object$draws)[-1]
  draws <- object$draws[object$draws$chain == chain, params]
  draws <- draws[seq(burnin+1, object$iter, by=thin), ]
  
  if(how == "samples"){
    if(is.element(type, c("lp-f", "lf-g", "lp-fg"))){
      for(y in names(object$family)){
        fX <- model.matrix(object$formula[[y]]$fixed, newdata)
        gX <- model.matrix(object$formula[[y]]$group, newdata)
        fX <- fX[, object$lfixnames[[y]]]
        gX <- gX[, object$lgrpnames[[y]]]
        
        beta_fix_name <- paste0("beta_",object$family[[y]],"_fix")
        beta_grp_name <- paste0("beta_",object$family[[y]])
        
        fnames <- params[grep(paste0(beta_fix_name, "_", y, "\\["), params)]
        beta_fix <- draws[, fnames]
        res_fix <- fX %*% t(beta_fix)
        
        samples <- array(0, c(dim(newdata)[1], 
                              ifelse(type=="lp-f", 1, length(gs)), 
                              dim(draws)[1]))
        if(is.element(type, c("lp-f", "lp-fg"))){
          for(j in 1:dim(samples)[2]){
            samples[, j, ] <- res_fix
          }
        }
        
        if(is.element(type, c("lp-g", "lp-fg"))){
          for(ig in 1:length(gs)){
            g <- gs[ig]
            gnames <- params[grep(paste0(beta_grp_name, "_", y, "\\(", g, "\\)"), params)]
            beta_grp <- draws[, gnames]
            
            samples[, ig, ] <- samples[, ig, ] + gX %*% t(beta_grp)
          }
        }
        
        res[[y]] <- list()
        res[[y]][["fit"]] <- apply(samples, c(1,2), fun_posterior)
        if(interval == "ET"){
          res[[y]][["lwr"]] <- apply(samples, c(1,2), quantile, probs = alpha_half)
          res[[y]][["upr"]] <- apply(samples, c(1,2), quantile, probs = 1-alpha_half)
        }
        if(interval == "HPD"){
          hpdis <- apply(samples, c(1,2), hdi, credMass = level)
          res[[y]][["lwr"]] <- hpdis[1,,1]
          res[[y]][["upr"]] <- hpdis[2,,1]
        }
        
        if(type == "lp-f"){
          res[[y]] <- as.data.frame(res[[y]])
        }
        if(is.element(type, c("lp-g", "lp-fg"))){
          colnames(res[[y]][["fit"]]) <- paste0("(",gs,")")
          if(is.element(interval, c("ET", "HPD"))){
            colnames(res[[y]][["lwr"]]) <- colnames(res[[y]][["upr"]]) <- paste0("(",gs,")")
          }
        }
      }
    }else{
      stop("Unknown type for predicting with how=samples.")
    }
  }
  
  if(how == "plugin"){
    if(is.element(type, c("lp-f", "lf-g", "lp-fg"))){
      for(y in names(object$family)){
        fX <- model.matrix(object$formula[[y]]$fixed, newdata)
        gX <- model.matrix(object$formula[[y]]$group, newdata)
        fX <- fX[, object$lfixnames[[y]]]
        gX <- gX[, object$lgrpnames[[y]]]
        
        beta_fix_name <- paste0("beta_",object$family[[y]],"_fix")
        beta_grp_name <- paste0("beta_",object$family[[y]])
        
        fnames <- params[grep(paste0(beta_fix_name, "_", y, "\\["), params)]
        beta_fix <- apply(draws[, fnames], 2, fun_posterior)
        res_fix <- fX %*% beta_fix
        
        res_grp <- matrix(0, nrow=dim(newdata)[1], ncol=ifelse(type=="lp-f", 1, length(gs)))
        
        for(ig in 1:length(gs)){
          g <- gs[ig]
          gnames <- params[grep(paste0(beta_fix_name, "_", y, "\\(", g, "\\)"), params)]
          beta_grp <- apply(draws[, gnames], 2, fun_posterior)
          res_grp[, ig] <- gX %*% beta_grp
        }
        colnames(res_grp) <- paste0("(", gs, ")")
        
        res[[y]] <- list()
        # fit
        if(type == "lp-f"){
          res[[y]][["fit"]] <- res_fix
        }
        if(type == "lp-g"){
          res[[y]][["fit"]] <- res_grp
        }
        if(type == "lp-fg"){
          res[[y]][["fit"]] <- res_fix + res_grp
        }
        # lwr and upr
        if(interval == "ET"){
          # todo
          res[[y]][["lwr"]] <- NULL
          res[[y]][["upr"]] <- NULL
          warning("Lower and upper bounds for how=plugin not implemented (yet?).")
        }
        if(interval == "HPD"){
          # todo
          res[[y]][["lwr"]] <- NULL
          res[[y]][["upr"]] <- NULL
          warning("Lower and upper bounds for how=plugin not implemented (yet?).")
        }
      }
    }else{
      stop("Unknown type for predicting with how=plugin.")
    }
  }
}