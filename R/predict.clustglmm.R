predict.clustglmm <- 
function(object, newdata, y=names(object$family)[1], k=NA,
         type=c("link", "response"), 
         method=c("samples", "plugin"), 
         what=c("fg", "f", "g"),
         fun_posterior=mean, 
         level=0.95, 
         interval=c("ET", "HPD"),
         chain=1, burnin=0, thin=1, gs=object$clusters[[chain]], ...){
  # we rather work with "matrix" format here
  if(object$howsave == "list"){
    object <- from_list_to_matrix(object)
  }

  type <- match.arg(type, c("link", "response"))
  method <- match.arg(method, c("samples", "plugin"))
  what <- match.arg(what, c("fg", "f", "g"))
  interval <- match.arg(interval, c("ET", "HPD"))
    
  alpha_half <- (1-level)/2
  params <- colnames(object$draws[[chain]])[-1]
  draws <- object$draws[[chain]][object$draws[[chain]]$m > burnin, params]
  draws <- draws[seq(1, dim(draws)[1], by=thin), ]
  if(object$family[y] == "ord"){
    if(missing(k) | (is.na(k)) | !is.element(k, 1:object$Kord[y])){
      warning(paste0("Outcome ", y, " is ordinal, but level k (for the choice of intercept) not given. Only values 1, ..., ", object$Kord[y], " are permitted. I continue with k=1."))
      k=1
    }
  }
  
  if(object$family[y] == "cat"){
    if(missing(k) | (is.na(k)) | !is.element(k, 1:object$Kcat[y])){
      warning(paste0("Outcome ", y, " is categorical, but level k not given. Only values 1, ..., ", object$Kcat[y], " are permitted. I continue with k=1."))
      k=1
    }
  }
  
  if(object$family[y] == "cat"){
    nametag <- paste0("\\[[0-9]+,",k,"\\]")
  }else{
    nametag <- "\\[[0-9]+\\]"
  }
  
  if(method == "samples"){
    if(type == "link"){
      fX <- model.matrix(object$formula[[y]]$fixed, newdata)
      gX <- model.matrix(object$formula[[y]]$group, newdata)
      fX <- fX[, object$lfixnames[[y]]]
      gX <- gX[, object$lgrpnames[[y]]]
      
      beta_fix_name <- paste0("beta_",object$family[[y]],"_fix")
      beta_grp_name <- paste0("beta_",object$family[[y]])
      
      fnames <- params[grep(paste0(beta_fix_name, "_", y, nametag), params)]
      beta_fix <- draws[, fnames]
      res_fix <- fX %*% t(beta_fix)
      
      samples <- array(0, c(dim(newdata)[1], 
                            ifelse(what=="f", 1, length(gs)), 
                            dim(draws)[1]))
      if(is.element(what, c("f", "fg"))){
        for(j in 1:dim(samples)[2]){
          samples[, j, ] <- res_fix
        }
      }
      
      if(is.element(what, c("g", "fg"))){
        for(ig in 1:length(gs)){
          g <- gs[ig]
          gnames <- params[grep(paste0(beta_grp_name, "_", y, "\\(", g, "\\)", nametag), params)]
          beta_grp <- draws[, gnames]
          
          samples[, ig, ] <- samples[, ig, ] + gX %*% t(beta_grp)
        }
      }
      
      if(object$family[y] == "ord"){
        if(object$varying["c_ord"]){
          for(ig in 1:length(gs)){
            g <- gs[ig]
            c_ord_name <- paste0("c_ord_", y, "(", g, ")[", k, "]")
            samples[, ig, ] <- samples[, ig, ] - draws[, c_ord_name]
          }
        }else{
          c_ord_name <- paste0("c_ord_", y, "[", k, "]")
          for(j in 1:dim(samples)[2]){
            samples[, j, ] <- samples[, j, ] - draws[, c_ord_name]
          }
        }
      }
      
      res <- list()
      res[["fit"]] <- apply(samples, c(1,2), fun_posterior)
      if(interval == "ET"){
        res[["lwr"]] <- apply(samples, c(1,2), quantile, probs = alpha_half)
        res[["upr"]] <- apply(samples, c(1,2), quantile, probs = 1-alpha_half)
      }
      if(interval == "HPD"){
        hpdis <- apply(samples, c(1,2), hdi, credMass = level)
        res[["lwr"]] <- hpdis[1,,1]
        res[["upr"]] <- hpdis[2,,1]
      }
      
      if(what == "f"){
        res <- as.data.frame(res)
      }
      if(is.element(what, c("g", "fg"))){
        colnames(res[["fit"]]) <- paste0("(",gs,")")
        if(is.element(interval, c("ET", "HPD"))){
          colnames(res[["lwr"]]) <- colnames(res[["upr"]]) <- paste0("(",gs,")")
        }
      }
    }else{
      if(type == "response"){
        stop("type 'response' not implemented yet.")
      }else{
        stop("Unknown type for predicting with method='samples'.")
      }
    }
  }
  
  if(method == "plugin"){
    if(type == "link"){
      fX <- model.matrix(object$formula[[y]]$fixed, newdata)
      gX <- model.matrix(object$formula[[y]]$group, newdata)
      fX <- fX[, object$lfixnames[[y]]]
      gX <- gX[, object$lgrpnames[[y]]]
      
      beta_fix_name <- paste0("beta_",object$family[[y]],"_fix")
      beta_grp_name <- paste0("beta_",object$family[[y]])
      
      fnames <- params[grep(paste0(beta_fix_name, "_", y, nametag), params)]
      beta_fix <- apply(draws[, fnames], 2, fun_posterior)
      res_fix <- fX %*% beta_fix
      
      res_grp <- matrix(0, nrow=dim(newdata)[1], ncol=ifelse(what=="f", 1, length(gs)))
      
      for(ig in 1:length(gs)){
        g <- gs[ig]
        gnames <- params[grep(paste0(beta_grp_name, "_", y, "\\(", g, "\\)", nametag), params)]
        beta_grp <- apply(draws[, gnames], 2, fun_posterior)
        res_grp[, ig] <- gX %*% beta_grp
      }
      colnames(res_grp) <- paste0("(", gs, ")")
      
      res <- list()
      # fit
      if(what == "f"){
        res[["fit"]] <- res_fix
      }
      if(what == "g"){
        res[["fit"]] <- res_grp
      }
      if(what == "fg"){
        res[["fit"]] <- matrix(rep(res_fix, dim(res_grp)[2]), ncol=dim(res_grp)[2]) + res_grp
      }
      if(object$family[y] == "ord"){
        if(object$varying["c_ord"]){
          if(what != "f"){
            for(ig in 1:length(gs)){
              g <- gs[ig]
              c_ord_name <- paste0("c_ord_", y, "(", g, ")[", k, "]")
              res[["fit"]][,ig] <- res[["fit"]][,ig] - fun_posterior(draws[, c_ord_name])
            }
          }
        }else{
          c_ord_name <- paste0("c_ord_", y, "\\[", k, "\\]")
          for(j in 1:dim(res[["fit"]])[2]){
            res[["fit"]][, j] <- res[["fit"]][, j] - fun_posterior(draws[, c_ord_name])
          }
        }
      }
      
      # lwr and upr
      if(interval == "ET"){
        # todo
        res[["lwr"]] <- NULL
        res[["upr"]] <- NULL
        warning("Lower and upper bounds for method='plugin' not implemented (yet?).")
      }
      if(interval == "HPD"){
        # todo
        res[["lwr"]] <- NULL
        res[["upr"]] <- NULL
        warning("Lower and upper bounds for method='plugin' not implemented (yet?).")
      }
    }else{
      if(type == "response"){
        stop("type 'response' not implemented yet.")
      }else{
        stop("Unknown type for predicting with method='samples'.")
      }
    }
  }
  
  return(res)
}
