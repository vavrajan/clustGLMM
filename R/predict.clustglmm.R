predict.clustglmm <- 
  function(object, newdata, y=names(object$family)[1],
           type=c("link", "response", "posterior_predictive_response"), 
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
  
  type <- match.arg(type, c("link", "response", "posterior_predictive_response"))
  method <- match.arg(method, c("samples", "plugin"))
  what <- match.arg(what, c("fg", "f", "g"))
  interval <- match.arg(interval, c("ET", "HPD"))
  
  alpha_half <- (1-level)/2
  params <- colnames(object$draws[[chain]])[-1]
  draws <- object$draws[[chain]][object$draws[[chain]]$m > burnin, params]
  draws <- draws[seq(1, dim(draws)[1], by=thin), ]
  settings <- object$settings[[chain]]
  
  maxK <- 1
  if(object$family[y] == "ord"){maxK <- object$Kord[y]}
  if(object$family[y] == "cat"){maxK <- object$Kcat[y]}
  
  res <- list()
  
  ### Estimates based on transforming every iteration
  if(method == "samples"){
    
    fX <- model.matrix(object$formula[[y]]$fixed, newdata)
    gX <- model.matrix(object$formula[[y]]$group, newdata)
    rX <- model.matrix(object$formula[[y]]$random, newdata)
    fX <- fX[, object$lfixnames[[y]]]
    gX <- gX[, object$lgrpnames[[y]]]
    rX <- rX[, object$lrannames[[y]]]
    ry <- grep(paste0("^",y,"[0-9]+"), unlist(object$lrannames))
    
    beta_fix_name <- paste0("beta_",object$family[y],"_fix")
    beta_grp_name <- paste0("beta_",object$family[y])
    
    # samples for linear predictor
    samples <- list()
    samples[["f"]] <- array(0, c(dim(newdata)[1], # newdata dimension
                                 1,               # group/cluster dimension
                                 dim(draws)[1],   # draws dimension
                                 maxK),           # level dimension (relevant only for Ord and Cat)
                            dimnames = list(1:dim(newdata)[1],
                                            1,
                                            1:dim(draws)[1],
                                            paste0("k=",1:maxK)))           
    samples[["g"]] <- samples[["fg"]] <- 
      array(0, c(dim(newdata)[1], # newdata dimension
                 length(gs),      # group/cluster dimension
                 dim(draws)[1],   # draws dimension
                 maxK),           # level dimension (relevant only for Ord and Cat)
            dimnames = list(1:dim(newdata)[1],
                            paste0("(",gs,")"),
                            1:dim(draws)[1],
                            paste0("k=",1:maxK))) 
    samples[["response"]] <-
      array(0, c(dim(newdata)[1], # newdata dimension
                 length(gs),      # group/cluster dimension
                 dim(draws)[1],   # draws dimension
                 ifelse(is.element(object$family[y], c("num", "poi")),
                        1, maxK+1)),           # level dimension (relevant only for Ord and Cat)
            dimnames = list(1:dim(newdata)[1],
                            paste0("(",gs,")"),
                            1:dim(draws)[1],
                            paste0("k=",1:ifelse(is.element(object$family[y], c("num", "poi")),
                                                 1, maxK+1)))) # 1 more dimension for probs
    
    samples[["o"]] <- array(0, dim = c(dim(newdata)[1], 1, 1))
    if(length(object$loffnames[[y]]) == 1){
      samples[["o"]][, 1, 1] <- newdata[, object$loffnames[[y]]]
    }
    
    for(k in 1:maxK){
      if(object$family[y] == "cat"){
        nametag <- paste0("\\[[0-9]+,",k,"\\]")
      }else{
        nametag <- "\\[[0-9]+\\]"
      }
      
      fnames <- params[grep(paste0(beta_fix_name, "_", y, nametag), params)]
      beta_fix <- draws[, fnames]
      samples[["f"]][, 1, , k] <- fX %*% t(beta_fix)
      
      for(ig in 1:length(gs)){
        g <- gs[ig]
        gnames <- params[grep(paste0(beta_grp_name, "_", y, "\\(", g, "\\)", nametag), params)]
        beta_grp <- draws[, gnames]
        samples[["g"]][, ig, , k] <- gX %*% t(beta_grp)
      }
      
      if(object$family[y] == "ord"){
        if(object$varying["c_ord"]){
          for(ig in 1:length(gs)){
            g <- gs[ig]
            c_ord_name <- paste0("c_ord_", y, "(", g, ")[", k, "]")
            samples[["g"]][, ig, , k] <- samples[["g"]][, ig, , k] - draws[, c_ord_name]
          }
        }else{
          c_ord_name <- paste0("c_ord_", y, "[", k, "]")
          samples[["f"]][, 1, , k] <- samples[["f"]][, 1, , k] - draws[, c_ord_name]
        }
      }
    }
    
    for(ig in 1:length(gs)){
      samples[["fg"]][, ig, , ] <- samples[["g"]][, ig, , ] + samples[["f"]][, 1, , ]
    }
    
    
    
    
    # Just linear predictor
    if(type == "link"){
      print("ahoj type 1")
      res[["fit"]] <- apply(samples[[what]], c(1,2,4), fun_posterior)[,,,drop=TRUE]
      if(interval == "ET"){
        res[["lwr"]] <- apply(samples[[what]], c(1,2,4), quantile, probs = alpha_half)[,,,drop=TRUE]
        res[["upr"]] <- apply(samples[[what]], c(1,2,4), quantile, probs = 1-alpha_half)[,,,drop=TRUE]
      }
      if(interval == "HPD"){
        hpdis <- apply(samples[[what]], c(1,2,4), hdi, credMass = level)
        res[["lwr"]] <- hpdis[1,,,,drop=TRUE]
        res[["upr"]] <- hpdis[2,,,,drop=TRUE]
      }
    }
    
    # On response level
    if(type == "response"){
      print("ahoj type 2")
      # Depending on type of outcome:
      if(object$family[y] == "num"){
        samples[["response"]] <- samples[["fg"]] # ? + samples[["o"]]
      }
      
      if(object$family[y] == "poi"){
        samples[["response"]] <- exp(samples[["fg"]]) # ? # + samples[["o"]]
      }
      
      if(object$family[y] == "bin"){
        exp_eta <- exp(samples[["fg"]]) # ? # + samples[["o"]]
        samples[["response"]][, , , 2] <- exp_eta / (1 + exp_eta)
        samples[["response"]][, , , 1] <- 1 - samples[["response"]][, , , 2]
      }
      
      if(object$family[y] == "ord"){
        for(m in 1:dim(draws)[1]){
          for(ig in 1:length(gs)){
            exp_eta <- exp(samples[["fg"]][, ig, m, ]) # ? # + samples[["o"]]
            pk <- cbind(rep(0, dim(newdata)[1]), 
                        exp_eta / (1 + exp_eta),
                        rep(1, dim(newdata)[1]))
            samples[["response"]][, ig, m, ] <- pk[,2:(maxK+2)] - pk[,1:(maxK+1)]
          }
        }
      }
      
      if(object$family[y] == "cat"){
        for(m in 1:dim(draws)[1]){
          for(ig in 1:length(gs)){
            exp_eta <- exp(cbind(rep(0,dim(newdata)[1]), # 0 for baseline
                                 samples[["fg"]][, ig, m, ]))  # ? # + samples[["o"]]
            softmax <- exp_eta / apply(exp_eta, 1, sum)
            samples[["response"]][, ig, m, ] <- softmax
          }
        }
      }
      
      res[["fit"]] <- apply(samples[["response"]], c(1,2,4), fun_posterior)[,,,drop=TRUE]
      if(interval == "ET"){
        res[["lwr"]] <- apply(samples[["response"]], c(1,2,4), quantile, probs = alpha_half)[,,,drop=TRUE]
        res[["upr"]] <- apply(samples[["response"]], c(1,2,4), quantile, probs = 1-alpha_half)[,,,drop=TRUE]
      }
      if(interval == "HPD"){
        hpdis <- apply(samples[["response"]], c(1,2,4), hdi, credMass = level)
        res[["lwr"]] <- hpdis[1,,,,drop=TRUE]
        res[["upr"]] <- hpdis[2,,,,drop=TRUE]
      }
      
    }
    
    # Posterio predictive distribution of the response
    if(type == "posterior_predictive_response"){
      print("ahoj type 3")
      samples[["y"]] <- array(dim = c(dim(newdata)[1], 1, dim(draws)[1], 1))
      samples[["aux"]] <- array(dim = c(dim(newdata)[1], 1, dim(draws)[1], 
                                        ifelse(is.element(object$family[y], c("num", "poi")),
                                               1, maxK+1)))
      
      for(m in 1:dim(draws)[1]){
        samples[["r"]] <- array(dim = c(dim(newdata)[1], 1, 1))
        samples[["eta"]] <- array(dim = c(dim(newdata)[1], 1, maxK))
        
        # get w
        w <- draws[m, paste0("w(",gs,")")]
        
        # get Sigma matrix
        if(settings["InvSigma","gspec"]){
          Sigma <- list()
          for(g in gs){
            Sigma[[g]] <- get_Sigma(draws, m, settings, g = g, howsave = "matrix")
          }
        }else{
          Sigma <- get_Sigma(draws, m, settings, howsave = "matrix")
        }
        
        ### Generate latent data - U and b
        ### !!! only clusters in gs are considered here !!!
        U <- as.numeric(sample(as.character(gs), dim(newdata)[1], replace = TRUE, prob = w))
        if(settings["InvSigma","gspec"]){
          b <- t(sapply(U, function(u){rmvnorm(1, 
                                               mean = rep(0, object$totnran),
                                               sigma = Sigma[[u]])}))
        }else{
          b <- rmvnorm(dim(newdata)[1],
                       mean = rep(0, object$totnran),
                       sigma = Sigma)
        }
        # only random effects for outcome y
        b <- b[, ry, drop=FALSE]
        samples[["r"]][, 1, 1] <- apply(b * rX, 1, sum)
        # TODO kspec_bi_cat!!!
        
        
        # reconstruct predictor eta
        for(ig in 1:length(gs)){
          for(k in 1:maxK){
            for(j in 1:dim(newdata)[1]){
              samples[["eta"]][j, 1, k] <- 
                samples[["fg"]][j, U[j], m, k] + samples[["r"]][j, 1, 1] + samples[["o"]][j, 1, 1]
            }
          }
        }
        
        
        ### Sampling outcomes
        if(object$family[y] == "num"){
          sd_num <- get_sd_num(draws, m, settings, y = y, howsave="matrix")
          if(settings["prec_num","gspec"]){
            sd <- sd_num[U]
          }else{
            sd <- sd_num
          }
          samples[["y"]][, 1, m, 1] <- rnorm(dim(newdata)[1], 
                                             mean = samples[["eta"]][, 1, 1], 
                                             sd = sd)
          samples[["aux"]][, 1, m, 1] <- samples[["eta"]][, 1, 1]
        }
        
        if(object$family[y] == "poi"){
          exp_eta <- exp(samples[["eta"]][1, , 1])
          samples[["y"]][, 1, m, 1] <- rpois(dim(newdata)[1], lambda = exp_eta)
          samples[["aux"]][, 1, m, 1] <- exp_eta
          
        }
        
        if(object$family[y] == "bin"){
          exp_eta <- exp(samples[["eta"]][, 1, 1])
          prob <- exp_eta / (1 + exp_eta)
          samples[["y"]][, 1, m, 1] <- rbinom(dim(newdata)[1], 1, prob = prob) 
          samples[["aux"]][, 1, m, 2] <- prob
          samples[["aux"]][, 1, m, 1] <- 1 - prob
        }
        
        if(object$family[y] == "ord"){
          exp_eta <- exp(samples[["eta"]][, 1, ])
          pk <- cbind(rep(0, dim(newdata)[1]), 
                      exp_eta / (1 + exp_eta),
                      rep(1, dim(newdata)[1]))
          qk <- pk[,2:(maxK+2)] - pk[,1:(maxK+1)]
          for(j in 1:dim(newdata)[1]){
            # auxU <- runif(1)
            # samples[["y"]][j, 1, m, 1] <- sum(pk[j,] < auxU)
            samples[["y"]][j, 1, m, 1] <- sample(1:(maxK+1), 1, prob = qk[j,])
          }
          samples[["aux"]][, 1, m, ] <- qk
        }
        
        if(object$family[y] == "cat"){
          exp_eta <- exp(cbind(rep(0,dim(newdata)[1]), # 0 for baseline
                               samples[["eta"]][, 1, ]))
          softmax <- exp_eta / apply(exp_eta, 1, sum)
          for(j in 1:dim(newdata)[1]){
            samples[["y"]][j, 1, m, 1] <- sample(1:(maxK+1), 1, prob = softmax[j,])
          }
          samples[["aux"]][, 1, m, ] <- softmax
        }
      }
      
      ### return results
      # summarized samples
      res[["fit"]] <- apply(samples[["y"]], c(1,2,4), fun_posterior)[,,,drop=TRUE]
      if(interval == "ET"){
        res[["lwr"]] <- apply(samples[["y"]], c(1,2,4), quantile, probs = alpha_half)[,,,drop=TRUE]
        res[["upr"]] <- apply(samples[["y"]], c(1,2,4), quantile, probs = 1-alpha_half)[,,,drop=TRUE]
      }
      if(interval == "HPD"){
        hpdis <- apply(samples[["y"]], c(1,2,4), hdi, credMass = level)
        res[["lwr"]] <- hpdis[1,,,,drop=TRUE]
        res[["upr"]] <- hpdis[2,,,,drop=TRUE]
      }
      # all Y samples
      res[["samples"]] <- samples[["y"]][,,,,drop=TRUE]
      # aux: for Num: eta
      #      for Poi: exp(eta)
      #      for Bin: prob
      #      for Ord: probs qk
      #      for Cat: probs softmax
      res[["aux"]] <- samples[["aux"]][,,,,drop=TRUE]
    }
  }
  # end of method "samples"
  
  
  
  ### Estimates based on transforming an estimate of the parameters
  if(method == "plugin"){
    fX <- model.matrix(object$formula[[y]]$fixed, newdata)
    gX <- model.matrix(object$formula[[y]]$group, newdata)
    fX <- fX[, object$lfixnames[[y]]]
    gX <- gX[, object$lgrpnames[[y]]]
    
    beta_fix_name <- paste0("beta_",object$family[[y]],"_fix")
    beta_grp_name <- paste0("beta_",object$family[[y]])
    
    # samples for linear predictor
    samples <- list()
    samples[["f"]] <- array(0, c(dim(newdata)[1], # newdata dimension
                                 1,               # group/cluster dimension
                                 maxK),           # level dimension (relevant only for Ord and Cat)
                            dimnames = list(1:dim(newdata)[1],
                                            1,
                                            paste0("k=",1:maxK)))           
    samples[["g"]] <- samples[["fg"]] <- 
      array(0, c(dim(newdata)[1], # newdata dimension
                 length(gs),      # group/cluster dimension
                 maxK),           # level dimension (relevant only for Ord and Cat)
            dimnames = list(1:dim(newdata)[1],
                            paste0("(",gs,")"),
                            paste0("k=",1:maxK)))
    samples[["response"]] <-
      array(0, c(dim(newdata)[1], # newdata dimension
                 length(gs),      # group/cluster dimension
                 ifelse(is.element(object$family[y], c("num", "poi")),
                        1, maxK+1)),           # level dimension (relevant only for Ord and Cat)
            dimnames = list(1:dim(newdata)[1],
                            paste0("(",gs,")"),
                            paste0("k=",1:ifelse(is.element(object$family[y], c("num", "poi")),
                                                 1, maxK+1))))
    
    for(k in 1:maxK){
      if(object$family[y] == "cat"){
        nametag <- paste0("\\[[0-9]+,",k,"\\]")
      }else{
        nametag <- "\\[[0-9]+\\]"
      }
      
      fnames <- params[grep(paste0(beta_fix_name, "_", y, nametag), params)]
      beta_fix <- apply(draws[, fnames], 2, fun_posterior)
      samples[["f"]][, 1, k] <- fX %*% beta_fix
      
      for(ig in 1:length(gs)){
        g <- gs[ig]
        gnames <- params[grep(paste0(beta_grp_name, "_", y, "\\(", g, "\\)", nametag), params)]
        beta_grp <- apply(draws[, gnames], 2, fun_posterior)
        samples[["g"]][, ig, k] <- gX %*% beta_grp
      }
      
      if(object$family[y] == "ord"){
        if(object$varying["c_ord"]){
          for(ig in 1:length(gs)){
            g <- gs[ig]
            c_ord_name <- paste0("c_ord_", y, "(", g, ")[", k, "]")
            samples[["g"]][, ig, k] <- samples[["g"]][, ig, k] - fun_posterior(draws[, c_ord_name])
          }
        }else{
          c_ord_name <- paste0("c_ord_", y, "\\[", k, "\\]")
          samples[["f"]][, 1, k] <- samples[["f"]][, 1, k] - fun_posterior(draws[, c_ord_name])
        }
      }
      
    }
    
    for(ig in 1:length(gs)){
      samples[["fg"]][, ig, ] <- samples[["g"]][, ig, ] + samples[["f"]][, 1, ]
    }
    
    
    res <- list()
    
    # Just linear predictor
    if(type == "link"){
      res[["fit"]] <- samples[[what]][,,,drop=TRUE]
      
      # lwr and upr
      if(interval == "ET"){
        # todo
        res[["lwr"]] <- NULL
        res[["upr"]] <- NULL
        warning("Lower and upper bounds for method='plugin' not implemented.")
      }
      if(interval == "HPD"){
        # todo
        res[["lwr"]] <- NULL
        res[["upr"]] <- NULL
        warning("Lower and upper bounds for method='plugin' not implemented.")
      }
    }
    
    # On outcome level
    if(type == "response"){
      # Depending on the type of outcome
      if(object$family[y] == "num"){
        samples[["response"]] <- samples[["fg"]] # ? + samples[["o"]]
      }
      
      if(object$family[y] == "poi"){
        samples[["response"]] <- exp(samples[["fg"]]) # ? # + samples[["o"]]
      }
      
      if(object$family[y] == "bin"){
        exp_eta <- exp(samples[["fg"]]) # ? # + samples[["o"]]
        samples[["response"]] <- exp_eta / (1 + exp_eta)
      }
      
      if(object$family[y] == "ord"){
        for(ig in 1:length(gs)){
          exp_eta <- exp(samples[["fg"]][, ig, ]) # ? # + samples[["o"]]
          pk <- cbind(rep(0, dim(newdata)[1]), 
                      exp_eta / (1 + exp_eta),
                      rep(1, dim(newdata)[1]))
          samples[["response"]][, ig, ] <- pk[,2:(maxK+2)] - pk[,1:(maxK+1)]
        }
      }
      
      if(object$family[y] == "cat"){
        for(ig in 1:length(gs)){
          exp_eta <- exp(cbind(rep(0,dim(newdata)[1]), # 0 for baseline
                               samples[["fg"]][, ig, ]))  # ? # + samples[["o"]]
          softmax <- exp_eta / apply(exp_eta, 1, sum)
          samples[["response"]][, ig, ] <- softmax
        }
      }
      
      # return as fit
      res[["fit"]] <- samples[["response"]][,,,drop=TRUE]
      
      # lwr and upr
      if(interval == "ET"){
        # todo
        res[["lwr"]] <- NULL
        res[["upr"]] <- NULL
        warning("Lower and upper bounds for method='plugin' not implemented.")
      }
      if(interval == "HPD"){
        # todo
        res[["lwr"]] <- NULL
        res[["upr"]] <- NULL
        warning("Lower and upper bounds for method='plugin' not implemented.")
      }
    }
    
    # Posterior predictive distribution
    if(type=="posterior_predictive_response"){
      stop("type='posterior_predictive_response' not implemented for method='plugin', use method='samples' instead.")
    }
    
    
  }
  
  return(res)
  
}