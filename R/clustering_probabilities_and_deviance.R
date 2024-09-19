clustering_probabilities_and_deviance <-
function(mcmc,      # MCMC generated parameters
         #   output of Metropolis_within_Gibbs_MBC_NumBinOrdCat
         id,     # name of the id variable in data
         data,   # outcomes, regressors and ids of subjects
         dynamic_prob = FALSE, # Should probabilities be calculated dynamically?
         #dev = T,
         start = 1,
         end = mcmc$iter,
         thin = 1,
         chains = 1:mcmc$nchains,
         tolerance = mcmc$tuning$double$tolerance, # tolerance when computing norm of the shift within Newton-Raphson
         maxiter = mcmc$tuning$integer$maxiter, # maximum iterations allowed during Newton-Raphson update of proposal distribution
         maxnrep = mcmc$tuning$integer$maxnrep,
         NGQ = 1 # number of points for Adaptive Gaussian Quadrature
         # default value 1 corresponds to Laplacian approximation
         # do not use values higher than 7, since there are NGQ^totnran summands
){
  ### Every auxiliary values are stored within mcmc output
  # Such as save, calc, varying, ...
  # message("\nR function pUig_dev succesfully called.\n")
  
  # Find roots and weights for Adaptive Gaussian Quadrature
  # require("gaussquad")
  # require("polynom")
  # require("orthopolynom")
  if(NGQ == 1){
    rules <- hermite.he.quadrature.rules(n = 2, normalized=FALSE)[[1]]
  }else{
    rules <- hermite.he.quadrature.rules(n = NGQ, normalized=FALSE)[[NGQ]]
  }
  
  # Iterations to be taken
  iterations <- seq(start, end, by = thin)
  niter <- length(iterations)
  
  # For probability calculation we need the following parameters:
  #   w, tau, beta, InvSigma
  # If those were not saved during the mcmc_MBC_GLMMs function,
  # then we cannot use this function to calculate it.
  if(!(mcmc$save["w"] & mcmc$save["c_ord"] & mcmc$save["tau_num"] 
       & mcmc$save["beta_num_fix"] & mcmc$save["beta_num"] 
       & mcmc$save["beta_poi_fix"] & mcmc$save["beta_poi"] 
       & mcmc$save["beta_bin_fix"] & mcmc$save["beta_bin"] 
       & mcmc$save["beta_ord_fix"] & mcmc$save["beta_ord"] 
       & mcmc$save["beta_cat_fix"] & mcmc$save["beta_cat"] 
       & mcmc$save["InvSigma"])){
    stop("Some of the parameters w, c_ord, tau_num, betas, InvSigma was not saved.
Change the  settings and use mcmc_MBC_GLMMs() again.")
  }
  
  # Calculation will be done by C function, to speed up calculation
  # parameters need to be given to C function in a similar way to Metropolis_within_Gibbs_MBC_NumBinOrdCat
  
  n <- length(unique(data[,id]))
  N <- dim(data)[1]
  unique_ids <- unique(data[,id])
  n <- length(unique_ids)
  numbered_unique_ids <- c(1:n)
  names(numbered_unique_ids) <- unique_ids
  nsubj <- table(data[,id])
  max_n_i <- max(nsubj)
  
  Ys <- c(mcmc$Nums, mcmc$Pois, mcmc$Bins, mcmc$Ords, mcmc$Cats)
  nY <- sapply(list(mcmc$Nums, mcmc$Pois, mcmc$Bins, mcmc$Ords, mcmc$Cats), length)
  names(nY) <- c("Nums", "Pois", "Bins", "Ords", "Cats")
  
  ydata <- data
  # Numeric
  if(nY["Nums"] > 0){
    for(y in mcmc$Nums){
      yval <- data[,y]
      if(!is.numeric(yval)){
        stop(paste0("Outcome ", y, " is not numeric!"))
      }
      if(sum(infs <- is.infinite(yval))>0){
        yval[infs] <- NA
        warning(paste0("Infinite values of outcome ", y, " have been converted to NAs."))
      }
      ydata[,y] <- yval # numeric values
    }
  }
  
  # Poisson
  if(nY["Pois"] > 0){
    for(y in mcmc$Pois){
      yval <- data[,y]
      if(!is.integer(yval)){
        stop(paste0("Outcome ", y, " is not integer-valued!"))
      }
      if(sum(infs <- is.infinite(yval))>0){
        yval[infs] <- NA
        warning(paste0("Infinite values of outcome ", y, " have been converted to NAs."))
      }
      ydata[,y] <- yval # integer values
    }
  }
  
  # Binary
  if(nY["Bins"] > 0){
    for(y in mcmc$Bins){
      if(is.factor(data[,y])){
        fy <- data[,y]
      }else{
        fy <- factor(data[,y])
      }
      if(nlevels(fy) != 2){
        stop(paste0("Binary outcome ", y, " does not have 2 levels."))
      }
      ydata[,y] <- as.numeric(fy)-1 # 0, 1 values
    }
  }
  
  # Ordinal
  if(nY["Ords"] > 0){
    for(y in mcmc$Ords){
      if(is.factor(data[,y])){
        fy <- data[,y]
      }else{
        fy <- factor(data[,y])
      }
      if(!is.ordered(fy)){
        fy <- factor(fy, ordered = TRUE)
        warning(paste0("Converting ordinal outcome ", y, " to ordered factor with levels: ",
                       paste(levels(fy), collapse = " < ")))
      }
      if(nlevels(fy) <= 2){
        stop(paste0("Ordinal outcome ", y, " has less than 3 levels."))
      }
      ydata[,y] <- as.numeric(fy)-1 # 0 < 1 < ... < Kord[y]
    }
  }
  
  # Categorical
  if(nY["Cats"] > 0){
    for(y in mcmc$Cats){
      if(is.factor(data[,y])){
        fy <- data[,y]
      }else{
        fy <- factor(data[,y])
      }
      if(!is.factor(fy)){
        fy <- factor(fy)
        warning(paste0("Converting outcome  ", y, " into categorical factor with levels: ",
                       paste(levels(fy), collapse = ", ")))
      }
      if(nlevels(fy) <= 2){
        stop(paste0("Categorical outcome ", y, " has less than 3 levels."))
      }
      # First category is assumed to be the first one here
      # Later for C function changed to the last one
      ydata[,y] <- as.numeric(fy)-1 # 0, 1, ..., Kcat[y]
    }
  }
  
  
  # Creating model matrix containing all needed columns for C
  finX <- data.frame(id = data[,id])
  Xcolnames <- numeric()
  fcols <- gcols <- rcols <- ocols <- list()
  for(y in Ys){
    fauxX <- model.matrix(mcmc$formula[[y]]$fixed, data)
    gauxX <- model.matrix(mcmc$formula[[y]]$group, data)
    rauxX <- model.matrix(mcmc$formula[[y]]$random, data)
    if(is.element(y, mcmc$Ords)){
      fcols[[y]] <- colnames(fauxX)[-1]
      gcols[[y]] <- colnames(gauxX)[-1]
    }else{
      fcols[[y]] <- colnames(fauxX)
      gcols[[y]] <- colnames(gauxX)
    }
    fcols[[y]] <- setdiff(fcols[[y]], gcols[[y]])
    rcols[[y]] <- colnames(rauxX)
    
    notf <- setdiff(fcols[[y]], Xcolnames)
    Xcolnames <- c(Xcolnames, notf)
    addX <- data.frame(fauxX[,notf])
    colnames(addX) <- notf
    finX <- cbind(finX, addX)
    
    notg <- setdiff(gcols[[y]], Xcolnames)
    Xcolnames <- c(Xcolnames, notg)
    addX <- data.frame(gauxX[,notg])
    colnames(addX) <- notg
    finX <- cbind(finX, addX)
    
    notr <- setdiff(rcols[[y]], Xcolnames)
    Xcolnames <- c(Xcolnames, notr)
    addX <- data.frame(rauxX[,notr])
    colnames(addX) <- notr
    finX <- cbind(finX, addX)
    
    ocols[[y]] <- mcmc$formula[[y]]$offset
    noff <- ifelse(mcmc$formula[[y]]$offset == "", 0, 1)
    if(noff > 0){
      noto <- setdiff(ocols[[y]], Xcolnames)
      Xcolnames <- c(Xcolnames, noto)
      addX <- data.frame(data[,noto])
      colnames(addX) <- noto
      finX <- cbind(finX, addX)
    }
    
  }
  
  params <- c("beta_num_fix", "beta_num", "tau_num", 
              "beta_poi_fix", "beta_poi", 
              "beta_bin_fix", "beta_bin", 
              "beta_ord_fix", "beta_ord", "c_ord", 
              "beta_cat_fix", "beta_cat", 
              "InvSigma", 
              "w", "pUig_int", "deviance", "dev_i")
  ydepparams <- c(paste0("beta_", c("num", "poi", "bin", "ord", "cat"), "_fix"),
                  paste0("beta_", c("num", "poi", "bin", "ord", "cat")),
                  "c_ord")
  
  
  settings <- mcmc$settings
  rownames(settings)[rownames(settings)=="loglik"] <- "deviance"
  rownames(settings)[rownames(settings)=="pUig"] <- "pUig_int"
  settings["deviance", "save"] <- T
  settings["pUig_int", "save"] <- T
  settings["pUig_int",  c("d1", "dims", "dimswithG")] = c(ifelse(dynamic_prob,N,n), 
                                                          ifelse(dynamic_prob,N,n)*mcmc$G, 
                                                          ifelse(dynamic_prob,N,n)*mcmc$G)
  addset <- matrix(c(1,0,mcmc$G,0,0,0,0,0,0,0,0,niter,
              n, 0, 1, 0, 0, 1, 1, 
              n, n),
              nrow = 1)
  rownames(addset) <- "dev_i"
  colnames(addset) <- colnames(settings)
  settings <- rbind(settings, addset)
  
  settings <- settings[params,]
  settings$iter <- niter
  
  
  ### Preparation of parameters for Gibbs sampler in C
  ## Data
  # First column is going to be the id variable (0-th column)
  cId <- numbered_unique_ids[as.character(data[,id])] - 1
  # -1 is there for C which works better with number beggining with 0
  # other columns  (beginning with 1st column)
  cY <- numeric()
  for(y in Ys){
    if(is.element(y, mcmc$Cats)){
      # Change the first category to be the last, otherwise preserve the order
      auxy <- ydata[,y]
      auxy[auxy==0] <- mcmc$Kcat[y]+1
      auxy <- auxy - 1
      cY <- c(cY, auxy)
    }else{
      cY <- c(cY, ydata[,y])
    }
    #cY <- c(cY, as.numeric(as.character(data[,y])))
  }
  cisYna <- is.na(cY)
  cY[cisYna] <- 0 # NA is not sent to C function
  
  cX <- numeric()
  # take only needed columns
  Xcolnums <- 1:length(Xcolnames)
  names(Xcolnums) <- Xcolnames
  for(x in Xcolnames){
    cX <- c(cX, finX[,x])
  } # no need for id
  
  # formula
  cFormulaF <- cFormulaG <- cFormulaR <- cFormulaO <- numeric()
  for(y in Ys){
    cFormulaF <- c(cFormulaF, Xcolnums[mcmc$lfixnames[[y]]])
    cFormulaG <- c(cFormulaG, Xcolnums[mcmc$lgrpnames[[y]]])
    cFormulaR <- c(cFormulaR, Xcolnums[mcmc$lrannames[[y]]])
    cFormulaO <- c(cFormulaO, Xcolnums[mcmc$loffnames[[y]]])
  }
  cFormulaF <- cFormulaF - 1 # index in C
  cFormulaG <- cFormulaG - 1 # index in C
  cFormulaR <- cFormulaR - 1 # index in C
  cFormulaO <- cFormulaO - 1 # index in C
  cnfix <- as.numeric(mcmc$nfix) # number of FIXED  regressors with variables y
  cngrp <- as.numeric(mcmc$ngrp) # number of GROUP-SPECIFIC  regressors with variables y
  cnran <- as.numeric(mcmc$nran) # number of RANDOM regressors with variables y
  cnoff <- as.numeric(mcmc$noff) # number of OFFSET regressors with variables y
  
  # head(finX)
  
  # all the following three are going to be saved
  nparams <- sum(settings[c("deviance", "dev_i", "pUig_int"), "dimswithG"])
  if(mcmc$howsave == "matrix"){
    param_names <- list()
    output <- matrix(0, nrow = length(chains)*niter, ncol = 2+nparams)
    output[,1] = rep(chains, each=niter)
    output[,2] = rep(iterations, length(chains))
    cnames <- character(2+nparams)
    cnames[1] = "chain"
    cnames[2] = "m"
    index = 2
    for(p in c("deviance", "dev_i", "pUig_int")){
      aux = from_C_to_matrix(values = rep(0, settings[p,"dimswithG"]),
                             p = p, settings = settings, family = mcmc$family)
      param_names[[p]] <- colnames(aux)
      cnames[(index+1):(index+settings[p,"dimswithG"])] = colnames(aux)
      index = index + settings[p,"dimswithG"]
    }
    colnames(output) <- cnames
    output = as.data.frame(output)
  }
  
  if(mcmc$howsave == "list"){
    output <- list()
  }
  
  for(ch in chains){
    # arrays/fields where those variables are and will be stored
    cpUig_int <- double(niter*settings["pUig_int", "dimswithG"])
    cdeviance <- double(niter*settings["deviance", "dimswithG"])
    cdev_i <- double(niter*settings["dev_i", "dimswithG"])
    #if(dev){cdeviance <- double(niter*cdimswithK["deviance"])}else{cdeviance <- as.double(0)}
    
    if(mcmc$howsave == "list"){
      cbeta_num_fix <- from_list_to_C(mcmc$draws[[ch]], iterations, "beta_num_fix", settings, 
                                      yspecd1 = mcmc$yspecd1$beta_num_fix,
                                      family=mcmc$family)
      cbeta_num <- from_list_to_C(mcmc$draws[[ch]], iterations, "beta_num", settings, 
                                  yspecd1 = mcmc$yspecd1$beta_num,
                                  family=mcmc$family)
      ctau_num  <- from_list_to_C(mcmc$draws[[ch]], iterations, "tau_num", settings,
                                  family=mcmc$family)
      cbeta_poi_fix <- from_list_to_C(mcmc$draws[[ch]], iterations, "beta_poi_fix", settings, 
                                      yspecd1 = mcmc$yspecd1$beta_poi_fix,
                                      family=mcmc$family)
      cbeta_poi <- from_list_to_C(mcmc$draws[[ch]], iterations, "beta_poi", settings,
                                  yspecd1 = mcmc$yspecd1$beta_poi,
                                  family=mcmc$family)
      cbeta_bin_fix <- from_list_to_C(mcmc$draws[[ch]], iterations, "beta_bin_fix", settings, 
                                      yspecd1 = mcmc$yspecd1$beta_bin_fix,
                                      family=mcmc$family)
      cbeta_bin <- from_list_to_C(mcmc$draws[[ch]], iterations, "beta_bin", settings, 
                                  yspecd1 = mcmc$yspecd1$beta_bin,
                                  family=mcmc$family)
      cbeta_ord_fix <- from_list_to_C(mcmc$draws[[ch]], iterations, "beta_ord_fix", settings,
                                      yspecd1 = mcmc$yspecd1$beta_ord_fix,
                                      family=mcmc$family)
      cbeta_ord <- from_list_to_C(mcmc$draws[[ch]], iterations, "beta_ord", settings, 
                                  yspecd1 = mcmc$yspecd1$beta_ord,
                                  family=mcmc$family)
      cc_ord    <- from_list_to_C(mcmc$draws[[ch]], iterations, "c_ord", settings, 
                                  yspecd1 = mcmc$yspecd1$c_ord,
                                  family=mcmc$family)
      cbeta_cat_fix <- from_list_to_C(mcmc$draws[[ch]], iterations, "beta_cat_fix", settings, 
                                      yspecd1 = mcmc$yspecd1$beta_cat_fix, 
                                      yspecd2 = mcmc$yspecd2$beta_cat_fix,
                                      family=mcmc$family)
      cbeta_cat <- from_list_to_C(mcmc$draws[[ch]], iterations, "beta_cat", settings,
                                  yspecd1 = mcmc$yspecd1$beta_cat, 
                                  yspecd2 = mcmc$yspecd2$beta_cat,
                                  family=mcmc$family)
      cInvSigma <- from_list_to_C(mcmc$draws[[ch]], iterations, "InvSigma", settings,
                                  family=mcmc$family)
      cw        <- from_list_to_C(mcmc$draws[[ch]], iterations, "w", settings,
                                  family=mcmc$family)
    }
    if(mcmc$howsave == "matrix"){
      cbeta_num_fix <- from_matrix_to_C(mcmc$draws, chain = ch, iterations, "beta_num_fix", settings)
      cbeta_num <- from_matrix_to_C(mcmc$draws, chain = ch, iterations, "beta_num", settings)
      ctau_num  <- from_matrix_to_C(mcmc$draws, chain = ch, iterations, "tau_num", settings)
      cbeta_poi_fix <- from_matrix_to_C(mcmc$draws, chain = ch, iterations, "beta_poi_fix", settings)
      cbeta_poi <- from_matrix_to_C(mcmc$draws, chain = ch, iterations, "beta_poi", settings)
      cbeta_bin_fix <- from_matrix_to_C(mcmc$draws, chain = ch, iterations, "beta_bin_fix", settings)
      cbeta_bin <- from_matrix_to_C(mcmc$draws, chain = ch, iterations, "beta_bin", settings)
      cbeta_ord_fix <- from_matrix_to_C(mcmc$draws, chain = ch, iterations, "beta_ord_fix", settings)
      cbeta_ord <- from_matrix_to_C(mcmc$draws, chain = ch, iterations, "beta_ord", settings)
      cc_ord    <- from_matrix_to_C(mcmc$draws, chain = ch, iterations, "c_ord", settings)
      cbeta_cat_fix <- from_matrix_to_C(mcmc$draws, chain = ch, iterations, "beta_cat_fix", settings)
      cbeta_cat <- from_matrix_to_C(mcmc$draws, chain = ch, iterations, "beta_cat", settings)
      cInvSigma <- from_matrix_to_C(mcmc$draws, chain = ch, iterations, "InvSigma", settings)
      cw        <- from_matrix_to_C(mcmc$draws, chain = ch, iterations, "w", settings)
    }
    
    #dyn.load("Cfun/pUig_dev.dll")
    #dyn.unload("Cfun/pUig_dev.dll")
    
    print("About to start C computations.")
    # print(ch)
    # print(mcmc$G)
    # print(niter)
    # print(N)
    # print(n)
    # print(mcmc$nY)
    # print(cFormulaF)
    # print(cFormulaG)
    # print(cFormulaR)
    # print(cFormulaO)
    # print(cnfix)
    # print(cngrp)
    # print(cnran)
    # print(cnoff)
    # print(mcmc$totnran)
    # print(mcmc$Kord)
    # print(mcmc$Kcat)
    # print(settings$dims)
    # print(settings$dimswithG)
    # print(mcmc$tuning$integer$kspec_bi_cat)
    # print(cbeta_num_fix)
    # print(cbeta_num)
    # print(ctau_num)
    # print(cbeta_poi_fix)
    # print(cbeta_poi)
    # print(cbeta_bin_fix)
    # print(cbeta_bin)
    # print(cbeta_ord_fix)
    # print(cbeta_ord)
    # print(cc_ord)
    # print(cbeta_cat_fix)
    # print(cbeta_cat)
    # print(cInvSigma)
    # print(cw)
    # print(cpUig_int)
    # print(cdeviance)
    # print(cdev_i)
    # print(dynamic_prob)
    # print(NGQ)
    # print(rules$x)
    # print(rules$w)
    # print(tolerance)
    # print(maxiter)
    # print(maxnrep)
    
    #system.time(
    cmcmc <-
      .C(C_pUig_dev,
         Id        = as.integer(cId),
         Y         = as.double(cY),
         isYna     = as.integer(cisYna),
         X         = as.double(cX),
         varying   = as.integer(mcmc$varying),
         # parameters describing dimensions
         chain     = as.integer(ch), # number of the chain
         G         = as.integer(mcmc$G), # number of classes
         iter      = as.integer(niter), # total number of generated states
         N         = as.integer(N), # total number of observations
         n         = as.integer(n), # total number of subjects (different ids in the dataset)
         nY        = as.integer(mcmc$nY), # 3 numbers: counts of Nums, Ords and Bins variables
         FormulaF  = as.integer(cFormulaF), # numbers of columns of X that should be used for FIXED  effects of modelled responses
         FormulaG  = as.integer(cFormulaG), # numbers of columns of X that should be used for GROUP-SPECIFIC  effects of modelled responses
         FormulaR  = as.integer(cFormulaR), # numbers of columns of X that should be used for RANDOM effects of modelled responses
         FormulaO  = as.integer(cFormulaO), # numbers of columns of X that should be used for OFFSET effects of modelled responses
         nfix      = as.integer(cnfix), 
         ngrp      = as.integer(cngrp), 
         nran      = as.integer(cnran),
         noff      = as.integer(cnoff),
         totnran   = as.integer(mcmc$totnran),
         Kord      = as.integer(mcmc$Kord), # the counts of categories of ordinal variables (-1)
         Kcat      = as.integer(mcmc$Kcat), # the counts of categories of categorical variables (-1)
         dims      = as.integer(settings$dims), # the length of subarray that corresponds to one state (disected by various parameters)
         dimswithG = as.integer(settings$dimswithG), # the length of subarray that corresponds to one state
         kspec_bi_cat = as.integer(mcmc$tuning$integer$kspec_bi_cat),
         #predictor = double(N*sum(nY)),
         # the function should count totnran, totnfix and cumsum versions of nfix and nran
         # arrays to store generated states
         beta_num_fix = as.double(cbeta_num_fix),
         beta_num     = as.double(cbeta_num),
         tau_num      = as.double(ctau_num),
         beta_poi_fix = as.double(cbeta_poi_fix),
         beta_poi     = as.double(cbeta_poi),
         beta_bin_fix = as.double(cbeta_bin_fix),
         beta_bin     = as.double(cbeta_bin),
         beta_ord_fix = as.double(cbeta_ord_fix),
         beta_ord     = as.double(cbeta_ord),
         c_ord        = as.double(cc_ord),
         beta_cat_fix = as.double(cbeta_cat_fix),
         beta_cat     = as.double(cbeta_cat),
         InvSigma     = as.double(cInvSigma),
         w            = as.double(cw),
         pUig_int     = as.double(cpUig_int),
         deviance     = as.double(cdeviance),
         dev_i        = as.double(cdev_i),
         dynamic_prob = as.integer(dynamic_prob),
         NGQ          = as.integer(NGQ),
         roots        = as.double(rules$x),
         weights      = as.double(rules$w),
         tolerance    = as.double(tolerance),
         maxiter      = as.integer(maxiter),
         maxnrep      = as.integer(maxnrep)
      )

    if(mcmc$howsave == "list"){
      # results will be returned in structured list - the same way as original R function
      chain <- list()
      for(p in c("deviance", "dev_i", "pUig_int")){
        if(settings[p, "save"]){
          chain[[p]] <- from_C_to_list(cmcmc[[p]], p, settings, family=mcmc$family)
        }
      }
      chain[["m"]] <- iterations

      output[[ch]] <- chain
    }

    if(mcmc$howsave == "matrix"){
      # results will be returned in matrix - the same way as original R function
      for(p in c("deviance", "dev_i", "pUig_int")){
        aux = from_C_to_matrix(cmcmc[[p]], p, settings, family=mcmc$family)
        output[iterations + (ch-1)*niter, colnames(aux)] <- aux
      }
    }
    
  }  # end of for ch in 1:chains
  
  return(output)
  
}
