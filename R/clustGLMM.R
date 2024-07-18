#mcmc_MBC_GLMMs <- function(...) UseMethod("mbcglmms")
#mbcglmms.default <-
clustGLMM <- 
function(
  formula,   # list of formulas for each of the variables
  id,        # column name of IDs declaring which observations come from the same unit
  family,    # named vector of family types 
  data,      # a data.frame object containing outcomes, regressors and id
  G,         # total number of latent Groups=clusters
  varying,   # which parameters should be class-specific
                # order: tau_num, a_ord, InvSigma, InvQ, naY
  save,      # what parameters should be saved, in order:
              # "beta_num", "beta_num_grp",
              # "tau_num", "sd_num", "var_num",
              # "beta_poi", "beta_poi_grp",
              # "beta_bin", "beta_bin_grp",
              # "beta_ord", "beta_ord_grp",
              # "a_ord", "c_ord", "pi_ord",
              # "beta_cat", "beta_cat_grp",
              # "InvSigma", "Sigma", "sdSigma", "corSigma", "detInvSigma",
              # "InvQ", "Q", "detInvQ",
              # "b",
              # "w", "ng",
              # "loglik", "pUig", "U",
              # "Gplus", "e0"
  howsave,   # how should be the output MCMC chains returned?
                # matrix - in a matrix with special column names
                # list - in lists of lists of parameters...
  param,     # list of hyperparameters for prior distributions
  inits,     # list of initial values from which to start
                # if missing or NULL --> create your own
  tuning,    # list of tuning parameters
  iter = 1000,# length of sampled chain
  nchains = 1 # number of chains
){
  # TODO make these function parameters
  clusteringU_modeGplus = TRUE
  
  # require("MASS")
  # require("nnet")
  
  ### Missing values from input
  if (missing(varying)){
    varying <- default_varying()
  }
  # correction of possible mistakes 
  # - if deeper variables is class-specific then the dependent one must be as well
  # prior for InvSigma depends on Invq
  if(varying["InvQ"]){varying["InvSigma"] <- T}
  
  if(missing(save)){
    save <- default_save()
  }
  
  important <- c("beta_num_fix", "beta_num",
                 "tau_num",
                 "beta_poi_fix", "beta_poi",
                 "beta_bin_fix", "beta_bin",
                 "beta_ord_fix", "beta_ord",
                 "c_ord",
                 "beta_cat_fix", "beta_cat",
                 "InvSigma", 
                 "w")
  
  if(!all(save[important])){
    warning("Some of the parameters: tau_num, betas, c_ord, InvSigma will not be saved.
Therefore, you will not be able to calculate unconditional classification probabilities of new subjects.
Change the save settings if you want to calculate these probabilities afterwards.")
  }
  
  if(missing(howsave)){
    howsave = "matrix"
  }
  
  if(missing(param)){
    param <- default_param()
  } # end of if missing(param)
  
  
  if(missing(tuning)){
    tuning <- default_tuning()
  }
  
  ### Most important auxiliary variables
  mcmc <- list()
  N <- dim(data)[1]
  unique_ids <- unique(data[,id])
  n <- length(unique_ids)
  numbered_unique_ids <- c(1:n)
  names(numbered_unique_ids) <- unique_ids
  nsubj <- table(data[,id])
  
  ### Sorting family types
  Nums <- names(family)[is.element(family, c("num", "gaussian"))]
  Pois <- names(family)[is.element(family, c("poi", "poisson"))]
  Bins <- names(family)[is.element(family, c("bin", "bernoulli"))]
  Ords <- names(family)[is.element(family, c("ord", "cumulative"))]
  Cats <- names(family)[is.element(family, c("cat", "categorical"))]
  Ys <- c(Nums, Pois, Bins, Ords, Cats)
  nY <- sapply(list(Nums, Pois, Bins, Ords, Cats), length)
  names(nY) <- c("Nums", "Pois", "Bins", "Ords", "Cats")
  
  # sorted_family <- c(family[Nums], family[Pois], family[Bins], family[Ords], family[Cats])
  sorted_family <- c(rep("num", length(Nums)),
                     rep("poi", length(Pois)),
                     rep("bin", length(Bins)),
                     rep("ord", length(Ords)),
                     rep("cat", length(Cats)))
  names(sorted_family) <- Ys
  if(length(sorted_family) < length(family)){
    not_matched_outcomes <- setdiff(names(family), names(sorted_family))
    stop(paste0("Families of outcomes: ", paste(not_matched_outcomes, collapse = ", "), 
                "... do not match any allowed family type. 
                See help page for the list of accepted family types."))
  }
  if(!identical(names(family), names(sorted_family))){
    warning("family argument has been reordered by type: num, poi, bin, ord, cat. 
            This is also the ordering used for the random effects.")
  }
  
  ### Checking family types + transformations
  ydata <- data
  # Numeric
  if(nY["Nums"] > 0){
    for(y in Nums){
      yval <- data[,y]
      if(!is.numeric(yval)){
        stop(paste0("Outcome ", y, " is not numeric!"))
      }
      if(sum(infs <- is.infinite(yval))>0){
        stop(paste0("There are infinite values for outcome ", y, "!"))
      }
      ydata[,y] <- yval # numeric values
    }
  }
  
  # Poisson
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5){abs(x - as.integer(x)) < tol}
  if(nY["Pois"] > 0){
    for(y in Pois){
      yval <- data[,y]
      if(!all(is.wholenumber(yval))){
        stop(paste0("Outcome ", y, " is not integer-valued!"))
      }
      if(!all(yval >= 0)){
        stop(paste0("Outcome ", y, " has some negative values!"))
      }
      if(sum(infs <- is.infinite(yval))>0){
        stop(paste0("There are infinite values for outcome ", y, "!"))
      }
      ydata[,y] <- yval # integer values
    }
  }
  
  # Binary
  if(nY["Bins"] > 0){
    for(y in Bins){
      fy <- factor(data[,y])
      if(nlevels(fy) != 2){
        stop(paste0("Binary outcome ", y, " does not have 2 levels."))
      }
      ydata[,y] <- as.numeric(fy)-1 # 0, 1 values
    }
  }
  
  # Ordinal
  if(nY["Ords"] > 0){
    Kord <- numeric(length(Ords))
    names(Kord) <- Ords
    for(y in Ords){
      fy <- data[,y]
      if(!is.ordered(fy)){
        fy <- factor(fy, ordered = TRUE)
        warning(paste0("Converting ordinal outcome ", y, " to ordered factor with levels: ",
                       paste(levels(fy), collapse = " < ")))
      }
      if(nlevels(fy) <= 2){
        stop(paste0("Ordinal outcome ", y, " has less than 3 levels."))
      }
      ydata[,y] <- as.numeric(fy)-1 # 0 < 1 < ... < Kord[y]
      Kord[y] <- nlevels(fy)-1
    }
  }else{
    Kord <- NULL
  }
  
  # Categorical
  if(nY["Cats"] > 0){
    Kcat <- numeric(length(Cats))
    names(Kcat) <- Cats
    for(y in Cats){
      fy <- data[,y]
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
      Kcat[y] <- nlevels(fy)-1
    }
  }else{
    Kcat <- NULL
  }
  
  if(length(unique(Ys)) != length(Ys)){
    stop("Some outcome names are being repeated.")
  }
  
  
 
  nfix <- ngrp <- nran <- noff <- numeric(sum(nY))
  names(nfix) <- names(ngrp) <- names(nran) <- names(noff) <- Ys
  fixnames <- grpnames <- rannames <- offnames <- numeric()
  lfixnames <- lgrpnames <- lrannames <- loffnames <- list()
  # Creating model matrix containing all needed columns for C
  finX <- data.frame(id = data[,id])
  Xcolnames <- numeric()
  fcols <- gcols <- rcols <- ocols <- list()
  for(y in Ys){
    fauxX <- model.matrix(formula[[y]]$fixed, data)
    gauxX <- model.matrix(formula[[y]]$group, data)
    rauxX <- model.matrix(formula[[y]]$random, data)
    if(is.element(y, Ords)){
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
    
    ocols[[y]] <- formula[[y]]$offset
    loffnames[[y]] <- as.character(c())
    noff[y] <- ifelse(formula[[y]]$offset == "", 0, 1)
    if(noff[y] > 0){
      loffnames[[y]] <- ocols[[y]]
      noto <- setdiff(ocols[[y]], Xcolnames)
      Xcolnames <- c(Xcolnames, noto)
      addX <- data.frame(data[,noto])
      colnames(addX) <- noto
      finX <- cbind(finX, addX)
    }
    
    nfix[y] <- length(fcols[[y]])
    ngrp[y] <- length(gcols[[y]])
    nran[y] <- length(rcols[[y]])
    lfixnames[[y]] <- fcols[[y]]
    lgrpnames[[y]] <- gcols[[y]]
    lrannames[[y]] <- rcols[[y]]
    
    fixnames <- c(fixnames, lfixnames[[y]])
    grpnames <- c(grpnames, lgrpnames[[y]])
    rannames <- c(rannames, lrannames[[y]])
    offnames <- c(offnames, loffnames[[y]])
  }
  totnran <- sum(nran)+ifelse(tuning$integer$kspec_bi_cat, sum(nran[Cats]*Kcat)-sum(nran[Cats]), 0)
  totnfix <- sum(nfix)-sum(nfix[Cats])+sum(nfix[Cats]*Kcat)
  totngrp <- G * (sum(ngrp)-sum(ngrp[Cats])+sum(ngrp[Cats]*Kcat))
  
  # if(N != dim(X)[1]){stop("Y and X have different number of rows.")}
  
  
  # update param depending on the total dimension of random effects
  param$nu_0 <- param$nu_0 + totnran
  param$nu_1 <- param$nu_1 + totnran
  param$InvV <- diag(param$InvV, totnran)
  
  
  ### Initial values computation
  initsgiven = 1
  if(missing(inits)){
    message("Inits are missing, calculating...\n")
    initsgiven = 0
    inits <- list()
    for(ch in 1:nchains){
      inits[[ch]] <- list()
      
      # parameters connected to classification
      inits[[ch]]$w <- rep(1/G, G)
      inits[[ch]]$U <- sample(1:G, n, replace = T, prob = inits[[ch]]$w)
      inits[[ch]]$e0 <- param$ae / param$be
      #inits[[ch]]$U <- as.numeric(as.character(data$g[data$j==1]))
      Ng <- selection <- list()
      for(g in 1:G){
        Ng[[g]] <- unique_ids[inits[[ch]]$U == g]
        selection[[g]] <- is.element(data[,id], Ng[[g]])
      }
      inits[[ch]]$pUig <- matrix(rep(inits[[ch]]$w, n), ncol = G, nrow = n, byrow = T)
      
      # parameters connected to numerical variables --> use linear regression
      inits[[ch]]$beta_num <- inits[[ch]]$beta_num_fix <- inits[[ch]]$tau_num <- list()
      for(y in Nums){
        Xvars <- c(fcols[[y]], gcols[[y]])
        # Using all data
        nonna <- !is.na(ydata[,y])
        pomY <- ydata[nonna,y] - apply(finX[nonna, loffnames[[y]]], 1, sum)
        boldX <- as.matrix(finX[nonna,Xvars])
        coefs <- as.numeric(solve(t(boldX) %*% boldX, t(boldX) %*% pomY))
        names(coefs) <- Xvars
        MSe <- sum((pomY - boldX %*% matrix(coefs, ncol = 1))^2)/(length(pomY) - length(Xvars))
        init_tau_num <- 1/MSe
        inits[[ch]]$beta_num_fix[[y]] <- coefs[fcols[[y]]]
        
        if(varying["tau_num"]){
          inits[[ch]]$tau_num[[y]] <- list()
        }else{
          inits[[ch]]$tau_num[[y]] <- init_tau_num
        }
        
        inits[[ch]]$beta_num[[y]] <- list()
        for(g in 1:G){
          inds <- selection[[g]] & !is.na(ydata[,y])
          pomY <- ydata[inds,y] - apply(finX[inds, loffnames[[y]]], 1, sum)
          boldX <- as.matrix(finX[inds, Xvars])
          XtX <- t(boldX) %*% boldX
          coefs <- try(as.numeric(solve(t(boldX) %*% boldX, t(boldX) %*% pomY)), silent = T)
          if(is(coefs, "try-error")){
            nonna <- !is.na(ydata[,y])
            pomY <- ydata[nonna,y] - apply(finX[nonna, loffnames[[y]]], 1, sum)
            boldX <- as.matrix(finX[nonna,Xvars])
            coefs <- as.numeric(solve(t(boldX) %*% boldX, t(boldX) %*% pomY))
            names(coefs) <- Xvars
            MSe <- sum((pomY - boldX %*% matrix(coefs, ncol = 1))^2)/(length(pomY) - length(Xvars))
          }else{
            names(coefs) <- Xvars
            MSe <- sum((pomY - boldX %*% matrix(coefs, ncol = 1))^2)/(length(pomY) - length(Xvars))
          }
          if(varying["tau_num"]){
            inits[[ch]]$tau_num[[y]][[g]] <- 1/MSe
          }
          inits[[ch]]$beta_num[[y]][[g]] <- coefs[gcols[[y]]]
        } # end of for g in 1:G
      } # end of for y in Nums
      
      # parameters connected to binary variables --> use logistic regression
      inits[[ch]]$beta_poi <- inits[[ch]]$beta_poi_fix <- list()
      for(y in Pois){
        aux_form <- paste0(y, " ~ ",
                           formula[[y]]$fixed[2],
                           " + ", 
                           formula[[y]]$group[2]
                           )
        if(noff[y] > 0){
          aux_form <- paste0(aux_form, "+ offset(",
                             formula[[y]]$offset,
                             ")")
        }
        auxglm0 <- glm(aux_form, data = ydata, family = poisson(link = "log"))
        coefs <- auxglm0$coefficients
        inits[[ch]]$beta_poi_fix[[y]] <- coefs[fcols[[y]]]
        
        inits[[ch]]$beta_poi[[y]] <- list()
        for(g in 1:G){
          auxdata <- as.data.frame(ydata[selection[[g]], ])
          auxglm <- tryCatch(glm(aux_form, data = auxdata,
                                 family = poisson(link = "log")),
                             error=function(e) e, warning=function(w) w)
          if(is(auxglm,"warning") || is(auxglm,"error")){
            auxglm <- auxglm0
          }else{
            if(!auxglm$converged){
              auxglm <- auxglm0
            }
          }
          coefs <- auxglm$coefficients
          if(sum(is.na(coefs))>0){
            coefs <- auxglm0$coefficients
          }
          inits[[ch]]$beta_poi[[y]][[g]] <- coefs[gcols[[y]]]
        }
      } # end for y in Pois
      
      # parameters connected to binary variables --> use logistic regression
      inits[[ch]]$beta_bin <- inits[[ch]]$beta_bin_fix <- list()
      for(y in Bins){
        aux_form <- paste0(y, " ~ ",
                           formula[[y]]$fixed[2],
                           " + ", 
                           formula[[y]]$group[2]
        )
        if(noff[y] > 0){
          aux_form <- paste0(aux_form, "+ offset(",
                             formula[[y]]$offset,
                             ")")
        }
        auxglm0 <- glm(aux_form, data = ydata, family = binomial(link = "logit"))
        coefs <- auxglm0$coefficients
        inits[[ch]]$beta_bin_fix[[y]] <- coefs[fcols[[y]]]
        
        inits[[ch]]$beta_bin[[y]] <- list()
        for(g in 1:G){
          auxdata <- as.data.frame(ydata[selection[[g]],])
          auxglm <- tryCatch(glm(aux_form, data = auxdata,
                                 family = binomial(link = "logit")),
                             error=function(e) e, warning=function(w) w)
          if(is(auxglm,"warning") || is(auxglm,"error")){
            auxglm <- auxglm0
          }else{
            if((!auxglm$converged) || (sum(is.na(auxglm$coefficients[gcols[[y]]]))>0)){
              auxglm <- auxglm0
            }
          }
          coefs <- auxglm$coefficients
          if(sum(is.na(coefs))>0){
            coefs <- auxglm0$coefficients
          }
          inits[[ch]]$beta_bin[[y]][[g]] <- coefs[gcols[[y]]]
        }
      } # end for y in Bins
      
      # parameters connected to ordinal variables --> use ordered logistic regression
      # library("MASS")
      inits[[ch]]$beta_ord <- inits[[ch]]$beta_ord_fix <- inits[[ch]]$a_ord <- inits[[ch]]$c_ord <- inits[[ch]]$pi_ord <-list()
      for(y in Ords){
        aux_form <- paste0(y, " ~ ",
                           formula[[y]]$fixed[2],
                           " + ", 
                           formula[[y]]$group[2]
        )
        if(noff[y] > 0){
          aux_form <- paste0(aux_form, "+ offset(",
                             formula[[y]]$offset,
                             ")")
        }
        auxdata0 <- ydata
        auxdata0[,y] <- factor(auxdata0[,y])
        auxpolr0 <- polr(aux_form, data = auxdata0, Hess = TRUE)
        inits[[ch]]$beta_ord_fix[[y]] <- auxpolr0$coefficients[fcols[[y]]]
        
        if(varying["c_ord"]){
          inits[[ch]]$c_ord[[y]] <- inits[[ch]]$a_ord[[y]] <- inits[[ch]]$pi_ord[[y]] <- list()
        }else{
          inits[[ch]]$c_ord[[y]] <- auxpolr0$zeta
          logits <- c(0, exp(auxpolr0$zeta)/(1+exp(auxpolr0$zeta)))
          inits[[ch]]$a_ord[[y]] <- log((logits[2:length(logits)] - logits[1:(length(logits)-1)])/(1-logits[length(logits)]))
          inits[[ch]]$pi_ord[[y]] <- exp(c(inits[[ch]]$a_ord[[y]], 0))/(sum(exp(c(inits[[ch]]$a_ord[[y]], 0))))
        }
        
        inits[[ch]]$beta_ord[[y]] <- list()
        for(g in 1:G){
          if(length(unique(ydata[selection[[g]],y])) == Kord[y]+1){
            auxdata <- as.data.frame(auxdata0[selection[[g]],])
          }else{
            auxdata <- auxdata0
          }
          auxdata[,y] <- factor(auxdata[,y])
          auxpolr <- tryCatch(polr(aux_form, data = auxdata, Hess = TRUE), 
                              error=function(e) e, warning=function(w) w)
          dowithall <- F
          if(is(auxpolr,"warning") | is(auxpolr,"error")){
            dowithall <- T
          }else{
            logits <- exp(auxpolr$zeta)/(1+exp(auxpolr$zeta))
            sdlogits <- sd(logits)
            if(is.na(sdlogits)){
              dowithall <- T
            }else{
              if(auxpolr$convergence || sdlogits < 1e-14 || (sum(is.na(auxpolr$coefficients[gcols[[y]]]))>0)){
                dowithall <- T
              }
            }
          }
          if(dowithall){
            auxpolr <- auxpolr0
          }
        
          inits[[ch]]$beta_ord[[y]][[g]] <- auxpolr$coefficients[gcols[[y]]]
          
          if(varying["c_ord"]){
            inits[[ch]]$c_ord[[y]][[g]] <- auxpolr$zeta
            logits <- c(0, exp(auxpolr$zeta)/(1+exp(auxpolr$zeta)))
            inits[[ch]]$a_ord[[y]][[g]] <- log((logits[2:length(logits)] - logits[1:(length(logits)-1)])/(1-logits[length(logits)]))
            inits[[ch]]$pi_ord[[y]][[g]] <- exp(c(inits[[ch]]$a_ord[[y]][[g]], 0))/(sum(exp(c(inits[[ch]]$a_ord[[y]][[g]], 0))))
          }
        } # end for g in 1:G
      } # end for y in Ords
      
      # parameters connected to ordinal variables --> use multinomial logistic regression
      # library("nnet")
      inits[[ch]]$beta_cat <- inits[[ch]]$beta_cat_fix <- list()
      for(y in Cats){
        aux_form <- paste0(y, " ~ ",
                           formula[[y]]$fixed[2],
                           " + ", 
                           formula[[y]]$group[2]
        )
        if(noff[y] > 0){
          aux_form <- paste0(aux_form, "+ offset(",
                             formula[[y]]$offset,
                             ")")
        }
        auxdata0 <- ydata
        auxdata0[,y] <- factor(auxdata0[,y]) # the first category is the baseline
        auxmult0 <- multinom(as.formula(aux_form), data = auxdata0)
        sumauxmult0 <- summary(auxmult0)
        inits[[ch]]$beta_cat_fix[[y]] <- t(sumauxmult0$coefficients[,fcols[[y]]])
        
        inits[[ch]]$beta_cat[[y]] <- list()
        for(g in 1:G){
          if(length(unique(ydata[selection[[g]],y])) == Kcat[y]+1){
            auxdata <- as.data.frame(auxdata0[selection[[g]],])
          }else{
            auxdata <- auxdata0
          }
          auxdata[,y] <- factor(auxdata[,y]) # the first category is the baseline
          auxmult <- tryCatch(multinom(as.formula(aux_form), data = auxdata),
                              error=function(e) e, warning=function(w) w)
          if(is(auxmult,"warning") || is(auxmult,"error")){
            auxmult <- auxmult0
          }else{
            if(auxmult$convergence || (sum(is.na(auxmult$coefficients[gcols[[y]]]))>0)){
              auxmult <- auxmult0
            }
          }
          sumauxmult <- summary(auxmult)
          inits[[ch]]$beta_cat[[y]][[g]] <- t(sumauxmult$coefficients[,gcols[[y]]])
        }
        
      } # end for y in Cats
    
      # InvSigma, InvQ parameters
      # E InvSigma = E W(Q,nu_0) = nu_0 * Q
      # InvQ approx Inv(InvSigma/nu_0)
      if(varying["InvQ"]){
        # so must be InvSigma
        inits[[ch]]$InvQ <- inits[[ch]]$InvSigma <- list()
        for(g in 1:G){
          inits[[ch]]$InvSigma[[g]] <- diag(totnran) / rgamma(totnran, 1/param$init_sd_b^2)
          inits[[ch]]$InvQ[[g]] <- chol2inv(chol(inits[[ch]]$InvSigma[[g]]/param$nu_0))  
        } # end of for(k in 1:K)
      }else{
        if(varying["InvSigma"]){
          inits[[ch]]$InvSigma <- list()
          meanInvSigma <- matrix(0, ncol = totnran, nrow = totnran)
          for(g in 1:G){
            inits[[ch]]$InvSigma[[g]] <- diag(totnran) / rgamma(totnran, 1/param$init_sd_b^2)
            meanInvSigma <- meanInvSigma + inits[[ch]]$InvSigma[[g]]/G
          } 
          inits[[ch]]$InvQ <- chol2inv(chol(meanInvSigma/param$nu_0))
        }else{
          inits[[ch]]$InvSigma <- diag(totnran) / rgamma(totnran, 1/param$init_sd_b^2)
          inits[[ch]]$InvQ <- chol2inv(chol(inits[[ch]]$InvSigma/param$nu_0))
        } # end of if(varying["InvSigma"])
      } # end of else if(varying["InvQ"])
      
      # random effects b
      inits[[ch]]$b <- matrix(rnorm(n*totnran, mean = 0, sd = param$init_sd_b), 
                              nrow = n, ncol = totnran)
      #colnames(inits[[ch]]$b) <- rannames
      
      
      if(ch == nchains){message("Inits are ready.\n")}
    } # end for ch in 1:nchains
    
  } # end of missing(inits)
  
  
  ### Preparations for transfer from C to List output
  # creating settings matrix for all (calculable) parameters
  params <- names(save)
  ydepparams <- c(paste0("beta_", c("num", "poi", "bin", "ord", "cat"), "_fix"),
                  paste0("beta_", c("num", "poi", "bin", "ord", "cat")),
                  paste0(c("c", "a", "pi"), "_ord"))
  
  settings <- data.frame(save = save,
                         isspec = sapply(params, function(p){ifelse(is.na(varying[p]), F, varying[p])}),
                         G = G,
                         isy = F,
                         ydepd1 = F,
                         ydepd2 = F,
                         ynums = 0,
                         ypois = 0,
                         ybins = 0,
                         yords = 0,
                         ycats = 0,
                         iter = rep(iter, length(params)),
                         d1 = 0,
                         d2 = 0, 
                         BYROW = T,
                         sym = 0,
                         diag = 0,
                         diagval = 1,
                         D = 0)
  
  rownames(settings) <- params
  settings[ydepparams, "ydepd1"] <- 1
  settings[paste0("beta_", c("num", "poi", "bin", "ord", "cat")), "isspec"] <- T
  #settings
  
  # Individual changes
  settings["w", c("d1", "D")] = c(G, 1)
  settings["U", c("d1", "D")] = c(n, 1)
  settings["pUig", c("d1", "d2", "BYROW", "D")] = c(n, G, T, 2)
  
  # ydepparams have to be done separately
  # due to d1 varying on y
  settings["tau_num", c("isy", "ynums")] <- c(T, 1)
  settings["sd_num", c("isy", "ynums", "isspec")] <- c(T, 1, varying["tau_num"])
  settings["var_num", c("isy", "ynums", "isspec")] <- c(T, 1, varying["tau_num"])
  settings["beta_num_fix", c("isy", "ynums", "D")] <- c(T,1,1)
  settings["beta_num", c("isy", "ynums", "D")] <- c(T,1,1)
  
  settings["beta_poi_fix", c("isy", "ypois", "D")] <- c(T,1,1)
  settings["beta_poi", c("isy", "ypois", "D")] <- c(T,1,1)
  
  settings["beta_bin_fix", c("isy", "ybins", "D")] <- c(T,1,1)
  settings["beta_bin", c("isy", "ybins", "D")] <- c(T,1,1)
  
  settings["beta_ord_fix", c("isy", "yords", "D")] <- c(T,1,1)
  settings["beta_ord", c("isy", "yords", "D")] <- c(T,1,1)
  settings["c_ord", c("isy", "yords", "D")] <- c(T,1,1)
  settings["a_ord", c("isspec", "isy", "yords", "D")] <- c(varying["c_ord"],T,1,1)
  settings["pi_ord", c("isspec", "isy", "yords", "D")] <- c(varying["c_ord"],T,1,1)
  
  settings["beta_cat_fix", c("isy", "ydepd2", "ycats", "D", "BYROW")] <- c(T,T,1,2,F)
  settings["beta_cat", c("isy", "ydepd2", "ycats", "D", "BYROW")] <- c(T,T,1,2,F)
  
  yspecd1 <- list()
  yspecd1$a_ord <- yspecd1$c_ord <- Kord
  yspecd1$pi_ord <- Kord+1
  yspecd1$beta_num_fix <- nfix[Nums]
  yspecd1$beta_poi_fix <- nfix[Pois]
  yspecd1$beta_bin_fix <- nfix[Bins]
  yspecd1$beta_ord_fix <- nfix[Ords]
  yspecd1$beta_cat_fix <- nfix[Cats]
  yspecd1$beta_num <- ngrp[Nums]
  yspecd1$beta_poi <- ngrp[Pois]
  yspecd1$beta_bin <- ngrp[Bins]
  yspecd1$beta_ord <- ngrp[Ords]
  yspecd1$beta_cat <- ngrp[Cats]
  yspecd2 <- list()
  yspecd2$beta_cat_fix <- yspecd2$beta_cat <- Kcat
  
  # Sigma parameters
  settings["InvSigma", c("d1","d2", "sym", "diag", "D")] <- c(totnran,totnran,1,1,2)
  settings["Sigma", c("isspec","d1","d2", "sym", "diag", "D")] <- c(varying["InvSigma"],totnran,totnran,1,1,2)
  settings["sdSigma", c("isspec", "d1", "D")] <- c(varying["InvSigma"], totnran, 1)
  settings["corSigma", c("isspec","d1","d2", "sym", "diag", "diagval", "D")] <- c(varying["InvSigma"],totnran,totnran,1,0,1,2)
  settings["detInvSigma", c("isspec")] <- c(varying["InvSigma"])
  
  # Q parameters
  settings["InvQ", c("d1","d2", "sym", "diag", "D")] <- c(totnran,totnran,1,1,2)
  settings["Q", c("isspec","d1","d2", "sym", "diag", "D")] <- c(varying["InvQ"],totnran,totnran,1,1,2)
  settings["detInvQ", c("isspec")] <- c(varying["InvQ"])
  
  settings["b", c("d1", "d2", "BYROW", "D")] <- c(n, totnran, T, 2)
  
  # sparse clusters parameters
  settings["ng", c("d1", "D")] <- c(G, 1)
  
  ### Preparation of parameters for Gibbs sampler in C
  ## Data
  # First column is going to be the id variable (0-th column)
  cId <- numbered_unique_ids[as.character(data[,id])] - 1
  # -1 is there for C which works better with number beggining with 0
  # other columns  (beginning with 1st column)
  cY <- numeric()
  for(y in Ys){
    if(is.element(y, Cats)){
      # Change the first category to be the last, otherwise preserve the order
      auxy <- ydata[,y]
      auxy[auxy==0] <- Kcat[y]+1
      auxy <- auxy - 1
      cY <- c(cY, auxy)
    }else{
      cY <- c(cY, ydata[,y])
    }
    #cY <- c(cY, as.numeric(as.character(data[,y])))
  }
  cisYna <- is.na(cY)
  cY[cisYna] <- 0 # NA is not sent to C function
  
  settings["naY", c("d1", "D")] <- 
    c(sum(cisYna), 1)
  
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
    cFormulaF <- c(cFormulaF, Xcolnums[lfixnames[[y]]])
    cFormulaG <- c(cFormulaG, Xcolnums[lgrpnames[[y]]])
    cFormulaR <- c(cFormulaR, Xcolnums[lrannames[[y]]])
    cFormulaO <- c(cFormulaO, Xcolnums[loffnames[[y]]])
  }
  cFormulaF <- cFormulaF - 1 # index in C
  cFormulaG <- cFormulaG - 1 # index in C
  cFormulaR <- cFormulaR - 1 # index in C
  cFormulaO <- cFormulaO - 1 # index in C
  cnfix <- as.numeric(nfix) # number of FIXED  regressors with variables y
  cngrp <- as.numeric(ngrp) # number of GROUP-SPECIFIC FIXED  regressors with variables y
  cnran <- as.numeric(nran) # number of RANDOM regressors with variables y
  cnoff <- as.numeric(noff) # number of OFFSETS regressors with variables y
  # dims - in the following order:
  cdims <- 
    c(sum(nfix[Nums]), sum(ngrp[Nums]), #"beta_num_fix", "beta_num",
      nY["Nums"], nY["Nums"], nY["Nums"], #"tau_num", "sd_num", "var_num",
      sum(nfix[Pois]), sum(ngrp[Pois]), #"beta_poi_fix", "beta_poi",
      sum(nfix[Bins]), sum(ngrp[Bins]), #"beta_bin_fix", "beta_bin",
      sum(nfix[Ords]), sum(ngrp[Ords]), #"beta_ord_fix", "beta_ord",
      sum(Kord), sum(Kord), sum(Kord+1), #"a_ord", "c_ord", "pi_ord",
      sum(nfix[Cats]*Kcat), sum(ngrp[Cats]*Kcat), #"beta_cat_fix", "beta_cat",
      totnran*(totnran+1)/2, totnran*(totnran+1)/2, totnran, totnran*(totnran-1)/2, 1,
      #"InvSigma", "Sigma", "sdSigma", "corSigma", "detInvSigma",
      totnran*(totnran+1)/2, totnran*(totnran+1)/2, 1, #"InvQ", "Q", "detInvQ",
      n*totnran, #"b",
      G, G, #"w", "ng",
      1, n*G, n, #"loglik", "pUig", "U",
      1, 1,  #"Gplus", "e0"
      sum(cisYna)
    )
  cdimswithG <- cdims * ((settings$isspec==0)*1 + (settings$isspec==1)*G)
  
  names(cdims) <- names(cdimswithG) <- names(save)
  settings$dims <- cdims
  settings$dimswithG <- cdimswithG
  
  # cparam - created in specific order (user is not required to have the same order)
  cparam <- list(nu_0 = param$nu_0, # for I-W prior of Sigma
                 nu_1 = param$nu_1, # for I-W prior of Q
                 gamma_a=param$gamma_a, 
                 gamma_b=param$gamma_b, # gamma prior for taus
                 # standard deviations for beta parameters depending on type and group-specificity
                 sd_beta_num_fix = param$sd_beta_num_fix, 
                 sd_beta_num = param$sd_beta_num,
                 sd_beta_poi_fix = param$sd_beta_poi_fix,
                 sd_beta_poi = param$sd_beta_poi,
                 sd_beta_bin_fix = param$sd_beta_bin_fix,
                 sd_beta_bin = param$sd_beta_bin,
                 sd_beta_ord_fix = param$sd_beta_ord_fix,
                 sd_beta_ord = param$sd_beta_ord,
                 api_prior = param$api_prior,
                 sd_beta_cat_fix = param$sd_beta_cat_fix,
                 sd_beta_cat = param$sd_beta_cat,
                 # gamma prior hyperparameters for e_0
                 ae = param$ae, 
                 be = param$be)
  cparam$InvV <- param$InvV[upper.tri(param$InvV, diag = TRUE)]
  # order is nu_0, nu_1, delta, gamma_a, gamma_b, ae, be, (all [1])
  # InvV [totnran*(totnran+1)/2]
  # and sd_beta [1]
  cparam <- unlist(cparam) # will be delivered as vector
  
  # csave - what parameters should be returned by C function
  csave <- save
  csave["U"] <- TRUE # in order to perform U-based clustering
  
  # ctuning - tuning in a given order
  ctuning <- list()
  ctuning$integer <- list(
    freq_proposal_update = tuning$integer$freq_proposal_update,
    times_proposal = tuning$integer$times_proposal,
    maxiter = tuning$integer$maxiter,
    maxnrep = tuning$integer$maxnrep,
    kspec_bi_cat = tuning$integer$kspec_bi_cat
  )
  ctuning$double <- list(
    const_proposal_beta_poi_fix = tuning$double$const_proposal_beta_poi_fix,
    const_proposal_beta_poi = tuning$double$const_proposal_beta_poi,
    const_proposal_beta_bin_fix = tuning$double$const_proposal_beta_bin_fix,
    const_proposal_beta_bin = tuning$double$const_proposal_beta_bin,
    const_proposal_beta_ord_fix = tuning$double$const_proposal_beta_ord_fix,
    const_proposal_beta_ord = tuning$double$const_proposal_beta_ord,
    const_proposal_beta_cat_fix = tuning$double$const_proposal_beta_cat_fix,
    const_proposal_beta_cat = tuning$double$const_proposal_beta_cat,
    const_proposal_a_ord = tuning$double$const_proposal_a_ord,
    const_proposal_b = tuning$double$const_proposal_b,
    const_proposal_e0 = tuning$double$const_proposal_e0,
    tolerance = tuning$double$tolerance
  )
  
  # arrays/fields where those variables will (might) be stored
  cstore <- list()
  for(p in params){
    if(p=="U" | p=="Gplus" | p=="ng"){
      # integer variables
      if(settings[p,"save"]){
        cstore[[p]] <- integer(settings[p,"iter"] * settings[p,"dimswithG"])
      }else{
        cstore[[p]] <- as.integer(0)  
      }
    }else{
      # double variables
      if(settings[p,"save"]){
        cstore[[p]] <- double(settings[p,"iter"] * settings[p,"dimswithG"])
      }else{
        cstore[[p]] <- as.double(0)  
      }
    }
  }
  
  # list of last generated states
  last <- list()
  Usamples <- list()
  
  message("Settings ready, about to start for ch in 1:nchains.\n")
  
  nparams <- sum(settings[settings$save, "dimswithG"])
  mcmc$param_names <- list()
  cnames <- character(1+nparams)
  cnames[1] <- "chain"
  index <- 1
  for(p in params){
    if(settings[p,"save"] & settings[p,"dims"] > 0){
      if(settings[p,"ydepd1"] & settings[p,"ydepd2"]){
        aux = from_C_to_matrix(values = rep(0, settings[p,"dimswithG"]),
                               p = p,
                               settings = settings,
                               yspecd1 = yspecd1[[p]],
                               yspecd2 = yspecd2[[p]],
                               family = sorted_family)
      }
      if(settings[p,"ydepd1"] & !settings[p,"ydepd2"]){
        aux = from_C_to_matrix(values = rep(0, settings[p,"dimswithG"]),
                               p = p,
                               settings = settings,
                               yspecd1 = yspecd1[[p]],
                               family = sorted_family)
      }
      if(!settings[p,"ydepd1"] & settings[p,"ydepd2"]){
        aux = from_C_to_matrix(values = rep(0, settings[p,"dimswithG"]),
                               p = p,
                               settings = settings,
                               yspecd2 = yspecd2[[p]],
                               family = sorted_family)
      }
      if(!settings[p,"ydepd1"] & !settings[p,"ydepd2"]){
        aux = from_C_to_matrix(values = rep(0, settings[p,"dimswithG"]),
                               p = p,
                               settings = settings,
                               family = sorted_family)
      }
      mcmc$param_names[[p]] <- colnames(aux)
      cnames[(index+1):(index+settings[p,"dimswithG"])] <- colnames(aux)
      index <- index + settings[p,"dimswithG"]
    }else{
      mcmc$param_names[[p]] <- c()
    }
  }

  if(howsave == "matrix"){
    # mcmc$draws <- data.frame()
    mcmc$draws <- matrix(0, nrow = nchains*iter, ncol = 1+nparams)
    mcmc$draws[,1] = rep(1:nchains, each=iter)
    colnames(mcmc$draws) <- cnames
    mcmc$draws = as.data.frame(mcmc$draws)
  }
  
  if((howsave == "list") || (howsave == "cmcmc")){
    mcmc$draws <- list()
  }
  
  ### Now prepare inits for each chain
  # and generate chains
  for(ch in 1:nchains){
    initsch <- inits[[ch]]
    # cinits
    cinits <- list()
    cinits$w <- initsch$w
    cinits$e0 <- initsch$e0
    cinits$U <- initsch$U-1 
    cinits$pUig <- c(t(initsch$pUig))
    
    cinits$tau_num <- numeric()
    if(varying["tau_num"]){
      for(g in 1:G){
        for(y in Nums){
          cinits$tau_num <- c(cinits$tau_num, initsch$tau_num[[y]][[g]])
        }
      }
    }else{
      cinits$tau_num <- unlist(initsch$tau_num)
    }
    
    cinits$beta_num_fix <- unlist(initsch$beta_num_fix)
    cinits$beta_num <- numeric()
    for(g in 1:G){
      for(y in Nums){
        cinits$beta_num <- c(cinits$beta_num, initsch$beta_num[[y]][[g]])
      }
    }
    
    cinits$beta_poi_fix <- unlist(initsch$beta_poi_fix)
    cinits$beta_poi <- numeric()
    for(g in 1:G){
      for(y in Pois){
        cinits$beta_poi <- c(cinits$beta_poi, initsch$beta_poi[[y]][[g]])
      }
    }
    
    cinits$beta_bin_fix <- unlist(initsch$beta_bin_fix)
    cinits$beta_bin <- numeric()
    for(g in 1:G){
      for(y in Bins){
        cinits$beta_bin <- c(cinits$beta_bin, initsch$beta_bin[[y]][[g]])
      }
    }
    
    cinits$beta_ord_fix <- unlist(initsch$beta_ord_fix)
    cinits$beta_ord <- numeric()
    for(g in 1:G){
      for(y in Ords){
        cinits$beta_ord <- c(cinits$beta_ord, initsch$beta_ord[[y]][[g]])
      }
    }
    
    cinits$c_ord <- numeric()
    if(varying["c_ord"]){
      for(g in 1:G){
        for(y in Ords){
          cinits$c_ord <- c(cinits$c_ord, initsch$c_ord[[y]][[g]])
        }
      }
    }else{
      cinits$c_ord <- unlist(initsch$c_ord)
    }
    
    cinits$a_ord <- numeric()
    if(varying["c_ord"]){
      for(g in 1:G){
        for(y in Ords){
          cinits$a_ord <- c(cinits$a_ord, initsch$a_ord[[y]][[g]])
        }
      }
    }else{
      cinits$a_ord <- unlist(initsch$a_ord)
    }
    
    cinits$pi_ord <- numeric()
    if(varying["c_ord"]){
      for(g in 1:G){
        for(y in Ords){
          cinits$pi_ord <- c(cinits$pi_ord, initsch$pi_ord[[y]][[g]])
        }
      }
    }else{
      cinits$pi_ord <- unlist(initsch$pi_ord)
    }
    
    cinits$beta_cat_fix <- unlist(initsch$beta_cat_fix)
    cinits$beta_cat <- numeric()
    for(g in 1:G){
      for(y in Cats){
        cinits$beta_cat <- c(cinits$beta_cat, initsch$beta_cat[[y]][[g]])
      }
    }
    
    # first [[1]], then regressor 1, then rows for different k values
    
    if(varying["InvSigma"]){
      InvSigma2 <- list()
      for(g in 1:G){
        InvSigma2[[g]] <- initsch$InvSigma[[g]][upper.tri(initsch$InvSigma[[g]], diag = T)]
      }
      cinits$InvSigma <- unlist(InvSigma2)
    }else{
      cinits$InvSigma <- initsch$InvSigma[upper.tri(initsch$InvSigma, diag = T)]
    }
    if(varying["InvQ"]){
      InvQ2 <- list()
      for(g in 1:G){
        InvQ2[[g]] <- initsch$InvQ[[g]][upper.tri(initsch$InvQ[[g]], diag = T)]
      }
      cinits$InvQ <- unlist(InvQ2)
    }else{
      cinits$InvQ <- initsch$InvQ[upper.tri(initsch$InvQ, diag = T)]
    }
    cinits$b <- c(t(initsch$b))
    cinits$naY <- rep(0.0, cdimswithG["naY"]) # does not matter what --> it is updated as first
    
    cinits <- unlist(cinits) # will be delivered as a vector
    #which(is.na(cinits))
    #cinits_U <- initsch$U - 1 # need values 0, ..., G-1
    
    message(paste0("Inits for chain ", ch, " are ready.\n"))
    message("Triggering C function.\n")
    #summary(cinits)
    
    #dyn.load(paste0(ROOT, "Cfun/Metropolis_within_Gibbs_MBC_NumPoiBinOrdCat.dll"))
    #dyn.unload(paste0(ROOT, "Cfun/Metropolis_within_Gibbs_MBC_NumPoiBinOrdCat.dll"))
    #save.image(paste0(ROOT, "image.RData"))
    #load(paste0(ROOT, "image.RData"))
    
    #system.time(
    cmcmc <-
      .C(C_Metropolis_within_Gibbs_MBC_NumPoiBinOrdCat,
         Id        = as.integer(cId),
         Y         = as.double(cY),
         isYna     = as.integer(cisYna),
         X         = as.double(cX),
         varying   = as.integer(varying),
         save      = as.integer(unlist(csave)),
         vecparam  = as.double(cparam), # passed as vector of values 
         vecinits  = as.double(cinits), # passed as vector of double values (except U)
         veclast   = double(length(cinits)), # passed as vector of double values (except U)
         # parameters describing dimensions
         chain     = as.integer(ch), # number of the chain
         G         = as.integer(G), # number of classes
         iter      = as.integer(iter), # total number of generated states
         N         = as.integer(N), # total number of observations
         n         = as.integer(n), # total number of subjects (different ids in the dataset)
         nY        = as.integer(nY), # 4 numbers: counts of Nums,  Bins, Ords and Cats variables
         FormulaF  = as.integer(cFormulaF), # numbers of columns of X that should be used for FIXED  effects of modelled responses
         FormulaG  = as.integer(cFormulaG), # numbers of columns of X that should be used for GROUP-SPECIFIC  effects of modelled responses
         FormulaR  = as.integer(cFormulaR), # numbers of columns of X that should be used for RANDOM effects of modelled responses
         FormulaO  = as.integer(cFormulaO), # numbers of columns of X that should be used for OFFSET effects of modelled responses
         nfix      = as.integer(cnfix),
         ngrp      = as.integer(cngrp),
         nran      = as.integer(cnran),
         noff      = as.integer(cnoff),
         Kord      = as.integer(Kord), # the counts of categories of ordinal variables (-1)
         Kcat      = as.integer(Kcat), # the counts of categories of categorical variables (-1)
         dims      = as.integer(cdims), # the length of subarray that corresponds to one state (disected by various parameters)
         dimswithG = as.integer(cdimswithG), # the length of subarray that corresponds to one state 
         # (disected by various parameters, also multiplication by K incorporated when such parameters is class-specific)
         # arrays to store generated states
         beta_num_fix= cstore$beta_num_fix,
         beta_num    = cstore$beta_num,
         tau_num     = cstore$tau_num,
         sd_num      = cstore$sd_num,
         var_num     = cstore$var_num,
         beta_poi_fix= cstore$beta_poi_fix,
         beta_poi    = cstore$beta_poi,
         beta_bin_fix= cstore$beta_bin_fix,
         beta_bin    = cstore$beta_bin,
         beta_ord_fix= cstore$beta_ord_fix,
         beta_ord    = cstore$beta_ord,
         a_ord       = cstore$a_ord,
         c_ord       = cstore$c_ord,
         pi_ord      = cstore$pi_ord,
         beta_cat_fix= cstore$beta_cat_fix,
         beta_cat    = cstore$beta_cat,
         InvSigma    = cstore$InvSigma,
         Sigma       = cstore$Sigma,
         sdSigma     = cstore$sdSigma,
         corSigma    = cstore$corSigma,
         detInvSigma = cstore$detInvSigma,
         InvQ        = cstore$InvQ,
         Q           = cstore$Q,
         detInvQ     = cstore$detInvQ,
         b           = cstore$b,
         w           = cstore$w,
         ng          = cstore$ng,
         loglik      = cstore$loglik,
         pUig        = cstore$pUig,
         U           = cstore$U,
         Gplus       = cstore$Gplus,
         e0          = cstore$e0,
         naY         = cstore$naY,
         # Tuning parameters
         vectuningdouble  = as.double(unlist(ctuning$double)),
         vectuninginteger = as.integer(unlist(ctuning$integer))
      )
    #)
    message(paste0("Sampling of chain ", ch, " is completed.\n"))

    ### reconstruction of last state from cmcmc$last_U and cmcmc$veclast
    last[[ch]] <- list()
    # order of construction of inits: "w","U","pUig","tau_num",
    #  "beta_num","beta_bin","beta_ord","c_ord","a_ord","beta_cat"
    #  "InvSigma","InvQ","b"
    if(howsave != "cmcmc"){
      last[[ch]]$w <- cmcmc$veclast[1:(lastcopied <- cdims["w"])]
      last[[ch]]$e0 <- cmcmc$veclast[(lastcopied + 1):(lastcopied <- lastcopied + 1)]
      last[[ch]]$U <- 1+cmcmc$veclast[(lastcopied + 1):(lastcopied <- lastcopied + n)]
      last[[ch]]$pUig <- matrix(cmcmc$veclast[(lastcopied + 1):(lastcopied <- lastcopied + n*G)],
                                n, G, byrow = T)

     # tau
      if(varying["tau_num"]){
        last[[ch]]$tau_num <- list()
        for(g in 1:G){
          for(y in Nums){
            last[[ch]]$tau_num[[y]][[g]] <- cmcmc$veclast[(lastcopied + 1):(lastcopied <- lastcopied + 1)]
          }
        }
      }else{
        for(y in Nums){
          last[[ch]]$tau_num[[y]] <- cmcmc$veclast[(lastcopied + 1):(lastcopied <- lastcopied + 1)]
        }
      }
      
      # beta_num
      last[[ch]]$beta_num_fix <- list()
      for(y in Nums){
        if(nfix[y]>0){
          last[[ch]]$beta_num_fix[[y]] <- cmcmc$veclast[(lastcopied + 1):(lastcopied <- lastcopied + nfix[y])]
        }else{
          last[[ch]]$beta_num_fix[[y]] <- numeric()
        }
      }
      last[[ch]]$beta_num <- list()
      for(y in Nums){
        last[[ch]]$beta_num[[y]] <- list()
      }
      for(g in 1:G){
        for(y in Nums){
          if(ngrp[y]>0){
            last[[ch]]$beta_num[[y]][[g]] <- cmcmc$veclast[(lastcopied + 1):(lastcopied <- lastcopied + ngrp[y])]
          }else{
            last[[ch]]$beta_num[[y]][[g]] <- numeric()
          }
        }
      }
      
      # beta_poi
      last[[ch]]$beta_poi_fix <- list()
      for(y in Pois){
        if(nfix[y]>0){
          last[[ch]]$beta_poi_fix[[y]] <- cmcmc$veclast[(lastcopied + 1):(lastcopied <- lastcopied + nfix[y])]
        }else{
          last[[ch]]$beta_poi_fix[[y]] <- numeric()
        }
      }
      last[[ch]]$beta_poi <- list()
      for(y in Pois){
        last[[ch]]$beta_poi[[y]] <- list()
      }
      for(g in 1:G){
        for(y in Pois){
          if(ngrp[y]>0){
            last[[ch]]$beta_poi[[y]][[g]] <- cmcmc$veclast[(lastcopied + 1):(lastcopied <- lastcopied + ngrp[y])]
          }else{
            last[[ch]]$beta_poi[[y]][[g]] <- numeric()
          }
        }
      }
      
      # beta_bin
      last[[ch]]$beta_bin_fix <- list()
      for(y in Bins){
        if(nfix[y]>0){
          last[[ch]]$beta_bin_fix[[y]] <- cmcmc$veclast[(lastcopied + 1):(lastcopied <- lastcopied + nfix[y])]
        }else{
          last[[ch]]$beta_bin_fix[[y]] <- numeric()
        }
      }
      last[[ch]]$beta_bin <- list()
      for(y in Bins){
        last[[ch]]$beta_bin[[y]] <- list()
      }
      for(g in 1:G){
        for(y in Bins){
          if(ngrp[y]>0){
            last[[ch]]$beta_bin[[y]][[g]] <- cmcmc$veclast[(lastcopied + 1):(lastcopied <- lastcopied + ngrp[y])]
          }else{
            last[[ch]]$beta_bin[[y]][[g]] <- numeric()
          }
        }
      }
      
      # beta_ord
      last[[ch]]$beta_ord_fix <- list()
      for(y in Ords){
        if(nfix[y]>0){
          last[[ch]]$beta_ord_fix[[y]] <- cmcmc$veclast[(lastcopied + 1):(lastcopied <- lastcopied + nfix[y])]
        }else{
          last[[ch]]$beta_ord_fix[[y]] <- numeric()
        }
      }
      last[[ch]]$beta_ord <- list()
      for(y in Ords){
        last[[ch]]$beta_ord[[y]] <- list()
      }
      for(g in 1:G){
        for(y in Ords){
          if(ngrp[y]>0){
            last[[ch]]$beta_ord[[y]][[g]] <- cmcmc$veclast[(lastcopied + 1):(lastcopied <- lastcopied + ngrp[y])]
          }else{
            last[[ch]]$beta_ord[[y]][[g]] <- numeric()
          }
        }
      }
      
      # c_ord
      if(varying["c_ord"]){
        last[[ch]]$c_ord <- list()
        for(y in Ords){
          last[[ch]]$c_ord[[y]] <- list()
        }
        for(g in 1:G){
          for(y in Ords){
            last[[ch]]$c_ord[[y]][[g]] <- cmcmc$veclast[(lastcopied + 1):(lastcopied <- lastcopied + Kord[y])]
          }
        }
      }else{
        last[[ch]]$c_ord <- list()
        for(y in Ords){
          last[[ch]]$c_ord[[y]] <- cmcmc$veclast[(lastcopied + 1):(lastcopied <- lastcopied + Kord[y])]
        }
      }
      # a_ord
      if(varying["c_ord"]){
        last[[ch]]$a_ord <- list()
        for(y in Ords){
          last[[ch]]$a_ord[[y]] <- list()
        }
        for(g in 1:G){
          for(y in Ords){
            last[[ch]]$a_ord[[y]][[g]] <- cmcmc$veclast[(lastcopied + 1):(lastcopied <- lastcopied + Kord[y])]
          }
        }
      }else{
        last[[ch]]$a_ord <- list()
        for(y in Ords){
          last[[ch]]$a_ord[[y]] <- cmcmc$veclast[(lastcopied + 1):(lastcopied <- lastcopied + Kord[y])]
        }
      }
      # pi_ord
      if(varying["c_ord"]){
        last[[ch]]$pi_ord <- list()
        for(y in Ords){
          last[[ch]]$pi_ord[[y]] <- list()
        }
        for(g in 1:G){
          for(y in Ords){
            last[[ch]]$pi_ord[[y]][[g]] <- cmcmc$veclast[(lastcopied + 1):(lastcopied <- lastcopied + Kord[y]+1)]
          }
        }
      }else{
        last[[ch]]$pi_ord <- list()
        for(y in Ords){
          last[[ch]]$pi_ord[[y]] <- cmcmc$veclast[(lastcopied + 1):(lastcopied <- lastcopied + Kord[y]+1)]
        }
      }
      
      # beta_cat
      last[[ch]]$beta_cat_fix <- list()
      for(y in Cats){
        if(nfix[y]>0){
          last[[ch]]$beta_cat_fix[[y]] <- matrix(cmcmc$veclast[(lastcopied + 1):(lastcopied <- lastcopied + nfix[y]*Kcat[y])],
                                                 nrow = nfix[y])
        }else{
          last[[ch]]$beta_cat_fix[[y]] <- numeric()
        }
      }
      last[[ch]]$beta_cat <- list()
      for(y in Cats){
        last[[ch]]$beta_cat[[y]] <- list()
      }
      for(g in 1:G){
        for(y in Cats){
          if(ngrp[y]>0){
            last[[ch]]$beta_cat[[y]][[g]] <- matrix(cmcmc$veclast[(lastcopied + 1):(lastcopied <- lastcopied + ngrp[y]*Kcat[y])],
                                                    nrow = ngrp[y])
          }else{
            last[[ch]]$beta_cat[[y]][[g]] <- numeric()
          }
        }
      }
      # InvSigma
      if(varying["InvSigma"]){
        last[[ch]]$InvSigma <- list()
        for(g in 1:G){
          pommatrix <- matrix(0, totnran, totnran)
          pomvec <- cmcmc$veclast[(lastcopied + 1):(lastcopied <- lastcopied + cdims["InvSigma"])]
          pommatrix[upper.tri(pommatrix, diag = T)] <- pomvec
          pommatrix <- t(pommatrix)
          pommatrix[upper.tri(pommatrix, diag = T)] <- pomvec
          last[[ch]]$InvSigma[[g]] <- pommatrix
        }
      }else{
        pommatrix <- matrix(0, totnran, totnran)
        pomvec <- cmcmc$veclast[(lastcopied + 1):(lastcopied <- lastcopied + cdims["InvSigma"])]
        pommatrix[upper.tri(pommatrix, diag = T)] <- pomvec
        pommatrix <- t(pommatrix)
        pommatrix[upper.tri(pommatrix, diag = T)] <- pomvec
        last[[ch]]$InvSigma <- pommatrix
      }
      # InvQ
      if(varying["InvQ"]){
        last[[ch]]$InvQ <- list()
        for(g in 1:G){
          pommatrix <- matrix(0, totnran, totnran)
          pomvec <- cmcmc$veclast[(lastcopied + 1):(lastcopied <- lastcopied + cdims["InvQ"])]
          pommatrix[upper.tri(pommatrix, diag = T)] <- pomvec
          pommatrix <- t(pommatrix)
          pommatrix[upper.tri(pommatrix, diag = T)] <- pomvec
          last[[ch]]$InvQ[[g]] <- pommatrix
        }
      }else{
        pommatrix <- matrix(0, totnran, totnran)
        pomvec <- cmcmc$veclast[(lastcopied + 1):(lastcopied <- lastcopied + cdims["InvQ"])]
        pommatrix[upper.tri(pommatrix, diag = T)] <- pomvec
        pommatrix <- t(pommatrix)
        pommatrix[upper.tri(pommatrix, diag = T)] <- pomvec
        last[[ch]]$InvQ <- pommatrix
      }
      # b
      last[[ch]]$b <- matrix(cmcmc$veclast[(lastcopied + 1):(lastcopied <- lastcopied + cdimswithG["b"])],
                             n, totnran, byrow = T)
      if(cdims["naY"]>0){
        if(varying["naY"]){
          last[[ch]]$naY <- list()
          for(g in 1:G){
            last[[ch]]$naY[[g]] <- cmcmc$veclast[(lastcopied + 1):(lastcopied <- lastcopied + cdims["naY"])]
          }
        }else{
          last[[ch]]$naY <- cmcmc$veclast[(lastcopied + 1):(lastcopied <- lastcopied + cdims["naY"])]
        }
      }
    }
    message(paste0("Saving the last state of chain ", ch, " is completed.\n"))

    ### Saving results as structured list (converted from cmcmc)
    # uses previously created matrix settings
    if(howsave == "cmcmc"){
      mcmc$draws[[ch]] <- cmcmc
      Usamples[[ch]] <- 1 + cmcmc[["U"]]
    }
    
    if(howsave == "list"){
      # results will be returned in structured list - the same way as original R function
      chain <- list()
      for(p in params){
        if(p == "U"){
            cmcmc[[p]] = Usamples[[ch]] = 1 + cmcmc[[p]] # to make cluster allocations in [1, ..., G]
          }
        if(settings[p, "save"] & settings[p,"dims"] > 0){
          if(settings[p,"ydepd1"] & settings[p,"ydepd2"]){
            chain[[p]] <- from_C_to_list(values = cmcmc[[p]],
                                         p = p,
                                         settings = settings,
                                         yspecd1 = yspecd1[[p]],
                                         yspecd2 = yspecd2[[p]],
                                         family = sorted_family)
          }
          if(settings[p,"ydepd1"] & !settings[p,"ydepd2"]){
            chain[[p]] <- from_C_to_list(values = cmcmc[[p]],
                                         p = p,
                                         settings = settings,
                                         yspecd1 = yspecd1[[p]],
                                         family = sorted_family)
          }
          if(!settings[p,"ydepd1"] & settings[p,"ydepd2"]){
            chain[[p]] <- from_C_to_list(values = cmcmc[[p]],
                                         p = p,
                                         settings = settings,
                                         yspecd2 = yspecd2[[p]],
                                         family = sorted_family)
          }
          if(!settings[p,"ydepd1"] & !settings[p,"ydepd2"]){
            chain[[p]] <- from_C_to_list(values = cmcmc[[p]],
                                         p = p,
                                         settings = settings,
                                         family = sorted_family)
          }
        }
      }

      mcmc$draws[[ch]] <- chain
      message(paste0("Saving the chain ", ch, " into lists is completed.\n"))
    }

    if(howsave == "matrix"){
      # results will be returned in matrix (row = state number, col = variable)
      # column chain will distinguish from what chains it comes from
      iters = ((ch-1)*iter+1):(ch*iter)
      for(p in params){
        if(p == "U"){
          cmcmc[[p]] = Usamples[[ch]] = 1 + cmcmc[[p]] # to make cluster allocations in [1, ..., G]
        }
        if(settings[p,"save"] & settings[p,"dims"] > 0){
          if(settings[p,"ydepd1"] & settings[p,"ydepd2"]){
            aux = from_C_to_matrix(values = cmcmc[[p]],
                                   p = p,
                                   settings = settings,
                                   yspecd1 = yspecd1[[p]],
                                   yspecd2 = yspecd2[[p]],
                                   family = sorted_family)
          }
          if(settings[p,"ydepd1"] & !settings[p,"ydepd2"]){
            aux = from_C_to_matrix(values = cmcmc[[p]],
                                   p = p,
                                   settings = settings,
                                   yspecd1 = yspecd1[[p]],
                                   family = sorted_family)
          }
          if(!settings[p,"ydepd1"] & settings[p,"ydepd2"]){
            aux = from_C_to_matrix(values = cmcmc[[p]],
                                   p = p,
                                   settings = settings,
                                   yspecd2 = yspecd2[[p]],
                                   family = sorted_family)
          }
          if(!settings[p,"ydepd1"] & !settings[p,"ydepd2"]){
            aux = from_C_to_matrix(values = cmcmc[[p]],
                                   p = p,
                                   settings = settings,
                                   family = sorted_family)
          }
          mcmc$draws[iters, colnames(aux)] <- aux
        }
      }

      message(paste0("Saving the chain ", ch, " into one big matrix is completed.\n"))
    }
  } # end of chain
  message("All chains have been already sampled.\n")
  
  
  ## U-based clustering
  # the most common Gplus value
  if(howsave == "matrix"){
    modeGplus <- as.numeric(tapply(mcmc$draws$Gplus, mcmc$draws$chain, 
                                   function(x){names(table(x))[which.max(table(x))]}))
  }else{
    modeGplus <- as.numeric(lapply(mcmc$draws, 
                                   function(x){
                                     tab <- table(x$Gplus)
                                     names(tab)[which.max(tab)]}))
  }
  
  fun <- function(U){
    tab = table(U)
    if(clusteringU_modeGplus){
      tab <- tab[as.character(clusters[[ch]])]
    }
    c(as.numeric(names(which.max(tab))), max(tab, na.rm=TRUE)/iter)
  }
    
    
  clusters <- list()
  clusteringU <- matrix(-1, nrow = as.numeric(n), ncol = nchains)
  certaintyU <- matrix(-1, nrow = as.numeric(n), ncol = nchains)
  for(ch in 1:nchains){
    Usamples[[ch]] <- from_C_to_matrix(values = Usamples[[ch]],
                                       p = "U",
                                       settings = settings,
                                       family = sorted_family)
    TAB <- table(unlist(Usamples[[ch]]))
    TAB <- TAB[order(TAB, decreasing = T)]
    # Only the top modeGplus in frequency are considered
    clusters[[ch]] <- sort(as.numeric(names(TAB)[1:modeGplus[ch]]))
    
    ccU <- apply(Usamples[[ch]], 2, fun)
    clusteringU[,ch] <- ccU[1, ]
    certaintyU[,ch] <- ccU[2, ]
  }
  # units with frequency ratio smaller than threshold remain unclassified (group 0)
  # clusteringU[certaintyU < threshold] <- 0
  rownames(clusteringU) <- rownames(certaintyU) <- unique_ids
  
  
  mcmc$modeGplus <- modeGplus
  mcmc$clusters <- clusters
  mcmc$clustering <- clusteringU
  mcmc$certainty <- certaintyU
  
  
  mcmc$inits <- inits
  mcmc$last <- last
  mcmc$InitType <- ifelse(initsgiven, "given", "created")

  mcmc$call <- c("###------------------------------------------------------------###\n")
  mcmc$call <- paste0(mcmc$call, "### Model Based Clustering for Generalized Linear Mixed Models ###\n")
  mcmc$call <- paste0(mcmc$call, "###------------------------------------------------------------###\n")
  mcmc$call <- paste0(mcmc$call, "N = sample size: ", N, "\n")
  mcmc$call <- paste0(mcmc$call, "n = number of units: ", n, "\n")
  mcmc$call <- paste0(mcmc$call, "observations per unit (minimum): ", min(nsubj), " \n")
  mcmc$call <- paste0(mcmc$call, "observations per unit (median): ", median(nsubj), " \n")
  mcmc$call <- paste0(mcmc$call, "observations per unit (maximum): ", max(nsubj), " \n")
  mcmc$call <- paste0(mcmc$call, "G = maximal number of clusters: ", G, "\n")
  mcmc$call <- paste0(mcmc$call, "\nOverview of outcomes: \n")
  mcmc$call <- paste0(mcmc$call, "num: ", paste0(Nums, collapse = ", "), "\n")
  mcmc$call <- paste0(mcmc$call, "poi: ", paste0(Pois, collapse = ", "), "\n")
  mcmc$call <- paste0(mcmc$call, "bin: ", paste0(Bins, collapse = ", "), "\n")
  mcmc$call <- paste0(mcmc$call, "ord: ", paste0(Ords, collapse = ", "), "\n")
  mcmc$call <- paste0(mcmc$call, "cat: ", paste0(Cats, collapse = ", "), "\n")
  mcmc$call <- paste0(mcmc$call, "\nFixed effects invariant towards clustering: \n")
  for(y in Ys){
    mcmc$call <- paste0(mcmc$call, y, " ~ ", formula[[y]]$fixed[2], "\n")
  }
  mcmc$call <- paste0(mcmc$call, "\nRandom effects invariant towards clustering: \n")
  for(y in Ys){
    mcmc$call <- paste0(mcmc$call, y, " ~ ", formula[[y]]$random[2], "\n")
  }
  mcmc$call <- paste0(mcmc$call, "\nThe structure of covariance matrix Sigma for random effects: \n")
  for(y in Ys){
    for(l in length(lrannames[[y]])){
      paste0(mcmc$call, y, "_", lrannames[[y]][l], "\n")
    }
  }
  mcmc$call <- paste0(mcmc$call, "\nGroup-specific formulae for outcomes: \n")
  for(y in Ys){
    mcmc$call <- paste0(mcmc$call, y, " ~ ", formula[[y]]$group[2], "\n")
  }
  gparams <- rownames(settings)[settings$save & settings$isspec]
  gparams <- gparams[!grepl("beta", gparams)]
  mcmc$call <- paste0(mcmc$call, "\nOther saved group-specific parameters: ", 
                      paste0(gparams, collapse = ", "), "\n")
  mcmc$call <- paste0(mcmc$call, "\nMCMC sampling: \n")
  mcmc$call <- paste0(mcmc$call, "nchains = number of sampled chains: ", nchains, "\n")
  mcmc$call <- paste0(mcmc$call, "iter = length of the chain: ", iter, "\n")
  mcmc$call <- paste0(mcmc$call, ifelse(initsgiven, 
                                        "Initialized with given values.", 
                                        "Initialized by uniformely random partition."),
                      "\n")
  

  mcmc$iter <- iter
  mcmc$nchains <- nchains
  mcmc$n <- n
  mcmc$numbered_unique_ids <- numbered_unique_ids
  mcmc$Nums <- Nums
  mcmc$Pois <- Pois
  mcmc$Bins <- Bins
  mcmc$Ords <- Ords
  mcmc$Cats <- Cats
  mcmc$family <- sorted_family
  mcmc$formula <- formula
  mcmc$G <- G
  mcmc$varying <- varying
  mcmc$save <- save
  mcmc$howsave <- howsave
  mcmc$param <- param
  mcmc$nfix <- nfix
  mcmc$ngrp <- ngrp
  mcmc$nran <- nran
  mcmc$noff <- noff
  mcmc$totnran <- totnran
  mcmc$lfixnames <- lfixnames
  mcmc$lgrpnames <- lgrpnames
  mcmc$lrannames <- lrannames
  mcmc$loffnames <- loffnames
  mcmc$fixnames <- fixnames
  mcmc$grpnames <- grpnames
  mcmc$rannames <- rannames
  mcmc$offnames <- offnames
  mcmc$nY <- nY
  mcmc$Kord <- Kord
  mcmc$Kcat <- Kcat
  mcmc$settings <- settings
  mcmc$tuning <- tuning
  mcmc$yspecd1 <- yspecd1
  mcmc$yspecd2 <- yspecd2
  mcmc$isYna <- cisYna
  
  class(mcmc) <- "clustglmm"
  return(mcmc)
  
}
