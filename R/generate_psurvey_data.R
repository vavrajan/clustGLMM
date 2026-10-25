generate_psurvey_data <- function(n){
  
  # require(mvtnorm)
  # require(splines)
  
  ### Some functions to be used
  logit <- function(p){log(p/(1-p))}
  
  logit_inv <- function(x){
    if(x == -Inf){
      return(0)
    }else{
      if(x == Inf){
        return(1)
      }else{
        return(exp(x)/(1+exp(x)))
      }
    }
  }
  
  softmax <- function(x){
    return(as.numeric(exp(x)/sum(exp(x))))
  }
  
  ### Preparation
  G <- 3
  sd_num <- c(10, 10, 10)
  totnran <- 3
  corr <- matrix(c(1, 0.1, 0.5,
                   0.1, 1, 0.8, 
                   0.5, 0.8, 1), byrow = T, ncol = totnran)
  sdSigmas <- c(3, 0.5, 0.5)
  Sigma <- diag(sdSigmas) %*% corr %*% diag(sdSigmas)
  
  degree <- 2
  knots <- c(18, 20, 25, 30, 35, 40)
  inner <- knots[-c(1, length(knots))]
  bound <- knots[c(1, length(knots))]
  xgrid <- seq(18,40,by=0.1)
  splB <- bs(xgrid, knots = inner, Boundary.knots = bound, 
             degree = degree, intercept = FALSE)
  spldim <- dim(splB)[2]
  bsnames <- paste0("bs",1:spldim)
  
  beta_num <- list()
  beta_num[[1]] <- c(100, 1.5, 12, 26.5, 37, 43.5, 46)
  beta_num[[2]] <- c(90, -1.0, 19, 48.5, 68, 82, 83.5)
  beta_num[[3]] <- c(90, -1.0, 20, 52, 72, 87, 93)
  names(beta_num[[1]]) <- names(beta_num[[2]]) <- names(beta_num[[3]]) <- c(
    "Intercept", bsnames
  )
  
  
  beta_poi <- list()
  beta_poi[[1]] <- c(1.4, -0.2, -1.0, -0.6, -0.2, -0.2, 0.0)
  beta_poi[[2]] <- c(1.5, -0.7, -1.8, 0.2, -0.8, 0.5, 0.8)
  beta_poi[[3]] <- c(1.5, -0.8, -1.6, -0.6, -1.0, -0.6, -0.9)
  names(beta_poi[[1]]) <- names(beta_poi[[2]]) <- names(beta_poi[[3]]) <- c(
    "(Intercept)", bsnames
  )
  
  beta_ord <- list() # no intercept here
  beta_ord[[1]] <- c(0, 0, 0, 0, 0, 0)
  beta_ord[[2]] <- c(-1.0, -0.5, 0.0, 0.5, 1.0, 1.5)
  beta_ord[[3]] <- c(-0.5, 0.0, 0.5, 0.0, -1.0, -2.0)
  names(beta_ord[[1]]) <- names(beta_ord[[2]]) <- names(beta_ord[[3]]) <- c(
    bsnames
  )
  
  Kord <- 3
  c_ord <- list()
  c_ord[[1]] <- c(-Inf, -1, 1, Inf)
  c_ord[[2]] <- c(-Inf, -1.5, 0.5, Inf)
  c_ord[[3]] <- c(-Inf, -0.5, 1.5, Inf)
  
  data <- data.frame()
  # sample by subjects - covariate + group + random effects
  ni <- 8 + rpois(n, 1) 
  for(i in 1:n){
    datai <- data.frame(i = i, j=1:ni[i], 
                        age = seq(18,40, length.out=ni[i]) + c(runif(1, 0, 1),
                                                               runif(ni[i]-2, -1, 1),
                                                               runif(1, -1, 0)))
    datai <- datai[order(datai$age),]
    # datai <- transform(datai, x=age-18)
    # spline base
    splx <- bs(datai$age, knots = inner, Boundary.knots = bound, 
               degree = degree, intercept = FALSE)
    for(j in 1:spldim){
      datai[, paste0("bs",j)] <- splx[,j]
    }
    # latent group 
    # g = 1: non-uni -> income earlier, but low after 
    # g = 2: happy uni -> high satisfaction, more trips
    # g = 3: uni hard workers -> low satisfaction, less trips
    datai$g <- rep(sample(1:G, 1, prob = c(0.4, 0.3, 0.3)), ni[i])
    # random intercepts
    bi <- rmvnorm(1, mean = rep(0,totnran), sigma = Sigma)
    datai$b_num <- rep(bi[1], ni[i])
    datai$b_poi <- rep(bi[2], ni[i])
    datai$b_ord <- rep(bi[3], ni[i])
    
    data <- rbind(data, datai)
  }
  N <- dim(data)[1]
  rownames(data) <- 1:N
  
  # calculate predictors
  for(i in 1:N){
    g <- data$g[i]
    # Numeric
    data$pred_num[i] <- beta_num[[g]][1] + sum(beta_num[[g]][bsnames] * data[i,bsnames]) + data$b_num[i]
    # Poisson count
    data$pred_poi[i] <- beta_poi[[g]][1] + sum(beta_poi[[g]][bsnames] * data[i,bsnames]) + data$b_poi[i]
    # Ord
    data$pred_ord[i] <- sum(beta_ord[[g]][bsnames] * data[i,bsnames]) + data$b_ord[i]
    
    data$c_ord1[i] <- c_ord[[g]][2]
    
  }
  
  ## sample outcomes
  # Numeric
  data$income <- rnorm(N, mean = data$pred_num, sd = sd_num[data$g])
  
  # Poisson count
  data$ntrips <- rpois(N, lambda = exp(data$pred_poi))
  
  # Ordinal
  probs <- matrix(0, nrow = N, ncol = Kord)
  for(i in 1:N){
    g <- data$g[i]
    # response sampled by probabilities
    loginvs <- sapply(data$pred_ord[i]-c_ord[[g]], logit_inv)
    probs[i,] <- loginvs[1:Kord] - loginvs[2:(Kord+1)]
    data$satisfaction[i] <- sample(0:(Kord-1), 1, prob = probs[i,])
  }
  data$satisfaction = factor(data$satisfaction, ordered = TRUE)
  
  return(list(data = data, 
              sd_num = sd_num, 
              beta_num = beta_num, 
              beta_poi = beta_poi, 
              beta_ord = beta_ord, 
              c_ord = c_ord))
}

# # Source code to generate the dataset used in the manuscript
# set.seed(31415)
# n <- 500
# ldata <- generate_psurvey_data(n)
# psurvey <- ldata$data[,c("i", "j", "age", "income", "ntrips", "satisfaction")]
# psurvey_latent <- ldata$data[,c("i", "j", "g", "b_num", "b_poi", "b_ord")]
