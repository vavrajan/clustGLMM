plot_ACF <-
function(mcmc, scalar_samples, 
         what, gspec, dimspec, yspec,
         burnin = 0, thin = 1, lag.max = 30, iterations, chains,
         COL, move.width = 0.3, labcex = 1, whatLAB){
  opar <- par(mar = c(3.5,3.5,0.8,0.8))
  on.exit(opar)
  if(missing(chains)){chains = mcmc$chains}
  #if (missing(gspec)){gspec <- 1}
  if(!missing(yspec)){if(!is.element(yspec, names(mcmc$family))){stop("Unexisting response variable specified.")}}
  #if (missing(yspec)){yspec <- ifelse(is.null(mcmc$Nums), mcmc$Ords[1], mcmc$Nums[1])}
  #if (missing(yspec)){yspec <- ""}
  if(missing(scalar_samples)){
    scalar_samples <- get_scalar_samples(mcmc = mcmc, what = what,
                                         gspec = gspec, dimspec = dimspec, yspec = yspec,
                                         burnin = burnin, thin = thin, 
                                         iterations = iterations, chains = chains)
  }
  subiter <- max(unlist(lapply(scalar_samples, function(x){dim(x)[1]})), na.rm=TRUE)
  if(lag.max >= subiter){lag.max <- subiter-1}
  
  if(missing(whatLAB)){
    whatLAB <- what
    
    if(!missing(yspec)){
      whatLAB <- paste(whatLAB, yspec, sep = "_")
    }
    
    if(!missing(gspec)){
      whatLAB <- paste0(whatLAB, "(", gspec, ")")
    }
    
    if(missing(dimspec)){
      # no dimension supplied = no need for supplying -> 1-dim parameter
    }else{
      if(length(dimspec) == 1){
        whatLAB <- paste0(whatLAB,"[",dimspec[1],"]")
      }
      if(length(dimspec) == 2){
        whatLAB <- paste0(whatLAB,"[",dimspec[1],",",dimspec[2],"]")
      }
    }
  }
  
  if(missing(COL)){
    if(mcmc$nchains == 1){
      COL <- "grey"
    }else{
      COL <- diverge_hcl(mcmc$nchains, c = 80, l =70)
    }
  }
  
  plot(c(-0.2,1) ~ c(0,lag.max), 
       ylim = c(-0.2,1),
       ylab = "", xlab = "", type= "n")
  mtext(text = "Lag", side = 1, line = 2.5, cex = labcex)
  mtext(text = "ACF", side = 2, line = 2.5, cex = labcex)
  #mtext(whatLAB, side = 3, line = 0.5)
  pom <- seq(-move.width, move.width, length.out = length(chains)+1)
  posun <- (pom[1:length(chains)]+pom[2:(1+length(chains))])/2
  chcount <- 0
  for(ch in chains){
    chcount <- chcount+1
    pomACF <- acf(scalar_samples[[ch]]$y, lag.max = lag.max, plot = F, na.action = na.pass)$acf[,1,1]
    lines(pomACF ~ c(0:lag.max + posun[chcount]),
          col = COL[ch], lwd = 2)
    segments(x0 = c(0:lag.max)+posun[chcount], 
             y0 = 0, y1 = pomACF,
             col = COL[ch], lty = 1, lwd = 1)
    abline(h = 0, col = "grey", lty = 2)
  }
}
