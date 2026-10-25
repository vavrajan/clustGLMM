plot_acceptance <-
function(mcmc, what, 
         burnin = 0, thin = 1, 
         iterations, chains,
         dolegend = FALSE, legend.pos = "bottomright",
         COL, labcex = 1, whatLAB,
         adjust.ylim = 0){
  opar <- par(mar = c(3.5, 3.5, 2.0, 0.8))
  on.exit(par(opar))
  
  if(missing(chains)){chains = mcmc$chains}
  misiter = missing(iterations)
  
  cummeanaccept <- list()
  if(adjust.ylim > 0){
    YLIM <- c(100, 0)
  }else{
    YLIM <- c(0, 100)
  }
  XLIM <- c(Inf, -Inf)
  
  if(grepl("^b\\[", what)){
    what <- sub(",\\d", "", what)
  }else{
    what <- sub("\\[[0-9]+(,[0-9]+)?\\]", "", what)
  }
  if(missing(whatLAB)){whatLAB = what}
  
  what <- sub("^c_ord", "a_ord", what)
  what <- sub("^pi_ord", "a_ord", what)
  
  for(ch in chains){
    if(misiter){
      # iterations = seq(burnin+1, mcmc$iter, by = thin)
      # iterations <- mcmc$draws[[ch]]$m
      if(!grepl("\\(\\d\\)", what)){
        iterations <- mcmc$draws[[ch]]$m
      }else{
        # to avoid NA values created by post_processing()
        iterations <- mcmc$iterations[[ch]]
      }
      iterations <- iterations[iterations > burnin]
      iterations <- iterations[seq(1, length(iterations), by=thin)]
    }
    
    ind_iterations <- is.element(mcmc$draws[[ch]]$m, iterations) 
    
    cummeanaccept[[ch]] <- data.frame(chain = ch, m = iterations)
      
    if(!is.element(what, colnames(mcmc$accept[[ch]]))){
      stop(paste("Wrong specification, column ", what, 
                 " does not exist in mcmc$accept[[",ch,"]]!"))
    }
    auxcummean <- cumsum(mcmc$accept[[ch]][, what]) / 1:nrow(mcmc$accept[[ch]])
    cummeanaccept[[ch]]$y = 100*auxcummean[ind_iterations]
    if(adjust.ylim > 0){
      YLIM[1] <- min(YLIM[1], quantile(cummeanaccept[[ch]]$y, probs = adjust.ylim/2, na.rm = TRUE))
      YLIM[2] <- max(YLIM[2], quantile(cummeanaccept[[ch]]$y, probs = 1-adjust.ylim/2, na.rm = TRUE))
    }
    XLIM[1] <- min(XLIM[1], min(iterations))
    XLIM[2] <- max(XLIM[2], max(iterations))
  }
  
  if(missing(COL)){
    if(mcmc$nchains == 1){
      COL <- "grey"
    }else{
      COL <- diverge_hcl(mcmc$nchains, c = 80, l =70)
    }
  }
  
  plot(range(cummeanaccept[[ch]]$m), YLIM, 
       ylim = YLIM,
       ylab = "", xlab = "", type= "n")
  mtext(text = "Iteration", side = 1, line = 2.5, cex = labcex)
  mtext(text = "Acceptance [%]", side = 2, line = 2.5, cex = labcex)
  mtext(text = whatLAB, side = 3, line = 0.2, cex = labcex, font = 2)
  for(ch in chains){
    lines(y ~ m, data = cummeanaccept[[ch]], col = COL[ch], lty = 1, lwd = 2)
  }
  # legend
  if(dolegend){
    legend(legend.pos, legend = chains, ncol = length(chains),
           col = COL[chains], lty = 1, lwd = 2, bty = "n", title = "Chain", cex = 0.8)
  }
  
  
}