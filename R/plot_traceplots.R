plot_traceplots <-
function(mcmc, scalar_samples, 
         what, gspec, dimspec, yspec,
         burnin = 0, thin = 1, chains,
         COL, labcex = 1, whatLAB,
         restrictedtoGplus = FALSE){
  opar <- par(mar = c(3.5,3.5,0.8,0.8))
  on.exit(par(opar))
  if (missing(chains)){chains <- 1:mcmc$nchains}
  #if (missing(gspec)){gspec <- 1}
  if (!missing(yspec)){if(!is.element(yspec, c(mcmc$Nums, mcmc$Pois, mcmc$Ords, mcmc$Bins, mcmc$Cats))){stop("Unexisting response variable specified.")}}
  #if (missing(yspec)){yspec <- ifelse(is.null(mcmc$Nums), mcmc$Ords[1], mcmc$Nums[1])}
  #if (missing(yspec)){yspec <- ""}
  if (missing(scalar_samples)){
    scalar_samples <- get_scalar_samples(mcmc = mcmc, what = what, 
                                         gspec = gspec, dimspec = dimspec, yspec = yspec,
                                         burnin = burnin, thin = thin, chains = chains,
                                         restrictedtoGplus = restrictedtoGplus)
  }
  if (missing(whatLAB)){
    whatLAB <- what
    
    if(!missing(yspec)){
      whatLAB <- paste(whatLAB, yspec, sep = "_")
    }
    
    if(!missing(gspec)){
      if(restrictedtoGplus){
        if(gspec > min(unlist(mcmc$modeGplus)[as.character(chains)])){
          stop("You have asked for g such that may not be included in one of the chains")
        }
      }
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
  
  if(restrictedtoGplus){
    Max <- -Inf
    Min <- Inf
    for(ch in chains){
      scalar_samples[[ch]]$y[is.infinite(scalar_samples[[ch]]$y)] <- NA
      Max <- max(c(scalar_samples[[ch]]$y, Max), na.rm = T)
      Min <- min(c(scalar_samples[[ch]]$y, Min), na.rm = T)
    }
  }else{
    scalar_samples$y[is.infinite(scalar_samples$y)] <- NA
    Max <- max(scalar_samples$y, na.rm = T)
    Min <- min(scalar_samples$y, na.rm = T)
  }
  
  if(missing(COL)){
    if(mcmc$nchains == 1){
      COL <- "grey"
    }else{
      COL <- diverge_hcl(mcmc$nchains, c = 80, l =70)
    }
  }
  
  plot(c(Min, Max) ~ c((burnin+1), mcmc$iter), 
       ylim = c(Min-0.05*(Max-Min), Max+0.05*(Max-Min)),
       ylab = "", xlab = "", type= "n")
  mtext(text = "Iteration", side = 1, line = 2.5, cex = labcex)
  mtext(text = whatLAB, side = 2, line = 2.5, cex = labcex)
  for(ch in chains){
    if(restrictedtoGplus){
      yy <- scalar_samples[[ch]]$y
      xx <- scalar_samples[[ch]]$inds[scalar_samples[[ch]]$chain == ch]
    }else{
      yy <- scalar_samples$y[scalar_samples$chain == ch]
      xx <- c(seq(burnin+1, mcmc$iter, by = thin))
    }
    lines(yy ~ xx, col = COL[ch])
  }
  # legend
  #legend("top", legend = chains, ncol = length(chains),
  #       col = COL[chains], lty =1, bty = "n", title = "Chain", cex = 0.8)
}
