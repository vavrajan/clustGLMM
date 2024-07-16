plot_kerneldensity <-
function(mcmc, scalar_samples, 
         what, gspec, dimspec, yspec,
         burnin = 0, thin = 1, chains,
         dolegend = FALSE,
         COL, labcex = 1, whatLAB,
         restrictedtoGplus = FALSE){
  opar <- par(mar = c(3.5,3.5,0.8,0.8))
  on.exit(opar)
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
  if(missing(COL)){
    if(mcmc$nchains == 1){
      COL <- "grey"
    }else{
      COL <- diverge_hcl(mcmc$nchains, c = 80, l =70)
    }
  }
  
  kernest <- data.frame()
  for(ch in chains){
    if(restrictedtoGplus){
      dd <- density(scalar_samples[[ch]]$y, na.rm = T)
    }else{
      dd <- density(scalar_samples$y[scalar_samples$chain == ch], na.rm = T)
    }
    kernest <- rbind(kernest, data.frame(x=dd$x, y=dd$y, chain = ch))
  }
  plot(c(min(kernest$y),max(kernest$y)) ~ quantile(kernest$x, probs = c(0.03,0.97)), type = "n",
       xlab = "", ylab = "")
  mtext(text = whatLAB, side = 1, line = 2.5, cex = labcex)
  mtext(text = "Density", side = 2, line = 2.5, cex = labcex)
  for(ch in chains){
    lines(kernest$y[kernest$chain==ch] ~ kernest$x[kernest$chain == ch],
          col = COL[ch], lty = 1, lwd = 2)
  }
  #theme_set(theme_classic())
  #options(scipen = 10000)
  #require(scales)
  #g <- ggplot(scalar_samples, aes(y)) + 
  #    geom_density(aes(fill=factor(chain)), alpha=0.7) + 
  #    labs(x = whatLAB,
  #         y = "Density",
  #         fill = "Chain") + 
  #    scale_fill_manual(values = COL[chains]) 
  #g
  if(dolegend){
    legend("topright", legend = chains,
          col = COL[chains], lty =1, bty = "n", title = "Chain", cex = 0.8)
  }
  
}
