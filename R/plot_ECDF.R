plot_ECDF <-
function(mcmc, AllData, 
         what, gspec, dimspec, yspec,
         burnin = 0, thin = 1, ChainsToPlot,
         COL, labcex = 1, whatLAB, restrictedtoGplus = FALSE){
  opar <- par(mar = c(3.5,3.5,0.8,0.8))
  on.exit(opar)
  if (missing(ChainsToPlot)){ChainsToPlot <- 1:mcmc$chains}
  #if (missing(gspec)){gspec <- 1}
  if (!missing(yspec)){if(!is.element(yspec, c(mcmc$Nums, mcmc$Pois, mcmc$Ords, mcmc$Bins, mcmc$Cats))){stop("Unexisting response variable specified.")}}
  #if (missing(yspec)){yspec <- ifelse(is.null(mcmc$Nums), mcmc$Ords[1], mcmc$Nums[1])}
  #if (missing(yspec)){yspec <- ""}
  if (missing(AllData)){
    AllData <- FindAllData(mcmc = mcmc, what = what, 
                           gspec = gspec, dimspec = dimspec, yspec = yspec,
                           burnin = burnin, thin = thin, ChainsToPlot = ChainsToPlot, 
                           restrictedtoGplus = restrictedtoGplus)
  }
  if (missing(whatLAB)){
    whatLAB <- what
    
    if(!missing(yspec)){
      whatLAB <- paste(whatLAB, yspec, sep = "_")
    }
    
    if(!missing(gspec)){
      if(restrictedtoGplus){
        if(gspec > min(unlist(mcmc$modeGplus)[as.character(ChainsToPlot)])){
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
    for(ch in ChainsToPlot){
      Max <- max(c(AllData[[ch]]$y, Max), na.rm = T)
      Min <- min(c(AllData[[ch]]$y, Min), na.rm = T)
    }
  }else{
    Max <- max(AllData$y)
    Min <- min(AllData$y)
  }
  
  #if(missing(whatLAB)){
  #  if(is.element(what, names(mcmc$spec)) & mcmc$spec[what]){
  #    whatLAB <- paste0(what, 
  #                      ifelse(is.list(mcmc[[1]][[what]][[gspec]]), ifelse(yspec=="",paste0("_",yspec), ""), ""),
  #                      "(",gspec,")",
  #                      ifelse(is.na(dimspec[1]),"",
  #                             paste0("[",paste(dimspec,collapse=","),"]")))
  #  }else{
  #    whatLAB <- paste0(what, 
  #                      ifelse(is.list(mcmc[[1]][[what]]), ifelse(yspec=="",paste0("_",yspec), ""), ""),
  #                      ifelse(is.na(dimspec[1]),"",
  #                             paste0("[",paste(dimspec,collapse=","),"]")))
  #  }
  #}
  
  if(missing(COL)){
    if(mcmc$chains == 1){
      COL <- "grey"
    }else{
      COL <- diverge_hcl(mcmc$chains, c = 80, l =70)
    }
  }
  #if (missing(COL)){COL <- topo.colors(mcmc$chains, alpha = 1)}
  
  plot(x = 1, y = 1, type = "n", xlim = c(Min, Max), ylim = c(0,1),
       xlab = "", ylab = "")
  mtext(text = whatLAB, side = 1, line = 2.5, cex = labcex)
  mtext(text = "ECDF", side = 2, line = 2.5, cex = labcex)
  for(ch in ChainsToPlot){
    if(restrictedtoGplus){
      yy <- AllData[[ch]]$y
    }else{
      yy <- AllData$y[AllData$chain == ch]
    }
    plot(ecdf(yy), col = COL[ch], add = T, verticals = T, pch = NA)
  }
}
