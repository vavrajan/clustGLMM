plot_clusters <-
function(mcmc, AllData, 
         what, dimspec, yspec, whichg,
         burnin = 0, thin = 1, setparmfrow = TRUE,
         ChainsToPlot, doKern = TRUE, doECDF = TRUE,
         COL, labcex = 0.8, whatLAB, wherelegend = "topright",
         restrictedtoGplus = FALSE){
  if(setparmfrow){
    opar <- par(no.readonly = TRUE)
    on.exit(par(opar))
    par(mfrow = c(1,sum(c(doKern, doECDF))))
    }
  if(is.element(what,names(mcmc$varying))){
    if (!mcmc$varying[what]){stop("This parameter wasn't considered as cluster-specific.")}
  }
  if(missing(ChainsToPlot)){
    if(restrictedtoGplus){
      ChainsToPlot <- c(1)
    }else{
      ChainsToPlot <- 1:mcmc$chains
    }
  }else{
    if(restrictedtoGplus & length(ChainsToPlot) > 1){
      warning("It is recommended to use one chain only when restrictedtoGplus=T, only first element of ChainsToPlot will be used.")
      ChainsToPlot <- ChainsToPlot[1]
    }
  }
    
  #if (missing(yspec)){yspec <- ifelse(is.null(mcmc$Nums), mcmc$Ords[1], mcmc$Nums[1])}
  if(!missing(yspec)){
    if(!is.element(yspec, c(mcmc$Nums, mcmc$Pois, mcmc$Bins, mcmc$Ords, mcmc$Cats))){
      stop("Unexisting response variable specified.")
    }
  }
  G <- ifelse(restrictedtoGplus,
              mcmc$modeGplus[[ChainsToPlot]],
              mcmc$G)
  if(missing(whichg)){
    whichg <- 1:G
  }
  
  if(missing(AllData)){
    AllData <- data.frame()
    for(g in whichg){
      pomData <- FindAllData(mcmc = mcmc, what = what, 
                             gspec = g, dimspec = dimspec, yspec = yspec,
                             burnin = burnin, thin = thin, ChainsToPlot = ChainsToPlot, 
                             restrictedtoGplus = restrictedtoGplus)
      if(restrictedtoGplus){
        pomData <- pomData[[ChainsToPlot]]
      }
      colnames(pomData) <- c(paste0("y(",g,")"), "chain")
      if(is.element("chain", colnames(AllData))){AllData <- cbind(AllData, pomData)
      }else{AllData <- pomData}
    }
  }
  
  if (missing(whatLAB)){
    whatLAB <- what
    
    if(!missing(yspec)){
      whatLAB <- paste(whatLAB, yspec, sep = "_")
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
  if (missing(COL)){COL <- rainbow_hcl(G, c = 80, l =70)}
  
  
  # Density
  if(doKern){
    kernest <- data.frame()
    for(g in whichg){
      dd <- density(AllData[is.element(AllData$chain, ChainsToPlot), paste0("y(",g,")")], na.rm = T)
      kernest <- rbind(kernest, data.frame(x=dd$x, y=dd$y, cluster = g))
    }
    plot(c(min(kernest$y),max(kernest$y)) ~ quantile(kernest$x, probs = c(0.03,0.97)), type = "n",
         xlab = "", ylab = "")
    mtext(text = whatLAB, side = 1, line = 2.5, cex = labcex)
    mtext(text = "Density", side = 2, line = 2.5, cex = labcex)
    for(g in whichg){
      lines(kernest$y[kernest$cluster==g] ~ kernest$x[kernest$cluster==g],
            col = COL[g], lty = 1, lwd = 2)
    }
    legend(wherelegend, legend = whichg,
           col = COL[whichg], lty =1, bty = "n", title = "Cluster", cex = 0.8)
  }
  # ECDF
  if(doECDF){
    plot(x = 1, y = 1, type = "n", xlim = c(min(AllData[,grep("y",colnames(AllData))]), max(AllData[,grep("y",colnames(AllData))])), ylim = c(0,1),
         xlab = "", ylab = "")
    mtext(text = whatLAB, side = 1, line = 2.5, cex = labcex)
    mtext(text = "ECDF", side = 2, line = 2.5, cex = labcex)
    for(g in whichg){
      plot(ecdf(AllData[,paste0("y(", g,")")]), 
           col = COL[g], add = T, verticals = T, pch = NA)
    }
  }
}
