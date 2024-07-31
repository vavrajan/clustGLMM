plot_num_vs_x_grouped <- 
function(data, y = "y", x = "x", group = "g", id = "i", units, xbreaks,
         add_lowess = TRUE, legend_placement = "topright"){
  # opar <- par(no.readonly = TRUE)
  # on.exit(par(opar))
  
  if(missing(units)){
    units <- unique(data[,id])
  }
  
  gs = sort(unique(data[data[,group] != 0,group]))
  G = length(gs)
  # color definition
  has0 <- (sum(data[,group] == 0) > 1)
  COL <- rainbow_hcl(G, c = 80, l = 70)
  dCOL <- rainbow_hcl(G, c = 80, l = 50)
  names(COL) <- names(dCOL) <- ggs <- as.character(gs)
  if(has0){
    COL <- c("grey70", COL)
    dCOL <- c("grey40", dCOL)
    names(COL) <- names(dCOL) <- c("0", gs)
    ggs <- c(0, gs)
  }
  
  plot(x = c(0,1), y = c(0,1), type = "n", 
       xlim = range(data[,x], na.rm=TRUE),
       ylim = range(data[,y], na.rm=TRUE),
       xlab = x,
       ylab = y)
  for(i in units){
    datai <- data[data[[id]] == i, ]
    datai <- datai[order(datai[,x]),] # order by x value
    lines(datai[,x], datai[,y], col = COL[as.character(datai[1,group])], 
          lty = 1, lwd = 1)
  }
  
  # add lowess curves for each group
  if(add_lowess){
    for(g in ggs){
      gdata <- data[data[,group]==g,]
      gdata <- gdata[order(gdata[,x]),]
      lws <- loess.smooth(gdata[,x], gdata[,y])
      lines(lws, col = dCOL[as.character(g)], lwd = 5)
    }
  }
  
  legend(legend_placement, names(dCOL), col = dCOL, lty=1, bty="n", lwd = 2)
}
