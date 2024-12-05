plot.clustglmm <- function(x, what="ng", which="traceplots", ng_trace_chain_split=TRUE, ...){
  if(missing(what)){
    adepts <- c("ng", "Gplus", "w", "pUig_int", "deviance")
    availables <- intersect(adepts, rownames(x$settings[[x$chains[1]]][x$settings[[x$chains[1]]]$save, ]))
    what <- availables[1]
  }
  if(missing(which)){
    if(what == "pUig_int"){
      which="clusters_kerneldensity"
    }
  }
  
  if((what == "ng") & (which == "traceplots") & ng_trace_chain_split){
    plot_ng_trace_chain_split(x, ...)
  }else{
    switch(which,
           "ACF" = plot_ACF_param(x, what, ...),
           "ECDF" = plot_ECDF_param(x, what, ...),
           "kerneldensity" = plot_kerneldensity_param(x, what, ...),
           "traceplots" = plot_traceplots_param(x, what, ...),
           "clusters_ECDF" = plot_clusters_ECDF_param(x, what, ...),
           "clusters_kerneldensity" = plot_clusters_kerneldensity_param(x, what, ...)
           )
  }
}