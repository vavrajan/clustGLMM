plot.clustglmm <- function(x, what="ng", how="traceplots", ...){
  switch(how,
         "ACF" = plot_ACF_param(x, what, ...),
         "ECDF" = plot_ECDF_param(x, what, ...),
         "kerneldensity" = plot_kerneldensity_param(x, what, ...),
         "traceplots" = plot_traceplots_param(x, what, ...),
         "clusters_ECDF" = plot_clusters_ECDF_param(x, what, ...),
         "clustesr_kerneldensity" = plot_clusters_kerneldensity_param(x, what, ...)
         )
}