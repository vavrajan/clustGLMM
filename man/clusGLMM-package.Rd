\name{clusGLMM-package}
\alias{clusGLMM-package}
\docType{package}
\title{
\packageTitle{clusGLMM}
}
\description{
\packageDescription{clusGLMM}
}
\details{
The DESCRIPTION file:
\packageDESCRIPTION{clusGLMM}
\packageIndices{clusGLMM}
The key function of this package is \code{\link{clusGLMM}} which creates mcmc samples. To call this function the data, formulae and other auxiliary parameters need to be set up carefully before, so please follow the examples or tutorials.
Package also includes many other functions: one for generating longitudinal data of mixed-type (\code{\link{generate_longitudinal_mixed_type_data}}), functions to generate plots (\code{\link{plot_traceplots}}, \code{\link{plot_traceplots_param}}, \code{\link{plot_kerneldensity}}, \code{\link{plot_kerneldensity_param}}, \code{\link{plot_ACF}}, \code{\link{plot_ACF_param}}, \code{\link{plot_ECDF}}, \code{\link{plot_ECDF_param}}, \code{\link{plot_clusters}}, \code{\link{plot_clusters_ECDF_param}}, \code{\link{plot_clusters_kerneldensity_param}}, \code{\link{plot_all}}, \code{\link{plot_num_vs_x_grouped}}, \code{\link{plot_cat_vs_x_grouped}}), functions for post-processing the samples (\code{\link{PermuteClusterSpec}}, \code{\link{RestrictionToModeGplus}}) and many other functions for internal computations. 
}
\author{
\packageAuthor{clusGLMM}
Maintainer: \packageMaintainer{clusGLMM}
}
\references{
Vávra, J., Komárek, A., Grün, B., Malsiner-Walli, G. Clusterwise multivariate regression of mixed-type panel data. Stat Comput 34, 46 (2024). https://doi.org/10.1007/s11222-023-10304-5
}
\keyword{package}
