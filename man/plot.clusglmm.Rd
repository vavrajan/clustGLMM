\name{plot.clusglmm}
\alias{plot.clusglmm}
\title{
Plot the outputs of \code{clusGLMM}
}
\description{
Triggers plotting functions for selected parameter depending on the choice of the plot.
}
\usage{
plot(x, what="ng", how="traceplots", ...)
}
\arguments{
  \item{x}{
an object of class \code{"clusglmm"}, an output of \code{clusGLMM} function.
}
  \item{what}{
a string containing the name of the parameter. 
}
  \item{how}{
a string declaring the type of the plot to be created.
}
  \item{...}{
other arguments.
}
}
\details{
Depending on \code{how} the following functions are called:
\describe{
  \item{ACF}{creates auto-correlation function plot, calls \code{\link{plot_ACF_param}}.}
  \item{ECDF}{creates empirical cumulative distribution function plot, calls \code{\link{plot_ECDF_param}}.}
  \item{kerneldensity}{creates kernel density estimator plot, calls \code{\link{plot_kerneldensity_param}}.}
  \item{traceplots}{creates traceplots (samples vs iteration), calls \code{\link{plot_traceplots_param}}.}
  \item{clusters_ECDF}{creates empirical cumulative distribution function plot where different clusters are compared (\code{what} needs to be a cluster-specific parameter), calls \code{\link{plot_clusters_ECDF_param}}.}
  \item{clusters_kerneldensity}{creates kernel density estimator plot where different clusters are compared (\code{what} needs to be a cluster-specific parameter), calls \code{\link{plot_clusters_kerneldensity_param}}.}
}

These functions except also \code{burnin} (initial part of the chain to be ignored), \code{thin} (plot only every \code{thin}th sample) and \code{ChainsToPlot} (a vector of chain numbers to be plotted).
}
\value{
  No returned value, only creates plots.
}
\references{
Vávra, J., Komárek, A., Grün, B., Malsiner-Walli, G. Clusterwise multivariate regression of mixed-type panel data. Stat Comput 34, 46 (2024). https://doi.org/10.1007/s11222-023-10304-5
}
\author{
Jan Vávra
}
\seealso{
\code{\link{clusGLMM}}, \code{\link{plot_all}}
}
\examples{

###-----------------------------------------------------------------------------
###     Generate some longitudinal mixed-type data
###-----------------------------------------------------------------------------
corr <- matrix(c(1, -0.5, -0.5, -0.4, -0.4,
                 -0.5, 1, 0.3, 0.4, 0.5,
                 -0.5, 0.3, 1, 0.2, 0.3,
                 -0.4, 0.4, 0.2, 1, 0.1,
                 -0.4, 0.5, 0.3, 0.1, 1), byrow = TRUE, ncol = 5)
sdSigmas <- c(0.5, 0.5, 0.5, 0.5, 0.5)
Sigma <- diag(sdSigmas) %*% corr %*% diag(sdSigmas)

set.seed(31415)
n <- 500
ldata <- generate_longitudinal_mixed_type_data(n = n, n_i = 4,
                     G = 2, timepar = "cross", catpar = "fixed", 
                     catsubjfix = TRUE,
                     xin1 = 1, Sigma = Sigma,
                     Kord = 5, Kcat = 4, kspec_bi_cat = FALSE)
data <- ldata$data

###-----------------------------------------------------------------------------
###     Prepare inputs for clusGLMM
###-----------------------------------------------------------------------------
id <- "i"

family <- c("num", "poi", "bin", "ord", "cat")
names(family) <- Ys <- c("ynum", "ypoi", "ybin", "yord", "ycat")

## Create list of formulae for each outcome
formula <- list()
for(y in Ys){
  formula[[y]] <- list()
  formula[[y]]$fixed <- as.formula("~ f")
  formula[[y]]$group <- as.formula("~ x")
  formula[[y]]$random <- as.formula("~ 1")
  formula[[y]]$offset <- ""
}
formula
# intercept term will be group-specific 
# since it appears in $group which has a preference over $fixed

## All other inputs need to be prepared,
## unless the default values are desired (varying, save, howsave, param, tuning).
## See the tutorial for how these can be adjusted.

###-----------------------------------------------------------------------------
###     Final call of clusGLMM
###-----------------------------------------------------------------------------

iter <- 100    # the length of sampled chain
Gmax <- 5      # maximal number of mixture components to be considered
nchains <- 2   # number of chains to be sampled 

## No inits specified
## Saved in "matrix" format
mcmc <- clusGLMM(formula = formula, id = id, family = family, data = data, 
                 G = Gmax, iter = iter, nchains = nchains,
                 howsave = "matrix")

# traceplots by default                       
plot(mcmc)
plot(mcmc, "Gplus")
plot(mcmc, "beta_num", burnin=10, thin=2, chains=c(1))

# change the plot type
plot(mcmc, "w", "ACF")
plot(mcmc, "sdSigma", "ECDF")
plot(mcmc, "corSigma", "kerneldensity")
plot(mcmc, "beta_bin", "clusters_kerneldensity")
plot(mcmc, "beta_ord", "clusters_ECDF")

}
