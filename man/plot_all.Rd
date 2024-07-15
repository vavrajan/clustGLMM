\name{plot_all}
\alias{plot_all}
\alias{plot_ACF}
\alias{plot_ECDF}
\alias{plot_kerneldensity}
\alias{plot_traceplots}
\alias{plot_ACF_param}
\alias{plot_ECDF_param}
\alias{plot_kerneldensity_param}
\alias{plot_traceplots_param}
\alias{plot_clusters}
\alias{plot_clusters_ECDF_param}
\alias{plot_clusters_kerneldensity_param}
\title{
Plot mcmc samples from \code{clusGLMM}
}
\description{
There are 4 main plots to be created:
\describe{
  \item{\code{ACF}}{auto-correlation function,}
  \item{\code{ECDF}}{empirical cumulative distribution function,}
  \item{\code{kerneldensity}}{one-dimensional kernel density estimate,}
  \item{\code{traceplots}}{parameter samples against the iteration number.}
}
All of these for one specific parameter can be plotted by \code{plot_all}. These functions above work only for a univariate parameter. Plots for all elements of a multivariate parameter can be created using \code{..._param} functions. A special function is \code{plot_clusters} which allows to compare densities or ECDFs of univariate group-specific parameter across non-empty clusters.  
}
\usage{
plot_ACF(mcmc, AllData, what, gspec, dimspec, yspec,
         burnin = 0, thin = 1, lag.max = 30, ChainsToPlot, 
         COL, move.width = 0.3, labcex = 1, whatLAB, restrictedtoGplus = FALSE)
plot_ECDF(mcmc, AllData, what, gspec, dimspec, yspec,
          burnin = 0, thin = 1, ChainsToPlot, 
          COL, labcex = 1, whatLAB, restrictedtoGplus = FALSE)
plot_kerneldensity(mcmc, AllData, what, gspec, dimspec, yspec, 
                   burnin = 0, thin = 1, ChainsToPlot, 
                   dolegend = FALSE, COL, labcex = 1, whatLAB, restrictedtoGplus = FALSE)
plot_traceplots(mcmc, AllData, what, gspec, dimspec, yspec,
                burnin = 0, thin = 1, ChainsToPlot, 
                COL, labcex = 1, whatLAB, restrictedtoGplus = FALSE)
plot_all(mcmc, AllData, what, gspec, dimspec, yspec,
         burnin = 0, thin = 1, lag.max, ChainsToPlot, trueval, setparmfrow = TRUE, 
         COL, move.width = 0.3, whatLAB, restrictedtoGplus = FALSE)
plot_ACF_param(mcmc, what, burnin = 0, thin = 1, lag.max = 30, ChainsToPlot)
plot_ECDF_param(mcmc, what, burnin = 0, thin = 1, ChainsToPlot)
plot_kerneldensity_param(mcmc, what, burnin = 0, thin = 1, ChainsToPlot)
plot_traceplots_param(mcmc, what, burnin = 0, thin = 1, ChainsToPlot)
plot_clusters(mcmc, AllData, what, dimspec, yspec, whichg,
              burnin = 0, thin = 1, setparmfrow = TRUE, 
              ChainsToPlot, doKern = TRUE, doECDF = TRUE,
              COL, labcex = 0.8, whatLAB, wherelegend = "topright",
              restrictedtoGplus = FALSE)
plot_clusters_kerneldensity_param(mcmc, what, burnin = 0, thin = 1, ChainsToPlot)
plot_clusters_ECDF_param(mcmc, what, burnin = 0, thin = 1, ChainsToPlot)
}
\arguments{
  \item{mcmc}{
an output of \code{\link{clusGLMM}} accepted in both formats (list or matrix).
}
  \item{AllData}{
an output of \code{FindAllData}, overrides any given \code{burnin} or \code{thin}.
}
  \item{burnin}{
the length of burn-in period for each chain (number of not included initial values); default is 0.
}
  \item{thin}{
a thinning parameter, only every \code{thin}th value between \code{burnin} and \code{mcmc$iter} is gathered.
}
  \item{what}{
a name of the parameter to be extracted; see \code{names(mcmc$save)} for possible options. 
}
  \item{gspec}{
a group specifier (values in \code{1:mcmc$G}); if not group-specific parameter this \emph{should be missing}.
}
  \item{dimspec}{
a dimension specifier; not specifying this is possible only when the parameter is scalar.
}
  \item{yspec}{
an outcome specifier; if the parameter is not outcome-dependent then this \emph{should be missing}.
}
  \item{lag.max}{
maximal lag value to be plotted in ACF.
}
  \item{ChainsToPlot}{
a vector of chain numbers to be plotted.
}
  \item{trueval}{
a true value of the parameter to be plotted by vertical/horizontal line.
}
  \item{setparmfrow}{
whether \code{par(mfrow = c(2,2))} should be set before creating all 4 plots.
}
  \item{COL}{
colors for different chains; default values are set by \code{diverge_hcl(mcmc$chains, c = 80, l = 70)}.
}
  \item{move.width}{
ACF for different chains are jittered around the actual lag by \code{seq(-move.width, move.width, length.out = length(ChainsToPlot)+1)}
}
  \item{whatLAB}{
the label for the parameter; default corresponds to the label for that parameter that would be used in \emph{matrix} format of the output. 
}
  \item{labcex}{
\code{cex} value for label of plots. 
}
  \item{dolegend}{
a boolean declaring whether legend should be added to the plot.  
}
\item{restrictedtoGplus}{
a boolean declaring whether the values should be taken from post-processed data by \code{\link{RestrictionToModeGplus}}.
}
  \item{whichg}{
a set of cluster numbers to be compared in \code{plot_clusters}.  
}
  \item{doKern}{
a boolean declaring whether kernel density estimates should be plotted for comparing different clusters.
}
  \item{doECDF}{
a boolean declaring whether ECDFs should be plotted for comparing different clusters.
}
  \item{wherelegend}{
  a string that declares where the legend should be placed.
}
}
\details{
  \code{AllData} is a \code{data.frame} of two columns: \code{y} containing chosen samples and \code{chain} declaring the chain number. This could be created using \code{\link{FindAllData}}, see its help file for more details. If missing, \code{FindAllData} is called anyway with the given specifications.
  Iterations are selected by \code{seq(burnin+1, mcmc$iter, by = thin)}.
  Parameter is selected by filling/missing \code{what}, \code{gspec}, \code{yspec}, \code{dimspec}. More details here: \code{\link{FindAllData}}.
  
  Description of the plots:
  \describe{
    \item{\code{ACF}}{Autocorrelation function plots have the lag on x-axis [0, \code{lag.max}] and auto-correlations in interval (-1, 1) on y-axis. By default, \code{lag.max = min(30, mcmc$iter-1)}. Different chains have different colors and are slightly jittered for better visuality, see \code{move.width}. }
    \item{\code{ECDF}}{Cumulative distribution function plots have the range of sampled values on x-axis and the cumulative probabilities in [0, 1] on y-axis. Different colors depict different chains of non-decreasing ECDFs.}
    \item{\code{kerneldensity}}{Kernel density estimator plots have the range of sampled values on x-axis and the density value in [0, max(density)] on y-axis. Different colors depict different chains.}
    \item{\code{traceplot}}{Traceplots have the range of the considered iteration numbers on x-axis and the sampled values of the parameter on y-axis. Different colors depict different chains. Stationarity of these plots is desired for reliable posterior approximation, if the plots are not stationary then consider continuation in sampling. }
  }
  You can plot all 4 of these by \code{plot_all}. By default they are produced in 2 by 2 layout, but if you would prefer some different arrangement you can do that before calling the function and supressing \code{setparmfrow = FALSE}. 
  
  Using the functions above you can only plot a univariate model parameter which is specified by filling/missing \code{what}, \code{gspec}, \code{yspec}, \code{dimspec}. If you would like to plot results for all parameters that are associated with parameter \code{what} then \code{..._param} functions can do this for you automatically. Using function \code{\link{nice_nrow_ncol}} plots are arranged into a most effective grid. Moreover, parameters which are symmetric matrices (with or without diagonal) are arranged as expected and only the upper triangle is filled with plots. See examples below.
  
  Special function in this collection is \code{plot_clusters} which allows to plot kernel density estimates or empirical cumulative distribution function of the same parameters for different cluster labels (\code{whichg}). There are also \code{plot_clusters_kerneldensity_param} and \code{plot_clusters_ECDF_param} functions that create the plots automatically for a set of given parameter \code{what}, however, the plotted cluster labels only include those within the last sampled set of parameters \code{mcmc$last}. 
  
}
\value{
The functions do not return anything but produce one or more plots.
}
\references{
Vávra, J., Komárek, A., Grün, B., Malsiner-Walli, G. Clusterwise multivariate regression of mixed-type panel data. Stat Comput 34, 46 (2024). https://doi.org/10.1007/s11222-023-10304-5
}
\author{
Jan Vávra
}
\seealso{
\code{\link{clusGLMM}}, \code{\link{FindAllData}}, \code{\link{nice_nrow_ncol}}, \code{\link{RestrictionToModeGplus}}.
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
n <- 1000
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

Gmax <- 5     # maximal number of mixture components to be considered
chains <- 2   # number of chains to be sampled 

## No inits specified
mcmc <- clusGLMM(formula = formula, id = id, family = family, data = data, 
                 G = Gmax, iter = 1500, chains = chains)
                       
###-----------------------------------------------------------------------------
###     Creating plots
###-----------------------------------------------------------------------------

## Traceplots of number of non-empty clusters
# Plotting with defaults values:
# thin = 1, burnin = 0, ChainsToPlot = 1:mcmc$chains
opar <- par(mfrow = c(1, 1), mar = c(4, 4, 1, 1))
plot_traceplots(mcmc = mcmc, what = "Gplus")
# Omitting the burn-in-period
plot_traceplots(mcmc = mcmc, what = "Gplus", burnin = 500)
# Thinned plots
plot_traceplots(mcmc = mcmc, what = "Gplus", burnin = 500, thin = 5)

## First fixed beta coefficient for the first numerical outcome
plot_ACF(mcmc = mcmc, what = "beta_num_fix", yspec = mcmc$Nums[1], dimspec = 1)
plot_ECDF(mcmc = mcmc, what = "beta_num_fix", yspec = mcmc$Nums[1], dimspec = 1)
plot_kerneldensity(mcmc = mcmc, what = "beta_num_fix", yspec = mcmc$Nums[1],
                   dimspec = 1, dolegend = TRUE)
plot_traceplots(mcmc = mcmc, what = "beta_num_fix", yspec = mcmc$Nums[1],
                dimspec = 1)
plot_all(mcmc = mcmc, what = "beta_num_fix", yspec = mcmc$Nums[1], dimspec = 1)

## Correlation matrix between the random effects
par(mfrow = c(mcmc$totnran-1, mcmc$totnran-1), mar = c(4, 4, 1, 1))
for (i in 1:(mcmc$totnran-1)) {
  for (j in 2:(mcmc$totnran)) {
    if (i < j) {
      plot_traceplots(mcmc = mcmc, what = "corSigma", 
                      dimspec = c(i, j),
                      labcex = 0.6)
      abline(h = corr[i,j], col = "black", lty = 2)
    } else {
      plot(x = c(0, 1), y = c(0, 1), type = "n", xlab = "", ylab = "",
           xaxt = "n", yaxt = "n", bty = "n")
    }
  }
}
# or simply use
plot_traceplots_param(mcmc, "corSigma")

# plot_traceplots_param(mcmc, "Gplus", thin = 10)
# plot_traceplots_param(mcmc, "e0", thin = 10)
# plot_traceplots_param(mcmc, "ng", thin = 10)
# plot_traceplots_param(mcmc, "w", thin = 10)
# plot_traceplots_param(mcmc, "sdSigma", thin = 10)
# plot_traceplots_param(mcmc, "corSigma", thin = 10)
# plot_traceplots_param(mcmc, "Sigma", thin = 10)
# plot_traceplots_param(mcmc, "InvSigma", thin = 10)
# plot_traceplots_param(mcmc, "detInvSigma", thin = 10)
# plot_traceplots_param(mcmc, "tau_num", thin = 10)
# plot_traceplots_param(mcmc, "beta_num_fix", thin = 10)
# plot_traceplots_param(mcmc, "beta_num", thin = 10)
# plot_traceplots_param(mcmc, "beta_poi_fix", thin = 10)
# plot_traceplots_param(mcmc, "beta_poi", thin = 10)
# plot_traceplots_param(mcmc, "beta_bin_fix", thin = 10)
# plot_traceplots_param(mcmc, "beta_bin", thin = 10)
# plot_traceplots_param(mcmc, "beta_ord_fix", thin = 10)
# plot_traceplots_param(mcmc, "beta_ord", thin = 10)
# plot_traceplots_param(mcmc, "c_ord", thin = 10)
# plot_traceplots_param(mcmc, "beta_cat_fix", thin = 10)
# plot_traceplots_param(mcmc, "beta_cat", thin = 10)
# plot_traceplots_param(mcmc, "loglik", thin = 10)
# plot_traceplots_param(mcmc, "naY", thin = 10)

plot_kerneldensity_param(mcmc, "beta_num")
plot_ECDF_param(mcmc, "tau_num")
plot_ACF_param(mcmc, "w")

## Comparing different cluster-specific parameters
## We recommend to do it for chains separately
ch <- 1
TAB <- table(mcmc$last[[ch]]$U)
nonempty <- as.numeric(names(TAB)[which(TAB>0)]) 
# non-empty clusters decided based on last sampled values

neff <- mcmc$ngrp[[mcmc$Nums[1]]]
par(mfrow = c(length(mcmc$Nums), neff), mar = c(4, 4, 1, 1))

for (y in mcmc$Nums) {
  for (i in 1:neff) {
    plot_clusters(mcmc = mcmc, what = "beta_num", 
                  yspec = y, dimspec = i, whichg = nonempty,
                  setparmfrow = FALSE, ChainsToPlot = ch,
                  doECDF = FALSE)
    for (gg in 1:2) {
      abline(v = ldata$beta_num[[gg]][i], col = "black", lty = 2)
    }
  }
}
par(opar)
}
