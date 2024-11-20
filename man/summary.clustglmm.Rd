\name{summary.clustglmm}
\alias{summary.clustglmm}
\alias{print.summary.clustglmm}
\title{
Computing the summary statistics of outputs of \code{clustGLMM}
}
\description{
A method for computing and printing the summary statistics of class \code{clustglmm}, the outputs of \code{clustGLMM} function.
}
\usage{
summary(object, ...)
print(x, chains=c(1), which = c("mcmc", "inv_param", "clust_param"), ...)
}
\arguments{
  \item{object}{
an object of class \code{"clustglmm"}, an output of \code{clustGLMM} function.
}
  \item{x}{
an object of class \code{"summary.clustglmm"}, an output of \code{summary.clustGLMM} function.
}
  \item{chains}{
a vector of chain numbers to print summary results of.
}
  \item{...}{
other arguments (such as \code{burnin} or \code{thin} for \code{\link{as_coda}}).
}
}
\details{
  \code{print.summary} prints a summary for each selected chain in \code{chains}. First, a brief summary of clustering is provided. The raw output of \code{clustGLMM} is based on sampled allocation indicators, output of \code{clustering_probabilities_and_deviance} provides summary of clustering based on approximated clustering probabilities. Then, the summary statistics for the key model parameters are provided in several blocks in \code{library(coda)} style. The output of \code{summary} contains both tables (\code{"statistics"} or \code{"quantiles"}) but only a combination of these two is printed (all quantiles and Mean). 
}
\value{
  \code{summary} returns a list of:
  \item{call}{a string describing the model defined to be estimated,}
  \item{G}{maximal number of underlying clusters,}
  \item{modeGplus}{the most frequent number of non-empty clusters per each chain,}
  \item{clustering}{matrix of assigned clusters for each unit based on sampled indicators,}
  \item{certainty}{matrix of certainties (the highest clustering probability among G other) for each unit based on sampled indicators,}
  \item{clusters}{labels of non-empty clusters for each chain separately,}
  \item{coda_summary}{output of \code{coda} summary for each chain separately,}
  \item{output}{the type of output that is summarized: \code{mcmc} samples from \code{clustGLMM}, \code{probs} approximation of clustering probabilities and deviance,}
  \item{fnames}{vector of parameter names invariant to clustering,}
  \item{gnames}{cluster-specific parameter names only for cluster 1,}
  \item{lnames}{latent variable names (allocation indicators, random effects, missing values),}
  \item{cnames}{clustering related parameter names.}
}
\references{
Vávra, J., Komárek, A., Grün, B., Malsiner-Walli, G. Clusterwise multivariate regression of mixed-type panel data. Stat Comput 34, 46 (2024). https://doi.org/10.1007/s11222-023-10304-5
}
\author{
Jan Vávra
}
\seealso{
\code{\link{clustGLMM}}, \code{\link{clustering_probabilities_and_deviance}}, \code{\link{print.summary.clustglmm}}, \code{\link{summary.mcmc}}
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
###     Prepare inputs for clustGLMM
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
###     Final call of clustGLMM
###-----------------------------------------------------------------------------

iter <- 100    # the length of sampled chain
Gmax <- 5      # maximal number of mixture components to be considered
nchains <- 2   # number of chains to be sampled 

## No inits specified
## Saved in "matrix" format
mcmc <- clustGLMM(formula = formula, id = id, family = family, data = data, 
                  G = Gmax, iter = iter, nchains = nchains,
                  howsave = "matrix")

mcmc
s <- summary(mcmc)
print(s)
print(s, which=c("mcmc", "clustering", "inv_param", "clust_param", "latent_param", "g_param"))
}
