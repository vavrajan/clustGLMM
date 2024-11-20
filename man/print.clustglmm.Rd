\name{print.clustglmm}
\alias{print.clustglmm}
\title{
Printing the outputs of \code{clustGLMM}
}
\description{
A method for printing objects of class \code{clustglmm}, the outputs of \code{clustGLMM}.
}
\usage{
print(x, which = "call", ...)
}
\arguments{
  \item{x}{
an object of class \code{"clustglmm"}, an output of \code{clustGLMM} function.
}
  \item{which}{
a vector containing the following strings 
\describe{
\item{\code{"call"}}{to print the overview of the estimated model,}
\item{\code{"clustering"}}{to print the brief summary of the clustering.}
}
}
  \item{...}{
other arguments.
}
}
\details{
  First, if \code{"call"} is in \code{which}, it prints \code{x$call} that summarizes 
  what has been fitted including sample sizes, number of units, observations per unit, 
  chosen maximal number of clusters, overview of outcomes and selected types, 
  fixed-effects common to all groups, random-effects with the structure of \code{Sigma} matrix, 
  group-specific formulae for outcomes, list of other group-specific parameters and summary of the MCMC sampling. 

  Next, if \code{"clustering"} is in \code{which}, 
  it summarizes the clustering based on sampled indicators \code{U} for each chain separately. 
  First, \code{modeGplus}, the most frequent value of \code{Gplus} - number of non-empty components, is determined. 
  Then, all samples of cluster indicators \code{U} are gathered from all units to determine \code{modeGplus} most frequent cluster labels, stored in \code{clusters}. 
  To each unit \code{i} we assign the cluster label from \code{clusters} that maximizes the frequency of samples of \code{U[i]}, stored in \code{clustering}. 
  We also keep the relative frequency to the assigned cluster, stored in \code{certainty}. 
  By choosing reasonable \code{threshold}, e.g. \code{0.5}, you can make a unit unclassified (0) if the \code{certainty} is not high enough: \code{clustering[certainty < threshold] <- 0}. 
  If \code{G = 1}, then clustering summary is irrelevant. 
}
\value{
  No returned value, only prints.
}
\references{
Vávra, J., Komárek, A., Grün, B., Malsiner-Walli, G. Clusterwise multivariate regression of mixed-type panel data. Stat Comput 34, 46 (2024). https://doi.org/10.1007/s11222-023-10304-5
}
\author{
Jan Vávra
}
\seealso{
\code{\link{clustGLMM}}
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
print(mcmc, which = "clustering")
print(mcmc, which = c("call", "clustering"))

}
