\name{get_Sigma}
\alias{get_Sigma}
\alias{get_sd_num}
\title{
Extract the m-th Sample of Sigma or sd_num
}
\description{
These two parameters are saved under different formats or tranformation.
These functions return the m-th sample of the corresponding parameter even if it has to be transformed from another equivalent saved parameter.
}
\usage{
get_Sigma(draws, m, settings, g = NA, howsave = "matrix")  
get_sd_num(draws, m, settings, y = NA, howsave = "matrix") 
}
\arguments{
  \item{draws}{
samples of model parameters in either list or matrix structure.
}
  \item{m}{
an iteration to choose.
}
  \item{settings}{
a matrix of model settings (in output of \code{clustglmm}).
}
  \item{g}{
a cluster label.
}
  \item{y}{
a name of an outcome.
}
  \item{howsave}{
the format under which \code{draws} is saved.
}
}
\details{
From \code{draws}, only the \code{m}-th element is used, which slightly differs depending on \code{howsave}.
Functions are able to detect group-specificity and dimension from \code{settings}.

If \code{Sigma} is group-specific, then only matrix for gorup \code{g} is returned.
If \code{Sigma} is not saved among draws, then it is reconstructed from \code{InvSigma} and inversion.
If neither of \code{Sigma} and \code{InvSigma} is saved among draws, then it is reconstructed from standard deviations \code{sdSigma} or \code{corSigma}.
If none of the above, a warning is triggered and diagonal unit matrix is returned. 

On the other hand \code{get_sd_num} returns standard deviations for all clusters if cluster-specific. 
First, it tries to get it from \code{sd_num}, if saved.
If not it transforms either from \code{prec_num} or \code{var_num}. 
If none of the above, error is triggered. 
}
\value{
A single covariance matrix \code{Sigma} or standard deviation(s) of numeric outcome \code{y}.
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

iter <- 100    # length of the sampled chain
Gmax <- 5      # maximal number of mixture components to be considered
nchains <- 2   # number of chains to be sampled 

mcmc <- clustGLMM(formula = formula, id = id, family = family, data = data, 
                  G = Gmax, iter = iter, nchains = nchains)

Sigma <- get_Sigma(mcmc$draws[[1]], 1, mcmc$settings[[1]],
                   1, howsave = mcmc$howsave)
sd_num <- get_sd_num(mcmc$draws[[1]], 1, mcmc$settings[[1]],
                     "ynum", howsave = mcmc$howsave)                  

}

