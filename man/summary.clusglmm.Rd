\name{summary.clusglmm}
\alias{summary.clusglmm}
\title{
Computing the summary statistics of outputs of \code{clusGLMM}
}
\description{
A method for computing the summary statistics of class \code{clusglmm}, the outputs of \code{clusGLMM}.
}
\usage{
summary(object, ...)
}
\arguments{
  \item{object}{
an object of class \code{"clusglmm"}, an output of \code{clusGLMM} function.
}
  \item{...}{
other arguments.
}
}
\details{
  TODO describe what and how is computed for the summary.
}
\value{
  Returns a list of:
  \item{call}{A string describing the model defined to be estimated.}
  \item{clustering}{Summary of the clustering based on sampled indicators.}
}
\references{
Vávra, J., Komárek, A., Grün, B., Malsiner-Walli, G. Clusterwise multivariate regression of mixed-type panel data. Stat Comput 34, 46 (2024). https://doi.org/10.1007/s11222-023-10304-5
}
\author{
Jan Vávra
}
\seealso{
\code{\link{clusGLMM}}, \code{\link{summary.clusglmm}}
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
chains <- 2    # number of chains to be sampled 

## No inits specified
## Saved in "matrix" format
mcmc <- clusGLMM(formula = formula, id = id, family = family, data = data, 
                 G = Gmax, iter = iter, chains = chains,
                 howsave = "matrix")

mcmc
s <- summary.clusglmm(mcmc)
s
}
