\name{as_coda}
\alias{as_coda}
\alias{from_C_to_list}
\alias{from_C_to_matrix}
\alias{from_list_to_C}
\alias{from_list_to_coda}
\alias{from_list_to_matrix}
\alias{from_matrix_to_C}
\alias{from_matrix_to_coda}
\alias{from_matrix_to_list}
\title{
Change the output formats of \code{clustGLMM}
}
\description{
Functions are used internally in \code{clustGLMM} or \code{pUig_dev} to prepare data as an input of C function or to process data coming from C functions. Works for separate parameters, hence to transfer all the samples \code{for} cycles through chains and parameters are necessary. This is done in \code{from_list_to_matrix} and \code{from_matrix_to_list} which are able to switch between \code{howsave = "list"} and \code{howsave = "matrix"}.
}
\usage{
as_coda(mcmc, burnin=0, thin=1)
from_C_to_list(values, p, settings, yspecd1, yspecd2, family)
from_C_to_matrix(values, p, settings, yspecd1, yspecd2, family)
from_list_to_C(draws, iterations, p, settings, yspecd1, yspecd2, family)
from_list_to_coda(mcmc, burnin=0, thin=1)
from_list_to_matrix(mcmc)
from_matrix_to_C(draws, iterations, p, settings)
from_matrix_to_coda(mcmc, burnin=0, thin=1)
from_matrix_to_list(mcmc)
}
\arguments{
  \item{values}{
a vector of values to be stored.
}
  \item{p}{
a string containing the name of the parameter. 
}
  \item{settings}{
a matrix containing the information about model parameters.
}
  \item{yspecd1}{
a list of outcome-specific sizes along the first dimension.
}
  \item{yspecd2}{
a list of outcome-specific sizes along the second dimension.
}
  \item{family}{
a vector of distributional families named by the outcome name in \code{data}.
}
  \item{iterations}{
what iterations from a given chain to include.
}
  \item{draws}{
an element of \code{clustGLMM$draws} for a certain chain, either a list or a matrix.
}
  \item{mcmc}{
an output of \code{\link{clustGLMM}}.
}
  \item{burnin}{
the length of burn-in period (default 0), i.e. number of initial samples to ignore for inference.
}
  \item{thin}{
thinning rate (default 1); only every \code{thin}th value is used for inference (must be a positive integer). 
}
}
\details{
Most of the functions are used internally in \code{clustGLMM} or \code{clustering_probabilities_and_deviance} to process data coming from C functions. Functions that start or end in C only convert a specific parameter. The functions that change between \code{howsave = "list"} and \code{howsave = "matrix"} go through every chain and every saved parameter.
}
\value{
If \code{...to_list} then the output is a carefully structured list. Each chain \code{ch in 1:mcmc$chains} is saved in \code{mcmc$draws[[ch]]}. Then you select a parameter of your choice which may have even other sublists. For example, group-specific effects for numeric outcome would be accessed by
\code{mcmc$draws[[ch]]$beta_num[[g]][[mcmc$Nums[1]]]}
where \code{g} is the group.  

If \code{...to_matrix} then the output is again a list for each chain separately. But an item of that list is a matrix with carefully named columns. First column is \code{m} and contains the iteration number to which the row belongs to. Other columns are model parameters. For example, the first effect specific for group 2 for numeric outcome called "ynum" would be named \code{"beta_num_ynum(2)[1]"}. Note that the name of the variable has been pasted after the parameter name, then follows a group number in parentheses (only if parameter group-specific) and it ends with the dimensional coordinates in square brackets.

If \code{...to_C} then the output is a long vector where all the parameter values are stacked together in the order assumed in C functions. 

If \code{...to_coda} then the output is a \code{mcmc.list} object from \code{library(coda)}. \code{as_coda} is a wrapper that triggers \code{...to_coda} corresponding to \code{mcmc$howsave}. Then you can use environment for monitoring, investigating, plotting if you are familiar with it. However, a minor warning: the number of parameters is quite huge and can get overwhelming. Therefore, a set of plotting functions is provided with this package, see \code{\link{plot.clustglmm}}.
}
\references{
Vávra, J., Komárek, A., Grün, B., Malsiner-Walli, G. Clusterwise multivariate regression of mixed-type panel data. Stat Comput 34, 46 (2024). https://doi.org/10.1007/s11222-023-10304-5
}
\author{
Jan Vávra
}
\seealso{
\code{\link{clustGLMM}}, \code{\link{plot.clustglmm}}, \code{\link{clustering_probabilities_and_deviance}}
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
                       
###-----------------------------------------------------------------------------
###     Changing formats
###-----------------------------------------------------------------------------

## How to change from "matrix" to "list":
mcmc <- from_matrix_to_list(mcmc)

## How to change from "list" to "matrix":
mcmc <- from_list_to_matrix(mcmc)

## How to change from "matrix" to mcmc from 'coda' package
codamcmc <- from_matrix_to_coda(mcmc, burnin = 10, thin = 1)

## How to change from any howsave type to mcmc from 'coda' package
codamcmc <- as_coda(mcmc, burnin = 10, thin = 1)
}
