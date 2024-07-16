\name{clustering_probabilities_and_deviance}
\alias{clustering_probabilities_and_deviance}
\title{
Evaluate clustering probabilities and deviance of clusGLMM output
}
\description{
To evaluate clustering probabilities all latent variables need to be integrated out of log-pdf contributions. Jointly modelled random effects for mixed-type outcomes cannot be integrated out in a closed formula. Hence, we utilize Laplace's approximation or even Adaptive Gaussian Quadrature (AGQ) in order to approximate these integrals. Samples of model parameters from \code{clusGLMM} are sent to C environment for faster evaluation.  
}
\usage{
clustering_probabilities_and_deviance(mcmc, id, data, dynamic_prob = FALSE, 
  start = 1, end = mcmc$iter, thin = 1, chains = 1:mcmc$nchains, 
  tolerance = mcmc$tuning$double$tolerance, 
  maxiter = mcmc$tuning$integer$maxiter, 
  maxnrep = mcmc$tuning$integer$maxnrep, 
  NGQ = 1)
}
\arguments{
  \item{mcmc}{
an output of \code{\link{clusGLMM}}.
}
  \item{id}{
a string, name of the column containing unit ID. Moreover, new units not used in MCMC sampling can be used here.  
}
  \item{data}{
matrix of outcomes, regressors including id-variable \code{id}.
}
  \item{dynamic_prob}{
a boolean; should the clustering probabilities be computed dynamically? That is, if unit \code{i} is observed \code{ni} times, then compute \code{ni} different probabilities calculated using the first \code{j} observations for \code{j in 1:ni}.
}
  \item{start}{
an iteration number where to start.
}
  \item{end}{
an iteration number where to end.
}
  \item{thin}{
a thinning that should be used between \code{start} and \code{end}; e.g. if \code{thin = 10} only every 10th iteration will be approximated. 
}
  \item{chains}{
chain number(s) for which the approximation should be computed. 
}
  \item{tolerance}{
tolerance when computing norm of the shift within Newton-Raphson.
}
  \item{maxiter}{
maximum iterations allowed during Newton-Raphson.
}
  \item{maxnrep}{
maximum number of repetitions of Newton-Raphson from different starting points.
}
  \item{NGQ}{
number of points for Adaptive Gaussian Quadrature; if \code{NGQ = 1} then simple Laplace's approximation is used. Do not use values higher than 7, since there are \code{NGQ^totnran} summands.
}
}
\details{
An \code{mcmc} output can be used regardless of the format \code{howsave = "list"} or \code{howsave = "matrix"}. Output will be of the same format.

Matrix of outcomes and regressors \code{data} does not have to be the same as for MCMC sampling. Even completely new units are allowed here. However, keep in mind that IDs in \code{id} are aligned with numbers from 1 to number of unique IDs within the current version. For simplicity, the names of output variables are indexed by these numbers and not the original ID. 

The indexing of probabilities gets more complicated when \code{dynamic_prob = TRUE} because for one unit there will be as many probabilities as many observations belong to this unit. If unit \code{i} is observed 3 times, there will be 3 probabilities: one where only the first observation is used, second where first two observations are used and third where all observations are used. With more available data the classification towards one of the clusters should be more prominent. In the end, there will be as many probabilities as rows of \code{data}. However, it does not mean that the rows correspond to each other, e.g. when \code{data} is not ordered by unit IDs. The reason is that we use a \code{for} cycle that goes through all units and then for a specific unit through all its observations. 

The code can handle unsorted rows of \code{data} by ID. However, for one specific unit the data are assumed to be ordered in time. This is internally used as ordering for the dynamic probability calculation. However, we recommend to insert ordered matrices so that pairing rows with the probabilities under \code{dynamic_prob = TRUE} is easier. 

Computations are usually expensive and division of the workload into several jobs where different chains or even subsamples are processed simultaneously. This can be achieved by setting \code{start}, \code{end} and \code{thin} differently for different jobs. However, you need to divide this on your own just as gathering the results back together. To do this properly keep in mind that the used iteration numbers are created by \code{seq(start, end, by = thin)}.

To approximate \code{totnran}-dimensional integral over random effects we use Adaptive Gaussian Quadrature. First, Newton-Raphson method is used to find the Taylor expansion of the log-pdf as a function of random effects. If \code{NGQ = 1} then the integral is approximated by evaluation only in the central point found by Newton-Raphson, which is just simple Laplace's approximation. If \code{NGQ = 2, 3, ..., 7} then \code{hermite.he.quadrature.rules(n = NGQ, normalized=FALSE)} from \code{library(gaussquad)} is used to find the roots and weights of hermite quadrature rules. These are combined for every dimension of random effects, which yields \code{NGQ^totnran} points at which the log-pdf needs to be evaluated. So, think carefully whether higher approximation precision with high \code{NGQ} is worth the expensive evaluation. 

}
\value{
  If \code{howsave = "list"}, then it is a list of outputs for different chains (from \code{chains}). Each chain then contains the following quantities:
  \describe{
    \item{\code{deviance}}{a scalar containing the model deviance,}
    \item{\code{dev_i}}{a vector of individual unit contribution to the model deviance,}
    \item{\code{pUig_int}}{a \code{n} by \code{G} matrix containing clustering probability towards one of \code{G} clusters.}
  }
  If \code{howsave = "matrix"}, then the output is a matrix where each row is one sample iteration and columns are named by the output variables above. 
}
\references{
Pinheiro, J.C., Chao, E.C. Efficient Laplacian and adaptive Gaussian quadrature algorithms for multilevel generalized linear mixed models. J. Comput. Graph. Stat. 15(1), 58–81 (2006). https://doi.org/10.1198/106186006x96962

Vávra, J., Komárek, A., Grün, B., Malsiner-Walli, G. Clusterwise multivariate regression of mixed-type panel data. Stat Comput 34, 46 (2024). https://doi.org/10.1007/s11222-023-10304-5
}
\author{
Jan Vávra
}
\seealso{
\code{\link{clusGLMM}}
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

Gmax <- 5     # maximal number of mixture components to be considered
nchains <- 2  # number of chains to be sampled 

## No inits specified
mcmc <- clusGLMM(formula = formula, id = id, family = family, data = data, 
                 G = Gmax, iter = 100, nchains = nchains)

###-----------------------------------------------------------------------------
###     Sample completely new data
###-----------------------------------------------------------------------------
nnew <- 20
set.seed(27818)
ldatanew <- generate_longitudinal_mixed_type_data(n = nnew, n_i = 4,
                     G = 2, timepar = "cross", catpar = "fixed", 
                     catsubjfix = TRUE,
                     xin1 = 1, Sigma = Sigma,
                     Kord = 5, Kcat = 4, kspec_bi_cat = FALSE)
newdata <- ldatanew$data
id <- "i"

probs <- clustering_probabilities_and_deviance(mcmc, id, newdata,
                                               dynamic_prob = FALSE,
                                               start = 1,
                                               end = 100,
                                               thin = 1,
                                               chains = 1:nchains,
                                               tolerance = 1e-7, 
                                               maxiter = 100,
                                               maxnrep = 20, 
                                               NGQ = 1)

###-----------------------------------------------------------------------------
###     Final clustering
###-----------------------------------------------------------------------------
probmean <- ps <- list()
clustering <- matrix(-1, nrow = as.numeric(nnew), ncol = mcmc$nchains)
certainty <- matrix(-1, nrow = as.numeric(nnew), ncol = mcmc$nchains)
for(ch in 1:mcmc$nchains){
  ps[[ch]] <- probs[probs$chain == ch, grep("pUig", colnames(probs))]
  psum <- summary(coda::mcmc(ps[[ch]]))
  probmean[[ch]] <- matrix(psum$statistics[,"Mean"], ncol = mcmc$G,
                           nrow = nnew, byrow = TRUE)
  
  clustering[,ch] <- apply(probmean[[ch]], 1, which.max)
  certainty[,ch] <- apply(probmean[[ch]], 1, max)
}
# To be classified in the cluster with maximal clustering probability,
# it has to overcome the given threshold, otherwise remains unclassified.
threshold <- 0.5
clustering[certainty < threshold] <- 0

ch = 1
probmean[[ch]]
rowSums(probmean[[ch]])
newdata$clustering <- clustering[,ch]

newdata1 <- newdata[newdata$j==1,]
table(newdata1$clustering, newdata1$g)
# row 0: unclassified

}

