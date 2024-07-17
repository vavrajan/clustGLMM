\name{predict.clusglmm}
\alias{predict.clusglmm}
\title{
Estimate predictors from formulae
}
\description{
Predicts the linear combinations of fixed and group part of the model formulae of \code{clusGLMM} output.
}
\usage{
predict(object, newdata, 
        type=c("lp", "lp-f", "lp-g", "lp-fg"), how=c("samples", "plugin"),
        fun_posterior=mean, level=0.95, interval=c("none", "ET", "HPD"),
        chain=1, burnin=0, thin=1, gs=object$clusters[[chain]], ...)
}
\arguments{
  \item{object}{
an object of class \code{"clusglmm"}, an output of \code{clusGLMM} function.
}
  \item{newdata}{
a \code{data.frame} containing all columns needed to construct model matrices.
}
  \item{type}{
what should be predicted; one of \code{"lp-f"}, \code{"lp-g"}, \code{"lp"="lp-fg"} (default), see Details below.
}
  \item{how}{
the method of predicting; one of \code{"samples"} (default) or \code{"plugin"}, see Details below.
}
  \item{fun_posterior}{
a function summarizing samples into an estimator of location parameter of posterior distribution; e.g. \code{mean} or \code{median}.
}
  \item{level}{
the credibility level of credible intervals; 0.95 by default.
}
  \item{interval}{
the method of constructing credible intervals; one of \code{"none"}, \code{"ET"} (equal-tailed by default) or \code{"HPD"} (highest-posterior density), see Details below.
}
  \item{chain}{
chain number to be used to do the prediction, the default is the first \code{chain=1}.
}
  \item{burnin}{
the length of burn-in period (default 0), i.e. number of initial samples to ignore for inference.
}
  \item{thin}{
thinning rate (default 1); only every \code{thin}th value is used for inference (must be a positive integer). 
}
  \item{gs}{
cluster labels, for which prediction of cluster-specific should be made; all non-empty clusters by default.
}
  \item{...}{
other arguments.
}
}
\details{
Depending on \code{type} the following is computed:
\describe{
  \item{lp}{reconstructs the linear predictor for each of the considered outcomes. You can specify by adding "-f", "-g", "-fg" which part of the predictor should be reconstructed: the fixed part invariant to all clusters, the group-specific part for chosen clusters \code{gs} or sum of both.}
  \item{...}{other to be implemented?!?}
}

As we approximate the posterior distribution with samples, there may be several ways on obtaining the posterior estimate of desired quantity (\code{type}). Depending on \code{how} we proceed in the following way:
\describe{
  \item{samples}{treat the desired quantity as a parametric function, evaluate it with all (closely specified) samples in \code{object$draws} and estimate it by \code{fun_posterior}. This option allows for construction of credible intervals.}
  \item{plugin}{estimate the model parameters with \code{fun_posterior} first and then plug them into the function defining the quantity. Does not provide credible intervals (yet?).}
}

The computation of credible intervals is given by \code{interval}:
\describe{
  \item{none}{no lower and upper bound desired.}
  \item{ET}{equal-tailed credible interval determined out of samples from posterior distribution with \code{quantile} function with \code{probs = c((1-level)/2, 1 - (1-level)/2}.}
  \item{HPD}{highest-posterior density credible interval determined out of samples from posterior distribution with \code{\link{hdi}} function from \code{library(HDInterval)}.}
}
}
\value{
Returns a list named by all considered outcomes. Each contains another list with
\describe{
  \item{$fit}{the estimate of the posterior location.}
  \item{$lwr}{the lower bound of credible interval.}
  \item{$lwr}{the upper bound of credible interval.}
}
Depending on cluster-specificity of the task, the results are either one column matrix or a matrix with each column corresponding to one of the required clusters named by \code{paste0("(", gs, ")")}. 
}
\references{
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

iter <- 100    # the length of sampled chain
Gmax <- 5      # maximal number of mixture components to be considered
nchains <- 2   # number of chains to be sampled 

## No inits specified
## Saved in "matrix" format
mcmc <- clusGLMM(formula = formula, id = id, family = family, data = data, 
                 G = Gmax, iter = iter, nchains = nchains,
                 howsave = "matrix")
                 
predict(mcmc, newdata = data, type = "lp")

# predict the evolution in time for covariate f=0
newdata0 <- data.frame("x" = seq(0, 1, length.out=100), "f"=0)
predict(mcmc, newdata = newdata0, type = "lp")

# predict the evolution in time for covariate f=1
newdata1 <- data.frame("x" = seq(0, 1, length.out=100), "f"=1)
predict(mcmc, newdata = newdata1, type = "lp")

}
