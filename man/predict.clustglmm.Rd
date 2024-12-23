\name{predict.clustglmm}
\alias{predict.clustglmm}
\title{
Estimate predictors from formulae
}
\description{
Predicts the linear combinations of fixed and group part of the model formulae of \code{clustGLMM} output.
}
\usage{
\method{predict}{clustglmm}(object, newdata, y = names(object$family)[1],
    k = NA, type = c("lp", "response"), method = c("samples", "plugin"),
    what = c("fg", "f", "g"), fun_posterior = mean, level = 0.95,
    interval = c("none", "ET", "HPD"), chain = 1, burnin = 0, thin = 1,
    gs = object$clusters[[chain]], ...)
}
\arguments{
  \item{object}{
    an object of class \code{"clustglmm"}, an output of \code{clustGLMM} function.
}
  \item{newdata}{
    a \code{data.frame} containing all columns needed to construct model matrices.
}
  \item{y}{
    a string; outcome to be predicted.
}
  \item{k}{
    an integer; level of ordinal or categorical outcome to be predicted.
}
  \item{type}{
    what should be predicted; one of \code{"lp"} (default) or \code{"response"}
    (not implemented yet), see Details below.
}
  \item{method}{
    the method of predicting; one of \code{"samples"} (default) or \code{"plugin"}, see Details below.
}
  \item{what}{
    more detailed clarification to the selected \code{type}.
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
    chain number to be used to do the prediction, the default is the first \code{chain = 1}.
}
  \item{burnin}{
    the length of burn-in period (default 0), i.e., number of initial samples to ignore for inference.
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
  \item{lp}{reconstructs the linear predictor for given outcome
    \code{y}. You can specify by \code{what}: "f", "g", "fg" which part
    of the predictor should be reconstructed: the fixed part invariant
    to all clusters, the group-specific part for chosen clusters
    \code{gs} or sum of both (default).}
  \item{response}{not yet implemented.}
}

As we approximate the posterior distribution with samples, there may be
several ways on obtaining the posterior estimate of desired quantity
(\code{type}). Depending on \code{method} we proceed in the following
way:
\describe{
  \item{samples}{treat the desired quantity as a parametric function,
    evaluate it with all (closely specified) samples in
    \code{object$draws} and estimate it by \code{fun_posterior}. This
    option allows for construction of credible intervals.}
  \item{plugin}{estimate the model parameters with \code{fun_posterior}
    first and then plug them into the function defining the
    quantity. Does not provide credible intervals (yet?).}
}

The computation of credible intervals is given by \code{interval}:
\describe{
  \item{none}{no lower and upper bound desired.}
  \item{ET}{equal-tailed credible interval determined out of samples
    from posterior distribution with \code{quantile} function with
    \code{probs = c((1-level)/2, 1 - (1-level)/2}.}
  \item{HPD}{highest-posterior density credible interval determined out
    of samples from posterior distribution with \code{\link{hdi}}
    function from \pkg{HDInterval}.}
}
}
\value{
For given outcome \code{y} returns a list with
\describe{
  \item{$fit}{the estimate of the posterior location.}
  \item{$lwr}{the lower bound of credible interval.}
  \item{$lwr}{the upper bound of credible interval.}
}
Depending on cluster-specificity of the task, the results are either one
column matrix or a matrix with each column corresponding to one of the
required clusters named by \code{paste0("(", gs, ")")}.
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
predict(mcmc, newdata = data, type = "lp")

# predict the evolution in time for covariate f=0
newdata <- data.frame("x" = seq(0, 1, length.out=100), "f"=factor(0, levels=c(0,1)))

# numeric outcome
predict(mcmc, newdata, y="ynum", type="lp", what="fg")

# categorical outcome
predict(mcmc, newdata, y="ycat", k=1, type="lp", what="fg")

}
