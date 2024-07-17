\name{predict.clusglmm}
\alias{predict.clusglmm}
\title{
Estimate predictors from formulae
}
\description{
Predicts the linear combinations of fixed and group part of the model formulae of \code{clusGLMM} output.
}
\usage{
predict(object, newdata, type="f", posterior="E", chain=1, gs=object$clusters[[chain]], ...)
}
\arguments{
  \item{object}{
an object of class \code{"clusglmm"}, an output of \code{clusGLMM} function.
}
  \item{newdata}{
a \code{data.frame} containing all columns needed to construct model matrices.
}
  \item{type}{
one of \code{"f"}, \code{"g"}, \code{"fg"}, \code{"as_parametric_function"}, see Details below.
}
  \item{posterior}{
one of \code{"E"} \code{"median"} to choose the summary statistic for approximation of posterior mean.
}
  \item{chain}{
chain number to be used to do the prediction.
}
  \item{gs}{
cluster labels for which prediction of cluster-specific should be made.
}
  \item{...}{
other arguments.
}
}
\details{
Depending on \code{type} the following is computed:
\describe{
  \item{f}{computes only the fixed part of the predictor by using the estimated posterior mean.}
  \item{g}{computes only the group-specific part of the predictor by using the estimated posterior mean.}
  \item{fg}{computes the sum of fixed and group-specific part of the predictor by using the estimated posterior mean.}
  \item{as_parametric_function}{computes the sum of fixed and group-specific part with all samples of regression coefficients and then takes the corresponding posterior mean estimate.}
}
}
\value{
  Returns a list named by all considered outcomes. 
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
                 
predict(mcmc, newdata = data)

# predict the evolution in time for covariate f=0
newdata0 <- data.frame("x" = seq(0, 1, length.out=100), "f"=0)
predict(mcmc, newdata = newdata0)

# predict the evolution in time for covariate f=1
newdata1 <- data.frame("x" = seq(0, 1, length.out=100), "f"=1)
predict(mcmc, newdata = newdata1)

}
