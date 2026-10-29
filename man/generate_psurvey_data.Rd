\name{generate_psurvey_data}
\alias{generate_psurvey_data}
\title{
Generate Simulated Panel Survey Data
}
\description{
  A longitudinal dataset with explanatory variable age and dependent variables income (numeric), ntrips (count) and satisfaction (ordinal) is created. This function also allows to re-create the dataset \code{psurvey} included in the package.  
}
\usage{
generate_psurvey_data(n)
}
\arguments{
  \item{n}{
number of sample units.
}
}
\details{
The example below shows how the dataset \code{psurvey} included in the package has been generated.
}
\value{
A list containing:
\item{data }{a data frame with covariates, random effects, predictors and sampled outcomes,}
\item{sd_num }{a vector of length three containing the true values of the standard deviations for the numeric covariate,}
\item{beta_num }{a list of true values (named vector) of \eqn{\beta} parameters for numeric outcome for each cluster,}
\item{beta_poi }{a list of true values (named vector) of \eqn{\beta} parameters for count outcome for each cluster,}
\item{beta_ord }{a list of true values (named vector) of \eqn{\beta} parameters for ordinal outcome for each cluster,}
\item{c_ord }{a list of true values of intercepts specific to each level of the ordinal outcome for each cluster.}
}
\author{
Implemented by Jan Vávra to demonstrate the use of the package on artificial data in a vignette.
}
\examples{
### Generate psurvey data
set.seed(31415)
n <- 500
ldata <- generate_psurvey_data(n)
psurvey1 <- ldata$data[,c("i", "j", "age", "income", "ntrips", "satisfaction")]
psurvey_latent1 <- ldata$data[,c("i", "j", "g", "b_num", "b_poi", "b_ord")]

data("psurvey", package = "clustGLMM")
data("psurvey_latent", package = "clustGLMM")
all.equal(psurvey, psurvey1)
all.equal(psurvey_latent, psurvey_latent1)
}
