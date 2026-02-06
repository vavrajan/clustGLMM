\name{longitudinal_mixed_type_data}
\alias{longitudinal_mixed_type_data}
\docType{data}
\title{
A Simulated Dataset of Longitudinal Mixed-Type Data
}
\description{
An artificially created dataset by \code{\link{generate_longitudinal_mixed_type}} 
to demonstrate the use of \code{\link{clustGLMM}}. 
For details about the parameter setting see the example at the bottom of \code{\link{generate_longitudinal_mixed_type}}.
}
\usage{data("students")}
\format{
  A data frame with 400 rows and 32 columns.
  \describe{
    \item{\code{i}}{an id variable}
    \item{\code{x}}{continuous regressor}
    \item{\code{j}}{order of observation for given id}
    \item{\code{bs1},...,\code{bs5}}{B-spline basis for \code{x}}
    \item{\code{f}}{factor regressor}
    \item{\code{g}}{a true group label from 1, 2, 3}
    \item{\code{b_type}}{random intercept for the given \code{type}}
    \item{\code{pred_type}}{predictor for the given \code{type}, categorical has multiple predictors for different categories}
    \item{\code{c_ordL}}{ordered intercept for level \code{L}}
    \item{\code{ytype}}{an outcome of the given \code{type}}
  }
}
\examples{
data(longitudinal_mixed_type_data)
summary(longitudinal_mixed_type_data)
}
\keyword{longitudinal_mixed_type_data}
