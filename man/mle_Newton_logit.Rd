% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/logit.R
\name{mle_Newton_logit}
\alias{mle_Newton_logit}
\title{MLE for linear coefficients}
\usage{
mle_Newton_logit(design, outcome, max_iter = 1000)
}
\arguments{
\item{design}{matrix, the design matrix X for the predictors.}

\item{outcome}{vector, the output Y for the response.}
}
\value{
Estimate of the linear coefficients
}
\description{
Calculate the MLE of linear coefficients using pseudo-inverse for the logit model
}
