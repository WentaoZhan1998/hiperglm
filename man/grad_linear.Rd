% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/linear.R
\name{grad_linear}
\alias{grad_linear}
\title{Gradient for the likelihood under linear model}
\usage{
grad_linear(design, outcome, beta, noise_var)
}
\arguments{
\item{design}{matrix, the design matrix X for the predictors.}

\item{outcome}{vector, the output Y for the response.}

\item{beta}{vector, the linear coefficients}

\item{noise_var}{number, the variance of the noise}
}
\value{
Gradient for the likelihood under linear model
}
\description{
Calculate the gradient for the linear model
}
