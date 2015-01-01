\name{regfac.merge}
\alias{regfac.merge}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Utility Function for Adding Two Functions and Their Derivatives
}
\description{
Combining two log-density functions by adding the corresponding elements of their lists (function, gradient, Hessian). This can be useful, e.g. in combining the likelihood and the prior (in log domain) to form the posterior according to Bayes rule.
}
\usage{
regfac.merge(fgh1, fgh2, fgh = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fgh1}{First log-density list, containing elements \code{f}, \code{g} and \code{h}, corresponding to log-density function, its gradient vector, and its Hessian matrix.}
  \item{fgh2}{Second log-density list, containing elements \code{f}, \code{g} and \code{h}, corresponding to log-density function, its gradient vector, and its Hessian matrix.}
  \item{fgh}{Integer flag with possible values \code{0,1,2}, indicating the maximum order of derivative to be returned.}
}

\value{
If \code{fgh==0}, \code{fgh1+fgh2} is returned. Otherwise, a list is returned with elements \code{f}, \code{g}, and \code{h}, each of which is the sum of corresponding elements of \code{fgh1} and \code{fgh2} lists.
}

\author{
Alireza S. Mahani, Mansour T.A. Sharabiani
}

\examples{
# constructing the log-posterior for Bayesian logistic regression
# log-likelihood
loglike.logistic <- function(beta, X, y, fgh) {
  regfac.expand.1par(beta, X, y, fbase1.binomial.logit, fgh, n=1)
}
# log-prior
logprior.logistic <- function(beta, mu.beta, sd.beta, fgh) {
  f <- sum(dnorm(beta, mu.beta, sd.beta, log=TRUE))
  if (fgh==0) return (f)
  g <- -(beta-mu.beta)/sd.beta^2
  if (fgh==1) return (list(f=f, g=g))
  #h <- diag(rep(-1/sd.beta^2,length(beta)))
  h <- diag(-1/sd.beta^2)
  return (list(f=f, g=g, h=h))
}
# adding log-likelihood and log-prior according to Bayes rule
logpost.logistic <- function(beta, X, y, mu.beta, sd.beta, fgh) {
  ret.loglike <- loglike.logistic(beta, X, y, fgh)
  ret.logprior <- logprior.logistic(beta, mu.beta, sd.beta, fgh)
  regfac.merge(ret.loglike,ret.logprior, fgh=fgh)
}

}