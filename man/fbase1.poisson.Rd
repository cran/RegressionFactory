\name{fbase1.poisson.log}
\alias{fbase1.poisson.log}

\title{
Single-Parameter Base Log-likelihood Function for Poisson GLM
}
\description{
Vectorized, single-parameter base log-likelihood functions for Poisson GLM using log link function. The base function(s) can be supplied to the expander function \code{\link{regfac.expand.1par}} in order to obtain the full, high-dimensional log-likleihood and its derivatives.
}
\usage{
fbase1.poisson.log(u, y, fgh=2)
}
\arguments{
  \item{u}{Varying parameter of the base log-likelihood function. This parameter is intended to be projected onto a high-dimensional space using the familiar regression transformation of \code{u <- X\%*\%beta}. In the typical use-case where the caller is \code{regfac.expand.1par}, a vector of values are supplied, and return objects will have the same length as \code{u}.}
  \item{y}{Fixed slot of the base distribution, corresponding to the response variable in the regression model. For \code{Poisson} family, it must be a vector of non-negative integers.}
  \item{fgh}{Integer with possible values 0,1,2. If \code{fgh=0}, the function only calculates and returns the log-likelihood vector and no derivatives. If \code{fgh=1}, it returns the log-likelihood and its first derivative in a list. If \code{fgh=2}, it returns the log-likelihood, as well as its first and second derivatives in a list.}
}

\value{
If \code{fgh==0}, the function returns \code{y*u-exp(u)-lfactorial(y)} for \code{log}. If \code{fgh==1}, a list is returned with elements \code{f} and \code{g}, where the latter is a vector of length \code{length(u)}, with each element being the first derivative of the above expressions. If \code{fgh==2}, the list will include an element named \code{h}, consisting of the second derivatives of \code{f} with respect to \code{u}.
}

\note{
In all base log-likelihood functions, we have dropped any additive terms that are independent of the distribution parameter, e.g. constant terms or those terms that are dependent on the response variable only. This is done for computational efficiency. Therefore, these functions cannot be used to obtain the absolute values of log-likelihood functions but only in the context of optimization and/or sampling. Users can write thin wrappers around these functions to add the constant terms to the function value. (Derivatives do not need correction. For Poisson family, the \code{lfactorial(y)} term is dropped.)
}

\author{
Alireza S. Mahani, Mansour T.A. Sharabiani
}

\seealso{
\code{\link{regfac.expand.1par}}
}

\examples{
\dontrun{
library(sns)
library(MfUSampler)

# using the expander framework and base distributions to define
# log-likelihood function for Poisson regression
loglike.poisson <- function(beta, X, y, fgh) {
  regfac.expand.1par(beta, X, y, fbase1.poisson.log, fgh)
}

# generate data for Poisson regression
N <- 1000
K <- 5
X <- matrix(runif(N*K, min=-0.5, max=+0.5), ncol=K)
beta <- runif(K, min=-0.5, max=+0.5)
y <- rpois(N, lambda = exp(X\%*\%beta))

# obtaining glm coefficients for comparison
beta.glm <- glm(y~X-1, family="poisson")$coefficients

# mcmc sampling of log-likelihood
nsmp <- 100

# Slice Sampler (no derivatives needed)
beta.smp <- array(NA, dim=c(nsmp,K)) 
beta.tmp <- rep(0,K)
for (n in 1:nsmp) {
  beta.tmp <- MfU.Sample(beta.tmp
    , f=loglike.poisson, X=X, y=y, fgh=0)
  beta.smp[n,] <- beta.tmp
}
beta.slice <- colMeans(beta.smp[(nsmp/2+1):nsmp,])

# Adaptive Rejection Sampler
# (only first derivative needed)
beta.smp <- array(NA, dim=c(nsmp,K)) 
beta.tmp <- rep(0,K)
for (n in 1:nsmp) {
  beta.tmp <- MfU.Sample(beta.tmp, uni.sampler="ars"
    , f=function(beta, X, y, grad) {
        if (grad)
          loglike.poisson(beta, X, y, fgh=1)$g
        else
          loglike.poisson(beta, X, y, fgh=0)
      }
    , X=X, y=y)
  beta.smp[n,] <- beta.tmp
}
beta.ars <- colMeans(beta.smp[(nsmp/2+1):nsmp,])

# SNS (Stochastic Newton Sampler)
# (both first and second derivative needed)
beta.smp <- array(NA, dim=c(nsmp,K)) 
beta.tmp <- rep(0,K)
for (n in 1:nsmp) {
  beta.tmp <- sns(beta.tmp, fghEval=loglike.poisson, X=X, y=y, fgh=2, rnd = n>nsmp/4)
  beta.smp[n,] <- beta.tmp
}
beta.sns <- colMeans(beta.smp[(nsmp/2+1):nsmp,])

# compare results
cbind(beta.glm, beta.slice, beta.ars, beta.sns)
}
}
