\name{fbase1.geometric.logit}
\alias{fbase1.geometric.logit}

\title{
Single-Parameter Base Log-likelihood Function for Exponential GLM
}
\description{
Vectorized, single-parameter base log-likelihood functions for geometric GLM using logit link function. The base function(s) can be supplied to the expander function \code{\link{regfac.expand.1par}} in order to obtain the full, high-dimensional log-likleihood and its derivatives.
}
\usage{
fbase1.geometric.logit(u, y, fgh=2)
}
\arguments{
  \item{u}{Varying parameter of the base log-likelihood function. This parameter is intended to be projected onto a high-dimensional space using the familiar regression transformation of \code{u <- X\%*\%beta}. In the typical use-case where the caller is \code{regfac.expand.1par}, a vector of values are supplied, and return objects will have the same length as \code{u}.}
  \item{y}{Fixed slot of the base distribution, corresponding to the response variable in the regression model. For \code{Geometric} family, it must be a vector of non-negative integers.}
  \item{fgh}{Integer with possible values 0,1,2. If \code{fgh=0}, the function only calculates and returns the log-likelihood vector and no derivatives. If \code{fgh=1}, it returns the log-likelihood and its first derivative in a list. If \code{fgh=2}, it returns the log-likelihood, as well as its first and second derivatives in a list.}
}

\value{
If \code{fgh==0}, the function returns \code{-(y*u+(1+y)*log(1+exp(-u)))} for \code{log}. If \code{fgh==1}, a list is returned with elements \code{f} and \code{g}, where the latter is a vector of length \code{length(u)}, with each element being the first derivative of the above expressions. If \code{fgh==2}, the list will include an element named \code{h}, consisting of the second derivatives of \code{f} with respect to \code{u}.
}

\note{
The logit function must be applied to the probability parameter to give \code{X\%*\%beta}, which is in turn the inverse of the mean of the geometric distribution. For brevity, we still call the link function 'logit'.
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
# log-likelihood function for geometric regression
loglike.geometric <- function(beta, X, y, fgh) {
  regfac.expand.1par(beta, X, y, fbase1.geometric.logit, fgh)
}

# generate data for geometric regression
N <- 1000
K <- 5
X <- matrix(runif(N*K, min=-0.5, max=+0.5), ncol=K)
beta <- runif(K, min=-0.5, max=+0.5)
y <- rgeom(N, prob = 1/(1+exp(-X\%*\%beta)))

# mcmc sampling of log-likelihood
nsmp <- 100

# Slice Sampler
beta.smp <- array(NA, dim=c(nsmp,K)) 
beta.tmp <- rep(0,K)
for (n in 1:nsmp) {
  beta.tmp <- MfU.Sample(beta.tmp
    , f=loglike.geometric, X=X, y=y, fgh=0)
  beta.smp[n,] <- beta.tmp
}
beta.slice <- colMeans(beta.smp[(nsmp/2+1):nsmp,])

# Adaptive Rejection Sampler
beta.smp <- array(NA, dim=c(nsmp,K)) 
beta.tmp <- rep(0,K)
for (n in 1:nsmp) {
  beta.tmp <- MfU.Sample(beta.tmp, uni.sampler="ars"
   , f=function(beta, X, y, grad) {
     if (grad)
       loglike.geometric(beta, X, y, fgh=1)$g
     else
       loglike.geometric(beta, X, y, fgh=0)
   }
   , X=X, y=y)
  beta.smp[n,] <- beta.tmp
}
beta.ars <- colMeans(beta.smp[(nsmp/2+1):nsmp,])

# SNS (Stochastic Newton Sampler)
beta.smp <- array(NA, dim=c(nsmp,K)) 
beta.tmp <- rep(0,K)
for (n in 1:nsmp) {
  beta.tmp <- sns(beta.tmp, fghEval=loglike.geometric, X=X, y=y, fgh=2, rnd = n>nsmp/4)
  beta.smp[n,] <- beta.tmp
}
beta.sns <- colMeans(beta.smp[(nsmp/2+1):nsmp,])

# compare sample averages with actual values
cbind(beta, beta.sns, beta.slice, beta.ars)
}
}
