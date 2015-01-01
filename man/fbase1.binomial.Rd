\name{fbase1.binomial.logit}
\alias{fbase1.binomial.logit}
\alias{fbase1.binomial.probit}
\alias{fbase1.binomial.cauchit}
\alias{fbase1.binomial.cloglog}

\title{
Single-Parameter Base Log-likelihood Function(s) for Binomial GLM
}
\description{
Vectorized, single-parameter base log-likelihood functions for binomial GLM using various link functions. These base functions can be supplied to the expander function \code{\link{regfac.expand.1par}} in order to obtain the full, high-dimensional log-likleihood and its derivatives.
}
\usage{
fbase1.binomial.logit(u, y, fgh=2, n=1)
fbase1.binomial.probit(u, y, fgh=2, n=1)
fbase1.binomial.cauchit(u, y, fgh=2, n=1)
fbase1.binomial.cloglog(u, y, fgh=2, n=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{u}{Varying parameter of the base log-likelihood function. This parameter is intended to be projected onto a high-dimensional space using the familiar regression transformation of \code{u <- X\%*\%beta}. In the typical use-case where the caller is \code{regfac.expand.1par}, a vector of values are supplied, and return objects will have the same length as \code{u}.}
  \item{y}{Fixed slot of the base distribution, corresponding to the response variable in the regression model. For \code{binomial} family, it must be an integer vector with values between \code{0} and \code{n}.}
  \item{fgh}{Integer with possible values 0,1,2. If \code{fgh=0}, the function only calculates and returns the log-likelihood vector and no derivatives. If \code{fgh=1}, it returns the log-likelihood and its first derivative in a list. If \code{fgh=2}, it returns the log-likelihood, as well as its first and second derivatives in a list.}
  \item{n}{Number of trials in the binomial model. This parameter is assumed to be fixed, and must be supplied by the user. If \code{n==1}, the model is reduced to binary logit/probit/cauchit/cloglog regression.}
}

\value{
If \code{fgh==0}, the \code{logit} version returns \code{-(n*log(1+exp(-u))+(n-y)*u)}, the \code{probit} returns \code{y*log(pnorm(u))+(n-y)*log(1-pnorm(u))}, the \code{cauchit} returns \code{y*log(pcauchy(u))+(n-y)*log(1-pcauchy(u))}, and the \code{cloglog} returns \code{y*log(1-exp(-exp(u)))-(n-y)*exp(u)}                                                                                                                                                                                                                                                                                             . If \code{fgh==1}, a list is returned with elements \code{f} and \code{g}, where the latter is a vector of length \code{length(u)}, with each element being the first derivative of the above expressions. If \code{fgh==2}, the list will include an element named \code{h}, consisting of the second derivatives of \code{f} with respect to \code{u}.
}

\note{
In all base log-likelihood functions, we have dropped any additive terms that are independent of the distribution parameter, e.g. constant terms or those terms that are dependent on the response variable only. This is done for computational efficiency. Therefore, these functions cannot be used to obtain the absolute values of log-likelihood functions but only in the context of optimization and/or sampling. Users can write thin wrappers around these functions to add the constant terms to the function value. (Derivatives do not need correction. For binomial family, all factorial terms are ignored since they only depend on \code{n} and {y}.)
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

# using the expander framework and binomial base log-likelihood
# to define log-likelihood function for binary logit regression
loglike.logit <- function(beta, X, y, fgh) {
  regfac.expand.1par(beta, X, y, fbase1.binomial.logit, fgh, n=1)
}

# generate data for logistic regression
N <- 1000
K <- 5
X <- matrix(runif(N*K, min=-0.5, max=+0.5), ncol=K)
beta <- runif(K, min=-0.5, max=+0.5)
y <- 1*(runif(N) < 1.0/(1+exp(-X\%*\%beta)))

# obtaining glm coefficients for comparison
beta.glm <- glm(y~X-1, family="binomial")$coefficients

# mcmc sampling of log-likelihood
nsmp <- 100

# Slice Sampler (no derivatives needed)
beta.smp <- array(NA, dim=c(nsmp,K)) 
beta.tmp <- rep(0,K)
for (n in 1:nsmp) {
  beta.tmp <- MfU.Sample(beta.tmp
    , f=function(beta, X, y) loglike.logit(beta, X, y, fgh=0), X=X, y=y)
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
          loglike.logit(beta, X, y, fgh=1)$g
        else
          loglike.logit(beta, X, y, fgh=0)
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
  beta.tmp <- sns(beta.tmp, fghEval=loglike.logit, X=X, y=y, fgh=2)
  beta.smp[n,] <- beta.tmp
}
beta.sns <- colMeans(beta.smp[(nsmp/2+1):nsmp,])

# compare results
cbind(beta.glm, beta.slice, beta.ars, beta.sns)
}
}
