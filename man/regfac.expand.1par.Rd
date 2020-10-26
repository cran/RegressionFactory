\name{regfac.expand.1par}
\alias{regfac.expand.1par}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Expander Function for Single-Parameter Base Distributions
}
\description{
This function produces the full, high-dimensional gradient and Hessian from the base-distribution derivatives for linear transformations of the arguments of a single-parameter base distribution.
}
\usage{
regfac.expand.1par(beta, X, y, fbase1, fgh=2, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beta}{Vector of coefficients in the regression model.}
  \item{X}{Matrix of covariates in the regression model. Note that \code{ncol(X)} must be equal to \code{length(beta)}.}
  \item{y}{Vector of response variable in the regression model. Note that \code{length(y)} must be equal to \code{nrow(X)}.}
  \item{fbase1}{Base distribution function \code{fbase1(u, y, ...)} for the regression model. It must return a list with elements \code{f,g,h} corresponding to the function and its first and second derivatives relative to its first argument, \code{u}.}
  \item{fgh}{Integer with possible values 0,1,2. If \code{fgh=0}, the function only calculates and returns the log-likelihood function. If \code{fgh=1}, it returns the log-likelihood and its gradient vector. If \code{fgh=2}, it returns the log-likelihood, the gradient vector and the Hessian matrix.}
  \item{...}{Other parameters to be passed to \code{fbase1}.}
}

\value{
A list with elements \code{f,g,h} corresponding to the function, gradient vector, and Hessian matrix of the function \code{fbase(X\%*\%beta,y)}, i.e. the base function \code{fbase(u,y)} projected onto the high-dimensional space of \code{beta} through the linear transformation of its first argument (\code{u <- X\%*\%beta}).
}
\references{
Mahani, Alireza S. and Sharabiani, Mansour T.A. (2013)
\emph{Metropolis-Hastings Sampling Using Multivariate Gaussian Tangents}
\url{https://arxiv.org/pdf/1308.0657v1.pdf} 
}
\author{
Alireza S. Mahani, Mansour T.A. Sharabiani
}

\seealso{
\code{\link{regfac.expand.2par}}
}


\examples{
\dontrun{
library(sns)
# simulating logistic regression data
N <- 1000 # number of observations
K <- 10 # number of variables
X <- matrix(runif(N*K, min=-0.5, max=+0.5), ncol=K)
beta <- runif(K, min=-0.5, max=+0.5)
Xbeta <- X\%*\%beta
y <- 1*(runif(N)<1/(1+exp(-Xbeta)))
beta.est <- rep(0,K)
# run sns in non-stochastic mode, i.e. Newton-Raphson optimization
for (i in 1:10) {
  beta.est <- sns(beta.est, regfac.expand.1par, rnd=F, X=X, y=y
    , fbase1=fbase1.binomial.logit)
}
# use glm to estimate beta and compare
beta.est.glm <- glm(y~X-1, family="binomial")$coefficients
cbind(beta.est, beta.est.glm)
}
}

