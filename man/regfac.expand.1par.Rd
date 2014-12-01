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
regfac.expand.1par(beta, X, y, fbase)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beta}{Vector of coefficients in the regression model.}
  \item{X}{Matrix of covariates in the regression model. Note that \code{ncol(X)} must be equal to \code{length(beta)}.}
  \item{y}{Vector of response variable in the regression model. Note that \code{length(y)} must be equal to \code{nrow(X)}.}
  \item{fbase}{Base distribution function \code{fbase(u,y)} for the regression model. It must return a list with elements \code{f,g,h} corresponding to the function and its first and second derivatives relative to its first argument, \code{u}.}
}

\value{
A list with elements \code{f,g,h} corresponding to the function, gradient vector, and Hessian matrix of the function \code{fbase(X\%*\%beta,y)}, i.e. the base function \code{fbase(u,y)} projected onto the high-dimensional space of \code{beta} through the linear transformation of its first argument (\code{u <- X\%*\%beta}).
}
\references{
Mahani, Alireza S. and Sharabiani, Mansour T.A. (2013)
\emph{Metropolis-Hastings Sampling Using Multivariate Gaussian Tangents}
\url{http://arxiv.org/pdf/1308.0657v1.pdf} 
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
    , fbase=fbase.bern.ilogit)
}
# use glm to estimate beta and compare
beta.est.glm <- glm(y~X-1, family="binomial")$coefficients
cbind(beta.est, beta.est.glm)
}
}

