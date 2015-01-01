\name{regfac.expand.2par}
\alias{regfac.expand.2par}

\title{
Expander Function for Two-Parameter Base Distributions
}
\description{
This function produces the full, high-dimensional gradient and Hessian from the base-distribution derivatives for linear transformations of the arguments of a two-parameter base distribution.
}
\usage{
regfac.expand.2par(coeff, X, Z=matrix(1.0, nrow=nrow(X), ncol=1)
  , y, fbase2, fgh=2, block.diag=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{coeff}{Vector of coefficients in the regression model. The first \code{ncol(X)} elements correspond to the first parameter of the base distribution \code{fbase2(u, v, y, ...)}, and the next \code{ncol(Z)} elements corresponds to the second parameter of the base distribution \code{fbase2(u, v, y, ...)}.}
  \item{X}{Matrix of covariates corresponding to the first parameter of the base distribution \code{fbase2(u, v, y, ...)}.}
  \item{Z}{Matrix of covariates corresponding to the second parameter of the base distribution \code{fbase2(u, v, y, ...)}. Default is a single column of 1's, corresponding to an intercept-only model for the second parameter, i.e. assuming the second parameter is constant across all observations. Note that \code{nrow(Z)} must be equal to \code{nrow(X)}.}
  \item{y}{Vector of response variables. Note that \code{length(y)} must be equal to \code{nrow(X)}.}
  \item{fbase2}{Base distribution function \code{fbase2(u, v, y, ...)} for the regression model. It must return a list with elements \code{f,g,h} corresponding to the function and its first and second derivatives relative to its first two argument, \code{u,v}. The gradient must be a matrix of dimensions \code{nrow(X)}-by-2, where the first column is the gradient of the log-likelihood function with respect to its first parameter (fbase2_u), evaluated at each of the \code{nrow(X)} observations, and the second column is the gradient of the log-likelihood function with repsect to its second parameter (fbase2_v), also evaluated at each observation point. Similarly, the Hessian must be a matrix of dimensions \code{nrow(X)}-by-3, with elements being equal to fbase2_uu, fbase2_vv and fbase2_uv evaluated at each observation point (taking advantage of the Hessian being symmetric).}
  \item{fgh}{Integer with possible values 0,1,2. If \code{fgh=0}, the function only calculates and returns the log-likelihood function. If \code{fgh=1}, it returns the log-likelihood and its gradient vector. If \code{fgh=2}, it returns the log-likelihood, the gradient vector and the Hessian matrix.}
  \item{block.diag}{If \code{TRUE}, Hessian matrix is block-diagonalized by setting cross-terms between \code{beta} and \code{gamma} to zero. This can be useful if the full - i.e. non-block-diagonalized - Hessian is not negative definite, but block-diagonalization leads to definiteness. If \code{TRUE}, third element of the Hessian of \code{fbase} is not needed and thus it can be vector of length 2 instead of 3.}
  \item{...}{Other arguments to be passed to \code{fbase2}.}
}

\value{
A list with elements \code{f,g,h} corresponding to the function, gradient vector, and Hessian matrix of the function \code{fbase2(X\%*\%beta, Z\%*\%gamma, y, ...)}, where \code{beta=coeff[1:ncol(X)]} and \code{gamma=coeff[ncol(X)+1:ncol(Z)]}. (Derivatives are evaluated relative to \code{coeff}.) In other words, the base function \code{fbase2(u, v, y, ...)} is projected onto the high-dimensional space of \code{c(beta, gamma)} through the linear transformations of its first argument (\code{u <- X\%*\%beta}) and its second argument (\code{v <- Z\%*\%gamma}).
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
\code{\link{regfac.expand.1par}}
}

\examples{
\dontrun{
library(dglm)
library(sns)

# defining log-likelihood function
loglike.linreg <- function(coeff, X, y) {
  regfac.expand.2par(coeff = coeff, X = X, Z = X, y = y
    , fbase2 = fbase2.gaussian.identity.log, fgh = 2, block.diag = T)
}

# simulating data according to generative model
N <- 1000 # number of observations
K <- 5 # number of covariates
X <- matrix(runif(N*K, min=-0.5, max=+0.5), ncol=K)
beta <- runif(K, min=-0.5, max=+0.5)
gamma <- runif(K, min=-0.5, max=+0.5)
mean.vec <- X\%*\%beta
sd.vec <- exp(X\%*\%gamma)
y <- rnorm(N, mean.vec, sd.vec)

# estimation using dglm
est.dglm <- dglm(y~X-1, dformula = ~X-1, family = "gaussian", dlink = "log")
beta.dglm <- est.dglm$coefficients
gamma.dglm <- est.dglm$dispersion.fit$coefficients

# estimation using RegressionFactory
coeff.tmp <- rep(0, 2*K)
for (n in 1:10) {
  coeff.tmp <- sns(coeff.tmp, fghEval=loglike.linreg
    , X=X, y=y, rnd = F)
}
beta.regfac.vd <- coeff.tmp[1:K]
gamma.regfac.vd <- coeff.tmp[K+1:K]

# comparing dglm and RegressionFactory results
# neither beta's nor gamma's will match exactly
cbind(beta.dglm, beta.regfac.vd)
cbind(gamma.dglm, gamma.regfac.vd)
}
}

