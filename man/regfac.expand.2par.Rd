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
  , y, fbase, block.diag=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{coeff}{Vector of coefficients in the regression model. The first \code{ncol(X)} elements correspond to the first parameter of the base distribution \code{fbase(u,v,y)}, and the next \code{ncol(Z)} elements corresponds to the second parameter of the base distribution \code{fbase(u,v,y)}.}
  \item{X}{Matrix of covariates corresponding to the first parameter of the base distribution \code{fbase(u,v,y)}.}
  \item{Z}{Matrix of covariates corresponding to the first parameter of the base distribution \code{fbase(u,v,y)}. Default is a single column of 1's, corresponding to an intercept-only model for the second parameter, i.e. assuming the second parameter is constant across all observations. Note that \code{nrow(Z)} must be equal to \code{nrow(X)}.}
  \item{y}{Vector of response variables. Note that \code{length(y)} must be equal to \code{nrow(X)}.}
  \item{fbase}{Base distribution function \code{fbase(u,v,y)} for the regression model. It must return a list with elements \code{f,g,h} corresponding to the function and its first and second derivatives relative to its first two argument, \code{u,v}. The gradient must be a vector of length 2, and the Hessian must be a vector of length 3, with elements being equal to fbase_uu, fbase_vv and fbase_uv (taking advantage of the Hessian being symmetric).}
  \item{block.diag}{If \code{TRUE}, Hessian matrix is block-diagonalized by setting cross-terms between \code{beta} and \code{gamma} to zero. This can be useful if the full - i.e. non-block-diagonalized - Hessian is not definite, but block-diagonalization leads to definiteness. If \code{FALSE}, third element of the Hessian of \code{fbase} is not needed and thus it can be vector of length 2 instead of 3.}
}

\value{
A list with elements \code{f,g,h} corresponding to the function, gradient vector, and Hessian matrix of the function \code{fbase(X\%*\%beta,Z\%*\%gamma,y)}, where \code{beta=coeff[1:ncol(X)]} and \code{gamma=coeff[ncol(X)+1:ncol(Z)]}. (Derivatives are evaluated relative to \code{coeff}.) In other words, the base function \code{fbase(u,v,y)} is projected onto the high-dimensional space of \code{beta} through the linear transformations of its first argument (\code{u <- X\%*\%beta}) and its second argument (\code{v <- Z\%*\%gamma}).
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

%\examples{
%}

