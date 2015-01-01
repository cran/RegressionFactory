\name{fbase2.inverse.gaussian.log.log}
\alias{fbase2.inverse.gaussian.log.log}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Double-Parameter Base Log-likelihood Function for Inverse-Gaussian GLM
}
\description{
Vectorized, double-parameter base log-likelihood functions for Inverse-Gaussian GLM. The link functions map the mean and dispersion parameter to linear predictors. The base function can be supplied to the expander function \code{\link{regfac.expand.2par}} in order to obtain the full, high-dimensional log-likleihood and its derivatives.
}
\usage{
fbase2.inverse.gaussian.log.log(u, v, y, fgh = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{u}{First parameter of the base log-likelihood function (usually the result of applying a link function to distribution mean). This parameter is intended to be projected onto a high-dimensional space using the familiar regression transformation of \code{u <- X\%*\%beta}. In the typical use-case where the caller is \code{regfac.expand.2par}, a vector of values are supplied, and return objects will have the same length as \code{u} or \code{v}.}
  \item{v}{Second parameter of the base log-likelihood function (usually the result of applying a link function to distribution dispersion parameter). This parameter is intended to be projected onto a high-dimensional space using the familiar regression transformation of \code{v <- Z\%*\%gamma}. In the typical use-case where the caller is \code{regfac.expand.2par}, a vector of values are supplied, and return objects will have the same length as \code{u} or \code{v}.}
  \item{y}{Fixed slot of the base distribution, corresponding to the response variable in the regression model.}
  \item{fgh}{Integer with possible values 0,1,2. If \code{fgh=0}, the function only calculates and returns the log-likelihood vector and no derivatives. If \code{fgh=1}, it returns the log-likelihood and its first derivative in a list. If \code{fgh=2}, it returns the log-likelihood, as well as its first and second derivatives in a list.}
}

\value{
If \code{fgh==0}, the function returns \code{-v/2 - 0.5*exp(-v-2*u)*(y - exp(u))^2/y}. If \code{fgh==1}, a list is returned with elements \code{f} and \code{g}, where \code{f} is the same object as in \code{fgh==0} and \code{g} is a matrix of dimensions \code{length(u)}-by-2, with first column being the derivative of the above expression with respect to \code{u}, and the second column being the derivative of the above expression with respect to \code{v}. If \code{fgh==2}, the list will include an element named \code{h}, which is a matrix of dimensions \code{length(u)}-by-3, with the first column being the second derivative of \code{f} with respect to \code{u}, the second column being the second derivative of \code{f} with respect to \code{v}, and the third column is the cross-derivative term, i.e. the derivative of \code{f} with respect to \code{u} and \code{v}.
}

\author{
Alireza S. Mahani, Mansour T.A. Sharabiani
}

\seealso{
regfac.expand.2par
}

\examples{
\dontrun{
# we use this library for univariate slice sampling
# of multivariate distributions
library(MfUSampler)
library(dglm)

# simulating data according to assumed generative model
# we assume log link functions for both mean and dispersion
# (shape parameter is inverse of dispersion)
N <- 10000
K <- 5
X <- cbind(1,matrix(runif(N*(K-1), min=-0.5, max=+0.5), ncol=K-1))
beta <- runif(K, min=-0.5, max=+0.5)
gamma <- runif(K, min=-0.5, max=+0.5)
mean.vec <- exp(X \%*\% beta)
dispersion.vec <- exp(X \%*\% gamma)
y <- rinvgauss(N, mean = mean.vec, dispersion = dispersion.vec)

# model estimation using dglm package
reg.dglm <- dglm(y~X-1, dformula = ~X-1, family=inverse.gaussian(link="log"), dlink = "log")
beta.dglm <- reg.dglm$coefficients
gamma.dglm <- reg.dglm$dispersion.fit$coefficients

# model estimation using RegressionFactory
# (with univariate slice sampling)
# defining the log-likelihood using the expander framework
# assumng same covariates for both slots, hence we set Z=X
# slice sampler does not need derivatives, hence we set fgh=0
loglike.inverse.gaussian <- function(coeff, X, y) {
  regfac.expand.2par(coeff, X=X, Z=X, y=y, fbase2=fbase2.inverse.gaussian.log.log, fgh=0)
}
nsmp <- 100
coeff.smp <- array(NA, dim=c(nsmp, 2*K)) 
coeff.tmp <- rep(0.1, 2*K)
for (n in 1:nsmp) {
  coeff.tmp <- MfU.Sample(coeff.tmp, f=loglike.inverse.gaussian, X=X, y=y)
  coeff.smp[n,] <- coeff.tmp
}
beta.slice <- colMeans(coeff.smp[(nsmp/2+1):nsmp, 1:K])
gamma.slice <- colMeans(coeff.smp[(nsmp/2+1):nsmp, K+1:K])

# compare results
cbind(beta.dglm, beta.slice)
cbind(gamma.dglm, gamma.slice)

}
}

