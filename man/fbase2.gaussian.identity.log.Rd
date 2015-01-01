\name{fbase2.gaussian.identity.log}
\alias{fbase2.gaussian.identity.log}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Double-Parameter Base Log-likelihood Function for Gaussian GLM
}
\description{
Vectorized, double-parameter base log-likelihood functions for Gaussian GLM. The link functions map the mean and variance to linear predictors. The base function can be supplied to the expander function \code{\link{regfac.expand.2par}} in order to obtain the full, high-dimensional log-likleihood and its derivatives.
}
\usage{
fbase2.gaussian.identity.log(u, v, y, fgh = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{u}{First parameter of the base log-likelihood function (usually the mean). This parameter is intended to be projected onto a high-dimensional space using the familiar regression transformation of \code{u <- X\%*\%beta}. In the typical use-case where the caller is \code{regfac.expand.1par}, a vector of values are supplied, and return objects will have the same length as \code{u}.}
  \item{v}{Second parameter of the base log-likelihood function (usually the mean). This parameter is intended to be projected onto a high-dimensional space using the familiar regression transformation of \code{v <- Z\%*\%gamma}. In the typical use-case where the caller is \code{regfac.expand.1par}, a vector of values are supplied, and return objects will have the same length as \code{u}.}
  \item{y}{Fixed slot of the base distribution, corresponding to the response variable in the regression model. For \code{binomial} family, it must be an integer vector with values between \code{0} and \code{n}.}
  \item{fgh}{Integer with possible values 0,1,2. If \code{fgh=0}, the function only calculates and returns the log-likelihood vector and no derivatives. If \code{fgh=1}, it returns the log-likelihood and its first derivative in a list. If \code{fgh=2}, it returns the log-likelihood, as well as its first and second derivatives in a list.}
}

\value{
If \code{fgh==0}, the function returns \code{-0.5*(v + exp(-v)*(u-y)*(u-y))} (ignoring additive terms that are independent of \code{u,v}). It will therefore be of the same length as \code{u}. If \code{fgh==1}, a list is returned with elements \code{f} and \code{g}, where \code{f} is the same object as in \code{fgh==1} and \code{g} is a matrix of dimensions \code{length(u)}-by-2, with first column being the derivative of the above expression with respect to \code{u}, and the second column being the derivative of the above expression with respect to \code{v}. If \code{fgh==2}, the list will include an element named \code{h}, which is a matrix of dimensions \code{length(u)}-by-3, with the first column being the second derivative of \code{f} with respect to \code{u}, the second column being the second derivative of \code{f} with respect to \code{v}, and the third column is the cross-derivative term, i.e. the derivative of \code{f} with respect to \code{u} and \code{v}.
}

\author{
Alireza S. Mahani, Mansour T.A. Sharabiani
}

\seealso{
regfac.expand.2par
}

\examples{
\dontrun{
library(sns)
library(MfUSampler)
library(dglm)

# defining log-likelihood function
# vd==FALSE leads to constant-dispersion model (ordinary linear regression)
# while vd==TRUE produces varying-dispersion model
loglike.linreg <- function(coeff, X, y, fgh, block.diag = F, vd = F) {
  if (vd) regfac.expand.2par(coeff = coeff, X = X, Z = X, y = y
    , fbase2 = fbase2.gaussian.identity.log, fgh = fgh, block.diag = block.diag)
  else regfac.expand.2par(coeff = coeff, X = X, y = y
    , fbase2 = fbase2.gaussian.identity.log, fgh = fgh, block.diag = block.diag)
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

# constant-dispersion model
# estimation using glm
est.glm <- lm(y~X-1)
beta.glm <- est.glm$coefficients
sigma.glm <- summary(est.glm)$sigma
# estimation using RegressionFactory
# (we set rnd=F in sns to allow for better comparison with glm)
nsmp <- 20
coeff.smp <- array(NA, dim=c(nsmp, K+1)) 
coeff.tmp <- rep(0, K+1)
for (n in 1:nsmp) {
  coeff.tmp <- sns(coeff.tmp, fghEval=loglike.linreg
    , X=X, y=y, fgh=2, block.diag = F, vd = F, rnd = F)
  coeff.smp[n,] <- coeff.tmp
}
beta.regfac.cd <- colMeans(coeff.smp[(nsmp/2+1):nsmp, 1:K])
sigma.regfac.cd <- sqrt(exp(mean(coeff.smp[(nsmp/2+1):nsmp, K+1])))
# comparing glm and RegressionFactory results
# beta's must match exactly between glm and RegressionFactory
cbind(beta, beta.glm, beta.regfac.cd)
# sigma's won't match exactly
cbind(mean(sd.vec), sigma.glm, sigma.regfac.cd)

# varying-dispersion model
# estimation using dglm
est.dglm <- dglm(y~X-1, dformula = ~X-1, family = "gaussian", dlink = "log")
beta.dglm <- est.dglm$coefficients
gamma.dglm <- est.dglm$dispersion.fit$coefficients
# estimation using RegressionFactory
coeff.smp <- array(NA, dim=c(nsmp, 2*K)) 
coeff.tmp <- rep(0, 2*K)
for (n in 1:nsmp) {
  coeff.tmp <- sns(coeff.tmp, fghEval=loglike.linreg
    , X=X, y=y, fgh=2, block.diag = F, vd = T, rnd = F)
  coeff.smp[n,] <- coeff.tmp
}
beta.regfac.vd <- colMeans(coeff.smp[(nsmp/2+1):nsmp, 1:K])
gamma.regfac.vd <- colMeans(coeff.smp[(nsmp/2+1):nsmp, K+1:K])
# comparing dglm and RegressionFactory results
# neither beta's nor gamma's will match exactly
cbind(beta, beta.dglm, beta.regfac.vd)
cbind(gamma, gamma.dglm, gamma.regfac.vd)

}
}
