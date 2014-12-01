\name{fbase.bern.ilogit}
\alias{fbase.bern.ilogit}

\title{
Single-Parameter Base Distribution for Logistic Regression
}
\description{
Single-Parameter Base Distribution for Logistic Regression (Bernoulli distribution + inverse logit transformation), to be used with expander function \code{\link{regfac.expand.1par}} to obtain full, high-dimensional derivatives.
}
\usage{
fbase.bern.ilogit(u,y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{u}{Parameter slot of the base distribution. This slot is intended to be projected onto a high-dimensional space using the familiar regression transformation of \code{u <- X\%*\%beta}.}
  \item{y}{The fixed slot of the base distribution, corresponding to the response variable.}
}

\value{
The function value is \code{-(log(1+exp(-u))+(1-y)*u)}. (Derivatives are calculated by differentiating with respect to \code{u}.) This corresponds to the log-likelihood of observing the binary variable \code{y} if it follows the Bernoulli distribution with parameter \code{1/(1+exp(-u))}.
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
