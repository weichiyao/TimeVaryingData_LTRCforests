\name{Polr}
\alias{Polr}
\title{
  Ordered Categorical Regression
}
\description{
  Some regression models for ordered categorical responses
}
\usage{
Polr(formula, data, subset, weights, offset, cluster, na.action = na.omit, 
     method = c("logistic", "probit", "loglog", "cloglog"), ...)
}
\arguments{
  \item{formula}{an object of class \code{"formula"}: a symbolic description 
          of the model structure to be
          fitted.  The details of model specification are given under
          \code{\link{tram}} and in the package vignette.
  }
  \item{data}{an optional data frame, list or environment (or object
          coercible by \code{as.data.frame} to a data frame) containing the
          variables in the model.  If not found in \code{data}, the
          variables are taken from \code{environment(formula)}.
  }
  \item{subset}{an optional vector specifying a subset of observations to be
          used in the fitting process.
  }
  \item{weights}{an optional vector of weights to be used in the fitting
          process.  Should be \code{NULL} or a numeric vector. If present,
          the weighted log-likelihood is maximised.
  }
  \item{offset}{this can be used to specify an _a priori_ known component to
          be included in the linear predictor during fitting.  This
          should be \code{NULL} or a numeric vector of length equal to the
          number of cases.
  }
  \item{cluster}{optional factor with a cluster ID employed for computing
          clustered covariances.
  }
  \item{na.action}{a function which indicates what should happen when the data
          contain \code{NA}s.  The default is set by the \code{na.action} setting
          of \code{options}, and is \code{na.fail} if that is unset.
  }
  \item{method}{a character describing the link function.
  }
  \item{\dots}{additional arguments to \code{\link{tram}}.}
}
\details{
Models for ordered categorical responses reusing the interface of
\code{\link[MASS]{polr}}. Allows for stratification, censoring and
trunction.

The model is defined with a negative shift term, thus \code{exp(coef())}
is the multiplicative change of the odds ratio (conditional odds for 
reference divided by conditional odds of treatment or for a one unit
increase in a numeric variable). Large values of the
linear predictor correspond to large values of the conditional 
expectation response (but this relationship is nonlinear).

}
\value{
An object of class \code{Polr}, with corresponding \code{coef},
\code{vcov}, \code{logLik}, \code{estfun}, \code{summary}, 
\code{print}, \code{plot} and \code{predict} methods.
}
\references{
  Torsten Hothorn, Lisa Moest, Peter Buehlmann (2018), Most Likely
  Transformations, \emph{Scandinavian Journal of Statistics}, \bold{45}(1),
  110--134, \doi{10.1111/sjos.12291}.
}
\examples{

  data("wine", package = "ordinal")

  library("MASS")
  polr(rating ~ temp + contact, data = wine)

  Polr(rating ~ temp + contact, data = wine)

}
\keyword{models}
\keyword{regression}
