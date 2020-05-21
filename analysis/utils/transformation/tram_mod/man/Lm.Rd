\name{Lm}
\alias{Lm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Normal Linear Model
}
\description{
  Normal linear model with benefits
}
\usage{
Lm(formula, data, subset, weights, offset, cluster, na.action = na.omit, ...)
}
%- maybe also 'usage' for other objects documented here.
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
          contain \code{NA}s.  The default is set to \code{na.omit}.
  }
  \item{\dots}{additional arguments to \code{\link{tram}}.}
}
\details{
A normal linear model with simulaneous estimation of regression coefficients
and scale parameter(s). This function also allows for stratum-specific
intercepts and variances as well as censoring and truncation in the
response.

Note that the scale of the parameters is different from what is reported by
\code{\link[stats]{lm}}; the discrepancies are explained in the package
vignette.

The model is defined with a negative shift term. Large values of the
linear predictor correspond to large values of the conditional 
expectation response.

}
\value{
An object of class \code{BoxCox}, with corresponding \code{coef},
\code{vcov}, \code{logLik}, \code{estfun}, \code{summary}, 
\code{print}, \code{plot} and \code{predict} methods.
}
\references{
  Torsten Hothorn, Lisa Moest, Peter Buehlmann (2018), Most Likely
  Transformations, \emph{Scandinavian Journal of Statistics}, \bold{45}(1),
  110--134, \doi{10.1111/sjos.12291}.
}
\examples{

  data("BostonHousing2", package = "mlbench")

  lm(cmedv ~ crim + zn + indus + chas + nox + rm + age + dis + 
             rad + tax + ptratio + b + lstat, data = BostonHousing2)

  Lm(cmedv ~ chas + crim + zn + indus + nox + 
             rm + age + dis + rad + tax + ptratio + b + lstat, 
             data = BostonHousing2)
}
\keyword{models}
\keyword{regression}
