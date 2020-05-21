\name{BoxCox}
\alias{BoxCox}
\title{
  (Similar to) Box-Cox Models
}
\description{
  Non-normal linear regression inspired by Box-Cox models
}
\usage{
BoxCox(formula, data, subset, weights, offset, cluster, na.action = na.omit, ...)
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
          contain \code{NA}s.  The default is set to \code{na.omit}.
  }
  \item{\dots}{additional arguments to \code{\link{tram}}.}
}
\details{
A normal model for transformed responses, where the transformation is
estimated from the data simultaneously with the regression coefficients.
This is similar to a Box-Cox transformation, but the technical details
differ. Examples can be found in the package vignette.

The model is defined with a negative shift term. Large values of the
linear predictor correspond to large values of the conditional 
expectation response (but this relationship is potentially nonlinear).

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

  BoxCox(cmedv ~ chas + crim + zn + indus + nox + 
                 rm + age + dis + rad + tax + ptratio + b + lstat, 
                 data = BostonHousing2)
}
\keyword{models}
\keyword{regression}
\keyword{smooth}
