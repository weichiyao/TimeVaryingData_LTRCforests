\name{Colr}
\alias{Colr}
\title{
  Continuous Outcome Logistic Regression
}
\description{
  A proportional-odds model for continuous variables
}
\usage{
Colr(formula, data, subset, weights, offset, cluster, na.action = na.omit, ...)
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
          contain \code{NA}s.  The default is set by the \code{na.action} setting
          of \code{options}, and is \code{na.fail} if that is unset.
  }
  \item{\dots}{additional arguments to \code{\link{tram}}.}
}
\details{
Simultanous estimation of all possible binary logistic models obtained 
by dichotomisation of a continuous response. The regression coefficients can
be constant allowing for an interpretation as log-odds ratios.

The model is defined with a positive shift term, thus \code{exp(coef())} is
the multiplicative change of the odds ratio (conditional odds of treatment
or for a one unit increase in a numeric variable divided by conditional odds
of reference).  Large values of the linear predictor correspond to small
values of the conditional expectation response (but this relationship is
nonlinear).

}
\value{
An object of class \code{Colr}, with corresponding \code{coef},
\code{vcov}, \code{logLik}, \code{estfun}, \code{summary}, 
\code{print}, \code{plot} and \code{predict} methods.
}
\references{

  Tina Lohse, Sabine Rohrmann, David Faeh and Torsten Hothorn (2017),
  Continuous Outcome Logistic Regression for Analyzing Body Mass Index
  Distributions, \emph{F1000Research}, \bold{6}(1933),
  \doi{10.12688/f1000research.12934.1}.

  Torsten Hothorn, Lisa Moest, Peter Buehlmann (2018), Most Likely
  Transformations, \emph{Scandinavian Journal of Statistics}, \bold{45}(1),
  110--134, \doi{10.1111/sjos.12291}.

}
\examples{

  data("BostonHousing2", package = "mlbench")

  lm(cmedv ~ crim + zn + indus + chas + nox + rm + age + dis + 
             rad + tax + ptratio + b + lstat, data = BostonHousing2)

  Colr(cmedv ~ chas + crim + zn + indus + nox + 
               rm + age + dis + rad + tax + ptratio + b + lstat, 
               data = BostonHousing2)
}
\keyword{models}
\keyword{regression}
\keyword{smooth}
