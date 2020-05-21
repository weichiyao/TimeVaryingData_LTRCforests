\name{Coxph}
\alias{Coxph}
\title{
  Cox Proportional Hazards Model
}
\description{
  Cox model with fully parameterised baseline hazard function
}
\usage{
Coxph(formula, data, subset, weights, offset, cluster, na.action = na.omit, ...)
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
The original implementation of Cox models via the partial likelihood,
treating the baseline hazard function as a nuisance parameter, is available
in \code{\link[survival]{coxph}}. This function allows simultaneous
estimation of the log-hazard ratios and the log-cumulative baseline hazard,
the latter parameterised by a Bernstein polynomial. The model can be fitted
under stratification (time-varying coefficients), all types of random
censoring and trunction. An early reference to this parameterisation is
McLain and Ghosh (2013).

The responses is bounded (\code{bounds = c(0, Inf)}) when specified as a
\code{Surv} object. Otherwise, \code{bounds} can be specified via
\code{\dots}.

Parameters are log-hazard ratios comparing treatment (or a one unit increase
in a numeric variable) with a reference.

}
\value{
An object of class \code{Coxph}, with corresponding \code{coef},
\code{vcov}, \code{logLik}, \code{estfun}, \code{summary}, 
\code{print}, \code{plot} and \code{predict} methods.
}
\references{

  Alexander C. McLain and Sujit K. Ghosh (2013).
  Efficient Sieve Maximum Likelihood Estimation of
  Time-Transformation Models, \emph{Journal of Statistical Theory and Practice},
  \bold{7}(2), 285--303, \doi{10.1080/15598608.2013.772835}.

  Torsten Hothorn, Lisa Moest, Peter Buehlmann (2018), Most Likely
  Transformations, \emph{Scandinavian Journal of Statistics}, \bold{45}(1),
  110--134, \doi{10.1111/sjos.12291}.

}
\examples{

  data("GBSG2", package = "TH.data")

  library("survival")
  coxph(Surv(time, cens) ~ horTh, data = GBSG2)

  Coxph(Surv(time, cens) ~ horTh, data = GBSG2)

}
\keyword{models}
\keyword{regression}
\keyword{survival}
\keyword{smooth}
