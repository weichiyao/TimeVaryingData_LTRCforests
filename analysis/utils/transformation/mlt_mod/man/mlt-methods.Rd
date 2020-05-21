\name{mlt-methods}
\alias{coef<-}
\alias{coef<-.mlt}
\alias{coef.mlt}
\alias{weights.mlt}
\alias{logLik.mlt}
\alias{vcov.mlt}
\alias{Hessian}
\alias{Gradient}
\alias{Hessian.mlt}
\alias{Gradient.mlt}
\alias{estfun.mlt}
\alias{mkgrid.mlt}
\alias{bounds.mlt}
\alias{variable.names.mlt}
\alias{update.mlt_fit}
\alias{as.mlt}
\alias{as.mlt.mlt}
\title{
  Methods for mlt Objects
}
\description{
  Methods for objects of class mlt
}
\usage{
\method{coef}{mlt}(object, fixed = TRUE, ...)
coef(object) <- value
\method{weights}{mlt}(object, ...)
\method{logLik}{mlt}(object, parm = coef(object, fixed = FALSE), w = NULL, newdata, ...)
\method{vcov}{mlt}(object, parm = coef(object, fixed = FALSE), complete = FALSE, ...)
Hessian(object, ...)
\method{Hessian}{mlt}(object, parm = coef(object, fixed = FALSE), ...)
Gradient(object, ...)
\method{Gradient}{mlt}(object, parm = coef(object, fixed = FALSE), ...)
\method{estfun}{mlt}(object, parm = coef(object, fixed = FALSE),
       w = NULL, newdata, ...)
\method{mkgrid}{mlt}(object, n, ...)
\method{bounds}{mlt}(object)
\method{variable.names}{mlt}(object, ...)
\method{update}{mlt_fit}(object, weights, subset = NULL, theta, ...)
\method{as.mlt}{mlt}(object)
}
\arguments{
  \item{object}{a fitted conditional transformation model as returned by \code{\link{mlt}}}
  \item{fixed}{a logical indicating if only estimated coefficients (\code{fixed = FALSE}) 
               should be returned}
  \item{value}{coefficients to be assigned to the model}
  \item{parm}{model parameters}
  \item{w}{model weights}
  \item{weights}{model weights}
  \item{newdata}{an optional data frame of new observations. Allows
                 evaluation of the log-likelihood for a given
                 model \code{object} on these new observations. The
                 parameters \code{parm} and \code{w} are ignored in this situation.}
  \item{n}{number of grid points}
  \item{subset}{an optional integer vector indicating the subset of
                observations to be used for fitting.}
  \item{theta}{optional starting values for the model parameters}
  \item{complete}{currently ignored}
  \item{\dots}{additional arguments}
}
\details{
  \code{coef} can be used to get and set model parameters, \code{weights} and 
  \code{logLik} extract weights and evaluate the log-likelihood (also for
  parameters other than the maximum likelihood estimate). \code{Hessian}
  returns the Hessian and \code{vcov} the inverse thereof. \code{Gradient}
  gives the gradient (sum of the score contributions) 
  and \code{estfun} the score contribution by each observation. \code{mkgrid}
  generates a grid of all variables (as returned by \code{variable.names}) in the model.
  \code{update} allows refitting the model with alternative weights and potentially
  different starting values. \code{bounds} gets bounds for bounded variables in the model.
}
