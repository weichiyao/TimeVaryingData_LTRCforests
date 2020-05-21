\name{confband}
\alias{confband}
\alias{confband.mlt}
\title{
  Confidence Bands
}
\description{
  Confidence bands for transformation, distribution, survivor or
  cumulative hazard functions
}
\usage{
confband(object, newdata, level = 0.95, ...)
\method{confband}{mlt}(object, newdata, level = 0.95, 
       type = c("trafo", "distribution", "survivor", "cumhazard"), 
       K = 20, cheat = K, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object of class \code{\link{mlt}}}
  \item{newdata}{a data frame of observations}
  \item{level}{the confidence level}
  \item{type}{the function to compute the confidence band for}
  \item{K}{number of grid points the function is evaluated at}
  \item{cheat}{number of grid points the function is evaluated at when
               using the quantile obtained for \code{K} grid points}
  \item{\dots}{additional arguments to \code{\link[multcomp]{confint.glht}}}
}
\details{
  The function is evaluated at \code{K} grid points and simultaneous
  confidence intervals are then interpolated in order to construct the band.

  A smoother band can be obtained by setting \code{cheat} to something larger 
  than \code{K}: The quantile is obtained for \code{K} grid points but
  the number of evaluated grid points \code{cheat} can be much larger at no 
  additional cost. Technically, the nominal level is not maintained in
  this case but the deviation will be small for reasonably large \code{K}.
}
\value{
  For each row in \code{newdata} the function and corresponding confidence
  band evaluated at the \code{K} (or \code{cheat}) grid points is returned.
}
