\name{R}
\alias{R}
\alias{R.Surv}
\alias{R.factor}
\alias{R.ordered}
\alias{R.numeric}
\alias{R.integer}
\alias{R.list}
\alias{R.response}
\alias{as.Surv}
\alias{as.Surv.response}
\title{
    Response Variable
}
\description{
    Represent a possibly censored or truncated response variable
}
\usage{
R(object, ...)
\method{R}{numeric}(object = NA, cleft = NA, cright = NA,
   tleft = NA, tright = NA, tol = sqrt(.Machine$double.eps), ...)
\method{R}{ordered}(object, cleft = NA, cright = NA, ...)
\method{R}{integer}(object, cleft = NA, cright = NA, bounds = c(0L, Inf), ...)
\method{R}{factor}(object, ...)
\method{R}{Surv}(object, ...)
as.Surv(object)
\method{as.Surv}{response}(object)
}
\arguments{
  \item{object}{
  A vector of (conceptually) exact measurements or an object of class
  \code{response} (for \code{as.Surv}) or a list.}
  \item{cleft}{
  A vector of left borders of censored measurements}
  \item{cright}{
  A vector of right borders of censored measurements}
  \item{tleft}{
  A vector of left truncations}
  \item{tright}{
  A vector of right truncations}
  \item{tol}{
  Tolerance for checking if \code{cleft} < \code{cright}}
  \item{bounds}{
  Range of possible values for integers}
  \item{\dots}{
  other arguments, ignored except for \code{tleft} and \code{tright} to
  \code{R.ordered} and \code{R.integer}}
}
\details{
  \code{R} is basically an extention of \code{\link[survival]{Surv}} for the
  representation of arbitrarily censored or truncated measurements at any scale.

  \code{R} applied to a list calls \code{R} for each of the list elements
  and returns a joint object.
}
\examples{

 ### ordered factor
 R(gl(3, 3, labels = LETTERS[1:3]))

}
