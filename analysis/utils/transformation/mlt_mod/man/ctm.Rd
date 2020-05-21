\name{ctm}
\alias{ctm}
\title{
  Conditional Transformation Models
}
\description{
  Specification of conditional transformation models
}
\usage{
ctm(response, interacting = NULL, shifting = NULL, data = NULL, 
    todistr = c("Normal", "Logistic", "MinExtrVal", "MaxExtrVal"), 
    sumconstr = inherits(interacting, c("formula", "formula_basis")), ...)
}
\arguments{
  \item{response}{a basis function, ie, an object of class \code{basis}}
  \item{interacting}{a basis function, ie, an object of class \code{basis}}
  \item{shifting}{a basis function, ie, an object of class \code{basis}}
  \item{data}{either a \code{data.frame} containing the model variables
              or a formal description of these variables in an object of class \code{vars}}
  \item{todistr}{a character vector describing the distribution to be transformed}
  \item{sumconstr}{a logical indicating if sum constraints shall be applied}
  \item{\dots}{arguments to \code{as.basis} when \code{shifting} is a formula}
}
\details{
  This function only specifies the model which can then be fitted using
  \code{\link{mlt}}. The shift term is positive by default.

  Possible choices of the distributions the model transforms to (the inverse
  link functions) include the 
  standard normal (\code{"Normal"}), the standard logistic
  (\code{"Logistic"}), the standard minimum extreme value
  (\code{"MinExtrVal"}, also known as Gompertz distribution), and the
  standard maximum extreme value (\code{"MaxExtrVal"}, also known as Gumbel
  distribution) distributions.

}
\value{
  An object of class \code{ctm}.
}
\references{

  Torsten Hothorn, Lisa Moest, Peter Buehlmann (2018), Most Likely
  Transformations, \emph{Scandinavian Journal of Statistics}, \bold{45}(1),
  110--134, \doi{10.1111/sjos.12291}.

}
