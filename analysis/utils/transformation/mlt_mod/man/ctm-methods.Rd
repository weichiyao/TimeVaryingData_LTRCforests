\name{ctm-methods}
\alias{variable.names.ctm}
\alias{coef<-.ctm}
\alias{coef.ctm}
\title{
  Methods for ctm Objects
}
\description{
  Methods for objects of class ctm
}
\usage{
\method{variable.names}{ctm}(object, 
              which = c("all", "response", "interacting", "shifting"), 
              ...)
\method{coef}{ctm}(object, ...)
}
\arguments{
  \item{object}{an unfitted conditional transformation model as returned by \code{\link{ctm}}}
  \item{which}{a character specifying which names shall be returned}
  \item{\dots}{additional arguments}
}
\details{
\code{coef} can be used to get and set model parameters.
}