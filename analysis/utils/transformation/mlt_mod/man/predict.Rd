\name{plot-predict-simulate}
\alias{predict.ctm}
\alias{predict.mlt}
\alias{simulate.ctm}
\alias{simulate.mlt}
\alias{plot.ctm}
\alias{plot.mlt}
\title{
  Plots, Predictions and Samples from mlt Objects
}
\description{
  Plot, predict and sample from objects of class mlt
}
\usage{
\method{plot}{ctm}(x, newdata, type = c("distribution", "survivor", "density", 
     "logdensity", "hazard", "loghazard", "cumhazard", "quantile", "trafo"),
     q = NULL, prob = 1:(K - 1) / K, K = 50, col = rgb(.1, .1, .1, .1), lty = 1, 
     add = FALSE, ...)
\method{plot}{mlt}(x, ...)
\method{predict}{ctm}(object, newdata, type = c("trafo", "distribution", 
         "survivor", "density", "logdensity", "hazard", "loghazard", "cumhazard", 
         "quantile"), terms = c("bresponse", "binteracting", "bshifting"), 
         q = NULL, prob = NULL, K = 50, interpolate = TRUE, ...)
\method{predict}{mlt}(object, newdata = object$data, ...)
\method{simulate}{ctm}(object, nsim = 1, seed = NULL, newdata, K = 50, q = NULL,
         interpolate = TRUE, bysim = TRUE, ...)
\method{simulate}{mlt}(object, nsim = 1, seed = NULL, newdata = object$data, bysim = TRUE, ...)
}
\arguments{
  \item{object}{a fitted conditional transformation model as returned by \code{\link{mlt}}
                or an unfitted conditional transformation model as returned by \code{\link{ctm}}}
  \item{x}{a fitted conditional transformation model as returned by \code{\link{mlt}}}
  \item{newdata}{an optional data frame of observations}
  \item{type}{type of prediction or plot to generate}
  \item{q}{quantiles at which to evaluate the model}
  \item{prob}{probabilities for the evaluation of the quantile function (\code{type = "quantile"})}
  \item{terms}{terms to evaluate for the predictions, corresponds to the argument
     \code{response}, \code{interacting} and \code{shifting} in \code{\link{ctm}}}
  \item{K}{number of grid points to generate (in the absence of \code{q})}
  \item{col}{color for the lines to plot}
  \item{lty}{line type for the lines to plot}
  \item{add}{logical indicating if a new plot shall be generated (the default)}
  \item{interpolate}{logical indicating if quantiles shall be interpolated linearily}
  \item{nsim}{number of samples to generate}
  \item{seed}{optional seed for the random number generator}
  \item{bysim}{logical, if \code{TRUE} a list with \code{nsim} elements is returned,
               each element is of length \code{nrow(newdata)} and 
               contains one sample from the conditional distribution for each
               row of \code{newdata}. If \code{FALSE}, a list of length \code{nrow(newdata)}
               is returned, its ith element of length \code{nsim} contains \code{nsim} samples
               from the conditional distribution given \code{newdata[i,]}.} 
  \item{\dots}{additional arguments}
}
\details{
  \code{plot} evaluates the transformation function over a grid of \code{q} values
  for all observations in \code{newdata} and plots these functions (according to 
  \code{type}). \code{predict} evaluates the transformation function over a grid
  of \code{q} values for all observations in \code{newdata} and returns the
  result as a matrix (where _columns_ correspond to _rows_ in \code{newdata}).
  Note that the \code{predict} method for \code{ctm} objects requires all
  model coefficients to be specified in this unfitted model.
  \code{simulate} draws samples from \code{object} by numerical inversion of the
  quantile function.
}
