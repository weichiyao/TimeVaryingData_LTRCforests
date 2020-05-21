\name{tram}
\alias{tram}
\alias{tram_data}
\title{
Stratified Linear Transformation Models
}
\description{
Likelihood-inference for stratified linear transformation models
}
\usage{
tram(formula, data, subset, weights, offset, cluster, na.action = na.omit, 
     distribution = c("Normal", "Logistic", "MinExtrVal"), 
     transformation = c("discrete", "linear", "logarithmic", "smooth"), 
     LRtest = TRUE, prob = c(0.1, 0.9), support = NULL, 
     bounds = NULL, add = c(0, 0), order = 6, 
     negative = TRUE, scale = TRUE, extrapolate = FALSE, 
     log_first = FALSE, model_only = FALSE, ...)
tram_data(formula, data, subset, weights, offset, cluster, na.action = na.omit) 
}
\arguments{
  \item{formula}{an object of class \code{"formula"}: a symbolic description 
          of the model structure to be
          fitted.  The details of model specification are given under
          Details and in the package vignette.
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
  \item{distribution}{character specifying how the transformation function
                      is mapped into probabilities. Available choices
                      include the cumulative distribution functions of the standard normal, the
                      standard logistic and the standard minimum extreme value distribution.
  }
  \item{transformation}{character specifying the complexity of the
       response-transformation. For discrete responses, one parameter is assigned
       to each level (except the last one), for continuous responses linear,
       log-linear and smooth (parameterised as a Bernstein polynomial) function
       are implemented.
  }
  \item{LRtest}{logical specifying if a likelihood-ratio test for the null
                of all coefficients in the linear predictor being zero shall be performed.
  }
  \item{prob}{two probabilities giving quantiles of the response defining the support of a smooth
      Bernstein polynomial (if \code{transformation = "smooth"}).
  }
  \item{support}{a vector of two elements; the support of a smooth
      Bernstein polynomial (if \code{transformation = "smooth"}).
  }
  \item{bounds}{an interval defining the bounds of a real sample space.
  }
  \item{add}{add these values to the support before generating a grid via
             \code{\link[variables]{mkgrid}}.
  }
  \item{order}{integer >= 1 defining the order of the Bernstein polynomial 
               (if \code{transformation = "smooth"}).
  }
  \item{negative}{logical defining the sign of the linear predictor.
  }
  \item{scale}{logical defining if variables in the linear predictor shall
               be scaled. Scaling is internally used for model estimation,
               rescaled coefficients are reported in model output.
  }
  \item{extrapolate}{logical defining the behaviour of the Bernstein transformation
                     function outside \code{support}. The default
                     \code{FALSE} is to extrapolate linearily without requiring the second
                     derivative of the transformation function to be zero at \code{support}.
                     If \code{TRUE}, this additional constraint is respected.
  }
  \item{log_first}{logical; if \code{TRUE}, a Bernstein polynomial is
                   defined on the log-scale.}
  \item{model_only}{logical, if \code{TRUE} the unfitted model is returned.
  }
  \item{\dots}{additional arguments.}
}
\details{
The model formula is of the form \code{y | s ~ x} where \code{y} is an at
least ordered response variable, \code{s} are the variables defining strata
and \code{x} defines the linear predictor. \code{y ~ x} defines a model
without strata (but response-varying intercept function) and \code{y | s ~
0} sets-up response-varying coefficients for all variables in \code{s}.

The two functions \code{tram} and \code{tram_data} are not intended 
to be called directly by users. Instead,
functions \code{\link{Coxph}} (Cox proportional hazards models), 
\code{\link{Survreg}} (parametric survival models), 
\code{\link{Polr}} (models for ordered categorical responses), 
\code{\link{Lm}} (normal linear models),
\code{\link{BoxCox}} (non-normal linear models) or 
\code{\link{Colr}} (continuous outcome logistic regression) allow
direct access to the corresponding models.

The model class and the specific models implemented in \pkg{tram} are 
explained in the package vignette of package \pkg{tram}.
The underlying theory of most likely transformations 
is presented in Hothorn et al. (2018), computational
and modelling aspects in more complex situations 
are discussed by Hothorn (2018).

}
\value{
An object of class \code{tram} inheriting from \code{mlt}.
}
\references{
  Torsten Hothorn, Lisa Moest, Peter Buehlmann (2018), Most Likely
  Transformations, \emph{Scandinavian Journal of Statistics}, \bold{45}(1),
  110--134, \doi{10.1111/sjos.12291}.

  Torsten Hothorn (2018), Most Likely Transformations: The mlt Package,
  \emph{Journal of Statistical Software}, forthcoming. URL: 
  \url{https://cran.r-project.org/package=mlt.docreg}

}
\keyword{models}
\keyword{regression}
\keyword{smooth}
\keyword{survival}
