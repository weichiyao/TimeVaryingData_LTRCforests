\name{mltoptim}
\alias{mltoptim}
\title{
  Control Optimisation
}
\description{
  Define optimisers and their control parameters
}
\usage{
mltoptim(auglag = list(maxtry = 5, kkt2.check = FALSE), 
         spg = list(maxit = 10000, quiet = TRUE, checkGrad = FALSE), 
         trace = FALSE)
}
\arguments{
  \item{auglag}{
  A list with control parameters for the \code{\link[alabama]{auglag}} optimiser.
  \code{maxtry} is the number of times the algorithm is started on random starting
  values in case it failed with the precomputed ones.
}
  \item{spg}{
  A list with control parameters for the \code{\link[BB]{BBoptim}} optimiser (calling
   \code{\link[BB]{spg}} internally).
}
  \item{trace}{
  A logical switching trace reports by the optimisers off.
}
}
\details{
  This function sets-up functions to be called in \code{\link{mlt}} internally.
}
\value{
  A list of functions with arguments \code{theta} (starting values), \code{f} (log-likelihood),
  \code{g} (scores), \code{ui} and \code{ci} (linear inequality constraints).
  Adding further such functions is a way to add more optimisers to \code{\link{mlt}}.
  The first one in this list converging  defines the resulting model.
}
\keyword{list}
