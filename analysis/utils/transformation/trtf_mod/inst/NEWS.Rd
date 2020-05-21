
\name{NEWS}
\title{NEWS file for the trtf package}

\section{Changes in Version 0.3-5 (2019-02-01)}{
  \itemize{
    \item{Add R code reproducing empirical results for "Survival Forests under Test".}
  }
}

\section{Changes in Version 0.3-4 (2019-01-15)}{
  \itemize{
    \item{Respect \code{scale} argument from \code{object}.}
    \item{\code{traforest} broke in extremely rare cases due to incorrect
          dimensions of the score matrix.}
  }
}

\section{Changes in Version 0.3-3 (2018-09-18)}{
  \itemize{
    \item{All reparameterisation for max-type tests; this might
          give better (more powerful) variable selection.}
    \item{Allow fixed coefficients in trees and forests; this is
          useful for log-rank splitting, for example.}
  }
}


\section{Changes in Version 0.3-2 (2018-06-22)}{
  \itemize{
    \item{Document \code{predict} and \code{logLik} methods.}
    \item{Add \code{applyfun} and \code{cores} arguments to \code{predict}
          method for transformation forests (allowing parallel computation
          of predictions.}
  }
}

\section{Changes in Version 0.3-1 (2018-04-20)}{
  \itemize{
    \item{Some improvements in internal functions.}
    \item{Simulations with 500 (not 50) noise variables.}
    \item{Simplify likelihood evaluations based on new capabilities in
          \pkg{mlt}.}
    \item{New \code{update} argument for transformation forests.}
    \item{Update outputs (\code{scale = TRUE} in \pkg{partykit} 1.2-1).}
  }
}

\section{Changes in Version 0.3-0 (2018-01-04)}{
  \itemize{
    \item{The \pkg{trtf} package was first 
          published on CRAN.}
  }
}

