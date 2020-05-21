
\name{NEWS}
\title{NEWS file for the tram package}

\section{Changes in Version 0.2-5 (2019-02-25)}{
  \itemize{
    \item{\code{tram} objects are now in line with \code{multcomp::glht},
          allowing post-hoc multiple comparisons.}
  }
}

\section{Changes in Version 0.2-4 (2019-01-22)}{
  \itemize{
    \item{\code{log_first = TRUE} is now accepted by \code{tram} and all 
          model functions calling \code{tram}. This is still experimental.}
    \item{Silence regression tests on Windows.}
  }
}


\section{Changes in Version 0.2-3 (2018-09-17)}{
  \itemize{
    \item{New \code{extrapolate = TRUE} option for all models constraints
          the second derivate of the Bernstein transformation function to be
          zero at \code{support}. This helps to get rid of possible
          discontinuous densities at \code{support}.}
    \item{Add \code{Lehmann} for fitting shift transformation models those
          linear predictor can be interpreted as log-Lehmann-alternative.}
    \item{Add \code{method = "loglog"} to \code{Polr}.}
    \item{New \code{update} method.}
  }
}

\section{Changes in Version 0.2-2 (2018-06-22)}{
  \itemize{
    \item{Use default \code{bounds = c(0, Inf)} for \code{Surv} objects (for example
          in \code{Coxph}.}
    \item{Update to \pkg{mlt} 1.0-2.}
  }
}


\section{Changes in Version 0.2-1 (2018-04-20)}{
  \itemize{
    \item{Likelihood ratio tests with fixed parameters gave an error.}
    \item{Remove strata from linear predictor (experimental).}
    \item{Try to implement a better guess for the \code{add} argument
          later defining the axes corresponding to the response
          in plots.}
    \item{All models accept \code{add} and \code{bounds} arguments.}
  }
}

\section{Changes in Version 0.2-0 (2018-02-25)}{
  \itemize{
    \item{\pkg{tram} first published on CRAN.}
  }
}
