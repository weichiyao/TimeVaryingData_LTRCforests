
\name{NEWS}
\title{NEWS file for the mlt package}

\section{Changes in Version 1.0-4 (2018-10-03)}{
  \itemize{
    \item{Arguments \code{K} and \code{cheat} where
          ignored by \code{confband} when \code{newdata}
          had multiple rows.}
    \item{Computation of starting values more robust now.}
  }
}

\section{Changes in Version 1.0-3 (2018-09-17)}{
  \itemize{
    \item{Add \code{lty} argument to \code{plot.ctm}.}
    \item{\code{update} needs free coefficients only.}
    \item{Internal interface changes.}
  }
}


\section{Changes in Version 1.0-2 (2018-06-22)}{
  \itemize{
    \item{Make sure transformation functions outside \code{bounds} are minus
          or plus \code{Inf}.}
    \item{Initial guestimates for ordered responses were incorrect
          and may potentially have led to nonsense results.}
    \item{Some smaller improvements in computation of log-likelihoods and
          scores with respect to accuracy and speed.}
    \item{\code{print} respects \code{options(digits)}.}
  }
}

\section{Changes in Version 1.0-1 (2018-04-20)}{
  \itemize{
    \item{\code{estfun, parm = coef(object, fixed = TRUE))} evaluates
          scores for all model parameters, including fixed ones.}
    \item{\code{logLik(..., newdata, w)} ignored weights \code{w} when
          \code{newdata} was given. Same problem was also fixed for
          \code{estfun}.}
  }
}


\section{Changes in Version 1.0-0 (2018-03-07)}{
  \itemize{
    \item{A paper describing version 1.0-0 of the \pkg{mlt}, \pkg{basefun}, and
          \pkg{variables} packages was accepted for publication in the Journal of
          Statistical Software 2018-03-05.}
    \item{Documentation updates.}
  }
}

\section{Changes in Version 0.2-2 (2018-02-23)}{
  \itemize{
    \item{Use \pkg{coneprog} for getting the starting values.}
    \item{\code{logLik} and \code{estfun} accept matrices as \code{parm}
          argument for the evalution of log-likelihoods and scores
          with subject-specific parameters (for example in transformation
          trees or forests and boosting procedures.}
  }
}

\section{Changes in Version 0.2-1 (2017-12-12)}{
  \itemize{
    \item{\code{q} is forwarded to \code{qmlt} by \code{predict.ctm} now.}
    \item{\code{p} is now \code{prob} in \code{qmlt} and thus \code{predict.ctm}.}
    \item{Update citation.}
  }
}


\section{Changes in Version 0.2-0 (2017-06-20)}{
  \itemize{
    \item{Most Likely Transformations will be published in
          the Scandinavian Journal of Statistics.}
    \item{Import package \pkg{alabama}.}
    \item{\code{as.Surv(R(Surv(...)))} returns \code{Surv(...)},
          useful for converting output by \code{simulate}
          to \code{Surv} objects.}
  }
}

\section{Changes in Version 0.1-3 (2017-04-04)}{
  \itemize{
    \item{Add \code{subset} argument to \code{update} (for faster transformation 
          trees and forests).}
    \item{Sum over score contributions with positive weight only when
          evaluating the gradient.}
    \item{Having all response observations being interval-censored is
          allowed again (too heavy checking was in place).}
    \item{Don't try to numerically check KKT conditions automatically.}
    \item{Check for unused arguments in dots where necessary.}
    \item{Make sure the score doesn't get too large (avoid division by near
          zero probabilities).}
    \item{Improve \code{survfit} to compute non-parametric unconditional
          probabilities for obtaining starting values in the presence of
          censoring and truncation.}
 }
}

\section{Changes in Version 0.1-2 (2017-02-14)}{
  \itemize{
    \item{\code{logLik} with \code{newdata} argument ignored \code{parm} 
          and \code{weights} arguments.}
    \item{\code{estfun} now also has a \code{newdata} argument.}
    \item{Correct axes labelling when plotting quantile functions.}
    \item{make sure names are correct in \code{coef(model, fixed = FALSE)}.}
    \item{check if any exact or interval-censored response with non-zero 
          weight exists before trying to fit the model.}
    \item{make checks a little more robust against huge diffs under Windows.}
 }
}

\section{Changes in Version 0.1-1 (2016-11-29)}{
  \itemize{
    \item{Fix two bugs in computation of log-likelihood for possibly
          left-truncated responses such as \code{Surv(start, time, status)}.}
 }
}

\section{Changes in Version 0.1-0 (2016-11-10)}{
  \itemize{
    \item{Add augmented lagrangian minimization (\code{auglag()} from package \pkg{alabama}).}
    \item{Make optimiation procedure more general and adaptive, allow users to change
          defaults and even add their own optimiser.}
 }
}


\section{Changes in Version 0.0-34 (2016-09-13)}{
  \itemize{
    \item{fix bug when calling \code{survfit} for computing initial probabilities.}
    \item{add \code{bysim} argument to \code{simulate}.}
    \item{make sure \code{checkGrad} is respected by \code{update}.}
    \item{\code{predict} computes \code{q} with \code{K} elements 
          if not given (as \code{plot} always did).}
 }
}

\section{Changes in Version 0.0-33 (2016-05-25)}{
  \itemize{
    \item{Make sure \code{times} are ordered before calling
          \code{survival::summary.survfit}.}
  }
}

\section{Changes in Version 0.0-32 (2016-05-24)}{
  \itemize{
    \item{Introduce \code{as.mlt} generic.}
    \item{Introduce a \code{coef} slot in \code{ctm} objects
          and a corresponding \code{coef<-} and \code{coef}
          method for setting and extracting coefficients
          to and from unfitted conditional transformation models.}
    \item{\code{predict}, \code{simulate} and \code{plot}
          work on \code{ctm} objects (with meaningful coefficients)
          now.}
  }
}

\section{Changes in Version 0.0-31 (2016-03-11)}{
  \itemize{
    \item{Some small improvements wrt run time and memory consumption.}
    \item{Use \code{theta = coef(object)} as default starting parameters in \code{update()}.}
    \item{\code{logLik} has a new \code{newdata} argument.}
    \item{\code{simulate} has a new \code{q} argument.}
  }
}

\section{Changes in Version 0.0-30 (2016-02-19)}{
  \itemize{
    \item{The \pkg{mlt} package was first 
          published on CRAN.}
  }
}

