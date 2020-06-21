#' Constructs forest methods for left-truncated and right-censored (LTRC) survival data
#'
#' Constructs a LTRC conditional inference forest (LTRCCF) or
#' a LTRC random survival forest (LTRCRSF) for left-truncated and right-censored data,
#' it also allows for (left-truncated) right-censored survival data with
#' time-varying covariates.
#' The main functions of this package are \code{\link{ltrccf}} and \code{\link{ltrcrsf}}.
#'
#' \subsection{Problem setup and existing methods}{
#' Continuous-time survival data with time-varying covariates are common in practice.
#' Methods like the Cox proportional hazards model rely on restrictive assumptions such as
#' proportional hazards and a log-linear relationship between the hazard function and
#' covariates. Furthermore, because these methods are often parametric, nonlinear effects
#' of variables must be modeled by transformations or expanding the design matrix to
#' include specialized basis functions for more complex data structures in real world
#' applications. The functions \code{\link[LTRCtrees]{LTRCIT}} and
#' \code{\link[LTRCtrees]{LTRCART}} provide a conditional inference tree method and a survival tree method for
#' left-truncated right-censored survival data, which also allows for right-censored
#' survival data with time-varying covariates. Tree estimators are nonparametric and as such often exhibit
#' low bias and high variance. Ensemble methods like bagging and random forest can
#' reduce variance while preserving low bias.}
#'
#' \subsection{LTRC forests}{
#' This package implements \code{\link{ltrccf}} and \code{\link{ltrcrsf}}.
#' \code{\link{ltrccf}} extends the conditional inference forest
#' (see \code{\link[partykit]{cforest}}) to LTRC survival data.
#' It uses LTRC conditional inference survival trees
#' (see \code{\link[LTRCtrees]{LTRCIT}}) as base learners.
#' \code{\link{ltrcrsf}} extends the random survival forest
#' (see \code{\link[randomForestSRC]{rfsrc}}) to left-truncated right-censored survival data.
#' It uses LTRC survival trees with Poisson splitting rule
#' (see \code{\link[LTRCtrees]{LTRCART}}) as base learners.
#' The main functions \code{\link{ltrccf}} and \code{\link{ltrcrsf}} fit a corresponding
#' LTRC forest for LTRC data, with parameter
#' \code{mtry} tuned by \code{\link{tune.ltrccf}} or \code{\link{tune.ltrcrsf}}. This tuning
#' procedure relies on the evaluation of the out-of-bag errors, which is performed by the
#' function \code{\link{sbrier_ltrc}}. \code{\link{print}}
#' prints summary output for \code{ltrccf} objects and \code{ltrcrsf} objects.
#' \code{\link{predictProb}}
#' constructs survival function estimates for \code{ltrccf} objects and \code{ltrcrsf} objects.
#'
#' For (left-truncated) right-censored survival data with time-varying covariates,
#' one can first reformat the data structure to one with LTRC observations,
#' where the multiple records of a subject become a list of pseudo-subjects and
#' are treated independently. This procedure is usually referred to as the
#' Andersen-Gill method (Andersen and Gill, 1982). Then LTRC forest methods
#' can be applied on this reformatted dataset.
#' }
#'
#' Overall, the methods in this package can handle all combinations of left truncation,
#' right censoring, time-invariant covariates, and time-varying covariates.
#' If one is in the traditional case with right censored data
#' and time-invariant covariates, however, then it is recommended to use
#' the functions \code{\link[partykit]{cforest}} and \code{\link[randomForestSRC]{rfsrc}}
#' directly to construct conditional inference forests and random survival forests,
#' respectively.
#'
#' @references Andersen, P. and Gill, R. (1982). Cox’s regression model for counting
#' processes, a large sample study. \emph{Annals of Statistics}, \strong{10}, 1100-1120.
#' @seealso \code{\link{ltrccf}}, \code{\link{ltrcrsf}},
#' \code{\link{predictProb}}, \code{\link{print}},
#' \code{\link{tune.ltrccf}}, \code{\link{tune.ltrcrsf}}, \code{\link{sbrier_ltrc}}
#' @docType package
#' @name LTRCforests-package
NULL
