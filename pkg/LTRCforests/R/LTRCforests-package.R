#' Constructs forest methods for left-truncated and right-censored (LTRC) survival data
#'
#' Constructs a LTRC conditional inference forest (LTRCCIF) or
#' a LTRC relative risk forest (LTRCRRF) for left-truncated and right-censored data,
#' it also allows for (left-truncated) right-censored survival data with
#' time-varying covariates (Yao et al. 2022).
#' The main functions of this package are \code{\link{ltrccif}} and \code{\link{ltrcrrf}}.
#'
#' \subsection{Problem setup and existing methods}{
#' Continuous-time survival data with time-varying covariates are common in practice.
#' Methods like the Cox proportional hazards model rely on restrictive assumptions such as
#' proportional hazards and a log-linear relationship between the hazard function and
#' covariates. Furthermore, because these methods are often parametric, nonlinear effects
#' of variables must be modeled by transformations or expanding the design matrix to
#' include specialized basis functions for more complex data structures in real world
#' applications. The functions \code{\link[LTRCtrees]{LTRCIT}} and
#' \code{\link[LTRCtrees]{LTRCART}} provide a conditional inference tree method and 
#' a relative risk tree method for
#' left-truncated right-censored survival data, which also allows for right-censored
#' survival data with time-varying covariates. Tree estimators are nonparametric and 
#' as such often exhibit
#' low bias and high variance. Ensemble methods like bagging and random forest can
#' reduce variance while preserving low bias. 
#' The most popular survival forest methods, including conditional inference forest
#' (see \code{\link[partykit]{cforest}}), relative risk forest, and random survival 
#' forest method (see \code{\link[randomForestSRC]{rfsrc}}) can only be applied to 
#' right-censored survival data with time-invariant covariates.}
#'
#' \subsection{LTRC forests}{
#' This package implements \code{\link{ltrccif}} and \code{\link{ltrcrrf}}.
#' \code{\link{ltrccif}} extends the conditional inference forest
#' (see \code{\link[partykit]{cforest}}) to LTRC survival data.
#' It uses LTRC conditional inference survival trees
#' as base learners.
#' \code{\link{ltrcrrf}} extends the relative risk forest
#' (Ishwaran et al. 2004) to left-truncated right-censored survival data.
#' It uses LTRC risk relative tree
#' as base learners.
#' The main functions \code{\link{ltrccif}} and \code{\link{ltrcrrf}} 
#' fit a corresponding LTRC forest for LTRC data, with parameter
#' \code{mtry} tuned by \code{\link{tune.ltrccif}} or \code{\link{tune.ltrcrrf}}. This tuning
#' procedure relies on the evaluation of the out-of-bag errors, which is performed by the
#' function \code{\link{sbrier_ltrc}}. \code{\link{print}}
#' prints summary output for \code{ltrccif} objects and \code{ltrcrrf} objects.
#' \code{\link{predictProb}}
#' constructs survival function estimates for \code{ltrccif} objects and \code{ltrcrrf} objects.
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
#' If one is in the traditional case with right-censored data
#' and time-invariant covariates, however, then it is recommended to use
#' the functions \code{\link[partykit]{cforest}} and \code{\link[randomForestSRC]{rfsrc}}
#' directly to construct conditional inference forests and random survival forests,
#' respectively.
#'
#' @references Yao, W., Frydman, H., Larocque, D. and Simonoff, J. S. (2022). 
#' Ensemble methods for survival function estimation with time-varying covariates.
#' \emph{Statistical Methods in Medical Research}, \strong{31}(11):2217-2236.
#' @references Andersen, P. and Gill, R. (1982). Cox’s regression model for counting
#' processes, a large sample study. \emph{Annals of Statistics}, \strong{10}:1100-1120.
#' @references Ishwaran, H., Blackstone, E. H., Pothier, C., and Lauer, M. S. (2004).
#' Relative risk forests for exercise heart rate recovery as a predictor of mortality.
#' \emph{Journal of the American StatisticalAssociation}, \strong{99}(1):591–600.
#' @references Fu, W. and Simonoff, J. S. (2016). Survival trees for left-truncated and 
#' right-censored data, with application to time-varying covariate data. 
#' \emph{Biostatistics}, \strong{18}(2):352–369.
#' @seealso \code{\link{ltrccif}}, \code{\link{ltrcrrf}},
#' \code{\link{predictProb}}, \code{\link{print}},
#' \code{\link{tune.ltrccif}}, \code{\link{tune.ltrcrrf}}, \code{\link{sbrier_ltrc}}.
#' @docType package
#' @name LTRCforests-package
NULL
