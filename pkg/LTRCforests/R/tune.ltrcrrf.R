#' Tune \code{mtry} to the optimal value with respect to out-of-bag error for a LTRCRRF model
#'
#' Starting with the default value of \code{mtry}, search for the optimal value
#' (with respect to out-of-bag error estimate) of \code{mtry} for \code{\link{ltrcrrf}}.
#'
#' @param formula a formula object, with the response being a \code{\link[survival]{Surv}}
#' object, with form
#'
#'
#' \code{Surv(tleft, tright, event)}.
#' @param data a a data frame containing \code{n} rows of
#' left-truncated right-censored observations.
#' @param mtryStart starting value of \code{mtry}; default is \code{sqrt(nvar)}.
#' @param id variable name of subject identifiers. If this is present, it will be
#' searched for in the \code{data} data frame. Each group of rows in \code{data}
#' with the same subject \code{id} represents the covariate path through time of
#' a single subject. If not specified, the algorithm then assumes \code{data}
#' contains left-truncated and right-censored survival data with time-invariant
#' covariates.
#' @param bootstrap bootstrap protocol.
#' (1) If \code{id} is present,
#' the choices are: \code{"by.sub"} (by default) which bootstraps subjects,
#' \code{"by.root"} which bootstraps pseudo-subjects.
#' Both can be with or without replacement (by default sampling is without
#' replacement; see the option \code{samptype} below).
#' (2) If \code{id} is not specified, the default is \code{"by.root"} which
#' bootstraps the \code{data} by sampling with or without replacement;
#' if \code{"by.node"} is choosen, \code{data} is bootstrapped with replacement
#' at each node while growing the tree.
#' Regardless of the presence of \code{id}, if \code{"none"} is chosen, the
#' \code{data} is not bootstrapped at all. If \code{"by.user"} is choosen,
#' the bootstrap specified by \code{samp} is used.
#' @param samptype choices are \code{swor} (sampling without replacement) and
#' \code{swr} (sampling with replacement). The default action here is sampling
#' without replacement.
#' @param samp Bootstrap specification when \code{bootstype = "by.user"}.
#' Array of dim \code{n x ntree} specifying how many times each record appears
#' in each bootstrap sample.
#' @param trace whether to print the progress of the search of the optimal value
#' of \code{mtry}.
#' \code{trace = TRUE} is set by default.
#' @param stepFactor at each iteration, \code{mtry} is inflated (or deflated)
#' by this value, used when \code{mtry} is not specified (see \code{\link{ltrcrrf}}).
#' The default value is \code{2}.
#' @param nodedepth maximum depth to which a tree should be grown. The default behaviour
#' is that this parameter is ignored.
#' @param nsplit an non-negative integer value for number of random splits to consider
#' for each candidate splitting variable. This significantly increases speed.
#' When zero or \code{NULL}, the algorithm uses much slower deterministic splitting
#' where all possible splits are considered. \code{nsplit = 10L} by default.
#' @param sampfrac a fraction, determining the proportion of subjects to draw
#' without replacement when \code{samptype = "swor"}. The default value is \code{0.632}.
#' To be more specific, if \code{id} is present, \code{0.632 * N} of subjects with their
#' pseudo-subject observations are drawn without replacement (\code{N} denotes the
#' number of subjects); otherwise, \code{0.632 * n} is the requested size
#' of the sample.
#' @param na.action action taken if the data contains \code{NA}’s. The default
#' \code{"na.omit"} removes the entire record if any of its entries is
#' \code{NA} (for x-variables this applies only to those specifically listed
#' in \code{formula}). See function \code{\link[randomForestSRC]{rfsrc}} for
#' other available options.
#' @param ntime an integer value used for survival to constrain ensemble calculations
#' to a grid of \code{ntime} time points. Alternatively if a vector of values
#' of length greater than one is supplied, it is assumed these are the time points
#' to be used to constrain the calculations (note that the constrained time points
#' used will be the observed event times closest to the user supplied time points).
#' If no value is specified, the default action is to use all observed event times.
#' @param trace whether to print the progress of the search of the optimal value of \code{mtry}
#' when \code{mtry} is not specified (see \code{\link{tuneLTRCCF}}). \code{trace = TRUE}
#' is set by default.
#' @param nodesizeTry forest average terminal node size used at the tuning step.
#' @param ntreeTry number of trees used at the tuning step.
#' @param trace whether to print the progress of the search. \code{trace = TRUE} is set by default.
#' @param plot whether to plot the out-of-bag error as a function of \code{mtry}.
#' \code{plot = FALSE} is set by default.
#' @param doBest whether to run a \code{\link{ltrcrrf}} object using the optimal \code{mtry} found.
#' \code{doBest = FALSE} is set by default.
#' @param time.eval a vector of time points, at which the estimated survival probabilities
#' are evaluated.
#' @param time.tau an optional vector, with the \emph{i}-th entry giving the upper time limit for the
#' computed survival probabilities for the \emph{i}-th data (i.e., only computes
#' survival probabilies at \code{time.eval[time.eval <= time.tau[i]]} for the \emph{i}-th
#' data of interest).
#' @return
#' If \code{doBest = FALSE} (default), this returns the optimal mtry value of those searched.
#' @return
#' If \code{doBest = TRUE}, this returns the \code{\link{ltrcrrf}} object produced with the optimal \code{mtry}.
#' @importFrom survival Surv
#' @importFrom graphics axis
#' @import stats
#' @import utils
#' @import prodlim
#' @seealso \code{\link{sbrier_ltrc}} for evaluation of model fit for the optimal value of \code{mtry}.
#' @examples
#' ### Example with data pbcsample
#' library(survival)
#' Formula = Surv(Start, Stop, Event) ~ age + alk.phos + ast + chol + edema
#' ## mtry tuned by the OOB procedure with stepFactor 3, number of trees built 10.
#' mtryT = tune.ltrcrrf(formula = Formula, data = pbcsample, stepFactor = 3,
#'                      ntreeTry = 10L)
#' @export

tune.ltrcrrf <- function(formula, data, id,
                         mtryStart = NULL, stepFactor = 2,
                         time.eval = NULL, time.tau = NULL,
                         ntreeTry = 100L,
                         bootstrap = c("by.sub","by.root","by.node","by.user","none"),
                         samptype = c("swor","swr"),
                         sampfrac = 0.632,
                         samp = NULL,
                         na.action = "na.omit",
                         trace = TRUE,
                         doBest = FALSE,
                         plot = FALSE,
                         ntime,
                         nsplit = 10L,
                         nodesizeTry = max(ceiling(sqrt(nrow(data))), 15),
                         nodedepth = NULL) {
  Call <- match.call()
  # Call[[1]] <- as.name('tuneltrcrrf')  #make nicer printout for the user
  # create a copy of the call that has only the arguments we want,
  #  and use it to call model.frame()
  indx <- match(c('formula', 'id'), names(Call), nomatch = 0)
  if (indx[1] == 0) stop("a formula argument is required")

  # temp <- Call[c(1, indx)]
  # temp[[1L]] <- quote(stats::model.frame)
  # mf <- eval.parent(temp)

  # Terms <- terms(formula)
  # ord <- attr(Terms, 'order')
  # if (length(ord) & any(ord != 1))
  #   stop("Interaction terms are not valid for this function")

  # pull y-variable names
  yvar.names <- all.vars(formula(paste(as.character(formula)[2], "~ .")), max.names = 1e7)
  yvar.names <- yvar.names[-length(yvar.names)]

  if (length(yvar.names) == 4){
    yvar.names = yvar.names[2:4]
  }
  n <- nrow(data)

  bootstrap <- match.arg(bootstrap)
  samptype <- match.arg(samptype)
  # na.action <- match.arg(na.action)

  # right-censored time from all observations
  Rtimes <- data[, yvar.names[2]]

  # extract the x-variable names
  xvar.names <- attr(terms(formula), 'term.labels')
  nvar <- length(xvar.names)

  if (is.null(mtryStart)){
    mtryStart <- ceiling(sqrt(nvar))
  }

  if (missing(ntime)){
    ##unique death times
    ntime = sort(unique(Rtimes[data[, yvar.names[3]] == 1]))
  }

  # this is a must, otherwise id cannot be passed to the next level in tune.ltrccf
  if (indx[2] == 0){
    ## If id is not present, then we add one more variable
    # mf$`(id)` <- 1:nrow(mf) ## No relabel, due to missing value problem do not need this for tuning output
    data$id <- 1:n # this is a must, otherwise id cannot be passed to the next level
  } else {
    ## If id is present, then we rename the column to be id
    names(data)[names(data) == deparse(substitute(id))] <- "id" # this is a must, otherwise id cannot be passed to the next level
  }

  data <- data[, c("id", yvar.names, xvar.names)]

  if (na.action == "na.omit") {
    takeid = which(complete.cases(data) == 1)
  } else if (na.action == "na.impute") {
    takeid = 1:n
  } else {
    stop("na.action can only be either 'na.omit' or 'na.pass'.")
  }

  id.sub <- unique(data$id[takeid])
  n.seu <- length(takeid)
  ## number of subjects
  n.sub <- length(id.sub)

  Rtimes <- Rtimes[takeid]
  ## This is to determine time.tau and time.eval, so is different from ltrccf.R
  if (n.seu == n.sub){ # time-invariant LTRC data
    # it includes the case 1) when id = NULL, which is that id.seu is not specified
    #                      2) when id is specified, but indeed LTRC time-invariant
    if (is.null(time.eval)){
      # estimated survival probabilities will be calculated at (a subset of) time.eval
      time.eval <- c(0, sort(unique(Rtimes)))
    }
    if (bootstrap == "by.sub") bootstrap = "by.root"
  } else { # time-varying subject data
    if (is.null(time.eval)){
      # estimated survival probabilities will be calculated at (a subset of) time.eval
      time.eval <- c(0, sort(unique(Rtimes)), seq(max(Rtimes), 1.5 * max(Rtimes), length.out = 50)[-1])
    }
    if (is.null(time.tau)){
      # For i-th data, estimated survival probabilities only calculated up time.tau[i]
      time.tau <- sapply(1:n.sub, function(ii){
        1.5 * max(Rtimes[data$id[takeid] == id.sub[ii]])
      })
    }
  }

  # integrated Brier score of out-of-bag samples for a mtry value at test
  errorOOB_mtry <- function(eformula, edata, id, emtryTest,
                            etpnt, etau,
                            entreeTry, enodesizeTry, enodedepth,
                            ensplit,
                            ebootstrap,
                            esamptype,
                            esampfrac,
                            esamp,
                            ena.action,
                            entime){
    rrfOOB <- ltrcrrf(formula = eformula, data = edata, id = id,
                      mtry = emtryTest,
                      ntree = entreeTry,
                      nodesize = enodesizeTry,
                      nodedepth = enodedepth,
                      nsplit = ensplit,
                      bootstrap = ebootstrap,
                      samptype = esamptype,
                      sampfrac = esampfrac,
                      samp = esamp,
                      na.action = ena.action,
                      ntime = entime)
    predOOB <- predictProb(object = rrfOOB, time.eval = etpnt, time.tau = etau, OOB = TRUE)
    errorOOB <- sbrier_ltrc(obj = predOOB$survival.obj, id = predOOB$survival.id,
                            pred = predOOB, type = "IBS")
    rm(rrfOOB)
    rm(predOOB)
    return(errorOOB)
  }
  # errorOld
  errorOld <- errorOOB_mtry(eformula = formula, edata = data, id = id,
                            emtryTest = mtryStart,
                            etpnt = time.eval, etau = time.tau,
                            entreeTry = ntreeTry,
                            enodesizeTry = nodesizeTry,
                            enodedepth = nodedepth,
                            ensplit = nsplit,
                            ebootstrap = bootstrap,
                            esamptype = samptype,
                            esampfrac = sampfrac,
                            esamp = samp,
                            ena.action = na.action,
                            entime = ntime)

  if (errorOld < 0) stop("Initial setting gave 0 error and no room for improvement.")
  if (trace) {
    cat("mtry = ", mtryStart, " OOB Brier score = ",
        errorOld, "\n")
  }

  oobError <- list()
  oobError[[1]] <- errorOld
  names(oobError)[1] <- mtryStart

  for (direction in c("left", "right")) {
    if (trace) cat("Searching", direction, "...\n")
    mtryCur <- mtryStart
    while (mtryCur != nvar) {
      mtryOld <- mtryCur
      mtryCur <- if (direction == "left") {
        max(1, ceiling(mtryCur / stepFactor))
      } else {
        min(nvar, floor(mtryCur * stepFactor))
      }
      if (mtryCur == mtryOld) break
      errorCur <- errorOOB_mtry(eformula = formula, edata = data, id = id,
                                emtryTest = mtryCur,
                                etpnt = time.eval, etau = time.tau,
                                entreeTry = ntreeTry,
                                enodesizeTry = nodesizeTry,
                                enodedepth = nodedepth,
                                ensplit = nsplit,
                                ebootstrap = bootstrap,
                                esamptype = samptype,
                                esampfrac = sampfrac,
                                esamp = samp,
                                ena.action = na.action,
                                entime = ntime)

      if (trace) {
        cat("mtry = ", mtryCur, "\tOOB error = ", errorCur, "\n")
      }
      oobError[[as.character(mtryCur)]] <- errorCur
      errorOld <- errorCur
    }
  }
  mtry <- sort(as.numeric(names(oobError)))
  res_all <- unlist(oobError[as.character(mtry)])
  res_all <- cbind(mtry = mtry, OOBError = res_all)
  res <- res_all[which.min(res_all[, 2]), 1]

  if (plot) {
    res = res_all
    plot(res_all, xlab = expression(m[try]), ylab = "OOB Error", type = "o", log = "x", xaxt = "n")
    axis(1, at = res_all[, "mtry"])
  }

  if (doBest) {
    res <- ltrcrrf(formula = formula, data = data, id = id,
                   mtry = res,
                   ntree = ntreeTry,
                   nodesize = nodesizeTry,
                   nodedepth = nodedepth,
                   nsplit = nsplit,
                   ntime = ntime,
                   bootstrap = bootstrap,
                   samptype = samptype,
                   sampfrac = sampfrac,
                   samp = samp,
                   na.action = na.action)
  }

  return(res)
}


