#' Tune \code{mtry} to the optimal value with respect to out-of-bag error for a LTRCCF model
#'
#' Starting with the default value of \code{mtry}, search for the optimal value
#' (with respect to Out-of-Bag error estimate) of \code{mtry} for \code{\link{ltrccf}}.
#'
#' @param formula a formula object, with the response being a \code{\link[survival]{Surv}}
#' object, with form
#'
#'
#' \code{Surv(tleft, tright, event)}.
#' @param data a data frame containing \code{n} rows of
#' left-truncated right-censored observations.
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
#' replacement; see the option \code{perturb} below);
#' (2) If \code{id} is not specified, it bootstraps the \code{data} by
#' sampling with or without replacement.
#' Regardless of the presence of \code{id}, if \code{"none"} is chosen, the
#' \code{data} is not bootstrapped at all. If \code{"by.user"} is choosen,
#' the bootstrap specified by \code{samp} is used.
#' @param samp Bootstrap specification when \code{bootstype = "by.user"}.
#' Array of dim \code{n x ntree} specifying how many times each record appears
#' inbag in the bootstrap for each tree.
#' @param mtryStart starting value of \code{mtry}; default is \code{sqrt(nvar)}.
#' @param stepFactor at each iteration, \code{mtry} is inflated (or deflated)
#' by this value. The default value is \code{2}.
#' @param na.action a function which indicates what should happen when the data contain
#' missing values.
#' @param samptype choices are \code{swor} (sampling without replacement) and
#' \code{swr} (sampling with replacement). The default action here is sampling
#' without replacement.
#' @param sampfrac a fraction, determining the proportion of subjects to draw
#' without replacement when \code{samptype = "swor"}. The default value is \code{0.632}.
#' To be more specific, if \code{id} is present, \code{0.632 * N} of subjects with their
#' pseudo-subject observations are drawn without replacement (\code{N} denotes the
#' number of subjects); otherwise, \code{0.632 * n} is the requested size
#' of the sample.
#' @param applyfun an optional \code{lapply}-style function with arguments
#' \code{function(X, FUN, ...)}.
#' It is used for computing the variable selection criterion. The default is to use the
#' basic \code{lapply} function unless the \code{cores} argument is specified (see below).
#' See \code{\link[partykit]{ctree_control}}.
#' @param cores numeric. If set to an integer the \code{applyfun} is set to
#' \code{\link[parallel]{mclapply}} with the desired number of cores.
#' See \code{\link[partykit]{ctree_control}}.
#' @param trace whether to print the progress of the search of the optimal value of \code{mtry}
#' when \code{mtry} is not specified (see \code{\link{tuneLTRCCF}}). \code{trace = TRUE}
#' is set by default.
#' @param control a list with control parameters, see \code{\link[partykit]{cforest}}.
#' The default values correspond to those of the default values used by \code{\link{ltrccf}}.
#' @param ntreeTry number of trees used at the tuning step.
#' @param trace whether to print the progress of the search. \code{trace = TRUE} is set by default.
#' @param plot whether to plot the out-of-bag error as a function of \code{mtry}.
#' @param doBest whether to run a \code{\link{ltrccf}} object using the optimal \code{mtry} found.
#' @param time.eval a vector of time points, at which the estimated survival probabilities
#' are evaluated.
#' @param time.tau an optional vector, with the \emph{i}-th entry giving the upper time limit for the
#' computed survival probabilities for the \emph{i}-th data (i.e., only compute
#' survival probabilies at \code{time.eval[time.eval <= time.tau[i]]} for the \emph{i}-th
#' data of interest).
#' @keywords mtry, out-of-bag errors, brier score
#' @return
#' If \code{doBest = FALSE} (default), this returns the optimal mtry value of those searched.
#' @return
#' If \code{doBest = TRUE}, this returns the \code{\link{ltrccf}} object produced with the optimal \code{mtry}.
#' @importFrom survival Surv
#' @importFrom graphics axis
#' @seealso \code{\link{sbrier_ltrc}} for evaluation of model fit when searching
#' for the optimal value of \code{mtry}.
#' @examples
#' ### Example with data pbcsample
#' Formula = Surv(Start, Stop, Event) ~ age + alk.phos + ast + chol + edema
#' ## mtry tuned by the OOB procedure with stepFactor 3, number of trees built 50.
#' mtryT = tune.ltrccf(formula = Formula, data = pbcsample, id = ID, stepFactor = 3,
#'                     ntreeTry = 50L, plot = TRUE)
#'
#'
#' @export

tune.ltrccf <- function(formula, data, id,
                        mtryStart = NULL, stepFactor = 2, ntreeTry = 100L,
                        bootstrap = c("by.sub", "by.root", "none", "by.user"),
                        samptype = c("swor","swr"),
                        sampfrac = 0.632,
                        samp = NULL,
                        na.action = na.pass,
                        time.eval = NULL, time.tau = NULL,
                        trace = TRUE, plot = FALSE, doBest = FALSE,
                        applyfun = NULL, cores = NULL,
                        control = partykit::ctree_control(teststat = "quad", testtype = "Univ",
                                                          mincriterion = 0, saveinfo = FALSE,
                                                          minsplit = max(ceiling(sqrt(nrow(data))), 20),
                                                          minbucket = max(ceiling(sqrt(nrow(data))), 7),
                                                          minprob = 0.01)) {

  # package version dependency
  # if (packageVersion("partykit") < "1.2.7") {
  #   stop("partykit >= 1.2.7 needed for this function.", call. = FALSE)
  # }

  Call <- match.call()
  Call[[1]] <- as.name('tuneLTRCCF')  #make nicer printout for the user
  # create a copy of the call that has only the arguments we want,
  #  and use it to call model.frame()
  indx <- match(c('formula', 'data', 'id'), names(Call), nomatch = 0)
  if (indx[1] == 0) stop("a formula argument is required")

  temp <- Call[c(1, indx)]
  temp[[1L]] <- quote(stats::model.frame)
  mf <- eval.parent(temp)

  Terms <- terms(formula)
  ord <- attr(Terms, 'order')
  if (length(ord) & any(ord != 1))
    stop("Interaction terms are not valid for this function")

  n <- nrow(mf)
  y <- model.extract(mf, 'response')
  id <- model.extract(mf, 'id')

  ## if not specified, the first one will be used as default
  bootstrap <- match.arg(bootstrap)
  samptype <- match.arg(samptype)

  # right-censored time from all observations
  Rtimes <- y[, 2L]

  # extract the x-variable names
  xvar.names <- attr(Terms, 'term.labels')
  nvar <- length(xvar.names)

  if (is.null(mtryStart)){
    mtryStart <- ceiling(sqrt(nvar))
  }

  if (length(id) == length(unique(id))){ # time-invariant LTRC data
    # it includes the case 1) when id = NULL, which is that id.seu is not specified
    #                      2) when id is specified, but indeed LTRC time-invariant
    id <- 1:n
    if (is.null(time.eval)){
      # estimated survival probabilities will be calculated at (a subset of) time.eval
      time.eval <- c(0, sort(unique(Rtimes)))
    }
  } else { # time-varying subject data
    id.sub <- unique(id)
    ## number of subjects
    n.sub <- length(id.sub)

    if (is.null(time.eval)){
      # estimated survival probabilities will be calculated at (a subset of) time.eval
      time.eval <- c(0, sort(unique(Rtimes)), seq(max(Rtimes), 1.5 * max(Rtimes), length.out = 50)[-1])
    }
    if (is.null(time.tau)){
      # For i-th data, estimated survival probabilities only calculated up time.tau[i]
      time.tau <- sapply(1:n.sub, function(ii){
        1.5 * max(Rtimes[id == id.sub[ii]])
      })
    }
  }
  # id is a vector, as requred by sbrier_ltrc function
  data$id <- id # this is a must, otherwise id cannot be passed to the next level
  # integrated Brier score of out-of-bag samples for a mtry value at test
  errorOOB_mtry <- function(eformula, edata, id,
                            emtryTest,
                            etpnt, etau,
                            entreeTry, econtrol,
                            ebootstrap,
                            esamptype,
                            esampfrac,
                            esamp,
                            ena.action, eapplyfun, ecores){
    cfOOB <- ltrccf(formula = eformula, data = edata, id = id,
                    mtry = emtryTest,
                    ntree = entreeTry,
                    control = econtrol,
                    bootstrap = ebootstrap,
                    samptype = esamptype,
                    sampfrac = esampfrac,
                    samp = esamp,
                    na.action = ena.action,
                    applyfun = eapplyfun,
                    cores = ecores)
    predOOB <- predict(object = cfOOB, time.eval = etpnt, time.tau = etau, OOB = TRUE)
    errorOOB <- sbrier_ltrc(obj = predOOB$survival.obj, id = id, pred = predOOB, type = "IBS")
    rm(cfOOB)
    rm(predOOB)
    gc()
    return(errorOOB)
  }

  # # errorOld
  errorOld <- errorOOB_mtry(eformula = formula, edata = data, id = id,
                            emtryTest = mtryStart,
                            etpnt = time.eval,
                            etau = time.tau,
                            entreeTry = ntreeTry,
                            econtrol = control,
                            ebootstrap = bootstrap,
                            esamptype = samptype,
                            esampfrac = sampfrac,
                            esamp = samp,
                            ena.action = na.action,
                            eapplyfun = applyfun,
                            ecores = cores)
  if (errorOld < 0) stop("Initial setting gave 0 error and no room for improvement.")
  if (trace) {
    cat("mtry =", mtryStart, " OOB Brier score =",
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
                                etpnt = time.eval,
                                etau = time.tau,
                                entreeTry = ntreeTry,
                                econtrol = control,
                                ebootstrap = bootstrap,
                                esamptype = samptype,
                                esampfrac = sampfrac,
                                esamp = samp,
                                ena.action = na.action,
                                eapplyfun = applyfun,
                                ecores = cores)

      if (trace) {
        cat("mtry =", mtryCur, "\tOOB error =", errorCur, "\n")
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
    res <- res_all
    plot(res_all, xlab = expression(m[try]), ylab = "OOB Error", type = "o", log = "x", xaxt = "n")
    axis(1, at=res_all[,"mtry"])
  }

  if (doBest)
    res <- ltrccf(formula = formula, data = data, id = id,
                  mtry = res, ntree = ntreeTry,
                  control = control,
                  bootstrap = bootstrap,
                  samptype = samptype,
                  sampfrac = sampfrac,
                  samp = samp,
                  na.action = na.action,
                  applyfun = applyfun,
                  cores = cores)

  return(res)
}
