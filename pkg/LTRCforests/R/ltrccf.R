.logrank_trafo2 <- function(x2){
  if (sum(x2[, 3] == 1) == 0) {
    result <- x2[,3]
  } else {
    unique.times <- unique(x2[,2][which(x2[, 3] == 1)])
    D <- rep(NA, length(unique.times))
    R <- rep(NA, length(unique.times))

    for(j in 1:length(unique.times)){
      D[j] = sum(unique.times[j] == x2[, 2])
    }

    for(k in 1:length(unique.times) ){
      value <- unique.times[k]
      R[k] <- sum(apply(x2[, 1:2], 1, function(interval){interval[1] < value & value <= interval[2]}))
    }

    Ratio <- D / R

    Ratio <- Ratio[order(unique.times)]
    Nelson.Aalen <- cumsum(Ratio)
    Event.time <- unique.times[order(unique.times)]
    Left <- sapply(x2[, 1], function(t){if(t < min(Event.time)) return(0) else return(Nelson.Aalen[max(which(Event.time <= t))])})
    Right <- sapply(x2[, 2], function(t){if(t < min(Event.time)) return(0) else return(Nelson.Aalen[max(which(Event.time <= t))])})

    result<- x2[, 3] - (Right - Left)
  }
  return(as.double(result))
}

#' Fit a LTRC conditional inference forest
#'
#' An implementation of the random forest and bagging ensemble algorithms utilizing
#' LTRC conditional inference trees \code{\link{LTRCIT}} as base learners for
#' left-truncated right-censored survival data with time-invariant covariates.
#' It also allows for (left-truncated) right-censored survival data with
#' time-varying covariates.
#'
#'
#' This function extends the conditional inference survival forest algorithm in
#' \code{\link[partykit]{cforest}} to fit left-truncated and right-censored data,
#' which allow for time-varying covariates. The traditional survival forests in
#' \code{\link[partykit]{cforest}} only applies for right-censored data with
#' time-invariant invariant covariates.
#'
#' @param formula a formula object, with the response being a \code{\link[survival]{Surv}}
#' object, with form
#'
#'
#' \code{Surv(tleft, tright, event)}.
#' @param data a data frame containing \code{n} rows of
#' left-truncated right-censored observations.
#' For time-varying data, this should be
#' a data frame containing pseudo-subject observations based on the Andersen-Gill
#' reformulation.
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
#' (2) If \code{id} is not specified, it
#' bootstraps the \code{data} by sampling with or without replacement.
#' Regardless of the presence of \code{id}, if \code{"none"} is chosen, the
#' \code{data} is not bootstrapped at all. If \code{"by.user"} is choosen,
#' the bootstrap specified by \code{samp} is used.
#' @param samp Bootstrap specification when \code{bootstype = "by.user"}.
#' Array of dim \code{n x ntree} specifying how many times each record appears
#' inbag in the bootstrap for each tree.
#' @param na.action a function which indicates what should happen when the data contain
#' missing values.
#' @param mtry number of input variables randomly sampled as candidates at each node for
#' random forest like algorithms. The default \code{mtry} is tuned by \code{\link{tune.ltrccf}}.
#' @param ntree an integer, the number of the trees to grow for the forest.
#' \code{ntree = 100L} is set by default.
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
#' @param trace whether to print the progress of the search of the optimal value of
#' \code{mtry}, when \code{mtry} is not specified (see \code{\link{tune.ltrccf}}).
#' \code{trace = TRUE} is set by default.
#' @param stepFactor at each iteration, \code{mtry} is inflated (or deflated)
#' by this value, used when \code{mtry} is not specified (see \code{\link{tune.ltrccf}}).
#' The default value is \code{2}.
#' @param control a list of control parameters, see \code{\link[partykit]{ctree_control}}.
#' \code{control} parameters \code{minsplit}, \code{minbucket} have been adjusted from the
#' \code{\link[partykit]{cforest}} defaults. Other default values correspond to those of the
#' default values used by \code{\link[partykit]{ctree_control}}.
#' @keywords Ensemble method, conditional inference forest, left-truncated right-censored data,
#' time-varying covariate data
#' @return An object belongs to the class \code{ltrccf}, as a subclass of
#' \code{\link[partykit]{cforest}}.
#' @import partykit
#' @import survival
#' @import stats
#' @import utils
#' @import prodlim
#' @seealso \code{\link{predict.ltrccf}} for prediction and \code{\link{tune.ltrccf}}
#' for \code{mtry} tuning.
#' @references Andersen, P. and Gill, R. (1982). Cox's regression model for counting
#' processes, a large sample study. \emph{Annals of Statistics}, \strong{10}, 1100-1120.
#' @examples
#' #### Example with time-varying data pbcsample
#' Formula = Surv(Start, Stop, Event) ~ age + alk.phos + ast + chol + edema
#' ## Fit an LTRCCF on the time-varying data, with mtry specified
#' LTRCCFobj = ltrccf(formula = Formula, data = pbcsample, id = ID, mtry = 3, ntree = 50L)
#'
#' ## Fit an LTRCCF on the time-invariant data, with mtry tuned with stepFactor = 3.
#' LTRCCFobj = ltrccf(formula = Formula, data = pbcsample, ntree = 50L, stepFactor = 3)
#' @export

ltrccf <- function(formula, data, id,
                   mtry = NULL, ntree = 100L,
                   bootstrap = c("by.sub","by.root","by.user","none"),
                   samptype = c("swor","swr"),
                   sampfrac = 0.632,
                   samp = NULL,
                   trace = TRUE, stepFactor = 2,
                   na.action = na.pass,
                   applyfun = NULL, cores = NULL,
                   control = partykit::ctree_control(teststat = "quad", testtype = "Univ",
                                                     minsplit = max(ceiling(sqrt(nrow(data))), 20),
                                                     minbucket = max(ceiling(sqrt(nrow(data))), 7),
                                                     minprob = 0.01,
                                                     mincriterion = 0, saveinfo = FALSE)){
  # package version dependency
  if (packageVersion("partykit") < "1.2.7") {
    stop("partykit >= 1.2.7 needed for this function.", call. = FALSE)
  }

  requireNamespace("inum")
  Call <- match.call()
  Call[[1]] <- as.name('LTRCCF')  #make nicer printout for the user
  # create a copy of the call that has only the arguments we want,
  #  and use it to call model.frame()
  indx <- match(c('formula', 'data', 'id'),
                names(Call), nomatch=0)
  if (indx[1] == 0) stop("a formula argument is required")
  Call$formula <- eval(formula)
  # print(id)
  temp <- Call[c(1, indx)]
  temp[[1L]] <- quote(stats::model.frame)
  mf <- eval.parent(temp)

  Terms <- terms(formula)
  ord <- attr(Terms, 'order')
  if (length(ord) & any(ord !=1))
    stop("Interaction terms are not valid for this function")

  n <- nrow(mf)
  y <- model.extract(mf, 'response')

  if (!is.Surv(y)) stop("Response must be a survival object")
  if (!attr(y, "type") == "counting") stop("The Surv object must be of type 'counting'.")

  ## if not specified, the first one will be used as default
  bootstrap <- match.arg(bootstrap)
  samptype <- match.arg(samptype)

  id <- model.extract(mf, 'id')
  if (is.null(id)){ # If present, do not want to relabel them
    id <- 1:n
    mf$`(id)` <- id
  }
  ## bootstrap case
  if (length(id) == length(unique(id))){ # time-invariant LTRC data
    # it includes the case 1) when id = NULL, which is that id is not specified
    #                      2) when id is specified, but indeed LTRC time-invariant
    if (bootstrap == "by.sub") bootstrap <- "by.root"
  } else { # time-varying subject data
    id.sub <- unique(id)
    ## number of subjects
    n.sub <- length(id.sub)
  }

  if (samptype == "swor"){
    perturb = list(replace = FALSE, fraction = sampfrac)
  } else if (samptype == "swr"){
    perturb = list(replace = TRUE)
  } else {
    stop("samptype must set to be either 'swor' or 'swr'\n")
  }
  if (bootstrap == "by.sub"){
    size <- n.sub
    if (!perturb$replace) size <- floor(n.sub * perturb$fraction)
    samp <- replicate(ntree,
                      sample(id.sub, size = size,
                             replace = perturb$replace),
                      simplify = FALSE) # a list of length ntree
    samp <- lapply(samp, function(y) unlist(sapply(y, function(x) which(id %in% x), simplify = FALSE)))
    samp <- sapply(samp, function(x) as.integer(tabulate(x, nbins = length(id)))) # n x ntree
  } else if (bootstrap == "none"){
    samp <- matrix(1, nrow = n, ncol = ntree)
  } else if (bootstrap == "by.user") {
    if (is.null(samp)) {
      stop("samp must not be NULL when bootstrapping by user")
    }
    if (is.matrix(samp)){
      if (!is.matrix(samp)) stop("samp must be a matrx")
      if (any(!is.finite(samp))) stop("samp must be finite")
      if (any(samp < 0)) stop("samp must be non-negative")
      if (all(dim(samp) != c(n, ntree))) stop("dimension of samp must be n x ntree")
      samp <- as.matrix(samp)  # transform into matrix
    }
  } else if (bootstrap == "by.root"){
    samp <- rep(1, n)
  } else {
    stop("Wrong bootstrap is given! ")
  }

  if (is.null(mtry)){
    data$id = id # this is a must, otherwise id cannot be passed to the next level
    mtry <- tune.ltrccf(formula = formula, data = data, id = id,
                        control = control, ntreeTry = ntree,
                        bootstrap = "by.user",
                        samptype = samptype,
                        sampfrac = sampfrac,
                        samp = samp,
                        na.action = na.action,
                        stepFactor = stepFactor,
                        applyfun = applyfun,
                        cores = cores,
                        trace = trace)
    print(sprintf("mtry is tuned to be %1.0f", mtry))
  }

  h2 <- function(y, x, start = NULL, weights, offset, estfun = TRUE, object = FALSE, ...) {
    if (all(is.na(weights)) == 1) weights <- rep(1, NROW(y))
    s <- .logrank_trafo2(y[weights > 0, , drop = FALSE])
    r <- rep(0, length(weights))
    r[weights > 0] <- s
    list(estfun = matrix(as.double(r), ncol = 1), converged = TRUE)
  }

  ret <- partykit::cforest(formula, data,
                           weights = samp,
                           perturb = perturb,
                           ytrafo = h2,
                           control = control,
                           na.action = na.action,
                           mtry = mtry,
                           ntree = ntree,
                           applyfun = applyfun,
                           cores = cores)
  ret$formulaLTRC <- formula
  ret$info$call <- Call
  ret$info$bootstrap <- bootstrap
  ret$info$samptype <- samptype
  ret$info$sampfrac <- sampfrac
  ret$data <- mf
  class(ret) <- c("ltrccf", "grow", class(ret))
  ret
}

