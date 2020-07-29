#' Fit a LTRC random survival forest with Poisson splitting rule
#'
#' An implementation of the random forest algorithms utilizing LTRC \code{rpart}
#' trees \code{\link[LTRCtrees]{LTRCART}} as base learners for left-truncated right-censored
#' survival data with time-invariant covariates. It also allows for (left-truncated)
#' right-censored survival data with time-varying covariates.
#'
#' This function extends the random survival forest algorithm in
#' \code{\link[randomForestSRC]{rfsrc}} to fit left-truncated and right-censored data,
#' which allows for time-varying covariates.
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
#' @param na.action a function which indicates what should happen when the data contain
#' missing values.
#' @param mtry number of input variables randomly sampled as candidates at each node for
#' random forest like algorithms. The default \code{mtry} is tuned by
#' \code{\link{tune.ltrcrsf}}.
#' @param ntree an integer, the number of the trees to grow for the forest.
#' \code{ntree = 100L} is set by default.
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
#' Regardless of the presence of \code{id}, if \code{"none"} is chosen,
#' \code{data} is not bootstrapped at all, and is used in
#' every individual tree. If \code{"by.user"} is choosen,
#' the bootstrap specified by \code{samp} is used.
#' @param samptype choices are \code{swor} (sampling without replacement) and
#' \code{swr} (sampling with replacement). The default action here is sampling
#' without replacement.
#' @param samp Bootstrap specification when \code{bootstype = "by.user"}.
#' Array of dim \code{n x ntree} specifying how many times each record appears
#' in each bootstrap sample.
#' @param trace whether to print the progress of the search of the optimal value
#' of \code{mtry} if \code{mtry} is not specified (see \code{\link{tune.ltrcrsf}}).
#' \code{trace = TRUE} is set by default.
#' @param stepFactor at each iteration, \code{mtry} is inflated (or deflated)
#' by this value, used when \code{mtry} is not specified (see \code{\link{tune.ltrcrsf}}).
#' The default value is \code{2}.
#' @param nodesize an integer, forest average terminal node size. It has been
#' adjusted from the \code{\link[randomForestSRC]{rfsrc}} default \code{15} for survival families.
#' @param nodedepth maximum depth to which a tree should be grown. The default behaviour
#' is that this parameter is ignored.
#' @param nsplit an non-negative integer value for number of random splits to consider
#' for each candidate splitting variable. This significantly increases speed.
#' When zero or \code{NULL}, the algorithm uses much slower deterministic splitting where
#' all possible splits are considered. \code{nsplit = 10L} by default.
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
#' Further demails can be found in \code{\link[randomForestSRC]{rfsrc}}.
#' @return An object belongs to the class \code{ltrcrsf},
#' as a subclass of \code{\link[randomForestSRC]{rfsrc}}.
#' @import survival
#' @import stats
#' @import utils
#' @import prodlim
#' @importFrom survival Surv
#' @seealso \code{\link{predictProb}} for prediction and \code{\link{tune.ltrcrsf}}
#' for \code{mtry} tuning.
#' @references Andersen, P. and Gill, R. (1982). Cox’s regression model for counting
#' processes, a large sample study. \emph{Annals of Statistics}, \strong{10}, 1100-1120.
#' @examples
#' #### Example with time-varying data pbcsample
#' library(survival)
#' Formula = Surv(Start, Stop, Event) ~ age + alk.phos + ast + chol + edema
#' # Built a LTRCRSF forest (based on bootstrapping subjects without replacement)
#' # on the time-varying data by specifying id:
#' LTRCRSFobj = ltrcrsf(formula = Formula, data = pbcsample, id = ID, stepFactor = 3,
#'                      ntree = 20L)
#' LTRCRSFobj
#'
#'
#' @export
ltrcrsf <- function(formula, data, id, ntree = 100L, mtry = NULL,
                    nodesize = max(ceiling(sqrt(nrow(data))),15),
                    bootstrap = c("by.sub","by.root","by.node","by.user","none"),
                    samptype = c("swor", "swr"),
                    sampfrac = 0.632,
                    samp = NULL,
                    na.action = "na.omit",
                    stepFactor = 2,
                    trace = TRUE,
                    nodedepth = NULL,
                    nsplit = 10L,
                    ntime){
  Call <- match.call()
  Call[[1]] <- as.name('ltrcrsf')  #make nicer printout for the user
  # create a copy of the call that has only the arguments we want,
  #  and use it to call model.frame()
  indx <- match(c('formula', 'data', 'id'),
                names(Call), nomatch = 0)
  if (indx[1]==0) stop("a formula argument is required")
  Call$formula <- eval(formula)

  temp <- Call[c(1, indx)]
  temp[[1L]] <- quote(stats::model.frame)
  mf <- eval.parent(temp)
  y <- model.extract(mf, 'response')

  if (!is.Surv(y)) stop("Response must be a survival object")
  if (!attr(y, "type") == "counting") stop("The Surv object must be of type 'counting'.")
  rm(y)
  rm(mf)

  ## if not specified, the first one will be used as default
  bootstrap <- match.arg(bootstrap)
  samptype <- match.arg(samptype)
  # na.action <- match.arg(na.action)

  n <- nrow(data)
  # pull y-variable names
  yvar.names <- all.vars(formula(paste(as.character(formula)[2], "~ .")), max.names = 1e7)
  yvar.names <- yvar.names[-length(yvar.names)]
  if (length(yvar.names) == 4){
    yvar.names = yvar.names[2:4]
  }
  # try an example:
  # Formula = Surv(pbcsample$Start, pbcsample$Stop, pbcsample$Event, type="counting") ~ pbcsample$age + pbcsample$alk.phos + pbcsample$ast + pbcsample$chol + pbcsample$edema
  # and see the difference from:
  # yvar.names <- as.character(Formula[[2]])[2:4]

  # extract the x-variable names
  xvar.names <- attr(terms(formula), 'term.labels')
  rm(temp)

  # this is a must, otherwise id cannot be passed to the next level in tune.ltrccf
  if (indx[3] == 0){
    ## If id is not present, then we add one more variable
    data$id <- 1:n
  } else {
    ## If id is present, then we rename the column to be id
    names(data)[names(data) == deparse(substitute(id))] <- "id"
  }
  data <- data[, c("id", yvar.names, xvar.names), drop = FALSE]

  ## Transformation for LTRC data
  Status <- data[, yvar.names[3]]
  Times <- data[, yvar.names[2]]

  if (sum(Status) == 0) stop("All observations are right-censored with event = 0!")
  ##unique death times
  unique.times <- sort(unique(Times[Status == 1]))

  if (missing(ntime)){
    ntime = unique.times
  }

  y <- survival::Surv(data[, yvar.names[1]], data[, yvar.names[2]], data[, yvar.names[3]])
  temp <- survival::coxph(y ~ 1)
  cumhaz.table <- survival::basehaz(temp)

  # Check if Inf hazard exists
  if(sum(is.infinite(cumhaz.table$hazard))!=0){
    cumhaz.table <- cumhaz.table[cumhaz.table$hazard != Inf,] ## subset(cumhaz.table, hazard != Inf)
  }

  cumhaz.table2 <- cumhaz.table[cumhaz.table$time %in% unique.times,]

  cumhaz.times <- c(0, cumhaz.table2$time[-length(cumhaz.table2$time)], max(Times))
  cumhaz <- c(0, cumhaz.table2$hazard)

  Start.cumhaz <- stats::approx(cumhaz.times, cumhaz, y[, 1L])$y
  End.cumhaz <- stats::approx(cumhaz.times, cumhaz, y[, 2L])$y

  data$Newtime <- End.cumhaz - Start.cumhaz
  Formula = formula(paste(c(paste("Surv(Newtime,", yvar.names[3],")",sep = ""), formula[[3]]), collapse = "~"))
  rm(y)

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
  ## bootstrap case
  if (n.seu == n.sub){ # time-invariant LTRC data
    # it includes the case 1) when id = NULL, which is that id is not specified
    #                      2) when id is specified, but indeed LTRC time-invariant
    if (bootstrap == "by.sub") bootstrap = "by.root"
  }

  sampsize <- if (samptype == "swor") function(x){x * sampfrac} else function(x){x}

  bootstrap.org <- bootstrap
  if (bootstrap == "by.sub"){
    bootstrap <- "by.user"
    # dim n x ntree
    samp <- matrix(0, nrow = n.seu, ncol = ntree)
    if (samptype == "swr"){
      for (b in 1:ntree){
        while (sum(samp[, b]) < n.seu) {
          idx <- sample(id.sub, size = 1)
          inidx <- which(data$id[takeid] == idx)
          samp[inidx, b] = samp[inidx, b] + 1
        }
        if (sum(samp[, b]) > n.seu){
          seqn <- which(samp[, b] != 0)
          add <- length(seqn) - (sum(samp[, b]) - n.seu) + 1
          idx <- sample(add, size = 1)
          samp[seqn[idx:(idx + sum(samp[, b]) - n.seu - 1)],b] <- samp[seqn[idx:(idx + sum(samp[, b]) - n.seu - 1)], b] - 1
        }
      }
    } else if (samptype == "swor"){
      for (b in 1:ntree){
        nP <- floor(n.seu * sampfrac)
        idS <- sample(id.sub, size = n.sub, replace = FALSE)
        k <- 0
        while (sum(samp[, b]) < nP) {
          k <- k + 1
          inidx <- which(data$id[takeid] == idS[k])
          samp[inidx, b] <- samp[inidx, b] + 1
        }
        if (sum(samp[, b]) > nP){
          seqn <- which(samp[, b] != 0)
          add <- length(seqn) - (sum(samp[, b]) - nP) + 1
          idx <- sample(add, size = 1)
          samp[seqn[idx:(idx + sum(samp[, b]) - nP - 1)], b] <- 0
        }
      }
    } else {
      stop("Wrong samptype is given!")
    }
  }

  if (is.null(mtry)){
    # data$id = id # this is a must, otherwise id cannot be passed to the next level
    mtry <- tune.ltrcrsf(formula = formula,
                         data = data,
                         id = id,
                         bootstrap = bootstrap,
                         samptype = samptype,
                         sampfrac = sampfrac,
                         samp = samp,
                         ntreeTry = ntree,
                         nodesizeTry = nodesize,
                         nodedepth = nodedepth,
                         nsplit = nsplit,
                         na.action = na.action,
                         ntime = ntime,
                         stepFactor = stepFactor,
                         trace = trace,
                         doBest = FALSE)
    print(sprintf("mtry is tuned to be %1.0f", mtry))
  }
  ## Use randomSurvivalForest package
  forest.fit <- ltrcrfsrc(formula = Formula,
                          data = data,
                          ntree = ntree,
                          mtry = mtry,
                          nodesize = nodesize,
                          nodedepth = nodedepth,
                          splitrule = "custom1",
                          nsplit = nsplit,
                          bootstrap = bootstrap,
                          samptype = samptype,
                          sampsize = sampsize,
                          samp = samp,
                          membership = TRUE,
                          na.action = na.action,
                          ntime = ntime)
  forest.fit$yvarLTRC.names <- c("id", yvar.names)
  forest.fit$yvarLTRC = data[takeid, c("id", yvar.names), drop = FALSE]
  forest.fit$err.rate = NULL
  forest.fit$survival = NULL
  forest.fit$survival.oob = NULL
  forest.fit$chf = NULL
  forest.fit$chf.oob = NULL
  forest.fit$forest$bootstrapLTRC <- bootstrap.org
  forest.fit$forest$samptypeLTRC <- samptype
  forest.fit$forest$sampfracLTRC <- sampfrac
  forest.fit$callLTRC = Call
  forest.fit$splitruleLTRC = "Poisson"
  forest.fit$formulaLTRC <- formula
  class(forest.fit) <- "ltrcrsf"
  return(forest.fit)
}


