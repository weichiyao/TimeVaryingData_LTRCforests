#' Compute a Survival Curve from a LTRCCF model
#'
#' Constructs a monotone nonincreasing estimated survival curve from a LTRCCF model
#' for any given (left-truncated) right-censored survival data with time-varying covariates.
#' It can also compute survival function estimates for left-truncated right-censored data
#' with time-invariant covariates.
#'
#' @param object an object as returned by \code{\link{ltrccf}}.
#' @param newdata.id optional variable name of subject identifiers for \code{newdata}.
#' If this is present, it will be searched for in the \code{newdata} data frame.
#' Each group of rows in \code{newdata} with the same subject \code{id} represents
#' the covariate path through time of a single subject, and the result will
#' contain one curve per subject. If it is not specified, then an estimated survival
#' curve is returned for each row of \code{newdata}.
#' @param newdata an optional data frame containing the test data
#' (with the names of the variables the same as those in \code{data} from \code{object}).
#' @param OOB a logical specifying whether out-of-bag predictions are desired
#'
#'
#' (only if \code{newdata = NULL}).
#' @param time.eval a vector of time points, at which the estimated survival probabilities
#' will be computed.
#' @param time.tau an optional vector, with the \emph{i}-th entry giving the upper time limit for the
#' computed survival probabilities for the \emph{i}-th data of interest (i.e., only compute
#' survival probabilies at \code{time.eval[time.eval <= time.tau[i]]} for the \emph{i}-th
#' data of interest). If \code{OOB = TRUE}, the length of \code{time.tau} is equal to the length of
#' \code{data} used to train the \code{object};
#' If \code{OOB = FALSE}, the length of \code{time.tau} is equal to the length
#' of \code{newdata}, or equal to the length of \code{data} if \code{newdata} is not given.
#' The default \code{NULL} is simply to set all entries of \code{time.tau} equal to the maximum
#' value of \code{time.eval}, so that all estimated survival probabilities are computed at the
#' same \code{time.eval}.
#' @return A list containing:
#'    \item{survival.id}{subject identifiers.}
#'    \item{survival.obj}{an object of class \code{\link[survival]{Surv}}.}
#'    \item{survival.probs}{the estimated survival probabilities for each data of interest.
#'    It is a list if the length of the estimated values differs from one to another;
#'    otherwise, it is a matrix with the number of columns equal to the number of the data
#'    of interest, number of rows equal to the number of the time points at which the estimated
#'    survival probabilities are computed.}
#'    \item{survival.tau}{the input value \code{time.tau}.}
#'    \item{survival.times}{the input value \code{time.eval}. }
#' @import partykit
#' @import survival
#' @seealso \code{\link{sbrier_ltrc}} for evaluation of model fit
#' @examples
#' #### Example with data pbcsample
#' library(survival)
#' Formula <- Surv(Start, Stop, Event) ~ age + alk.phos + ast + chol + edema
#' ## Fit an LTRC conditional inference forest on time-varying data
#' LTRCCFobj <- ltrccf(formula = Formula, data = pbcsample, id = ID,
#'                     mtry = 3, ntree = 50L)
#'
#'
#' ## Construct an estimated survival estimate for the second subject
#' tpnt <- seq(0, max(pbcsample$Stop), length.out = 500)
#' newData <- pbcsample[pbcsample$ID == 2, ]
#' Pred <- predict(object = LTRCCFobj, newdata = newData, newdata.id = ID,
#'                 time.eval = tpnt)
#' ## Since time.tau = NULL, Pred$survival.probs is in the matrix format, with dimensions:
#' dim(Pred$survival.probs) # length(time.eval) x nrow(newdata)
#' ## Plot the estimated survival curve
#' plot(Pred$survival.times, Pred$survival.probs, type = "l", col = "red",
#'      xlab = "Time", ylab = "Survival probabilities")
#'
#'
#' ## When time.tau is specified and some entries are different from the others
#' Pred2 = predict(object = LTRCCFobj, newdata = pbcsample, newdata.id = ID,
#'                 time.eval = tpnt, time.tau = seq(100, 400,
#'                                                  length.out = length(unique(pbcsample$ID))))
#' ## Then Pred2$survival.probs is a list:
#' class(Pred2$survival.probs)
#' ## Plot the estimated survival curve for the subject with id = 18
#' plot(Pred2$survival.times[Pred2$survival.times <= Pred2$survival.tau[18]],
#'      Pred2$survival.probs[[18]], type = "l", col = "red",
#'      xlab = "Time", ylab = "Survival probabilities")
#'
#'
#' @export
predict.ltrccf <- function(object, newdata = NULL, newdata.id, OOB = FALSE,
                          time.eval, time.tau = NULL){

  pred <- partykit::predict.cforest(object = object, newdata = newdata, OOB = OOB, type = "prob",
                                    FUN = .pred_Surv_nolog)
  xvar.names <- attr(object$terms,"term.labels")
  yvar.names <- as.character(object$formulaLTRC[[2]])[2:4]
  Rname <- yvar.names[2]
  idname <- "id"

  # missing values can be present in the prediction
  if (is.null(newdata) || OOB){
    # first column: Surv(tleft,tright,event), second column: (id)
    newdata <- as.data.frame(as.matrix(object$data[, c(1, ncol(object$data)), drop = FALSE]))
    names(newdata) = c(yvar.names, idname)
  } else {
    if (missing(newdata.id)){
      newdata$id <- 1:nrow(newdata)
    } else {
      names(newdata)[names(newdata) == deparse(substitute(newdata.id))] <- idname
    }
    newdata <- as.data.frame(newdata[, c(yvar.names, idname)])
  }

  rm(object)
  N <- length(unique(newdata[, "id"])) # number of subjects
  if (is.null(time.tau)){
    time.tau <- rep(max(time.eval), N)
  } else {
    if (N != length(time.tau)) stop("time.tau should be a vector of length equaling to number of SUBJECT observation! \n
                                     In the time-varying case, check whether newdata.id has been correctly specified!")
  }

  Shat <- sapply(1:N, function(Ni) .shatfunc(Ni, data = newdata, pred = pred, tpnt = time.eval, tau = time.tau))
  obj <- Surv(newdata[, yvar.names[1]],
              newdata[, yvar.names[2]],
              newdata[, yvar.names[3]])
  RES <- list(survival.probs = Shat,
              survival.times = time.eval,
              survival.tau = time.tau,
              survival.obj = obj,
              survival.id = newdata$id)
  rm(newdata)
  rm(Shat)
  rm(time.eval)
  rm(time.tau)
  return(RES)
}

.shatfunc <- function(Ni, data, pred, tpnt, tau){
  ## This function is to compute the estimated survival probability of the Ni-th subject
  id.seu <- data[, "id"] # id
  id.sub <- unique(id.seu)

  ## the i-th data
  TestData <- data[id.seu == id.sub[Ni], ]

  TestT <- c(TestData[1, 1], TestData[, 2])
  TestTIntN <- nrow(TestData)

  tpnt <- tpnt[tpnt <= tau[Ni]]
  tlen <- length(tpnt)

  ## Compute the estimated survival probability of the Ni-th subject
  Shat_temp <- matrix(0, nrow = 1, ncol = tlen)

  r.ID <- findInterval(tpnt, TestT)
  # r.ID[tpnt >= TestData[, 2][TestTIntN]] <- TestTIntN
  r.ID[r.ID >= TestTIntN] <- TestTIntN

  jall <- unique(r.ID[r.ID > 0])
  nj <- length(jall)

  ## Deal with left-truncation
  Shat_temp[1, r.ID == 0] = 1
  if(nj == 1){
    ## Get the index of the Pred to compute Shat
    II <- which(id.seu == id.sub[Ni])[jall[nj]]
    Shat_i = ipred::getsurv(pred[[II]], tpnt[r.ID == jall[nj]])
    Shat_temp[1, r.ID == jall[nj]] <- Shat_i / Shat_i[1]
  } else if (nj > 1) {
    # c(1, S_{1}(R_{1}), ..., S_{nj}(R_{nj}))
    ShatR_temp <- matrix(0, nrow = 1, ncol = nj + 1)
    ShatR_temp[1, 1] <- 1
    ## Get the position of L_1, ..., L_n.
    ## Note that, L_2=R_1, ..., L_{j+1} = R_{j}, ..., L_n = R_{n-1}
    r.IDmax <- c(min(which(r.ID == jall[1])), sapply(jall[-nj], function(j){
      max(which(r.ID == j)) + 1
    }))

    # S_1(L_1), S_2(L_2), S_3(L_3), ..., S_{nj}(L_{nj})
    qL = rep(0, nj)
    for (j in 1:nj){
      ## Get the index of the Pred to compute Shat
      II <- which(id.seu == id.sub[Ni])[1] + jall[j] - 1
      Shat_j = ipred::getsurv(pred[[II]], tpnt[r.ID == jall[j]])

      qL[j] <- Shat_j[1]
      # S_{j}(R_{j}), j=1,...nj-1
      jR = ipred::getsurv(pred[[II]], TestT[j + 1])
      ShatR_temp[1, j + 1] = jR / qL[j]
      Shat_temp[1, r.ID == jall[j]] <- Shat_j / qL[j]
    }

    ql0 <- which(qL == 0)
    if(length(ql0) > 0){
      maxqlnot0 <- max(which(qL > 0))
      for(j in ql0){
        if (j < maxqlnot0) {
          ShatR_temp[1, j + 1] <- 1
          Shat_temp[1, r.ID == jall[j]] <- 1
        } else{
          ShatR_temp[1, j + 1] <- 0
          Shat_temp[1, r.ID == jall[j]] <- 0
        }
      }
    }
    m <- cumprod(ShatR_temp[1, 1:nj])
    for (j in 1:nj){
      Shat_temp[1, r.ID == jall[j]] <- Shat_temp[1, r.ID == jall[j]] * m[j]
    }
  }

  return(Shat_temp[1, ])
  rm(Shat_temp)
  rm(ShatR_temp)
  rm(id.seu)
  rm(id.sub)
}

.pred_Surv_nolog <- function(y, w) {
  if (length(y) == 0) return(NA)
  survfit(y ~ 1, weights = w, subset = w > 0, conf.type = "none", se.fit = FALSE)
}

