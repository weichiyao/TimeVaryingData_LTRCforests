#' Model fit evaluation for LTRC forests.
#'
#' Compute the (integrated) Brier score to evaluate the model fit for
#' (left-truncated) right-censored survival data with time-varying covariates,
#' as well as left-truncated right-censored data with time-invariant covariates.
#'
#' @param obj an object of class \code{\link[survival]{Surv}}, formed on
#' left-truncated right-censored observations (which are pseudo-subject
#' observations in the time-varying case).
#' @param id an optional vector as subject identifiers for \code{obj}.
#' @param pred a list. This should contain 1) either a matrix
#' or a list of survival probabilies named \code{survival.probs}; 2) a sequence
#' of time points \code{survival.times}; 3) a vector of upper time limits
#' \code{survival.tau}.
#' Please see the values returned by \code{\link{pred.ltrccf}} and
#' \code{\link{pred.ltrcrsf}}.
#' @param type a character string denoting the type of scores returned. If \code{type = "IBS"},
#' the integrated Brier score up to the last time point in \code{pred$surv.times} that is
#' not larger than the minimmum value of \code{pred$surv.tau} is returned.
#' If \code{type = "BS"}, the Brier score at every time point in \code{pred$surv.times} up to
#' the minimum value of \code{pred$surv.tau} is returned.
#' @keywords Brier score, integrated Brier score
#' @return
#' If \code{type = "IBS"}, this returns the integrated Brier score.
#' @return
#' If \code{type = "BS"}, this returns \code{BScore}, the Brier scores
#' and \code{Time}, the time points at which the scores are computed.
#' @import partykit
#' @import ipred
#' @import prodlim
#' @importFrom survival Surv
#' @import stats
#' @import utils
#' @examples
#' ### Example with dataset pbcsample
#' Formula = Surv(Start, Stop, Event) ~ age + alk.phos + ast + chol + edema
#' ## Fit an LTRC conditional inference forest on time-varying data
#' LTRCCFobj = ltrccf(formula = Formula, data = pbcsample, id = ID, mtry = 3, ntree = 50L)
#'
#' # Time points
#' tpnt = seq(0, 6000, by = 200)
#' # Set different upper time limits for each of the 20 subjects
#' tau = seq(4001, 6200, length.out = 20)
#' ## Obstain estimation at time points tpnt
#' Predobj = predict(object = LTRCCFobj, time.eval = tpnt, time.tau = tau)
#'
#' ## Compute the integrated Brier score:
#' pbcobj = Surv(pbcsample$Start, pbcsample$Stop, pbcsample$Event)
#' IBS = sbrier_ltrc(obj = pbcobj, id = pbcsample$ID, pred = Predobj, type = "IBS")
#'
#' ## Compute the Brier score at each value of tpnt
#' BS = sbrier_ltrc(obj = pbcobj, id = pbcsample$ID, pred = Predobj, type = "BS")
#' ## Plot the Brier scores
#' plot(BS$Time, BS$BScore, pch = 20, xlab = "Time", ylab = "Brier score", col = 2)
#' ## As one can see, the Brier scores are returned at all tpnt up to 4000,
#' ## this is because the algorithm set the last evaluation time point
#' ## to be 4000 based on the value of time.eval and time.tau
#' ## (max(tpnt[tpnt <= min(tau)]) == 4000).
#' @export

sbrier_ltrc <- function(obj, id = NULL, pred, type = c("IBS","BS")){
  if(!inherits(obj, "Surv"))
    stop("obj is not of class Surv")

  # check for left-trunctation and right-censoring
  if (attr(obj, "type") != "counting")
    stop("only dataset with left-truncated and right-censored (pseudo-subject) observations allowed")

  n <- nrow(obj)

  if (!is.null(id)){
    if (n != length(id)) stop("The length of id is different from the Surv object!")
  } else { # id is NULL => LTRC observations with time-invariant covariates
    id <- 1:n # label
  }

  id.sub = unique(id)
  n.sub = length(id.sub)

  obj <- as.data.frame(as.matrix(obj))

  if (n == n.sub){# ltrc data with time-invariant covariates
    data_sbrier = obj
    data_sbrier$id = 1:n
  } else {# ltrc data with time-varying covariates
    data_sbrier <- data.frame(matrix(0, nrow = n.sub, ncol = 3))
    names(data_sbrier) <- c("start", "stop", "status")
    data_sbrier$id <- id.sub
    for (ii in 1:n.sub){
      data_sbrier[ii, ]$start = min(obj[id == id.sub[ii], ]$start)
      data_sbrier[ii, ]$stop = max(obj[id == id.sub[ii], ]$stop)
      data_sbrier[ii, ]$status = sum(obj[id == id.sub[ii], ]$status)
    }
  }

  if (type[1] == "IBS"){
    ret <- sapply(1:n.sub, function(Ni) ibsfunc(Ni = Ni, data_sbrier = data_sbrier, pred = pred))
    ret <- mean(ret)
    names(ret) = "Integrated Brier score"
  } else if (type[1] == "BS"){
    # Brier score will be evaluated up to the last time point where survival probabilities of all data are computed.
    tpnt <- pred$survival.times[pred$survival.times <= min(pred$survival.tau)]
    bsres <- sapply(1:n.sub, function(Ni) bsfunc(Ni = Ni, data_sbrier = data_sbrier, pred = pred, tpnt = tpnt))
    bsres <- rowMeans(bsres)
    ret <- data.frame(matrix(0, ncol = 2, nrow = length(tpnt)))
    colnames(ret) <- c("Time", "BScore")
    ret$Time <- tpnt
    ret$BScore <- bsres
  } else {
    stop("type can only be 'IBS' or 'BS'")
  }
  return(ret)
}

ibsfunc <- function(Ni, data_sbrier, pred){
  id_uniq <- unique(data_sbrier$id)
  tpnt = pred$survival.times[pred$survival.times <= pred$survival.tau[Ni]]
  tlen = length(tpnt)
  ## Get the estimated survival probabilities
  if (class(pred$survival.probs)[1] == "matrix"){
    Shat = pred$survival.probs[1:tlen, Ni]
  } else if(class(pred$survival.probs)[1] == "list"){
    Shat = pred$survival.probs[[Ni]][1:tlen]
  }

  ######================ reverse Kaplan-Meier: estimate censoring distribution ====== ########
  # deal with ties
  hatcdist <- prodlim(Surv(start, stop, status) ~ 1, data = data_sbrier, reverse = TRUE)

  Ttildei <- data_sbrier[data_sbrier$id == id_uniq[Ni], ]$stop
  ### conditional survival for Observed value < t, G(Obs)
  csurv_obs <- predict(hatcdist, times = Ttildei, type = "surv")
  csurv_obs[csurv_obs == 0] <- Inf

  # conditional survival for Observed value > t, G(t)
  csurv_t <- predict(hatcdist, times = tpnt[tpnt < Ttildei], type = "surv")
  csurv_t[is.na(csurv_t)] <- min(csurv_t, na.rm = TRUE)
  csurv_t[csurv_t == 0] <- Inf

  ## c(G^{-1}(t), G^{-1}(Obs))
  csurv <- c(1/csurv_t, rep(1 / csurv_obs, sum(tpnt >= Ttildei)))

  ######================ indicator ================#################
  Indicator_t <- as.integer(tpnt < Ttildei)
  Indicator_t[Indicator_t == 0] = as.integer(data_sbrier[data_sbrier$id == id_uniq[Ni],]$status == 1)

  ######================ Brier score =================#################
  fibs_itg = (as.integer(tpnt < Ttildei) - Shat) ^ 2 * csurv * Indicator_t
  ibs = diff(tpnt) %*% (fibs_itg[-length(fibs_itg)] + fibs_itg[-1]) / 2
  ibs = ibs / diff(range(tpnt))
  ibs
}

bsfunc <- function(Ni, data_sbrier, pred, tpnt){
  id_uniq <- unique(data_sbrier$id)
  tlen = length(tpnt)
  ## Get the estimated survival probabilities
  if (class(pred$survival.probs)[1] == "matrix"){
    Shat = pred$survival.probs[1:tlen, Ni]
  } else if(class(pred$survival.probs)[1] == "list"){
    Shat = pred$survival.probs[[Ni]][1:tlen]
  }

  ######================ reverse Kaplan-Meier: estimate censoring distribution ================###########
  # deal with ties
  hatcdist <- prodlim(Surv(start, stop, status) ~ 1, data = data_sbrier, reverse = TRUE)

  Ttildei <- data_sbrier[data_sbrier$id == id_uniq[Ni], ]$stop
  ### conditional survival for Observed value < t, G(Obs)
  csurv_obs <- predict(hatcdist, times = Ttildei, type = "surv")
  csurv_obs[csurv_obs == 0] <- Inf

  # conditional survival for Observed value > t, G(t)
  csurv_t <- predict(hatcdist, times = tpnt[tpnt < Ttildei], type = "surv")
  csurv_t[is.na(csurv_t)] <- min(csurv_t, na.rm = TRUE)
  csurv_t[csurv_t == 0] <- Inf

  ## c(G^{-1}(t), G^{-1}(Obs))
  csurv <- c(1 / csurv_t, rep(1 / csurv_obs, sum(tpnt >= Ttildei)))

  ######================ indicator ================#################
  Indicator_t <- as.integer(tpnt < Ttildei)
  Indicator_t[Indicator_t == 0] = as.integer(data_sbrier[data_sbrier$id == id_uniq[Ni], ]$status == 1)

  ######================ Brier score =================#################
  fibs_itg = (as.integer(tpnt < Ttildei) - Shat) ^ 2 * csurv * Indicator_t
  fibs_itg
}


