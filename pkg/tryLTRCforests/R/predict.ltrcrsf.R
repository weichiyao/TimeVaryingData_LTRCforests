#' Compute a Survival Curve from a LTRCRSF model
#'
#' Constructs a monotone nonincreasing estimated survival curve from a LTRCRSF model
#' for any given (left-truncated) right-censored survival data with time-varying covariates.
#' It can also conputes survival function estimates for left-truncated right-censored data
#' with time-invariant covariate.
#'
#' @param object an object as returned by \code{\link{ltrcrsf}}.
#' @param newdata.id optional variable name of subject identifiers for \code{newdata}.
#' If this is present, it will be search for in the \code{newdata} data frame.
#' Each group of rows in \code{newdata} with the same subject \code{id} represents
#' the covariate path through time of a single subject, and the result will
#' contain one curve per subject. If it is not specified, then predictions are
#' returned for each row of \code{newdata}.
#' @param newdata an optional data frame containing the test data
#' (with the names of the variables the same as those in \code{data} from \code{object}).
#' @param OOB a logical specifying whether out-of-bag predictions are desired
#'
#'
#' (only if \code{newdata = NULL}).
#' @param time.eval a vector of time points, at which the estimated survival probabilities
#' are computed.
#' @param time.tau an optional vector, with the \emph{i}-th entry giving the upper time limit for the
#' computed survival probabilities for the \emph{i}-th data of interest (i.e., only compute
#' survival probabilies at \code{time.eval[time.eval <= time.tau[i]]} for the \emph{i}-th
#' data of interest). If \code{OOB = TRUE}, the length of \code{time.tau} is equal to the length of
#' \code{data} used to train the \code{object};
#' If \code{OOB = FALSE}, the length of \code{time.tau} is equal to the length
#' of \code{newdata}, or equal to the length of \code{data} if \code{newdata} is not given.
#' The default \code{NULL} is simply to set all entries of \code{time.tau} equal to the maximum
#' value of \code{time.eval}, therefore all estimated survival are computed at the
#' same \code{time.eval}.
#' @return A list containing:
#'    \item{survival.obj}{an object of class \code{\link[survival]{Surv}}.}
#'    \item{survival.probs}{the estimated survival probabilities for each data of interest.
#'    It is a list if the length of the estimated values differs from one to another;
#'    otherwise, it is a matrix with the number of columns equal to the number of the data
#'    of interest, number of rows equal to the number of the time points at which the estimated
#'    survival probabilities are computed.}
#'    \item{survival.tau}{the input value \code{time.tau}.}
#'    \item{survival.times}{the input value \code{time.eval}.}
#' @import ipred
#' @import survival
#' @seealso \code{\link{sbrier_ltrc}} for evaluation of model fit
#' @examples
#' #### Example with time-varying data pbcsample
#' ## View the prebuilt object
#' LTRCRSFobj
#'
#' ## Construct an estimated survival estimate for the second subject
#' tpnt <- seq(0, max(pbcsample$Stop), length.out = 500)
#' newData <- pbcsample[pbcsample$ID == 2,]
#' Pred <- predict(object = LTRCRSFobj, newdata = newData, newdata.id = ID, time.eval = tpnt)
#'
#' ## Since time.tau = NULL, Pred$survival.probs is in the matrix format, with dimensions:
#' dim(Pred$survival.probs) # length(time.eval) x nrow(newdata)
#' ## Plot the estimated survival curve
#' plot(Pred$survival.times, Pred$survival.probs, type = "l", col = "red",
#'      xlab = "Time", ylab = "Survival probabilities")
#'
#' ## When time.tau is specified and some entries are different from the others
#' Pred2 = predict.ltrcrsf(object = LTRCRSFobj, newdata = pbcsample, newdata.id = ID,
#'                         time.eval = tpnt, time.tau = seq(100, 400, length.out =
#'                                                            length(unique(pbcsample$ID))))
#'
#' ## Then Pred2$survival.probs is a list:
#' class(Pred2$survival.probs)
#' ## Plot the estimated survival curve for the subject with id = 20
#' plot(Pred2$survival.times[Pred2$survival.times <= Pred2$survival.tau[20]],
#'      Pred2$survival.probs[[20]], type = "l", col = "red",
#'      xlab = "Time", ylab = "Survival probabilities")
#' @export

predict.ltrcrsf <- function(object, newdata = NULL, newdata.id, OOB = FALSE,
                           time.eval, time.tau = NULL){
  ntree <- object$ntree
  formula <- object$formulaLTRC
  formula[[3]] <- 1

  wt <- object$inbag # of size Ndata x ntree

  yvar.names <- object$yvarLTRC.names
  Rname <- yvar.names[2]
  traindata <- object$yvarLTRC
  traindata$id <- object$id
  if (OOB){
    # relabel the traindata
    traindata$I <- 1:nrow(traindata) # make sure partial/baseline has done this

    id_uniq <- unique(traindata$id) # id of subjects
    n_uniq <- length(id_uniq) # number of subjects
    node_all <- object$membership # of size Ndata*ntree
    rm(object)
    if (is.null(time.tau)){
      time.tau <- rep(max(time.eval), n_uniq)
    }
    pred <- sapply(1:n_uniq, function(wi){
      newi <- traindata[traindata$id == id_uniq[wi], , drop = FALSE]
      n_newi <- nrow(newi)

      ## up to tau_i
      tpnt <- time.eval[time.eval <= time.tau[wi]]
      tlen <- length(tpnt)

      r.ID <- findInterval(tpnt, c(0, newi[, Rname]))
      r.ID[tpnt >= newi[, Rname][n_newi]] = n_newi

      jall <- unique(r.ID[r.ID > 0])
      nj <- length(jall)
      if (nj == 1){
        survival <- matrix(0, nrow = 1, ncol = tlen)
        ## deal with left-truncation
        survival[1, r.ID == 0] <- 1

        ## find out which trees does not contain the I[jall[j]]-th wi-th data
        id_tree_wi_j <- which(wt[newi$I[jall[nj]], ] == 0)

        for (ti in 1:length(id_tree_wi_j)){
          ## In each tree of id in idTree_wi, it falls into terminal id_node_witi_j
          id_node_witi_j <- node_all[newi$I[jall[nj]], id_tree_wi_j[ti]]
          ## id of samples that fall into the same node
          id_samenode_witi_j <- which(node_all[,id_tree_wi_j[ti]] == id_node_witi_j)
          ## Pick out those appearing in the bootstrapped samples
          id_inbag_j <- which(wt[, id_tree_wi_j[ti]] == 1)
          id_buildtree_j <- id_samenode_witi_j[id_samenode_witi_j %in% id_inbag_j]
          ## Build the survival tree
          KM <- survival::survfit(formula = formula, data = traindata[id_buildtree_j, ])
          ## Get survival probabilities
          Shat_ti <- ipred::getsurv(KM, tpnt[r.ID == jall[nj]])

          survival[1, r.ID == jall[nj]] <- survival[1, r.ID == jall[nj]] + Shat_ti / Shat_ti[1]
        }

        survival[1, r.ID == jall[nj]] <- survival[1, r.ID == jall[nj]] / length(id_tree_wi_j)
        RES <- survival[1, ]
      } else if (nj > 1) {
        # on [0, L_1), [L_1,R_1), [L_2,R_2), ..., [L_n,R_n]
        survival <- matrix(0, nrow = 1, ncol = tlen)
        # deal with left-truncation
        survival[1, r.ID == 0] <- 1
        # c(1, S_{1}(R_{1})/S_{1}(L_{1}),...,S_{n-1}(R_{n-1})/S_{n-1}(L_{n-1})),S_n(R_n)/S_n(L_n))
        survivalR <- matrix(0, nrow = 1, ncol = nj + 1)
        for (j in 1:nj){
          ## find out which trees does not contain the I[jall[j]]-th wi-th data
          id_tree_wi_j <- which(wt[newi$I[jall[j]], ] == 0)

          for (ti in 1:length(id_tree_wi_j)){
            ## In each tree of id in idTree_wi, it falls into terminal id_node_witi_j
            id_node_witi_j <- node_all[newi$I[jall[j]], id_tree_wi_j[ti]]
            ## id of samples that fall into the same node
            id_samenode_witi_j <- which(node_all[, id_tree_wi_j[ti]] == id_node_witi_j)
            ## Pick out those appearing in the bootstrapped samples
            id_inbag_j <- which(wt[,id_tree_wi_j[ti]] == 1)
            id_buildtree_j <- id_samenode_witi_j[id_samenode_witi_j %in% id_inbag_j]
            ## Build the survival tree
            KM <- survival::survfit(formula = formula, data = traindata[id_buildtree_j, ])
            Shat_ti <- ipred::getsurv(KM, tpnt[r.ID == jall[j]])

            if (Shat_ti[1] == 0){# jL = Shat_ti[1]
              survival[1, r.ID == jall[j]] <- survival[1, r.ID == jall[j]] + 1
              survivalR[1, j + 1] <- survivalR[1, j + 1] + 1
            } else {
              survival[1, r.ID == jall[j]] <- survival[1, r.ID == jall[j]] + Shat_ti / Shat_ti[1]
              survivalR[1, j + 1] <- survivalR[1, j + 1] + ipred::getsurv(KM, newi[, Rname][j]) / Shat_ti[1]
            }

          }
          survival[1, r.ID == jall[j]] <- survival[1, r.ID == jall[j]] / length(id_tree_wi_j)
          survivalR[1, j + 1] <- survivalR[1, j + 1] / length(id_tree_wi_j)
        }

        survivalR[1, 1] <- 1

        m <- cumprod(survivalR[1, 1:nj])
        for (j in 2:nj){
          survival[1, r.ID == jall[j]] <- m[j] * survival[1, r.ID == jall[j]]
        }
        RES <- survival[1, ]

      }
      return(RES)
    })
    obj <- Surv(traindata[, yvar.names[1]],
                traindata[, yvar.names[2]],
                traindata[, yvar.names[3]])
    rm(traindata)
  } else {
    nIDxdata <- object$membership # of size Ndata*ntree

    if (is.null(newdata)){
      newdata <- traindata
      nIDxnewdata <- nIDxdata
    } else {
      if (!is.data.frame(newdata)) stop("newdata must be a dataframe")
      x.IDs <- match(object$xvar.names, names(newdata))
      class(object) = class(object)[2:4] # for the prediction function in rfsrc to work
      nIDxnewdata <- predict.ltrcrfsrc(object, newdata = newdata[, x.IDs], membership = TRUE)$membership # of size Newdata*ntree

      if (missing(newdata.id)){
        newdata$id <- 1:nrow(newdata)
      } else {
        names(newdata)[names(newdata) == deparse(substitute(newdata.id))] <- "id"
      }

    }
    rm(object)
    obj.IDs <- match(c("id", yvar.names), names(newdata))
    newdata = newdata[, obj.IDs]

    # label the newdata
    newdata$I <- 1:nrow(newdata)

    id_uniq <- unique(newdata$id) # id of subjects
    n_uniq <- length(id_uniq) # number of subjects

    if (is.null(time.tau)){
      time.tau <- rep(max(time.eval), n_uniq)
    } else {
      if (n_uniq != length(time.tau)) stop("time.tau should be a vector of length equaling to number of SUBJECT observation! In the time-varying case, check whether newdata.id has been correctly specified!")
    }
    pred <- sapply(1:n_uniq, function(i){
      newi <- newdata[newdata$id == id_uniq[i], , drop = FALSE]
      n_newi <- nrow(newi)

      ## up to tau_i
      tpnt = time.eval[time.eval <= time.tau[i]]
      tlen <- length(tpnt)

      r.ID <- findInterval(tpnt, c(0, newi[, Rname]))
      r.ID[tpnt >= newi[, Rname][n_newi]] <- n_newi
      jall <- unique(r.ID[r.ID > 0])
      nj <- length(jall)

      if (nj == 1){
        # on [L_1,R_1), [L_2,R_2), ..., [L_n,R_n]
        survival <- matrix(0, nrow = 1, ncol = tlen)
        # deal with left truncation
        survival[1, r.ID == 0] <- 1
        for (b in 1:ntree){
          # observations that fall in the same terminal nodes as the new observation in b-th bootstrapped samples
          IDnew <- which(nIDxdata[,b] == nIDxnewdata[newi$I[jall[nj]], b])
          # ID of observations in the b-th bootstrap samples
          rw <- which(wt[, b] == 1)
          IDnew <- IDnew[IDnew %in% rw]
          KM <- survival::survfit(formula = formula, data = traindata[IDnew, ])
          Shat_b <- ipred::getsurv(KM, tpnt[r.ID == jall[nj]])
          survival[1, r.ID == jall[nj]] <- survival[1, r.ID == jall[nj]] + Shat_b / Shat_b[1]
        }

        survival[1, r.ID == jall[nj]] = survival[1, r.ID == jall[nj]] / ntree
        RES = survival[1, ]
      } else if (nj > 1) {
        # on [0, L_1), [L_1,R_1), [L_2,R_2), ..., [L_n,R_n]
        survival <- matrix(0, nrow = 1, ncol = tlen)

        # c(1, S_{1}(R_{1})/S_{1}(L_{1}),...,S_{n-1}(R_{n-1})/S_{n-1}(L_{n-1}))
        survivalR <- matrix(0,nrow = 1, ncol = nj)
        for (b in 1:ntree){
          # on [0, L_1), [L_1,R_1), [L_2,R_2), ..., [L_n,R_n]
          # deal with left truncation ==> all 1 to start, so that Shat_b[, r.ID == 0] == 1
          Shat_b <- matrix(1, nrow = 1, ncol = tlen)

          ShatR_b <- matrix(1, nrow = 1, ncol = nj + 1)
          # S_1(L_1), S_2(L_2), S_3(L_3), ..., S_{nj}(L_{nj})
          qL = rep(0, nj)
          for (j in 1:nj){
            # observations that fall in the same terminal nodes as the new observation in b-th bootstrapped samples
            IDnew <- which(nIDxdata[, b] == nIDxnewdata[newi$I[jall[j]], b])
            # ID of observations in the b-th bootstrap samples
            rw <- which(wt[, b] == 1)
            IDnew <- IDnew[IDnew %in% rw]
            KM <- survival::survfit(formula = formula, data = traindata[IDnew, ])
            Shat_bj <- ipred::getsurv(KM, tpnt[r.ID == jall[j]])

            qL[j] <- Shat_bj[1]
            # S_{j-1}(R_{j-1})
            jR <- ipred::getsurv(KM, newi[, Rname][j])
            # S_{j-1}(R_{j-1})/S_{j-1}(L_{j-1})
            ShatR_b[1, j + 1] = jR / qL[j]
            Shat_b[1, r.ID == jall[j]] <- Shat_bj / qL[j]
          }

          ql0 <- which(qL == 0)
          if(length(ql0) > 0){
            maxqlnot0 <- max(which(qL > 0))
            for(j in ql0){
              if (j < maxqlnot0) {
                ShatR_b[1, j + 1] <- 1
                Shat_b[1, r.ID == jall[j]] <- 1
              } else{
                ShatR_b[1, j + 1] <- 0
                Shat_b[1, r.ID == jall[j]] <- 0
              }
            }
          }

          survival <- survival + Shat_b
          survivalR <- survivalR + ShatR_b[1, 1:nj]
        }

        survival <- survival / ntree
        survivalR <- survivalR / ntree

        m <- cumprod(survivalR[1, ])
        ## construct the survival curve with the averaged ratio
        for (j in 2:nj){
          survival[1, r.ID == jall[j]] = m[j] * survival[1, r.ID == jall[j]]
        }
        RES <- survival[1, ]
      }
      return(RES)
    })

    rm(traindata)
    obj <- Surv(newdata[, yvar.names[1]],
                newdata[, yvar.names[2]],
                newdata[, yvar.names[3]])
  }


  RES = list(survival.probs = pred,
             survival.times = time.eval,
             survival.tau = time.tau,
             survival.obj = obj)
  return(RES)
}
