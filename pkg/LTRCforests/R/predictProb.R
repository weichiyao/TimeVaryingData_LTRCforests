#' Compute a Survival Curve from a LTRCCIF model or a LTRCRRF model
#'
#' Constructs a monotone nonincreasing estimated survival curve from a LTRCCIF model or a 
#' LTRCRRF model for any given (left-truncated) right-censored survival data with time-varying 
#' covariates.
#' It can also compute survival function estimates for left-truncated right-censored data
#' with time-invariant covariates.
#'
#' @param object an object as returned by \code{\link{ltrccif}} or by \code{\link{ltrcrrf}}.
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
#' computed survival probabilities for the \emph{i}-th data of interest (i.e., only computes
#' survival probabilies at \code{time.eval[time.eval <= time.tau[i]]} for the \emph{i}-th
#' data of interest). If \code{OOB = TRUE}, the length of \code{time.tau} is equal to the size of
#' \code{data} used to train the \code{object};
#' If \code{OOB = FALSE}, the length of \code{time.tau} is equal to the size
#' of \code{newdata}, or equal to the size of \code{data} if \code{newdata} is not given.
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
#' @import prodlim
#' @aliases predictProb.ltrccif, predictProb.ltrcrrf
#' @seealso \code{\link{sbrier_ltrc}} for evaluation of model fit
#' @examples
#' #### Example with data pbcsample
#' library(survival)
#' Formula <- Surv(Start, Stop, Event) ~ age + alk.phos + ast + chol + edema
#' ## Fit an LTRC conditional inference forest on time-varying data
#' LTRCCIFobj <- ltrccif(formula = Formula, data = pbcsample, id = ID,
#'                       mtry = 3, ntree = 50L)
#'
#'
#' ## Construct an estimated survival estimate for the second subject
#' tpnt <- seq(0, max(pbcsample$Stop), length.out = 50)
#' newData <- pbcsample[pbcsample$ID == 2, ]
#' Pred <- predictProb(object = LTRCCIFobj, newdata = newData, newdata.id = ID,
#'                     time.eval = tpnt)
#' ## Since time.tau = NULL, Pred$survival.probs is in the matrix format, with dimensions:
#' dim(Pred$survival.probs) # length(time.eval) x nrow(newdata)
#' ## Plot the estimated survival curve
#' plot(Pred$survival.times, Pred$survival.probs, type = "l", col = "red",
#'      xlab = "Time", ylab = "Survival probabilities")
#'
#'
#'
#' @export
#' 
predictProb <- function(object, newdata = NULL, newdata.id, OOB = FALSE,
                        time.eval, time.tau = NULL){
  UseMethod("predictProb", object)
}
#' @export
predictProb.ltrccif <- function(object, newdata = NULL, newdata.id, OOB = FALSE,
                                time.eval, time.tau = NULL){
  
  pred <- partykit::predict.cforest(object = object, newdata = newdata, OOB = OOB, type = "prob",
                                    FUN = .pred_Surv_nolog)
  xvar.names <- attr(object$terms,"term.labels")
  yvar.names <- as.character(object$formulaLTRC[[2]])[2:4]
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
  
  ################ Changes at July 29th
  tpntL <- c(TestT, tpnt)
  torder <- order(tpntL)
  tpntLod <- tpntL[torder]
  tlen <- length(tpntLod)
  
  ## Compute the estimated survival probability of the Ni-th subject
  Shat_temp <- matrix(0, nrow = 1, ncol = tlen)
  
  r.ID <- findInterval(tpntLod, TestT)
  r.ID[r.ID > TestTIntN] <- TestTIntN
  
  jall <- unique(r.ID[r.ID > 0])
  nj <- length(jall)
  
  ## Deal with left-truncation
  Shat_temp[1, r.ID == 0] <- 1
  if(nj == 1){
    ## Get the index of the Pred to compute Shat
    II <- which(id.seu == id.sub[Ni])[jall[nj]]
    Shat_i = ipred::getsurv(pred[[II]], tpntLod[r.ID == jall[nj]])
    Shat_temp[1, r.ID == jall[nj]] <- Shat_i / Shat_i[1]
  } else if (nj > 1) {
    # c(1, S_{1}(R_{1}), ..., S_{nj}(R_{nj}))
    ShatR_temp <- matrix(0, nrow = 1, ncol = nj + 1)
    ShatR_temp[1, 1] <- 1
    
    # S_1(L_1), S_2(L_2), S_3(L_3), ..., S_{nj}(L_{nj})
    qL = rep(0, nj)
    for (j in 1:nj){
      ## Get the index of the Pred to compute Shat
      II <- which(id.seu == id.sub[Ni])[1] + jall[j] - 1
      Shat_j = ipred::getsurv(pred[[II]], tpntLod[r.ID == jall[j]])
      
      qL[j] <- Shat_j[1]
      # S_{j}(R_{j}), j=1,...nj-1
      jR = ipred::getsurv(pred[[II]], TestT[j + 1])
      ShatR_temp[1, j + 1] = jR / qL[j]
      Shat_temp[1, r.ID == jall[j]] <- Shat_j / qL[j]
    }
    
    ql0 <- which(qL == 0)
    if (length(ql0) > 0){
      if (any(qL > 0)){
        maxqlnot0 <- max(which(qL > 0))
        
        ql0lmax <- ql0[ql0 < maxqlnot0]
        ql0mmax <- ql0[ql0 >= maxqlnot0]
        ShatR_temp[1, ql0lmax + 1] <- 1
        Shat_temp[1, r.ID %in% jall[ql0lmax]] <- 1
        ShatR_temp[1, ql0mmax + 1] <- 0
        Shat_temp[1, r.ID %in% jall[ql0mmax]] <- 0
        # for(j in ql0){
        #   if (j < maxqlnot0) {
        #     ShatR_temp[1, j + 1] <- 1
        #     Shat_temp[1, r.ID == jall[j]] <- 1
        #   } else{
        #     ShatR_temp[1, j + 1] <- 0
        #     Shat_temp[1, r.ID == jall[j]] <- 0
        #   }
        # }
      } else {
        ShatR_b[1, 2:(nj + 1)] <- 0
        Shat_temp[1, r.ID %in% jall] <- 0
      }
    }
    m <- cumprod(ShatR_temp[1, 1:nj])
    for (j in 1:nj){
      Shat_temp[1, r.ID == jall[j]] <- Shat_temp[1, r.ID == jall[j]] * m[j]
    }
  }
  
  # since: tpntLod[torder == 1] == TestData[1, 1]
  return(Shat_temp[1, -match(TestT, tpntLod)])
  rm(Shat_temp)
  rm(ShatR_temp)
  rm(id.seu)
  rm(id.sub)
}

.pred_Surv_nolog <- function(y, w) {
  if (length(y) == 0) return(NA)
  survfit(y ~ 1, weights = w, subset = w > 0, conf.type = "none", se.fit = FALSE)
}

#' @export
predictProb.ltrcrrf <- function(object, newdata = NULL, newdata.id, OOB = FALSE,
                                time.eval, time.tau = NULL){
  ntree <- object$ntree
  formula <- object$formulaLTRC
  formula[[3]] <- 1
  
  wt <- object$inbag # of size Ndata x ntree
  
  yvar.names <- object$yvarLTRC.names[2:4]
  traindata <- object$yvarLTRC
  
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
      
      ################ Changes at July 29th
      newiIntT <- c(newi[1, yvar.names[1]], newi[, yvar.names[2]])
      tpntL <- c(newiIntT, tpnt)
      torder <- order(tpntL)
      tpntLod <- tpntL[torder]
      tlen <- length(tpntLod)
      
      r.ID <- findInterval(tpntLod, newiIntT)
      r.ID[r.ID > n_newi] <- n_newi
      
      jall <- unique(r.ID[r.ID > 0])
      nj <- length(jall)
      if (nj == 1){
        survival <- matrix(0, nrow = 1, ncol = tlen)
        ## deal with left-truncation
        survival[1, r.ID == 0] <- 1
        
        ## find out which trees does not contain the I[jall[j]]-th wi-th data
        id_tree_wi_j <- which(wt[newi$I[jall[1]], ] == 0)
        ## add the if-else at Sept 16th, for id_tree_wi_j == integer(0)
        if (length(id_tree_wi_j) > 0){
          Shat_ti <- matrix(0, ncol = length(tpntLod[r.ID == jall[1]]), nrow = length(id_tree_wi_j))
          for (ti in 1:length(id_tree_wi_j)){
            ## In each tree of id in idTree_wi, it falls into terminal id_node_witi_j
            id_node_witi_j <- node_all[newi$I[jall[1]], id_tree_wi_j[ti]]
            ## id of samples that fall into the same node
            id_samenode_witi_j <- which(node_all[, id_tree_wi_j[ti]] == id_node_witi_j)
            ## Pick out those appearing in the bootstrapped samples
            id_inbag_j <- which(wt[, id_tree_wi_j[ti]] > 0)
            id_buildtree_j <- id_samenode_witi_j[id_samenode_witi_j %in% id_inbag_j]
            ## Build the survival tree
            traindata$KMwt <- 0
            traindata$KMwt[id_buildtree_j] = wt[id_buildtree_j, id_tree_wi_j[ti]] / sum(wt[id_buildtree_j, id_tree_wi_j[ti]])
            KMwt <- traindata$KMwt
            KM <- survival::survfit(formula = formula, data = traindata, se.fit = FALSE,
                                    weights = KMwt, subset = KMwt > 0, conf.type = "none")
            ## Get survival probabilities
            ## Changed at July 29th
            Shat_ti[ti, ] <- ipred::getsurv(KM, tpntLod[r.ID == jall[1]])
            
          }
          ## Changed at July 29th
          rowid.nz <- which(Shat_ti[, 1] != 0)
          Shat_ti[rowid.nz, ] <- Shat_ti[rowid.nz, ] / Shat_ti[rowid.nz, 1]
          # if (length(rowid.nz) > 0) {
          #   Shat_ti[rowid.nz, ] <- sweep(Shat_ti[rowid.nz, ], 1, Shat_ti[rowid.nz, 1], "/")
          # }
          survival[1, r.ID == jall[1]] <- apply(Shat_ti, 2, mean)
        } else {
          ## Changed at Sept 16th, for id_tree_wi_j == integer(0)
          KM <- survival::survfit(formula = formula, data = traindata, se.fit = FALSE,
                                  conf.type = "none")
          survival[1, r.ID == jall[1]] <- ipred::getsurv(KM, tpntLod[r.ID == jall[1]])
        }
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
            id_inbag_j <- which(wt[, id_tree_wi_j[ti]] > 0)
            id_buildtree_j <- id_samenode_witi_j[id_samenode_witi_j %in% id_inbag_j]
            ## Build the survival tree
            traindata$KMwt <- 0
            traindata$KMwt[id_buildtree_j] = wt[id_buildtree_j, id_tree_wi_j[ti]] / sum(wt[id_buildtree_j, id_tree_wi_j[ti]])
            KMwt <- traindata$KMwt
            KM <- survival::survfit(formula = formula, data = traindata, se.fit = FALSE,
                                    weights = KMwt, subset = KMwt > 0, conf.type = "none")
            Shat_ti <- ipred::getsurv(KM, tpntLod[r.ID == jall[j]])
            
            if (Shat_ti[1] == 0){# jL = Shat_ti[1]
              survival[1, r.ID == jall[j]] <- survival[1, r.ID == jall[j]] + 1
              survivalR[1, j + 1] <- survivalR[1, j + 1] + 1
            } else {
              survival[1, r.ID == jall[j]] <- survival[1, r.ID == jall[j]] + Shat_ti / Shat_ti[1]
              survivalR[1, j + 1] <- survivalR[1, j + 1] + ipred::getsurv(KM, newi[, yvar.names[2]][j]) / Shat_ti[1]
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
       
      }
      
      RES <- survival[1, -match(newiIntT, tpntLod)]
      return(RES)
    })
    obj <- Surv(traindata[, yvar.names[1]],
                traindata[, yvar.names[2]],
                traindata[, yvar.names[3]])
    id <- traindata$id
    rm(traindata)
  } else {
    nIDxdata <- object$membership # of size Ndata*ntree
    
    if (is.null(newdata)){
      newdata <- traindata
      nIDxnewdata <- nIDxdata
    } else {
      if (!is.data.frame(newdata)) stop("newdata must be a dataframe")
      x.IDs <- match(object$xvar.names, names(newdata))
      nIDxnewdata <- predict.ltrcrfsrc(object, newdata = newdata[, x.IDs], membership = TRUE)$membership # of size Newdata*ntree
      if (missing(newdata.id)){
        newdata$id <- 1:nrow(newdata)
      } else {
        names(newdata)[names(newdata) == deparse(substitute(newdata.id))] <- "id"
      }
      if (any(is.na.data.frame(newdata[, x.IDs]))) stop("newdata contains missing values in the covariates!")
    }
    
    rm(object)
    obj.IDs <- match(c("id", yvar.names), names(newdata), nomatch = 0)
    if (any(obj.IDs == 0)) stop("newdata has to be with variables as in formula (time1, time2, event)")
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
      ################ Changes at July 29th
      newiIntT <- c(newi[1, yvar.names[1]], newi[, yvar.names[2]])
      
      tpntL <- c(newiIntT, tpnt)
      torder <- order(tpntL)
      # torder == 1 corresponds with TestData[1, 1]: tpntLod[torder == 1] == TestData[1, 1]
      tpntLod <- tpntL[torder]
      tlen <- length(tpntLod)
      
      ################ Changes at July 29th
      # deal with left truncation in the training data
      r.ID <- findInterval(tpntLod, newiIntT)
      r.ID[r.ID > n_newi] <- n_newi
      jall <- unique(r.ID[r.ID > 0])
      nj <- length(jall)
      
      if (nj == 1){
        # on [L_1,R_1), [L_2,R_2), ..., [L_n,R_n]
        survival <- matrix(0, nrow = 1, ncol = tlen)
        # deal with left truncation
        survival[1, r.ID == 0] <- 1
        
        nlenb <- length(tpntLod[r.ID == jall[1]])
        Shat_b <- matrix(0, nrow = ntree, ncol = nlenb)
        for (b in 1:ntree){
          # observations that fall in the same terminal nodes as the new observation in b-th bootstrapped samples
          IDnew <- which(nIDxdata[, b] == nIDxnewdata[newi$I[jall[1]], b])
          # ID of observations in the b-th bootstrap samples
          rw <- which(wt[, b] > 0)
          IDnew <- IDnew[IDnew %in% rw]
          ## For each tree, we need to reset KMwt to be zero,
          ## otherwise subset = KMwt > 0 is not correct
          traindata$KMwt <- 0
          traindata$KMwt[IDnew] = wt[IDnew, b] / sum(wt[IDnew, b])
          KMwt <- traindata$KMwt
          KM <- survival::survfit(formula = formula, data = traindata, se.fit = FALSE,
                                  weights = KMwt, subset = KMwt > 0, conf.type = "none")
          Shat_b[b, ] <- ipred::getsurv(KM, tpntLod[r.ID == jall[1]]) # jall[nj] = 1
          # NaN problem if Shat_b[1] = 0
          # survival[1, r.ID == jall[1]] <- survival[1, r.ID == jall[1]] + Shat_b / Shat_b[1]
        }
        rowid.nz <- which(Shat_b[, 1] != 0)
        Shat_b[rowid.nz, ] <- Shat_b[rowid.nz, ] / Shat_b[rowid.nz, 1]
        # if (length(rowid.nz) > 0){
        #   Shat_b[rowid.nz, ] <- sweep(Shat_b[rowid.nz, ], 1, Shat_b[rowid.nz, 1], "/")
        # }
        survival[1, r.ID == jall[1]] <- apply(Shat_b, 2, mean)
        # Shat_b = apply(Shat_b, 2, mean)
        # survival[1, r.ID == jall[1]] <- Shat_b / Shat_b[1]
      } else if (nj > 1) {
        # on [0, L_1), [L_1,R_1), [L_2,R_2), ..., [L_n,R_n]
        survival <- matrix(0, nrow = 1, ncol = tlen)
        
        # c(1, S_{1}(R_{1})/S_{1}(L_{1}),...,S_{n-1}(R_{n-1})/S_{n-1}(L_{n-1}))
        survivalR <- matrix(0, nrow = 1, ncol = nj)
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
            rw <- which(wt[, b] > 0)
            IDnew <- IDnew[IDnew %in% rw]
            
            ## For each tree, we need to reset KMwt to be zero,
            ## otherwise subset = KMwt > 0 is not correct
            traindata$KMwt <- 0
            traindata$KMwt[IDnew] = wt[IDnew, b] / sum(wt[IDnew, b])
            KMwt <- traindata$KMwt
            KM <- survival::survfit(formula = formula, data = traindata, se.fit = FALSE,
                                    weights = KMwt, subset = KMwt > 0, conf.type = "none")
            Shat_bj <- ipred::getsurv(KM, tpntLod[r.ID == jall[j]])
            
            qL[j] <- Shat_bj[1]
            # S_{j-1}(R_{j-1})
            jR <- ipred::getsurv(KM, newi[, yvar.names[2]][j])
            # S_{j-1}(R_{j-1})/S_{j-1}(L_{j-1})
            ShatR_b[1, j + 1] = jR / qL[j]
            Shat_b[1, r.ID == jall[j]] <- Shat_bj / qL[j]
          }
          
          ql0 <- which(qL == 0)
          if (length(ql0) > 0){
            if (any(qL > 0)){
              maxqlnot0 <- max(which(qL > 0))

              ql0lmax <- ql0[ql0 < maxqlnot0]
              ql0mmax <- ql0[ql0 >= maxqlnot0]
              ShatR_b[1, ql0lmax + 1] <- 1
              Shat_b[1, r.ID %in% jall[ql0lmax]] <- 1
              ShatR_b[1, ql0mmax + 1] <- 0
              Shat_b[1, r.ID %in% jall[ql0mmax]] <- 0

            } else {
              ShatR_b[1, 2:(nj + 1)] <- 0
              Shat_b[1, r.ID %in% jall] <- 0
            }
          }
          
          
          survival <- survival + Shat_b
          survivalR <- survivalR + ShatR_b[1, 1:nj]
        }
        
        survival <- survival / ntree
        survivalR <- survivalR / ntree

        m <- cumprod(survivalR[1, ])
        # construct the survival curve with the averaged ratio
        for (j in 1:nj){
          survival[1, r.ID == jall[j]] = m[j] * survival[1, r.ID == jall[j]]
        }
        
      }
      RES <- survival[1, -match(newiIntT, tpntLod)]
      return(RES)
    })
    
    rm(traindata)
    
    obj <- Surv(newdata[, yvar.names[1]],
                newdata[, yvar.names[2]],
                newdata[, yvar.names[3]])
    id <- newdata$id
    
  }
  
  RES = list(survival.id = id,
             survival.probs = pred,
             survival.times = time.eval,
             survival.tau = time.tau,
             survival.obj = obj)
  return(RES)
}

predictProb.prodlim <- function(object, time.eval, ...) {
  predict(object, times = time.eval, type  = "surv")
}