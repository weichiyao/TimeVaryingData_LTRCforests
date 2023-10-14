###################################################################################
tsf_wy <- function(formula, data, ntree = 100L, mtry = NULL, log_first = TRUE, 
                   control = partykit::ctree_control(teststat = "quad", testtype = "Univ",
                                                     minsplit = max(ceiling(sqrt(nrow(data))),20),
                                                     minbucket = max(ceiling(sqrt(nrow(data))),7),
                                                     minprob = 0.01,
                                                     mincriterion = 0, saveinfo = FALSE),
                   bootstrap = "by.sub", stepFactor = 2, trace = TRUE){
  formula0 = formula
  formula0[[3]] = 1
  
  
  if (is.null(mtry)){
    mtry <- tune_tsf_wy(formula = formula, data = data, ntreeTry = ntree,
                        log_first = log_first, control = control, 
                        bootstrap = bootstrap, stepFactor = stepFactor,
                        trace = trace)
  } 
  
  rex <- Coxph(formula = formula0, data = data, log_first = log_first)
  if (bootstrap == "by.sub"){
    id_uniq <- unique(data$ID)
    n_uniq <- length(id_uniq)
    size = floor(n_uniq * 0.632)
    wt <- matrix(0, nrow = nrow(data), ncol = ntree)
    
    weights <- sapply(1:ntree, function(itree){
      wt = rep(0,nrow(data))
      idx <- sample(id_uniq, size = size, replace = FALSE)
      for (ii in 1:size) {
        inidx <- which(data$ID == idx[ii])
        wt[inidx] = 1
      }
      return(wt)
    })
    TSFfsub1 <- traforest(rex, formula = formula, data = data, ntree = ntree,
                          weights = weights, mtry = mtry, 
                          control = control)
  } else if (bootstrap == "by.root") {
    TSFfsub1 <- traforest(rex, formula = formula, data = data, ntree = ntree, 
                          mtry = mtry, 
                          control = control)
    
  } else {
    stop("wrong bootstrap is set")
  }
  
  TSFfsub1 
}
###################################################################################
predict_tsf_wy <- function(object, newdata=NULL, ensemble = TRUE, OOB = FALSE, varying = TRUE){
  if (OOB) {
    DATA = data.frame(as.matrix(object$data))
    ntree = length(object$weights)
    weights = matrix(unlist(object$weights),nrow = ntree,byrow = TRUE) ## row: number of trees
    if (ncol(DATA) > 2){
      names(DATA)[1:3] = c("start","stop","status")
    }  
    node_all <- predict(object, type = "node")
    
    rm(object)
    gc()
    
    pred = rep(list(0),nrow(DATA))
    for(wi in 1:nrow(DATA)){
      ## find out which trees does not contain wi-th data
      id_oobtree_wi = which(weights[,wi]==0)
      
      weights_wi = rep(0, nrow(DATA))
      for (ti in 1:length(id_oobtree_wi)){
        ## In each tree of id in id_oobtree_wi, it falls into terminal id_node_witi
        id_node_witi= node_all[[id_oobtree_wi[ti]]][wi]
        ## id of samples that fall into the same node
        id_samenode_witi = which(node_all[[id_oobtree_wi[ti]]]==id_node_witi)
        ## Pick out those appearing in the bootstrapped samples 
        id_inbag = which(weights[id_oobtree_wi[ti],]==1)
        id_buildtree = id_samenode_witi[id_samenode_witi %in% id_inbag]
        weights_wi[id_buildtree] = weights_wi[id_buildtree] + 1/length(id_buildtree) # with correct weighting
      }
      if (ncol(DATA) == 2){
        pred[[wi]] = survfit(formula = Surv(time,status)~1, data = DATA, 
                             weights = weights_wi, subset = weights_wi > 0)
      } else {
        pred[[wi]] = survfit(formula = Surv(start,stop,status)~1, data = DATA, 
                             weights = weights_wi, subset = weights_wi > 0)
      }
      
    }
  } else {
    if (is.null(newdata)){
      newdata=object$data
    } 
    nIDxdata <- predict(object, newdata = object$data, type = "node") # of size Ndata*ntree
    nIDxnewdata <- predict(object, newdata = newdata, type = "node") # of size Newdata*ntree
    if (ensemble) {
      ntree = length(object$weights)
      weights <- object$weights
      
      pred = lapply(1:nrow(newdata), function(newdata_i){
        same_nd <- rep(0,nrow(object$data))
        for (b in 1:ntree){
          ## ID of observations in the b-th bootstrap samples 
          rw = which(weights[[b]]==1)
          ## ID of observations that will fall in the same terminal node as the new observation in b-th tree
          tw <- which(nIDxdata[[b]] == nIDxnewdata[[b]][newdata_i])
          tw <- tw[tw%in%rw]
          same_nd[tw] <- same_nd[tw]+1/length(tw) # with correct weighting
        }
        if (ncol(as.matrix(object$data[,1])) == 2){
          survival::survfit(Surv(time,status)~1, data = data.frame(as.matrix(object$data[,1])), 
                            weights = same_nd, subset = same_nd > 0)
        }else {
          survival::survfit(Surv(start,stop,status)~1, data = data.frame(as.matrix(object$data[,1])), 
                            weights = same_nd, subset = same_nd > 0)
        }
      })
    } else {
      pred = lapply(1:nrow(newdata), function(newdata_i){
        same_nd <- rep(0,nrow(object$data))
        same_nd[nIDxdata == nIDxnewdata[newdata_i]] = 1
        same_nd = same_nd/sum(same_nd)
        if (ncol(as.matrix(object$data[,1])) == 2){
          survival::survfit(Surv(time,status)~1, data = data.frame(as.matrix(object$data[,1])), 
                            weights = same_nd, subset = same_nd > 0)
        }else {
          survival::survfit(Surv(start,stop,status)~1, data = data.frame(as.matrix(object$data[,1])), 
                            weights = same_nd, subset = same_nd > 0)
        }
      })
    }
  }
  
  return(pred)
}

###################################################################################
tune_tsf_wy <- function(formula, data, mtryStart = NULL, stepFactor = 1.5, ntreeTry = 100L, 
                        bootstrap = c("by.sub", "by.root"),
                        control = partykit::ctree_control(teststat = "quad", testtype = "Univ",
                                                          minsplit = max(ceiling(sqrt(nrow(data))),20),
                                                          minbucket = max(ceiling(sqrt(nrow(data))),7),
                                                          minprob = 0.01,
                                                          mincriterion = 0, saveinfo = FALSE),
                        log_first = TRUE, trace = TRUE, plot = FALSE) {
  
  # number of the variables
  allX <- substring(formula,1)[[3]]
  nameX <- strsplit(gsub("\\+", '', allX)," ")[[1]]
  nameX <- nameX[nameX!=""]
  nvar = length(nchar(nameX))
  if (is.null(mtryStart)){
    mtryStart = ceiling(sqrt(nvar))
  } 
  
  formula0 <- formula
  formula0[[3]] <- 1
  
  bootstrap = bootstrap[1]
  
  id_uniq = sort(unique(data$ID))
  n_uniq = length(id_uniq)
  Tpnt = c(0, sort(unique(data$Stop)), seq(max(data$Stop) + 5, 1.5*max(data$Stop), length.out = 50)[-1])
  Tau <- sapply(1:n_uniq, function(ii){
    1.5 * max(data[data$ID == id_uniq[ii], ]$Stop)
  })
  Sobj = Surv(data$Start, data$Stop, data$Event)
  # integrated Brier score of out-of-bag samples for a mtry value at test 
  errorOOB_mtry <- function(formula, data, mtryTest, ntreeTry, 
                            control, bootstrap, log_first, tpnt, tau,
                            sobj){
    tsfOOB <- tsf_wy(formula = formula, data = data, ntree = ntreeTry, mtry = mtryTest,
                     control = control, bootstrap = bootstrap, log_first = log_first)
    
    predOOB <- predict_tsf_wy(object = tsfOOB, OOB = TRUE)
    
    Shat = shat(data = data, pred = predOOB, tpnt = tpnt)
    predOOBnew = list(survival.probs = Shat, survival.times = tpnt, survival.tau = tau)
    rm(Shat)
    rm(predOOB)
    errorOOB <- sbrier_ltrc(obj = sobj, id = data$ID, pred = predOOBnew)
    rm(predOOBnew)
    gc()
    return(errorOOB)
  }
  
  
  # errorOld
  errorOld <- errorOOB_mtry(formula = formula, data = data, mtryTest = mtryStart, 
                            ntreeTry = ntreeTry,  
                            control = control, bootstrap = bootstrap, 
                            log_first = log_first, 
                            tpnt = Tpnt, tau = Tau,
                            sobj = Sobj)
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
      
      errorCur <- errorOOB_mtry(formula = formula, data = data, mtryTest = mtryCur, 
                                ntreeTry = ntreeTry,  
                                control = control, bootstrap = bootstrap, 
                                log_first = log_first, 
                                tpnt = Tpnt, tau = Tau,
                                sobj = Sobj)
      
      if (trace) {
        cat("mtry =",mtryCur, "\tOOB error =",errorCur, "\n")
      }
      oobError[[as.character(mtryCur)]] <- errorCur
      errorOld <- errorCur
    }
  }
  mtry <- sort(as.numeric(names(oobError)))
  res_all <- unlist(oobError[as.character(mtry)])
  res_all <- cbind(mtry = mtry, OOBError = res_all)
  res <- res_all[which.min(res_all[,2]), 1]
  
  if (plot) {
    res = res_all
    plot(res_all, xlab = expression(m[try]), ylab = "OOB Error", type = "o", log = "x", xaxt = "n")
    axis(1, at = res_all[, "mtry"])
  }
  
  return(res)
}
