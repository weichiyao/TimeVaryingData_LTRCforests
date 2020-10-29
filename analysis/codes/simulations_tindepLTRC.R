##########################################################################
## load the LTRCforests package
setwd("./TimeVaryingData_LTRCforests/pkg/LTRCforests")
devtools::load_all()

## load the modified transformation packages 
setwd("./TimeVaryingData_LTRCforests/analysis/utils/transformation/mlt_mod")
devtools::load_all()

setwd("./TimeVaryingData_LTRCforests/analysis/utils/transformation/tram_mod")
devtools::load_all()

setwd("./TimeVaryingData_LTRCforests/analysis/utils/transformation/trtf_mod")
devtools::load_all()

#################################################################
library(ipred)
library(partykit)
library(survival)
library(LTRCforests)
setwd("./TimeVaryingData_LTRCforests/analysis/")
source("./data/TimeindepLTRC_gnrt_PH.R")
source("./data/TimeindepLTRC_gnrt_nonPH.R")
source("./utils/Loss_funct_tindep.R")
predictTSF <- function(object, data, newdata=NULL, OOB=FALSE, tpnt){
  Formula0 = Surv(Start, Stop, Event) ~1
  ntree = length(object$weights)
  
  N = nrow(data)
  if (OOB) {
    weights = matrix(unlist(object$weights),nrow = ntree,byrow = TRUE) ## row: number of trees
    node_all <- predict(object, type = "node")
    rm(object)
    gc()
    
    pred = rep(list(0),N)
    for(wi in 1:N){
      ## find out which trees does not contain wi-th data
      id_oobtree_wi = which(weights[,wi]==0)
      
      weights_wi = rep(0, N)
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
      
      KM = survfit(formula = Formula0, data = data, 
                           weights = weights_wi, subset = weights_wi > 0)
      pred[[wi]] = getSurv(KM, times = tpnt)
    }
  } else {
    weights = matrix(unlist(object$weights),nrow = ntree,byrow = TRUE) ## row: number of trees
    nIDxdata <- predict(object, newdata = data, type = "node") # of size Ndata*ntree
    if (is.null(newdata)){
      nIDxnewdata = nIDxdata
      Nnew = N
    } else {
      nIDxnewdata <- predict(object, newdata = newdata, type = "node") # of size Newdata*ntree
      Nnew = nrow(newdata)
    }
    
    pred = rep(list(0),Nnew)
    for (n_i in 1:Nnew){
      same_nd = rep(0,N)
      for (b in 1:ntree){
        ## ID of observations in the b-th bootstrap samples 
        rw = which(weights[b,]==1)
        ## ID of observations that will fall in the same terminal node as the new observation in b-th tree
        tw <- which(nIDxdata[[b]] == nIDxnewdata[[b]][n_i])
        tw <- tw[tw%in%rw]
        same_nd[tw] <- same_nd[tw]+1/length(tw) # with correct weighting
      }
      KM =survival::survfit(formula=Formula0, data=data, weights=same_nd, subset=same_nd > 0)
      pred[[n_i]] = getSurv(KM, times = tpnt)
    }
  }
  RES = NULL
  RES$survival.probs = pred
  RES$survival.times = tpnt
  RES$survival.tau = rep(max(tpnt), N)
  return(RES)
  
}
#####################################################################################
Pred_funct <- function(N=1000, 
                       model = c("linear", "nonlinear", "interaction"), 
                       censor.rate = 1, 
                       ll, 
                       Nfold=10,
                       setting){

  L2mtry <- data.frame(matrix(0,nrow = 12,ncol = 3))
  names(L2mtry) <- c("cf","rrf","tsf")
  rownames(L2mtry) = c("20","10","5","3","2","1",
                       "20e","10e","5e","3e","2e","1e")
  
  
  OOBmtry <- data.frame(matrix(0,nrow = 6,ncol = 3))
  names(OOBmtry) <- c("cf","rrf","tsf")
  rownames(OOBmtry) = c("20","10","5","3","2","1")
  
  L2 <- data.frame(matrix(0,nrow = 2,ncol = 7))
  names(L2) = c("cx","cfD","cfT","rrfD","rrfT","tsfD","tsfT")
  
  mtryall = data.frame(matrix(0,nrow = 1,ncol = 3))
  names(mtryall) <- c("cf","rrf","tsf")
  
  ibsCVerr = data.frame(matrix(0,nrow = Nfold,ncol = 4))
  names(ibsCVerr) = c("cx","cf","rrf","tsf")
  
  RES = list(L2mtry = L2mtry, OOBmtry = OOBmtry, L2=L2, 
             mtryall=mtryall,
             ibsCVerr = ibsCVerr)
 
  
  if (RES$ibsCVerr[Nfold, 4] == 0){
    ###################### ---------------- Time-varying ----------------- ###############################
    set.seed(101)
    sampleID = sample(10000000,500)
    set.seed(sampleID[ll])
    if(model == "linear"){
      modelnum = 2
    } else if (model == "nonlinear"){
      modelnum = 3
    } else if (model == "interaction"){
      modelnum = 4
    } else{
      stop("wrong model specified")
    }
    if (setting == "PH"){
      RET <- TimeindepLTRC_gnrt_PH(N=N, model = modelnum, censor.rate=censor.rate)
    } else {
      RET <- TimeindepLTRC_gnrt_nonPH(N=N, model = model, censor.rate=censor.rate)
    }
    DATA <- RET$Data
    DATA <- DATA[, c("I","ID","X1","X2","X3","X4","X5","X6","X7","X8","X9","X10",
                     "X11","X12","X13","X14","X15","X16","X17","X18","X19","X20",
                     "Start","Stop","Event","Xi")]
    #### Time points to be evaluated
    Tpnt = sort(unique(round(c(0, DATA$Start, DATA$Stop), 2)))
    Info <- RET$Info
    rm(RET)
    gc()
    ###################### ---------------- noise variables ----------------- ###############################
    ntree = 100L
    Formula = Surv(Start,Stop,Event) ~ X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11+X12+X13+X14+X15+X16+X17+X18+X19+X20
    Formula_TD = Surv(Start,Stop,Event, type = "counting") ~ X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11+X12+X13+X14+X15+X16+X17+X18+X19+X20
  
    Formula_TD0 = Surv(Start,Stop,Event, type = "counting") ~ 1
    
    mtrypool = c(20,10,5,3,2,1)
    Test.obj <- Surv(DATA$Start, DATA$Stop, DATA$Event)
    print("Dataset has been created ...")
  }
  
  leftzero = which(RES$L2mtry$tsf[1:6] == 0)
  if (length(leftzero) > 0){
    for (jj in leftzero){
      #### Different mtry for LTRCTSF
      rex <- Coxph(formula = Formula_TD0, data = DATA, log_first = TRUE)
      
      modelT <- traforest(rex, formula = Formula_TD, data = DATA, ntree = ntree,
                          mtry = mtrypool[jj],
                          control = partykit::ctree_control(teststat = "quad", testtype = "Univ",
                                                            minsplit = max(ceiling(sqrt(nrow(DATA))),20),
                                                            minbucket = max(ceiling(sqrt(nrow(DATA))),7),
                                                            minprob = 0.01,
                                                            mincriterion = 0, saveinfo = FALSE))
      predT <- predictTSF(modelT, data=DATA, tpnt = Tpnt)
      RES$L2mtry$tsf[c(jj,jj+6)] = L2_tindep(KM=predT$survival.probs, Data=DATA, Info=Info, T.pnt=Tpnt)
      rm(predT)
      predOOB <- predictTSF(modelT, data=DATA, OOB = TRUE, tpnt = Tpnt)
      RES$OOBmtry$tsf[jj] <- sbrier_ltrc(Test.obj, pred = predOOB, type = "IBS")
      print(sprintf("TSF--mtry=%1.0f is done",mtrypool[jj]))
      rm(predOOB)
      rm(modelT)
      gc()
    }
  }
  if (RES$L2$tsfT[1] == 0){
    idxmin = which.min(RES$OOBmtry$tsf)
    RES$L2$tsfT = RES$L2mtry$tsf[c(idxmin,idxmin+6)]
    RES$mtryall$tsf = mtrypool[idxmin]
  }
  if (RES$L2$tsfD[1] == 0){
    rex <- Coxph(formula = Formula_TD0, data = DATA, log_first = TRUE)
    
    modelT <- traforest(rex, formula = Formula_TD, data = DATA, ntree = ntree,
                        mtry = ceiling(sqrt(varN)),
                        control = partykit::ctree_control(teststat = "quad", testtype = "Univ",
                                                          mincriterion = 0, saveinfo = FALSE,
                                                          minsplit = 20,
                                                          minbucket = 7,
                                                          minprob = 0.01))
    predT <- predictTSF(modelT, data=DATA, tpnt = Tpnt)
    rm(modelT)
    gc()
    RES$L2$tsfD = L2_tindep(KM = predT$survival.probs, Data=DATA, Info=Info, T.pnt=Tpnt)
    rm(predT)
    gc()
  }
  print("TSF is done ...")
  
  leftzero = which(RES$L2mtry$rrf[1:6] == 0)
  if (length(leftzero) > 0){
    for (jj in leftzero){
      #### Different mtry for LTRCRRF
      modelT = ltrcrrf(formula = Formula, data = DATA,  
                       mtry = mtrypool[jj], 
                       ntree = ntree)
      predT <- predictProb(object = modelT, newdata = DATA, time.eval = Tpnt)
      
      RES$L2mtry$rrf[c(jj,jj+6)] = L2_tindep(KM=predT$survival.probs, Data=DATA, Info=Info, T.pnt=Tpnt)
      rm(predT)
      predOOB = predictProb(object = modelT, time.eval = Tpnt, OOB=TRUE)
      rm(modelT)
      RES$OOBmtry$rrf[jj] = sbrier_ltrc(Test.obj, pred = predOOB, type = "IBS")
      print(sprintf("RRF--mtry=%1.0f is done",mtrypool[jj]))
      rm(predOOB)
      gc()
    }
  }
  if (RES$L2$rrfT[1] == 0){
    idxmin = which.min(RES$OOBmtry$rrf)
    RES$L2$rrfT = RES$L2mtry$rrf[c(idxmin,idxmin+6)]
    RES$mtryall$rrf = mtrypool[idxmin]

  }
  if (RES$L2$rrfD[1] == 0){
    ## Training
    modelT = ltrcrrf(formula = Formula, data = DATA,  
                     mtry = ceiling(sqrt(varN)), 
                     nodesize = 15, 
                     ntree = ntree)
    predT <- predictProb(object = modelT, newdata = DATA, time.eval = Tpnt)
    
    rm(modelT)
    RES$L2$rrfD = L2_tindep(KM=predT$survival.probs, Data=DATA, Info=Info, T.pnt=Tpnt)
    rm(predT)
    gc()
  }
  print("RRF is done ...")
  
  leftzero = which(RES$L2mtry$cf[1:6] == 0)
  if (length(leftzero) > 0){
    for (jj in leftzero){
      #### Different mtry for LTRCCIF
      modelT = ltrccif(formula = Formula, data = DATA, 
                      mtry = mtrypool[jj], 
                      ntree = ntree)
      predT <- predictProb(object = modelT, time.eval = Tpnt)
                           
      RES$L2mtry$cf[c(jj,jj+6)] = L2_tindep(KM=predT$survival.probs, Data=DATA, Info=Info, T.pnt=Tpnt)
      rm(predT)
      predOOB <- predictProb(object = modelT, time.eval = Tpnt, OOB = TRUE)
      RES$OOBmtry$cf[jj] = sbrier_ltrc(Test.obj, pred = predOOB, type = "IBS")
      print(sprintf("CF--mtry=%1.0f is done",mtrypool[jj]))
      rm(predOOB)
      rm(modelT)
      gc()
    }
  }
  if (RES$L2$cfT[1] == 0){
    idxmin = which.min(RES$OOBmtry$cf)
    RES$L2$cfT = RES$L2mtry$cf[c(idxmin,idxmin+6)]
    RES$mtryall$cf = mtrypool[idxmin]
  }
  if (RES$L2$cfD[1] == 0){
    modelT = ltrccif(formula = Formula, data = DATA, 
                    control = partykit::ctree_control(teststat = "quad", testtype = "Univ",
                                                      minsplit = 20,
                                                      minbucket = 7,
                                                      minprob = 0.01,
                                                      mincriterion = 0, saveinfo = FALSE),
                    mtry = ceiling(sqrt(varN)), 
                    ntree = ntree)
    predT <- predictProb(object = modelT, time.eval = Tpnt)
  
    rm(modelT)
    gc()
    RES$L2$cfD = L2_tindep(KM=predT$survival.probs, Data=DATA, Info=Info, T.pnt=Tpnt)
    rm(predT)
    gc()
  }
  print("CF is done ...")
  
  if (RES$L2$cx[1]==0){
    #### Training
    Coxfit <- coxph(formula = Formula, DATA)
    #### Prediction
    predT <- survfit(Coxfit, newdata = DATA)
    rm(Coxfit)
    RES$L2$cx = L2_tindep(KM=predT, Data=DATA, Info=Info, T.pnt=Tpnt)
    rm(predT)
    gc()
  }
  print("Cox is done ...")
  
  if (RES$ibsCVerr[Nfold,4]==0){
    id_uniq = unique(DATA$ID)
    set.seed(sampleID[ll])
    kfold = kfoldcv(k=Nfold,N=N)
    
    dataCV=rep(list(0),Nfold)
    id_uniq_spd=sample(id_uniq,N)
    for (kk in 1:Nfold){
      dataCV[[kk]] = id_uniq_spd[((kk-1)*kfold[kk]+1):(kk*kfold[kk])]
    }
    
    foldleft = which(RES$ibsCVerr[,4]==0)
    for (jj in foldleft){
      idxCV = which(DATA$ID %in% dataCV[[jj]])
      data_b = DATA[-idxCV,]
      newdata_b = DATA[idxCV,]
      Test.obj <- Surv(newdata_b$Start, newdata_b$Stop, newdata_b$Event)
      tpntt <- sort(unique(round(c(0, newdata_b$Start, newdata_b$Stop)), 3))
      ## Cox
      if (RES$ibsCVerr[jj,1]==0){
        Coxfit <- coxph(formula = Formula, data_b)
        #### Prediction
        predCox0 <- survfit(Coxfit, newdata = newdata_b)
        predT = NULL
        predT$survival.probs = rep(list(0),nrow(newdata_b))
        for( i in 1:nrow(newdata_b) ){
          predT$survival.probs[[i]] <- getSurv(predCox0[i], tpntt)
        }
        predT$survival.times <- tpntt
        predT$survival.tau <- rep(max(tpntt), nrow(newdata_b))
        RES$ibsCVerr[jj,1] = sbrier_ltrc(Test.obj, pred = predT, type = "IBS")
        rm(Coxfit)
        rm(predCox0)
        rm(predT)
      }
      ## cfT
      if (RES$ibsCVerr[jj,2]==0){
        mtryT = RES$mtryall$cf
        modelT = ltrccif(formula = Formula, data = data_b,  
                        mtry = mtryT, 
                        ntree = ntree)
        predT <- predictProb(object = modelT, newdata = newdata_b, time.eval = tpntt)
        
        RES$ibsCVerr[jj,2] = sbrier_ltrc(Test.obj, pred = predT, type = "IBS")
        rm(modelT)
        gc()
      }
      ## rrfT
      if (RES$ibsCVerr[jj,3]==0){
        mtryT = RES$mtryall$rrf
        modelT = ltrcrrf(formula = Formula, data = data_b,  
                         mtry = mtryT, 
                         ntree = ntree)
        predT <- predictProb(object = modelT, newdata = newdata_b, time.eval = tpntt)
        
        rm(modelT)
        RES$ibsCVerr[jj,3] = sbrier_ltrc(Test.obj, pred = predT, type = "IBS")
        rm(modelT)
        gc()
      }
      ## tsfT
      if (RES$ibsCVerr[jj,4]==0){
        mtryT = RES$mtryall$tsf
        rex <- Coxph(formula = Formula_TD0, data = data_b, log_first = TRUE)
        modelT <- traforest(rex, formula = Formula_TD, data = data_b, ntree = ntree,
                            mtry = mtryT,
                            control = partykit::ctree_control(teststat = "quad", testtype = "Univ",
                                                              minsplit = max(ceiling(sqrt(nrow(data_b))),20),
                                                              minbucket = max(ceiling(sqrt(nrow(data_b))),7),
                                                              minprob = 0.01,
                                                              mincriterion = 0, saveinfo = FALSE))
        predT <- predictTSF(modelT, data=data_b, newdata=newdata_b, tpnt = tpntt)
        rm(modelT)
        RES$ibsCVerr[jj,4] <- sbrier_ltrc(Test.obj, pred = predT, type = "IBS")
        rm(modelT)
        gc()
      }
      print(sprintf("IBSCV Round %1.0f is done",jj))
    }
    
  }
  return(RES)
}
############################################################################################
nndata = c(100,300,500,1000) 
Nfold = 10 # 10-fold IBS-based CV

setting = "nonPH" # or "PH"
model = "linear" # "nonlinear" or "interaction"
nn = 1 # nndata = c(50, 100, 300, 500) 

# ll-th iteration
ll = 1
RES <- Pred_funct(N = 20, model = model,
                  ll = ll, setting = setting, Nfold = Nfold)
##########################################################################################
