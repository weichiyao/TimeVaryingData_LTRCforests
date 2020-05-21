##########################################################################
## recompile the package randomForestSRC with the updated splitCustom.c
setwd("./TimeVaryingData_LTRCforests/analysis/utils/randomForestSRC/")
devtools::load_all()

library(trtf)
library(tram)
library(mlt)
library(ipred)
library(partykit)
library(survival)
library(prodlim)
library(pec)

setwd("./TimeVaryingData_LTRCforests/analysis/data/")
source("./Timefixed_gnrt_PH.R")
source("./Timefixed_gnrt_nonPH.R")
setwd("../utils/")
source("./Loss_funct_tfixed.R")
#####################################################################################
predictRSF <- function(object, data, newdata=NULL, tpnt, OOB=FALSE){
  Formula0 = Surv(Stop, Event) ~ 1
  weights <- object$inbag #of size Ndata x ntree
  ntree = ncol(weights)
  tlen = length(tpnt)
  N = nrow(data)
  
  if (OOB){
    node_all <- object$membership # of size Ndata*ntree
    
    predrsf <- sapply(1:N, function(wi){
      survival <- rep(0, tlen)
      ## find out which trees does not contain the wi-th data
      id_tree_wi_j = which(weights[wi,]==0)
      
      for (ti in 1:length(id_tree_wi_j)){
        ## In each tree of id in idTree_wi, it falls into terminal id_node_witi_j
        id_node_witi_j= node_all[wi,id_tree_wi_j[ti]]
        ## id of samples that fall into the same node
        id_samenode_witi_j = which(node_all[,id_tree_wi_j[ti]]==id_node_witi_j)
        ## Pick out those appearing in the bootstrapped samples 
        id_inbag_j = which(weights[,id_tree_wi_j[ti]]==1)
        id_buildtree_j = id_samenode_witi_j[id_samenode_witi_j %in% id_inbag_j]
        ## Build the survival tree
        KM <- survival::survfit(formula = Formula0, data = data[id_buildtree_j,])
        ## Get survival probabilities
        survival <- survival + ipred::getsurv(KM, tpnt)
      }
      
      survival = survival/length(id_tree_wi_j)
      return(survival)
    })
    
    ## cannot force the class to be "matrix"
  } else {
    nIDxdata <- object$membership # of size Ndata*ntree
    if (is.null(newdata)){
      nIDxnewdata = nIDxdata
      Nnew = N
    } else {
      nIDxnewdata <- predict(object, newdata = newdata, membership = TRUE)$membership # of size Newdata*ntree
      Nnew = nrow(newdata)
    }
    predrsf <- sapply(1:Nnew, function(j){
      pred <- rep(0,length(tpnt))
      for (b in 1:ntree){
        ## ID of observations in the b-th bootstrap samples 
        rw = which(weights[,b]==1)
        ## ID of observations that will fall in the same terminal node as the new observation in b-th tree
        tw <- which(nIDxdata[,b] == nIDxnewdata[j,b])
        #     
        tw <- tw[tw%in%rw]
        KM <- survival::survfit(formula = Formula0, data = data[tw,])
        pred <- getSurv(KM, tpnt) + pred
      }
      pred = pred/ntree
      return(pred)
    })
  }
  
  predrsf
}
predictTSF <- function(object, data, newdata=NULL, OOB=FALSE){
  Formula0 = Surv(Stop, Event) ~ 1
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
      
      pred[[wi]] = survfit(formula = Formula0, data = data, 
                           weights = weights_wi, subset = weights_wi > 0)
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
      pred[[n_i]]=survival::survfit(formula=Formula0, data=data, weights=same_nd, subset=same_nd > 0)
    }
  }
  
  return(pred)
  
}
#####################################################################################
Pred_CV <- function(N = 1000, Distribution = "WI", model = 2, 
                    censor.rate = 1, ll, Nfold=10, setting){
  
  L2mtry <- data.frame(matrix(0,nrow = 7, ncol = 3))
  names(L2mtry) <- c("cf","rsf","tsf")
  rownames(L2mtry) = c("20","10","5","3","2","1","opt")
  
  OOBmtry <- data.frame(matrix(0, nrow = 6, ncol = 3))
  names(OOBmtry) <- c("cf","rsf","tsf")
  rownames(OOBmtry) = c("20","10","5","3","2","1")
  
  L2 <- data.frame(matrix(0,nrow = 2, ncol = 7))
  names(L2) = c("cx","cfD","cfP","rsfD","rsfP","tsfD","tsfP")
  
  mtryall = data.frame(matrix(0, nrow = 1, ncol = 3))
  names(mtryall) <- c("cf", "rsf", "tsf")
  
  ibsCVerr = data.frame(matrix(0, nrow = Nfold, ncol = 4))
  names(ibsCVerr) = c("cx", "cf", "rsf", "tsf")
  
  RES = list(L2mtry = L2mtry, OOBmtry = OOBmtry, 
             L2 = L2, 
             mtryall = mtryall,
             ibsCVerr = ibsCVerr)

  if (RES$ibsCVerr[Nfold, 4]==0){
    ###################### ---------------- Time-varying ----------------- ###############################
    set.seed(101)
    sampleID = sample(10000000,500)
    set.seed(sampleID[ll])
    if (setting == "PH"){
      RET <- Timefixed_gnrt_PH(N = N, model = model, Distribution = Distribution, censor.rate = censor.rate)
    } else {
      RET <- Timefixed_gnrt_nonPH(N = N, model = model, censor.rate = censor.rate)
    }
    DATA <- RET$Data
    DATA <- DATA[,c("I","ID","X1","X2","X3","X4","X5","X6","X7","X8","X9","X10",
                    "X11","X12","X13","X14","X15","X16","X17","X18","X19","X20","Start","Stop","Event","Xi")]
    #### Time points to be evaluated
    Tpnt = c(0, sort(unique(DATA$Stop)))
    Info <- RET$Info
    rm(RET)
    gc()
    ###################### ---------------- noise variables ----------------- ###############################
    Formula = Surv(Stop,Event) ~ X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11+X12+X13+X14+X15+X16+X17+X18+X19+X20
    Formula0 = Surv(Stop,Event) ~ 1
    ntree = 100L
    ## pool of mtry searched in the OOB procedure 
    mtrypool = c(20, 10, 5, 3, 2, 1)
    mtryD <- ceiling(sqrt(20))
    Test.obj <- Surv(DATA$Stop, DATA$Event)
    print("Dataset has been created ...")
  }
  ###################### ---------------- L2 TSF ----------------- ###############################
  leftzero = which(RES$L2mtry$tsf[1:6] == 0)
  if (length(leftzero) > 0){
    for (jj in leftzero){
      #### Different mtry for TSF
      rex <- Coxph(formula = Formula0, data = DATA, log_first = TRUE)
      
      modelT <- traforest(rex, formula = Formula, data = DATA, ntree = ntree,
                          mtry = mtrypool[jj],
                          control = partykit::ctree_control(teststat = "quad", testtype = "Univ",
                                                            minsplit = max(ceiling(sqrt(nrow(DATA))),20),
                                                            minbucket = max(ceiling(sqrt(nrow(DATA))),7),
                                                            minprob = 0.01,
                                                            mincriterion = 0, saveinfo = FALSE))
      predT <- predictTSF(modelT, data = DATA)
      RES$L2mtry$tsf[jj] = L2_tfixed(KM = predT, Data = DATA, Info = Info, T.pnt = Tpnt)
      rm(predT)
      predOOB <- predictTSF(modelT, data = DATA, OOB = TRUE)
      RES$OOBmtry$tsf[jj] <- unname(ipred::sbrier(Test.obj, predOOB)[1])
      print(sprintf("TSF -- mtry = %1.0f is done",mtrypool[jj]))
      rm(predOOB)
      rm(modelT)
      gc()
    }
  }
  if (RES$L2$tsfT[1] == 0){
    idxmin <- which.min(RES$OOBmtry$tsf)
    RES$L2$tsfT <- RES$L2mtry$tsf[idxmin]
    RES$mtryall$tsf <- mtrypool[idxmin]
    RES$L2mtry$tsf[7] <- min(RES$OOBmtry$tsf[1:6])
  }
  
  if (RES$L2$tsfD[1] == 0){
    rex <- Coxph(formula = Formula0, data = DATA, log_first = TRUE)
    
    modelT <- traforest(rex, formula = Formula, data = DATA, ntree = ntree,
                        mtry = mtryD,
                        control = partykit::ctree_control(teststat = "quad", testtype = "Univ",
                                                          mincriterion = 0, saveinfo = FALSE,
                                                          minsplit = 20,
                                                          minbucket = 7,
                                                          minprob = 0.01))
    predT <- predictTSF(modelT, data = DATA)
    rm(modelT)
    gc()
    RES$L2$tsfD = L2_tfixed(KM = predT, Data = DATA, Info = Info, T.pnt = Tpnt)
    rm(predT)
    gc()
  }
  print("TSF is done ...")
  
  ###################### ---------------- L2 RSF ----------------- ###############################
  leftzero = which(RES$L2mtry$rsf[1:6] == 0)
  if (length(leftzero) > 0){
    for (jj in leftzero){
      #### Different mtry for rsf
      modelT <- rfsrc(formula = Formula, data = DATA, ntree = ntree,
                      mtry=mtrypool[jj],
                      nodesize = max(ceiling(sqrt(nrow(DATA))),15), 
                      splitrule = "custom1", 
                      membership = TRUE, forest = TRUE, samptype = "swor")
      predT = predictRSF(object = modelT, data = DATA, tpnt = Tpnt)
      RES$L2mtry$rsf[jj] = L2_tfixed(KM = predT, Data = DATA, Info = Info, T.pnt = Tpnt)
      rm(predT)
      predOOB = predictRSF(object = modelT,data = DATA, tpnt = Tpnt, OOB = TRUE)
      rm(modelT)
      RES$OOBmtry$rsf[jj] = unname(ipred::sbrier(Test.obj, predOOB[-1,])[1])
      print(sprintf("RSF -- mtry = %1.0f is done",mtrypool[jj]))
      rm(predOOB)
      gc()
    }
  }
  if (RES$L2$rsfT[1] == 0){
    idxmin <- which.min(RES$OOBmtry$rsf)
    RES$L2$rsfT <- RES$L2mtry$rsf[idxmin]
    RES$mtryall$rsf <- mtrypool[idxmin]
    RES$L2mtry$rsf[7] <- min(RES$OOBmtry$rsf[1:6])
  }
  if (RES$L2$rsfD[1] == 0){
    ## Training
    modelT <- rfsrc(formula = Formula, data = DATA, ntree = ntree,
                    mtry = mtryD,
                    nodesize = 15, splitrule = "custom1", 
                    membership = TRUE, forest = TRUE, samptype = "swor")
    
    predT = predictRSF(object = modelT,data = DATA, tpnt = Tpnt)
    rm(modelT)
    RES$L2$rsfD = L2_tfixed(KM = predT, Data = DATA, Info = Info, T.pnt = Tpnt)
    rm(predT)
    gc()
  }
  print("RSF is done ...")
  
  ###################### ---------------- L2 CF ----------------- ###############################
  leftzero = which(RES$L2mtry$cf[1:6] == 0)
  if (length(leftzero) > 0){
    for (jj in leftzero){
      #### Different mtry for CF
      modelT = cforest(formula = Formula, data = DATA,ntree=ntree, mtry=mtrypool[jj],
                       control = partykit::ctree_control(teststat = "quad", testtype = "Univ",
                                                         minsplit = max(ceiling(sqrt(nrow(DATA))),20),
                                                         minbucket = max(ceiling(sqrt(nrow(DATA))),7),
                                                         minprob = 0.01,
                                                         mincriterion = 0, saveinfo = FALSE))
      predT <- predict(object = modelT, type="prob")
      RES$L2mtry$cf[jj] = L2_tfixed(KM = predT, Data = DATA, Info = Info, T.pnt = Tpnt)
      rm(predT)
      predOOB <- predict(object = modelT, OOB = TRUE, type = "prob")
      RES$OOBmtry$cf[jj] = unname(ipred::sbrier(Test.obj, predOOB)[1])
      print(sprintf("CF -- mtry = %1.0f is done", mtrypool[jj]))
      rm(predOOB)
      rm(modelT)
      gc()
    }
  }
  if (RES$L2$cfT[1] == 0){
    idxmin = which.min(RES$OOBmtry$cf)
    RES$L2$cfT = RES$L2mtry$cf[c(idxmin,idxmin+6)]
    RES$mtryall$cf = mtrypool[idxmin]
    RES$L2mtry$cf[7] <- min(RES$OOBmtry$cf[1:6])
  }
  if (RES$L2$cfD[1] == 0){
    #### Training
    modelT <- cforest(formula = Formula, data = DATA, mtry = mtryD)
    #### Prediction
    predT <- predict(modelT, type = "prob")
    rm(modelT)
    gc()
    RES$L2$cfD = L2_tfixed(KM=predT, Data=DATA, Info=Info, T.pnt=Tpnt)
    rm(predT)
    gc()
  }
  print("CF is done ...")
  
  ###################### ---------------- L2 Cox ----------------- ###############################
  if (RES$L2$cx[1]==0){
    #### Training
    Coxfit <- coxph(formula = Formula, DATA)
    #### Prediction
    predT <- survfit(Coxfit, newdata = DATA)
    rm(Coxfit)
    RES$L2$cx = L2_tfixed(KM = predT, Data = DATA, Info = Info, T.pnt = Tpnt)
    rm(predT)
    gc()
  }
  print("Cox is done ...")
  
  ###################### ---------------- IBSCV ----------------- ###############################
  if (RES$ibsCVerr[Nfold, 4]==0){
    id_uniq = unique(DATA$ID)
    set.seed(sampleID[ll])
    kfold = kfoldcv(k = Nfold, N = N)
    
    dataCV = rep(list(0), Nfold)
    id_uniq_spd = sample(id_uniq, N)
    for (kk in 1:Nfold){
      dataCV[[kk]] = id_uniq_spd[((kk - 1) * kfold[kk] + 1):(kk * kfold[kk])]
    }
    
    foldleft = which(RES$ibsCVerr[, 4] == 0)
    for (jj in foldleft){
      idxCV = which(DATA$ID %in% dataCV[[jj]])
      data_b = DATA[-idxCV, ]
      newdata_b = DATA[idxCV, ]
      Test.obj <- Surv(newdata_b$Stop, newdata_b$Event)
      ## Cox
      if (RES$ibsCVerr[jj, 1]==0){
        Coxfit <- coxph(formula = Formula, data_b)
        #### Prediction
        predCox0 <- survfit(Coxfit, newdata = newdata_b)
        predT = rep(list(0),nrow(newdata_b))
        for( i in 1:nrow(newdata_b) ){
          predT[[i]] <- predCox0[i]
        }
        RES$ibsCVerr[jj, 1] = unname(ipred::sbrier(Test.obj, predT)[1])
        rm(Coxfit)
        rm(predT)
      }
      ## cfT
      if (RES$ibsCVerr[jj, 2]==0){
        mtryT = RES$mtryall$cf
        modelT = cforest(formula = Formula, data = data_b, ntree=ntree, mtry=mtryT,
                         control = partykit::ctree_control(teststat = "quad", testtype = "Univ",
                                                           minsplit = max(ceiling(sqrt(nrow(data_b))),20),
                                                           minbucket = max(ceiling(sqrt(nrow(data_b))),7),
                                                           minprob = 0.01,
                                                           mincriterion = 0, saveinfo = FALSE))
        predT <- predict(object=modelT,newdata =newdata_b,  type="prob")
        RES$ibsCVerr[jj, 2] = unname(ipred::sbrier(Test.obj, predT)[1])
        rm(modelT)
        gc()
      }
      ## rsf
      if (RES$ibsCVerr[jj,3]==0){
        mtryT = RES$mtryall$rsf
        modelT <- rfsrc(formula = Formula, data = data_b, ntree = ntree,
                        mtry=mtryT,
                        nodesize = max(ceiling(sqrt(nrow(data_b))),15), 
                        splitrule = "custom1", 
                        membership = TRUE, forest = TRUE, samptype = "swor")
        predT = predictRSF(object=modelT, data=data_b, newdata=newdata_b, tpnt = sort(unique(newdata_b$Stop)))
        rm(modelT)
        RES$ibsCVerr[jj, 3] = unname(ipred::sbrier(Test.obj, predT)[1])
        rm(modelT)
        gc()
      }
      ## tsf
      if (RES$ibsCVerr[jj,4]==0){
        mtryT = RES$mtryall$tsf
        rex <- Coxph(formula = Formula0, data = data_b, log_first = TRUE)
        modelT <- traforest(rex, formula = Formula, data = data_b, ntree = ntree,
                            mtry = mtryT,
                            control = partykit::ctree_control(teststat = "quad", testtype = "Univ",
                                                              minsplit = max(ceiling(sqrt(nrow(data_b))),20),
                                                              minbucket = max(ceiling(sqrt(nrow(data_b))),7),
                                                              minprob = 0.01,
                                                              mincriterion = 0, saveinfo = FALSE))
        predT <- predictTSF(modelT, data=data_b, newdata=newdata_b)
        rm(modelT)
        RES$ibsCVerr[jj, 4] <- unname(ipred::sbrier(Test.obj, predT)[1])
        rm(modelT)
        gc()
      }
      print(sprintf("IBSCV Round %1.0f is done",jj))
    }
    
  }
  return(RES)
}
#######################################################################################
nndata = c(50,100,300,500)

model = 2 # 2-linear; 3-nonlinear; 4-interaction
setting = "nonPH" # or "PH"
Distribution = "WI" # PH: "Exp", "WD", "WI", "Gtz"
nn = 1
Nfold = 10 # 10-fold IBS-based CV

## ll-th iteraction
ll = 1
LL2 <- Pred_CV(N = nndata[nn], Distribution = Distribution, model = model, 
               ll = ll, Nfold = Nfold, setting = setting)
print(LL2)

