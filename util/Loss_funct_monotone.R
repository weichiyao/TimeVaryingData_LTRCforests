############################################################################################################
# Created -- May 13th
############################################################################################################
shat_funct <- function(Ni, data, pred=NULL, tpnt, obj.roc = NULL){
  ## This function is to compute the estimated survival probability of the Ni-th subject
  id_uniq <- unique(data$ID)
  
  ## the i-th data 
  TestData = data[data$ID == id_uniq[Ni],]
  
  TestT = c(0,TestData$Stop)
  TestTIntN = nrow(TestData)
  
  ## tpnt can either be the original tpnt or can be the shorten one
  tlen = length(tpnt)
  if (is.null(obj.roc)) {
    if ("survfit.cox" %in% class(pred) || "survfitcox" %in% class(pred)){
      r.ID = findInterval(tpnt, TestT)
      r.ID[tpnt >= TestData$Stop[TestTIntN]] = TestTIntN
      jall <- unique(r.ID)
      nj <- length(jall)
      
      Shat_temp <- matrix(0,nrow = 1,ncol = tlen)
      ## Compute the estimated survival probability of the Ni-th subject
      if(nj == 1){
        ## Get the index of the Pred to compute Shat
        II = which(data$ID == id_uniq[Ni])[1] + jall[1] - 1
        Shat_temp[1,r.ID == jall[1]] <- getSurv(pred[II], tpnt[r.ID == jall[1]])
      } else {
        ShatR_temp <- matrix(1,nrow = 1,ncol = nj+1)
        ## Get the position of L_1, ..., L_n. Note that, L_2=R_1, ..., L_{j+1} = R_{j}, ..., L_n = R_{n-1}
        r.IDmax = c(1,sapply(jall[-nj], function(j){
          max(which(r.ID==j))+1
        }))
        for (j in 1:nj){
          ## Get the index of the Pred to compute Shat
          II = which(data$ID == id_uniq[Ni])[1] + jall[j] - 1
          KM = pred[II]
          Shat_temp[1,r.ID == jall[j]] <- getSurv(KM, tpnt[r.ID == jall[j]])
          # S_{j}(R_{j}), j=1,...nj-1
          jR = getSurv(KM, TestT[j+1]) 
          ShatR_temp[1,j+1] = jR
        }
        ShatR_temp=ShatR_temp[1,1:nj]
        
        # S_1(L_1), S_2(L_2), S_3(L_3), ..., S_{nj}(L_{nj})
        qL = Shat_temp[1,r.IDmax]
        
        ql0=which(qL==0)
        if(length(ql0)>0){
          maxqlnot0=max(which(qL>0))
          for(j in ql0){
            if (j<maxqlnot0) {
              Shat_temp[1,r.ID == jall[j]] = Shat_temp[1,r.ID == jall[j]]+.000000001
            }
          }
        }
        # Orignial non-monotone one:
        # Shat <- list(survprob = Shat_temp[1,])
        # S_1(L_1), S_2(L_2), S_3(L_3), ..., S_{nj}(L_{nj})
        qL = Shat_temp[1,r.IDmax]
        # S_1(L_1), S_1(R_1), S_2(R_2), ..., S_{nj-1}(R_{nj-1})
        qR = ShatR_temp
        m = cumprod(qR/qL)
        for (j in 1:nj){
          Shat_temp[1,r.ID==jall[j]] = Shat_temp[1,r.ID==jall[j]]*m[j]
        }
      }
      Shat = Shat_temp[1,]
    } else {
      r.ID = findInterval(tpnt, TestT)
      r.ID[tpnt >= TestData$Stop[TestTIntN]] = TestTIntN
      jall <- unique(r.ID)
      nj <- length(jall)
      
      Shat_temp <- matrix(0,nrow = 1,ncol = tlen)
      ## Compute the estimated survival probability of the Ni-th subject
      if(nj == 1){
        ## Get the index of the Pred to compute Shat
        II = which(data$ID == id_uniq[Ni])[1] + jall[1] - 1
        Shat_temp[1,r.ID == jall[1]] <- getSurv(pred[[II]], tpnt[r.ID == jall[1]])
      } else {
        ShatR_temp <- matrix(1,nrow = 1,ncol = nj+1)
        ## Get the position of L_1, ..., L_n. Note that, L_2=R_1, ..., L_{j+1} = R_{j}, ..., L_n = R_{n-1}
        r.IDmax = c(1,sapply(jall[-nj], function(j){
          max(which(r.ID==j))+1
        }))
        for (j in 1:nj){
          ## Get the index of the Pred to compute Shat
          II = which(data$ID == id_uniq[Ni])[1] + jall[j] - 1
          KM = pred[[II]]
          Shat_temp[1,r.ID == jall[j]] <- getSurv(KM, tpnt[r.ID == jall[j]])
          # S_{j}(R_{j}), j=1,...nj-1
          jR = getSurv(KM, TestT[j+1]) 
          ShatR_temp[1,j+1] = jR
        }
        ShatR_temp=ShatR_temp[1,1:nj]
        
        # S_1(L_1), S_2(L_2), S_3(L_3), ..., S_{nj}(L_{nj})
        qL = Shat_temp[1,r.IDmax]
        
        # ql0=which(qL==0)
        # if(length(ql0)>0){
        #   maxqlnot0=max(which(qL>0))
        #   for(j in ql0){
        #     if (j<maxqlnot0) {
        #       Shat_temp[1,r.ID == jall[j]] = Shat_temp[1,r.ID == jall[j]]+.000000001
        #     }
        #   }
        # }
        # Orignial non-monotone one:
        Shat <- Shat_temp[1,]
        # # S_1(L_1), S_2(L_2), S_3(L_3), ..., S_{nj}(L_{nj})
        # qL = Shat_temp[1,r.IDmax]
        # # S_1(L_1), S_1(R_1), S_2(R_2), ..., S_{nj-1}(R_{nj-1})
        # qR = ShatR_temp
        # 
        # m = cumprod(qR/qL)
        # for (j in 1:nj){
        #   Shat_temp[1,r.ID==jall[j]] = Shat_temp[1,r.ID==jall[j]]*m[j]
        # }
      }
      Shat = Shat_temp[1,]
    }
    
  } else {
    pred = predict(obj.roc, TestData)$pred
    Shat = getSurv(pred,tpnt)
  }
  
  Shat
}
############################################################################################################
shat <- function(data, pred=NULL, tpnt=NULL, obj.roc = NULL){
  if (is.null(tpnt)){
    tpnt = c(0,sort(unique(data$Stop)))
  }
  N = length(unique(data$ID))
  Shatt = sapply(1:N, function(Ni) shat_funct(Ni,data,pred,tpnt,obj.roc))
  return(Shatt)
}
############################################################################################################
Diff_funct <- function(Ni, data, fulldata=NULL, data_sbrier = NULL, info, pred=NULL, tpnt, obj.roc = NULL, IBS = TRUE){
  if (!is.null(obj.roc)){
    IBS = FALSE
  }
  id_uniq <- unique(data$ID)
  
  ## Only up to the last observed time
  maxT = max(data[data$ID==id_uniq[Ni],]$Stop)
  TTpnt = tpnt[tpnt<=maxT]
  
  ## Compute the estimated survival probabilities
  if (class(pred)[1] %in% c("rfsrc","matrix")){
    ShatA = pred[,Ni]
  } else {
    ShatA = shat_funct(Ni, data, pred, tpnt=tpnt, obj.roc)
  }
  
  ## Compute the true survival probability of the Ni-th subject
  St <- s_funct(Ni, fulldata, info, tpnt = TTpnt)
  
  f_itg = (ShatA[1:length(TTpnt)] - St)^2
  L2 = diff(TTpnt)%*%(f_itg[-length(f_itg)] + f_itg[-1])/2
  L2 = L2/diff(range(TTpnt)) 
  
  ibs = rep(0,3)
  if (IBS){
    ######================ reverse Kaplan-Meier: estimate censoring distribution 
    # deal with ties 
    hatcdist <- prodlim(Surv(times, cens) ~ 1, data = data_sbrier, reverse = TRUE)
    
    Ttildei <- data_sbrier[data_sbrier$ID == id_uniq[Ni],]$times
    if (Ttildei != maxT){
      stop("Check! data_sbrier gives different last observed time as data.")
    }
    ### conditional survival for Observed value < t, G(Obs)
    csurv_obs <- predict(hatcdist, times = Ttildei, type = "surv")
    csurv_obs[csurv_obs == 0] <- Inf
    
    # conditional survival for Observed value > t, G(t)
    csurv_t <- predict(hatcdist, times = tpnt[tpnt < Ttildei], type = "surv")
    csurv_t[is.na(csurv_t)] <- min(csurv_t, na.rm = TRUE)
    csurv_t[csurv_t == 0] <- Inf
    
    ## c(G^{-1}(t), G^{-1}(Obs))
    csurv <- c(1/csurv_t, rep(1/csurv_obs,sum(tpnt >= Ttildei)))
    
    ######================ indicator ================#################
    Indicator_t <- as.integer(tpnt<Ttildei)
    Indicator_t[Indicator_t==0] = as.integer(data_sbrier[data_sbrier$ID == id_uniq[Ni],]$cens == 1)
    
    ######================ Brier score =================#################
    ## up to last time point of tpnt
    fibs_itg = (as.integer(tpnt<Ttildei)-ShatA)^2*csurv*Indicator_t
    ibs[1] = diff(tpnt) %*% (fibs_itg[-length(fibs_itg)] + fibs_itg[-1]) / 2
    ibs[1] = ibs[1]/diff(range(tpnt))
    
    ttpnt = tpnt[tpnt<= 1.5*maxT]
    fibs_itg0 = fibs_itg[1:length(ttpnt)]
    ibs[2] = diff(ttpnt) %*% (fibs_itg0[-length(fibs_itg0)] + fibs_itg0[-1]) / 2
    ibs[2] = ibs[2]/diff(range(ttpnt))
    
    ttpnt = tpnt[tpnt<= max(data$Stop)]
    fibs_itg = fibs_itg[1:length(ttpnt)]
    ibs[3] = diff(ttpnt) %*% (fibs_itg[-length(fibs_itg)] + fibs_itg[-1]) / 2
    ibs[3] = ibs[3]/diff(range(ttpnt))
  }
  ret = c(L2,ibs)
  ret
}
###########################################################################################
Diff <- function(data, fulldata=NULL, info, pred=NULL, tpnt, obj.roc = NULL, IBS = TRUE){
  if (is.null(fulldata)){
    fulldata = data
  }
  id_uniq = unique(data$ID)
  N = length(id_uniq)
  
  data_sbrier <- data.frame(matrix(0,nrow = N, ncol = 3))
  names(data_sbrier) <- c("ID", "times", "cens")
  data_sbrier$ID <- id_uniq
  for (ii in 1:N){
    data_sbrier[ii,]$times = max(data[data$ID==id_uniq[ii],]$Stop)
    data_sbrier[ii,]$cens = sum(data[data$ID==id_uniq[ii],]$Event)
  }
  
  ret = sapply(1:N, function(Ni) Diff_funct(Ni, data, fulldata, data_sbrier, info, pred, tpnt, obj.roc, IBS))
  return(rowMeans(ret))
}
############################################################################################################
bs_funct <- function(Ni, data, data_sbrier, pred, tpnt, tl=FALSE){
  
  id_uniq <- unique(data$ID)
  ## Only up to the last observed time
  # maxT = max(data[data$ID==id_uniq[Ni],]$Stop)
  # TTpnt = tpnt[tpnt<=maxT]
  
  ## Compute the estimated survival probabilities
  if (class(pred)[1] %in% c("rfsrc","matrix")){
    ShatA = pred[1:length(tpnt),Ni] ## Therefore, do not need to adjust the length of input pred 
  } else {
    ShatA = shat_funct(Ni, data, pred, tpnt=tpnt)
  }
  
  ######================ reverse Kaplan-Meier: estimate censoring distribution 
  # deal with ties 
  hatcdist <- prodlim(Surv(times, cens) ~ 1, data = data_sbrier, reverse = TRUE)
  
  Ttildei <- data_sbrier[data_sbrier$ID == id_uniq[Ni],]$times
  ### conditional survival for Observed value < t, G(Obs)
  csurv_obs <- predict(hatcdist, times = Ttildei, type = "surv")
  csurv_obs[csurv_obs == 0] <- Inf
  
  # conditional survival for Observed value > t, G(t)
  csurv_t <- predict(hatcdist, times = tpnt[tpnt < Ttildei], type = "surv")
  csurv_t[is.na(csurv_t)] <- min(csurv_t, na.rm = TRUE)
  csurv_t[csurv_t == 0] <- Inf
  
  ## c(G^{-1}(t), G^{-1}(Obs))
  csurv <- c(1/csurv_t, rep(1/csurv_obs,sum(tpnt >= Ttildei)))
  
  ######================ indicator ================#################
  Indicator_t <- as.integer(tpnt<Ttildei)
  Indicator_t[Indicator_t==0] = as.integer(data_sbrier[data_sbrier$ID == id_uniq[Ni],]$cens == 1)
  
  ######================ Brier score =================#################
  fibs_itg = (as.integer(tpnt<Ttildei)-ShatA)^2*csurv*Indicator_t
  if (!tl){
    ibs = diff(tpnt) %*% (fibs_itg[-length(fibs_itg)] + fibs_itg[-1]) / 2
    ibs = ibs/diff(range(tpnt))
  } else {
    fibs_itg = fibs_itg[tpnt<= 1.5*Ttildei]
    ibs = diff(tpnt[tpnt<= 1.5*Ttildei]) %*% (fibs_itg[-length(fibs_itg)] + fibs_itg[-1]) / 2
    ibs = ibs/diff(range(tpnt[tpnt<= 1.5*Ttildei]))
  }
  ibs
}
############################################################################################################
bs <- function(data, pred, tpnt=NULL, tl=FALSE){
  if (is.null(tpnt)){
    tpnt = c(0,sort(unique(data$Stop)))
  }
  
  id_uniq = unique(data$ID)
  N = length(id_uniq)
  
  data_sbrier <- data.frame(matrix(0,nrow = N, ncol = 3))
  names(data_sbrier) <- c("ID", "times", "cens")
  data_sbrier$ID <- id_uniq
  for (ii in 1:N){
    data_sbrier[ii,]$times = max(data[data$ID==id_uniq[ii],]$Stop)
    data_sbrier[ii,]$cens = sum(data[data$ID==id_uniq[ii],]$Event)
  }
  
  ret = sapply(1:N, function(Ni) bs_funct(Ni, data, data_sbrier, pred, tpnt, tl))
  return(mean(ret))
}
############################################################################################################
l2_funct <- function(Ni, data, fulldata, info, pred, tpnt, obj.roc = NULL){
  
  id_uniq <- unique(data$ID)
  ## Only up to the last observed time
  maxT = max(data[data$ID==id_uniq[Ni],]$Stop)
  TTpnt = tpnt[tpnt<=maxT]
  
  ## Compute the estimated survival probabilities
  if (class(pred)[1] %in% c("rfsrc","matrix")){
    ShatA = pred[1:length(TTpnt),Ni]
  } else {
    ShatA = shat_funct(Ni, data, pred, tpnt=TTpnt, obj.roc)
  }
  
  ## Compute the true survival probability of the Ni-th subject
  St <- s_funct(Ni, fulldata, info, tpnt = TTpnt)
  
  f_itg = (ShatA - St)^2
  L2 = diff(TTpnt)%*%(f_itg[-length(f_itg)] + f_itg[-1])/2
  L2 = L2/diff(range(TTpnt)) 
  
  L2
}
############################################################################################################
l2 <- function(data, fulldata=NULL, info, pred=NULL, tpnt, obj.roc = NULL){
  if (is.null(fulldata)){
    fulldata = data
  }
  id_uniq = unique(data$ID)
  N = length(id_uniq)
  
  ret = sapply(1:N, function(Ni) l2_funct(Ni, data, fulldata, info, pred, tpnt, obj.roc))
  return(mean(ret))
}
############################################################################################################
bs_funct_tl <- function(Ni, data, data_sbrier, pred, tpnt){
  
  id_uniq <- unique(data$ID)
  # Only up to 1.5*last observed time
  maxT = max(data[data$ID==id_uniq[Ni],]$Stop)
  tpnt = tpnt[tpnt<=1.5*maxT]
  
  ## Compute the estimated survival probabilities
  if (class(pred)[1] %in% c("rfsrc","matrix")){
    ShatA = pred[[Ni]] ## for random forest, the class has been forced to be matrix
    if (length(pred[[Ni]]) != length(tpnt)) stop("Something is wrong that the length of prediction does not match")
  } else {
    ShatA = shat_funct(Ni, data, pred, tpnt=tpnt)
  }
  
  ######================ reverse Kaplan-Meier: estimate censoring distribution 
  # deal with ties 
  hatcdist <- prodlim(Surv(times, cens) ~ 1, data = data_sbrier, reverse = TRUE)
  
  Ttildei <- data_sbrier[data_sbrier$ID == id_uniq[Ni],]$times
  ### conditional survival for Observed value < t, G(Obs)
  csurv_obs <- predict(hatcdist, times = Ttildei, type = "surv")
  csurv_obs[csurv_obs == 0] <- Inf
  
  # conditional survival for Observed value > t, G(t)
  csurv_t <- predict(hatcdist, times = tpnt[tpnt < Ttildei], type = "surv")
  csurv_t[is.na(csurv_t)] <- min(csurv_t, na.rm = TRUE)
  csurv_t[csurv_t == 0] <- Inf
  
  ## c(G^{-1}(t), G^{-1}(Obs))
  csurv <- c(1/csurv_t, rep(1/csurv_obs,sum(tpnt >= Ttildei)))
  
  ######================ indicator ================#################
  Indicator_t <- as.integer(tpnt<Ttildei)
  Indicator_t[Indicator_t==0] = as.integer(data_sbrier[data_sbrier$ID == id_uniq[Ni],]$cens == 1)
  
  ######================ Brier score =================#################
  fibs_itg = (as.integer(tpnt<Ttildei)-ShatA)^2*csurv*Indicator_t
  
  ibs = diff(tpnt) %*% (fibs_itg[-length(fibs_itg)] + fibs_itg[-1]) / 2
  ibs = ibs/diff(range(tpnt))
  
  ibs
}
############################################################################################################
bs_tl <- function(data, pred, tpnt=NULL){
  if (is.null(tpnt)){
    tpnt = c(0,sort(unique(data$Stop)), seq(max(data$Stop)+5,1.5*max(data$Stop),by = 5))
  }
  
  id_uniq = unique(data$ID)
  N = length(id_uniq)
  
  data_sbrier <- data.frame(matrix(0,nrow = N, ncol = 3))
  names(data_sbrier) <- c("ID", "times", "cens")
  data_sbrier$ID <- id_uniq
  for (ii in 1:N){
    data_sbrier[ii,]$times = max(data[data$ID==id_uniq[ii],]$Stop)
    data_sbrier[ii,]$cens = sum(data[data$ID==id_uniq[ii],]$Event)
  }
  
  ret = sapply(1:N, function(Ni) bs_funct_tl(Ni, data, data_sbrier, pred, tpnt))
  return(mean(ret))
}
############################################################################################################
bs_funct_t <- function(Ni, data, data_sbrier, pred, tpnt, tl = TRUE){
  id_uniq <- unique(data$ID)
  # Only up to 1.5*last observed time
  maxT = max(data[data$ID==id_uniq[Ni],]$Stop)
  tlen = length(tpnt)
  if (tl) tpnt = tpnt[tpnt<=1.5*maxT]
  ## Compute the estimated survival probabilities
  if (class(pred)[1] %in% c("rfsrc","matrix")){
    ShatA = pred[1:length(tpnt),Ni] ## Therefore, do not need to adjust the length of input pred 
    #if (length(pred[[Ni]]) != length(tpnt)) stop("Something is wrong that the length of prediction does not match")
  } else {
    ShatA = shat_funct(Ni, data, pred, tpnt=tpnt)
  }
  
  ######================ reverse Kaplan-Meier: estimate censoring distribution 
  # deal with ties 
  hatcdist <- prodlim(Surv(times, cens) ~ 1, data = data_sbrier, reverse = TRUE)
  
  Ttildei <- data_sbrier[data_sbrier$ID == id_uniq[Ni],]$times
  ### conditional survival for Observed value < t, G(Obs)
  csurv_obs <- predict(hatcdist, times = Ttildei, type = "surv")
  csurv_obs[csurv_obs == 0] <- Inf
  
  # conditional survival for Observed value > t, G(t)
  csurv_t <- predict(hatcdist, times = tpnt[tpnt < Ttildei], type = "surv")
  csurv_t[is.na(csurv_t)] <- min(csurv_t, na.rm = TRUE)
  csurv_t[csurv_t == 0] <- Inf
  
  ## c(G^{-1}(t), G^{-1}(Obs))
  csurv <- c(1/csurv_t, rep(1/csurv_obs,sum(tpnt >= Ttildei)))
  
  ######================ indicator ================#################
  Indicator_t <- as.integer(tpnt<Ttildei)
  Indicator_t[Indicator_t==0] = as.integer(data_sbrier[data_sbrier$ID == id_uniq[Ni],]$cens == 1)
  
  ######================ Brier score =================#################
  fibs_itg = (as.integer(tpnt<Ttildei)-ShatA)^2*csurv*Indicator_t
  
  if (tl) fibs_itg = c(fibs_itg,rep(NA,tlen-length(fibs_itg)))
  
  fibs_itg
}
############################################################################################################
bs_t <- function(data, pred, tl = TRUE){
  if (tl) {
    tpnt = c(0,sort(unique(data$Stop)), seq(max(data$Stop)+5,1.5*max(data$Stop),by = 5))
  } else {
    tpnt = c(0,sort(unique(data$Stop)))
  }
  
  id_uniq = unique(data$ID)
  N = length(id_uniq)
  
  data_sbrier <- data.frame(matrix(0,nrow = N, ncol = 3))
  names(data_sbrier) <- c("ID", "times", "cens")
  data_sbrier$ID <- id_uniq
  for (ii in 1:N){
    data_sbrier[ii,]$times = max(data[data$ID==id_uniq[ii],]$Stop)
    data_sbrier[ii,]$cens = sum(data[data$ID==id_uniq[ii],]$Event)
  }
  
  ret = sapply(1:N, function(Ni) bs_funct_t(Ni, data, data_sbrier, pred, tpnt, tl=tl))
  
  return(rowMeans(ret,na.rm = TRUE))
}
