#####################################################################################################
Range_T_linr <- function(TALL, DIST, X, U, variation, snrhigh){  # TALL = c(0,TS)
  u = U
  x1 = X$X1
  x2 = X$X2
  x3 = X$X3
  x4 = X$X4
  x5 = X$X5
  x6 = X$X6
  if (variation) {
    x3 = 0.5
    x4 = 0.5
    x6 = 1
  }
  
  if (snrhigh){
    Beta <- c(5,5,-5,-5,2.5,-2.5,-4)
  } else {
    Beta <- c(1,1,-1,-1,0.5,-0.5,0)
  }
  tlen = length(TALL)
  
  if(DIST == "Exp"){
    Lambda = 0.003
    Alpha = 0
    V = 0
    
    R0 = exp(Beta[1]*x1 +Beta[2]*x2 + Beta[3]*x3 + Beta[4]*x4 + Beta[5]*x5 + Beta[6]*x6 + Beta[7])
    R = Lambda*R0[-tlen] * (TALL[-1] - TALL[-tlen])
    
    VEC = c(0,cumsum(R),Inf)
    
    R.ID <- findInterval(-log(u), VEC)
    TT = -log(u) - VEC[R.ID] 
    TT = TT/Lambda/R0[R.ID] + TALL[R.ID]
  }else if(DIST == "WD"){
    Lambda = 0.012
    V = 0.8
    Alpha = 0
    R0 = exp(Beta[1]*x1 +Beta[2]*x2 + Beta[3]*x3 + Beta[4]*x4 + Beta[5]*x5 + Beta[6]*x6+ Beta[7])
    R = Lambda*R0[-tlen]*(TALL[-1]^V - TALL[-tlen]^V)
    
    VEC = c(0,cumsum(R),Inf)
    
    R.ID <- findInterval(-log(u), VEC)
    
    TT = -log(u) - VEC[R.ID] 
    TT = (TT/Lambda/R0[R.ID] + TALL[R.ID]^V)^(1/V)
    
  }else if(DIST == "WI"){
    Lambda = 0.001
    V = 2
    Alpha = 0
    R0 = exp(Beta[1]*x1 +Beta[2]*x2 + Beta[3]*x3 + Beta[4]*x4 + Beta[5]*x5 + Beta[6]*x6+ Beta[7])
    R = Lambda*R0[-tlen] * (TALL[-1]^V - TALL[-tlen]^V)
    
    VEC = c(0,cumsum(R),Inf)
    
    R.ID <- findInterval(-log(u), VEC)
    
    TT = -log(u) - VEC[R.ID] 
    TT = (TT/Lambda/R0[R.ID] + TALL[R.ID]^V)^(1/V)
    
  }else if(DIST == "Gtz"){
    Alpha = 0.1
    Lambda = 0.008
    V = 0
    R0 = exp(Beta[1]*x1 +Beta[2]*x2 + Beta[3]*x3 + Beta[4]*x4 + Beta[5]*x5 + Beta[6]*x6+ Beta[7])
    R = Lambda*R0[-tlen]/Alpha*(exp(Alpha*TALL[-1])-exp(Alpha*TALL[-tlen]))
    
    VEC = c(0,cumsum(R),Inf)
    
    R.ID <- findInterval(-log(u), VEC)
    TT = Alpha*(-log(u)-VEC[R.ID])/Lambda/R0[R.ID]
    TT = log(TT + exp(Alpha*TALL[R.ID]))/Alpha
    
  }
  
  result = list(Time = TT, Row = R.ID, Lambda = Lambda, Beta = Beta, Alpha = Alpha, V = V, Xi = R0)
  return(result)
}


#####################======== Large number of pseudo-subjects ===============#####################
Timevarying_PH_linear_gnrt <- function(N = 200, Distribution = "WI", censor.rate = 1,
                                       partial = TRUE, variation = FALSE, snrhigh = FALSE){
  npseu = 11
  Data <- as.data.frame(matrix(NA,npseu*N,27))
  names(Data)<-c("I","ID","X1","X2","X3","X4","X5","X6","X7","X8","X9","X10", 
                 "X11","X12","X13","X14","X15","X16","X17","X18","X19","X20",
                 "Start","Stop","C","Event","Xi")
  Data$ID <- rep(1:N,each=npseu)
  ## time-invariant
  Data$X1 <- rep(sample(c(0,1),N,replace=TRUE),each=npseu)
  Data$X2 <- rep(runif(N,0,1),each=npseu)
  
  Data$X7 <- rep(runif(N,0,1),each=npseu)
  Data$X8 <- rep(runif(N,1,2),each=npseu)
  Data$X9 <- rep(sample(c(1:5),N,replace=TRUE),each=npseu)
  Data$X10 <- rep(runif(N,0,1),each=npseu)
  Data$X11 <- rep(sample(c(0,1),N,replace=TRUE),each=npseu)
  Data$X12 <- rep(sample(c(0,1,2),N,replace=TRUE),each=npseu) 
  ## time-varying
  Data$X3 <- rbinom(npseu*N,size = 1,prob = 0.5)
  Data$X4 <- runif(npseu*N,0,1)
  Data$X5 <- sample(c(1:5),npseu*N,replace = TRUE)
  Data$X6 <- 2
  Data$X14 <- sample(c(1:5),npseu*N,replace = TRUE)
  Data$X15 <- runif(npseu*N,0,1)
  Data$X17 <- runif(npseu*N,0,1)
  Data$X19 <- sample(c(0,1),npseu*N,replace=TRUE) 

  Count = 1
  tall <- rep(0, N)
  while(Count <= N){
    TS <- rep(0, npseu-1)
    if (Distribution == "Exp"){
      while(any(diff(sort(TS))<0.4)){
        TS <- sort(rtrunc(npseu-1, spec="beta", a=0.0001, b = 1, shape1 = 0.01, shape2 = 2))*900
      }
    } else if (Distribution == "WD"){
      while(any(diff(sort(TS))<0.4)){
        TS <- sort(rtrunc(npseu-1, spec="beta", a=0.0001, b = 1, shape1 = 0.001, shape2 = 2))*800
      }
    } else if (Distribution == "WI"){
      while(any(diff(sort(TS))<0.4)){
        TS <- sort(rtrunc(npseu-1, spec="beta", a=0.001, b = 1, shape1 = 0.05, shape2 = 5))*120
      }
    } else if (Distribution == "Gtz"){
      while(any(diff(sort(TS))<0.4)){
        TS <- sort(rtrunc(npseu-1, spec="beta", a=0.006, b = 5, shape1 = 0.05, shape2 = 5))*80
      }
    } else {
      stop("Wrong dististribution is given.")
    }   
    
    u = runif(1)
    k = runif(2)
    
    
    x61 <- sample(c(0,1,2),1)
    x62 <- sample(npseu,1)-1
    Data[Data$ID==Count,]$X6[1:(1+x62)] <- rep(x61*(x61!=2) + 2*(x61==2), x62+1)
    if (x61 == 0 && x62 != (npseu - 1)){
      x63 <- sample(npseu-(1+x62),1)
      Data[Data$ID==Count,]$X6[(1+x62+1):(1+x62+x63)] <- rep(1,x63)
    }
    
    x13r <- sample(npseu-1,1)
    Data[Data$ID==Count,]$X13 <- c(numeric(x13r), rep(1,npseu-x13r))
  
    x16r <- sample(npseu-1,1)
    x161 <- sample(c(0,1),1)
    Data[Data$ID==Count,]$X16 <- c(rep(x161, x16r), rep(x161+1, npseu-x16r))
    
    x18r <- sort(sample(npseu-1,2))
    Data[Data$ID==Count,]$X18 <- c(rep(0, x18r[1]), rep(1,x18r[2]-x18r[1]), rep(2,npseu-x18r[2]))
    
    Data[Data$ID==Count,]$X20 <- k[1] * c(0,TS) + k[2]
    
    Data[Data$ID==Count,]$Start <- c(0,TS)
    Data[Data$ID==Count,]$Stop <- c(TS,NA)
    
    RT <- Range_T_linr(TALL = c(0, TS), 
                       DIST = Distribution, 
                       X = Data[Data$ID == Count, ], 
                       U = u, 
                       variation = variation, 
                       snrhigh = snrhigh)
    t = RT$Time
    rID = RT$Row
    Data[Data$ID==Count,]$Xi = RT$Xi
    tall[Count] = t
    if(rID==1){
      Data[Data$ID==Count,][1,]$Stop = t
      Data[Data$ID==Count,][1,]$Event = 1
    }else{
      Data[Data$ID==Count,][1:(rID-1),]$Event=0
      Data[Data$ID==Count,][rID,]$Event = 1
      Data[Data$ID==Count,][rID,]$Stop = t
    }
    #print(Count)
    Count = Count + 1
    
  }## end of the while loop
  DATA <- Data[!is.na(Data$Event),]
  rm(Data)
  gc()
  ###================== Add Censoring =========================================
  DATA$C <- 0
  
  if(length(unique(DATA$ID)) != N){
    stop("ID length NOT equal to N")
  }
  RET <- NULL
  
  if (variation){ 
    if (snrhigh) {
      if(Distribution == "WI"){
        if (censor.rate == 0){
          Censor.time <- rep(Inf, N)
        } else if (censor.rate == 1){
          Censor.time = rexp(N, rate = 1/63) ## NEW
        } else if (censor.rate == 2) { 
          Censor.time = rexp(N,rate = 1/7.5) ## NEW
        } else {
          stop("Wrong censoring type")
        } 
      } else {
        stop("Wrong distribution")
      }
    } else {
      if(Distribution == "WI"){
        if (censor.rate == 0){
          Censor.time <- rep(Inf, N)
        } else if (censor.rate == 1){
          Censor.time = rexp(N, rate = 1/78)
        } else if (censor.rate == 2) {
          Censor.time = rexp(N,rate = 1/20) ## NEW
        } else {
          stop("Wrong censoring type")
        } 
      } else {
        stop("Wrong distribution")
      }
    }
    
  } else { # Basic scenario
    if (snrhigh) {
      if(Distribution == "WI"){
        if (censor.rate == 0) {
          Censor.time <- rep(Inf,N)
        } else if (censor.rate == 1){
          Censor.time = rexp(N,rate = 1/78) 
        } else if (censor.rate == 2){
          Censor.time = rexp(N,rate = 1/6) # NEW
        } else{
          stop("Wrong censoring type")
        }
      } else {
        stop("Wrong distribution")
      }
    } else { # SNR original
      if (censor.rate == 0){
        Censor.time <- rep(Inf,N)
      }else if(censor.rate == 1){
        if(Distribution == "Exp"){
          Censor.time = rexp(N,rate = 1/805)
        }else if(Distribution == "WD"){
          Censor.time = rexp(N,rate = 1/440)
        }else if(Distribution == "WI"){
          Censor.time = rexp(N,rate = 1/82)
        }else if(Distribution == "Gtz"){
          Censor.time = rexp(N,rate = 1/75)
        }else {
          stop("Wrong distribution")
        }
      }else if(censor.rate == 2){ 
        if(Distribution == "WI"){
          Censor.time = rexp(N,rate = 1/22)
        } else {
          stop("Wrong distribution")
        }
      }else{
        stop("Wrong censoring type")
      }
    }
  }
  
  
  for( j in 1:length(unique(DATA$ID)) ){
    Vec <- c(0,DATA[DATA$ID==j,]$Stop,Inf)
    ID <- findInterval(Censor.time[j], Vec)
    
    if( ID <= nrow(DATA[DATA$ID==j,]) ){
      DATA[DATA$ID==j,][ID,]$C = 1
      DATA[DATA$ID==j,][ID,]$Event = 0
      DATA[DATA$ID==j,][ID,]$Stop = Censor.time[j]
      if( ID != nrow(DATA[DATA$ID==j,]) ){
        DATA[DATA$ID==j,][(ID+1):nrow(DATA[DATA$ID==j,]),]$Event = NA
      }
    }
  }
  
  Data <- DATA[!is.na(DATA$Event),]
  rm(DATA)
  gc()
  if(length(unique(Data$ID))!=N){
    stop("ID length NOT equal to N")
  }
  Data$I <- 1:nrow(Data)
  RET$fullData <- Data
  RET$fullData$Start <- round(RET$fullData$Start, 3)
  RET$fullData$Stop <- round(RET$fullData$Stop, 3)
  RET$fullData$Stop[RET$fullData$Start == RET$fullData$Stop] = RET$fullData$Stop[RET$fullData$Start == RET$fullData$Stop] + 0.001
  RET$Info = list(Coeff=list(Lambda = RT$Lambda, Alpha = RT$Alpha, Beta = RT$Beta, V = RT$V),
                  Dist=Distribution,
                  Set = "PH")

  DATA = data.frame(matrix(0,nrow=N,ncol=ncol(Data)))
  names(DATA) = names(Data)
  for (ii in 1:N){
    DATA[ii,] = Data[Data$ID==ii,][1,]
    ni = nrow(Data[Data$ID==ii,])
    if (ni >1){
      DATA[ii,]$Stop = Data[Data$ID==ii,]$Stop[ni]
      DATA[ii,]$Event = Data[Data$ID==ii,]$Event[ni]
    }
  }
  DATA$I = 1:nrow(DATA)
  RET$baselineData = DATA
  rm(DATA)
  
  if (partial){
    IDrev = rev(unique(Data$ID))
    # partialInfo = data.frame(matrix(0,nrow = N,ncol = 7))
    # names(partialInfo) = c("n","n_unobv","percT_unobv","wXiabsDiff","n_mto","percT_unobv_mto","wXiabsDiff_mto")
    for (Count in IDrev) {
      pN = sum(Data$ID == Count)
      I_count = Data[Data$ID==Count,]$I
      ## at least 50% of the times are not observable
      n_unobv = floor(pN*0.5)
      # partialInfo[Count,1]=pN # number of pseudo-subject in the full dataset
      # partialInfo[Count,2]=n_unobv # number of unobserved in the partial dataset
      if (pN > 1){
        unobv = sort(sample(1:(pN-1),n_unobv))+1
        tLast = Data[Data$ID==Count,]$Stop[pN]
        EventLast = Data[Data$ID==Count,]$Event[pN]
        
        T_unobv = Data[I_count[unobv],]$Stop
        dt_unobv = Data[I_count[unobv],]$Stop - Data[I_count[unobv],]$Start
        Xi_unobv = Data[I_count[unobv],]$Xi
        # partialInfo[Count,3]= sum(dt_unobv)/tLast
        
        Data = Data[-I_count[unobv],]
        if (pN > 2){
          Data[Data$ID==Count,]$Stop[1:(pN-n_unobv-1)] = Data[Data$ID==Count,]$Start[2:(pN-n_unobv)]
        } 
        Data[Data$ID==Count,]$Stop[pN-n_unobv] = tLast
        Data[Data$ID==Count,]$Event[pN-n_unobv] = EventLast
        
        # Xi_new_rID = findInterval(T_unobv, c(0,Data[Data$ID==Count,]$Stop),
        #                           left.open = TRUE, rightmost.closed = TRUE)
        # partialInfo[Count,4]=sum(abs(Data[Data$ID==Count,]$Xi[Xi_new_rID]-Xi_unobv)/Xi_unobv *dt_unobv)/sum(dt_unobv)
      }    
    }
    Data$I = 1:nrow(Data)
    RET$partialData = Data
    RET$partialData$Start <- round(RET$partialData$Start, 3)
    RET$partialData$Stop <- round(RET$partialData$Stop, 3)
    # RET$partialInfo = round(colMeans(partialInfo),digits = 3)
    # RET$partialInfo[5] = sum(partialInfo$n_unobv==0)
    # RET$partialInfo[6] = round(mean(partialInfo$percT_unobv[partialInfo$percT_unobv!=0]),digits = 3)
    # RET$partialInfo[7] = round(mean(partialInfo$wXiabsDiff[partialInfo$wXiabsDiff!=0]),digits = 3)
  }
  
  rm(Data)
  gc()
  
  return(RET)
  # return(tall)
}
