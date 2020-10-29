#####################################################################################################
itct_term_WI <- function(beta,x1,x2,x3,x4,x5,x6){
  R1 = beta[9] * (x1 * x2 - log(x3 + x4) - x6 / x5) + beta[10]
  R2 = - beta[1] * x1 - beta[2] * x2 - beta[3] * x3 - beta[4] * x4 - beta[5] * x5 - beta[6] * x6 + beta[8]
  R3 = beta[11] * (cos(pi*(x1+x5))+sqrt(x6+x2)-x3) + beta[12]
  R4 = beta[1] * x1 + beta[2] * x2 + beta[3] * x3 + beta[4] * x4 + beta[5] * x5 + beta[6] * x6 + beta[7]
  R0 = (x4 >= 0.7) * ((x5 == 5)) * R1 + (x4 >= 0.7) * (1 - (x5 == 5)) * R3 +
    (x4 < 0.7) * ((x5 == 5)) * R2 + (x4 < 0.7) * (1 - (x5 == 5)) * R4
  return(R0)
}

itct_term_WI_variation <- function(beta,x1,x2,x3,x4,x5,x6){
  R1 = x1*x2 - 1/x5
  R2 = -beta[1]*x1-beta[2]*x2-beta[3]*0.5-beta[4]*0.5-beta[5]*x5-beta[6]
  R3 = cos(pi*(x1+x5))+sqrt(1+x2)-0.5
  R4 = beta[1]*x1+beta[2]*x2+beta[3]*0.5+beta[4]*0.5+beta[5]*x5+beta[6]
  R0 = (x2>=0.7)*((x5==5))*R1+(x2>=0.7)*(1-(x5==5))*R3+
    (x2<0.7)*((x5==5))*R2 + (x2<0.7)*(1-(x5==5))*R4
  return(R0)
}
#####################################################################################################
Comp_fstar_WI <- function(TALL, M, X, variation, snrhigh){  # TALL = c(0,TS)
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
    itct_term_WI <- itct_term_WI_variation
  }
  
  if (snrhigh){ 
    #         1 2 3 4  5 6  7       8 9  10 11      12 13       14
    Beta <- c(3,3,3,3,30,3,-64.5,64.5,5,0.6,5,1.101021, 5, -2.834679)
    # Beta <- c(5,5,5,5,10,5,-12,12,5,0.6,5,1.101021, 5, -2.834679)
  } else {
    Beta <- c(1,1,1,1,10,1,  0, 0,1,  0,1,       0, 1,         0)
  }
  tlen = length(TALL)
  if (M == "linear"){
    Fstar <- Beta[1] * x1 + Beta[2] * x2 + Beta[3] * x3 + Beta[4] * x4 + Beta[5] * x5 + Beta[6] * x6 + Beta[7]
  } else if (M == "nonlinear"){
    Fstar <- Beta[13] * cos(x1 + x3 + x5 + x6 + x4 + x2) + Beta[14]
  } else if (M == "interaction"){
    Fstar <- sapply(1:tlen, function(i) itct_term_WI(beta = Beta,
                                                     x1[i], x2[i], x3[i],
                                                     x4[i], x5[i], x6[i]))
  } else {
    stop("Wrong model type is given.")
  }
  
  return(Fstar)
}

Comp_xi_WI <- function(TALL, M, Fstar, U){  # TALL = c(0,TS)
  tlen = length(TALL)
  if (M == "linear"){
    Lambda = 0.001
  } else if (M == "nonlinear"){
    Lambda = 0.002
  } else if (M == "interaction"){
    Lambda = 0.001
  } else {
    stop("Wrong model is set.")
  }
  
  xi <- exp(Fstar)
  R = Lambda^xi[-tlen]*(TALL[-1]^xi[-tlen]-TALL[-tlen]^xi[-tlen])
  VEC = c(0,cumsum(R),Inf) ## H = cumsum(R) from Ht(1) to Ht(M-1)
  
  R.ID <- findInterval(-log(U), VEC)
  TT = -log(U) - VEC[R.ID] 
  TT = (TT/(Lambda^xi[R.ID]) + TALL[R.ID]^xi[R.ID])^(1/xi[R.ID])
  
  result = list(Time = TT, Row = R.ID, Coeff = Lambda, Xi = xi)
  return(result)
}

#####################======== Large number of pseudo-subjects ===============#####################
Timevarying_nonPH_gnrt <- function(N = 200, model = c("linear","nonlinear","interaction"), 
                                   Distribution = "WI",
                                   censor.rate = 1, partial = TRUE, variation = FALSE, snrhigh = FALSE){
  npseu = 11
  if (model == "linear"){
    Nadd = 100
  } else if (model == "interaction"){
    Nadd = 100
  } else if (model == "nonlinear"){
    if (N == 500){
      Nadd = 1200
    } else {
      Nadd = 500
    }
  } else {
    stop("Wrong model type is given")
  }
  N=N+Nadd
  
  Data <- as.data.frame(matrix(NA,npseu*N,12))
  names(Data)<-c("ID","X1","X2","X3","X4","X5","X6",
                 "Start","Stop","C","Event","Xi")
  Data$ID <- rep(1:N,each=npseu)
  ## time-invariant
  Data$X1 <- rep(sample(c(0,1),N,replace=TRUE),each=npseu)
  Data$X2 <- rep(runif(N,0,1),each=npseu)
  ## time-varying
  Data$X3 <- rbinom(npseu*N,size = 1,prob = 0.5)
  Data$X4 <- runif(npseu*N,0,1)
  Data$X5 <- sample(c(1:5),npseu*N,replace = TRUE)
  Data$X6 <- 2
  
  Fstar = rep(0, npseu*N)
  u <- runif(N)
  Count= 1
  while(Count <= N){
    TS <- rep(0, npseu-1)
    
    if (model == "linear") {
      while(any(diff(sort(TS)) < 0.2)){
        TS <- sort(rtrunc(npseu-1, spec="beta", a = 0.1, b = 3, shape1 = 0.2, shape2 = 4))*3000
      }
    } else if (model == "nonlinear"){
      while(any(diff(sort(TS)) < 2)){
        TS <- sort(rtrunc(npseu-1, spec="beta", b = 3, shape1 = 0.2, shape2 = 4))*3000
      }
    } else if (model == "interaction"){
      while(any(diff(sort(TS)) < 0.2)){
        TS <- sort(rtrunc(npseu-1, spec="beta", a = 0.1, b = 3, shape1 = 0.2, shape2 = 4))*3000
      }
    } else {
      stop("Wrong model type is given.")
    }
    
    x61 <- sample(c(0,1,2),1)
    x62 <- sample(npseu,1)-1
    Data[Data$ID==Count,]$X6[1:(1+x62)] <- rep(x61*(x61!=2) + 2*(x61==2), x62+1)
    if (x61 == 0 && x62 != (npseu - 1)){
      x63 <- sample(npseu-(1+x62),1)
      Data[Data$ID==Count,]$X6[(1+x62+1):(1+x62+x63)] <- rep(1,x63)
    }
    Data[Data$ID==Count,]$Start <- c(0,TS)
    Data[Data$ID==Count,]$Stop <- c(TS,NA)
    
    Fstar[Data$ID==Count] <- Comp_fstar_WI(TALL = c(0, TS), 
                                           M = model, 
                                           X = Data[Data$ID == Count, ], 
                                           variation = variation, 
                                           snrhigh = snrhigh)
    Count = Count + 1
  }## end of the while loop
  Fstar = (Fstar - min(Fstar)+0.005)/diff(range(Fstar))*3
  
  DATA <- NULL
  Count = 1
  CountR = 1
  N=N-Nadd
  tall <- rep(0,N)
  rIDall <- rep(0,N)
  while (CountR <= N){
    RT <- Comp_xi_WI(TALL=Data[Data$ID==Count,]$Start, M=model, Fstar=Fstar[Data$ID==Count], U=u[Count])
    Data[Data$ID==Count,]$Xi <- RT$Xi
    t = RT$Time
    rID = RT$Row
    checkseq <- c(Data[Data$ID==Count,][1:rID,]$Start,t)
    if (all(diff(checkseq)>0.2)){
      if(rID==1){
        Data[Data$ID==Count,][1,]$Stop = t
        Data[Data$ID==Count,][1,]$Event = 1
        Data[Data$ID==Count,][rID,]$Stop = t
      }else{
        Data[Data$ID==Count,][1:(rID-1),]$Event=0
        Data[Data$ID==Count,][rID,]$Event = 1
        Data[Data$ID==Count,][rID,]$Stop = t
      }
      DATA <- rbind(DATA,Data[Data$ID==Count,])
      DATA[DATA$ID==Count,]$ID = CountR
      tall[CountR] = t
      rIDall[CountR] = rID
      CountR = CountR+1
      # print(CountR)
    }
    Count = Count+1
  }
  rm(Data)
  rm(Fstar)
  gc()
  # sum((is.na(Data$Stop)==1)+(is.na(Data$Event)==0)==2)
  DATA <- DATA[!is.na(DATA$Event),]
  ###================== Add Censoring =========================================
  DATA$C <- 0
  RET <- NULL
  RET$evalData = DATA[,c("ID","X1","X2","X3","X4","X5","X6","Xi","Start","Stop","Event")]
  
  if (variation) { # only censor.rate = 1 is adjusted
    if (snrhigh) {
      if (model == "linear"){
        if(censor.rate == 0){
          Censor.time <- rep(Inf,N)
        }else if(censor.rate == 1){
          Censor.time = rexp(N,rate = 1/4202) # NEW
        }else if(censor.rate == 2){
          Censor.time = rexp(N,rate = 1/1242) # NEW
        } else{
          stop("Wrong censoring type")
        }
      }else if (model == "nonlinear"){
        if(censor.rate == 0){
          Censor.time <- rep(Inf,N)
        }else if(censor.rate == 1){
          Censor.time = rexp(N,rate = 1/2202) # NEW
        }else if(censor.rate == 2){
          Censor.time = rexp(N,rate = 1/642) # NEW
        } else{
          stop("Wrong censoring type.")
        }
      } else if (model == "interaction"){
        if(censor.rate == 0){
          Censor.time <- rep(Inf,N)
        }else if(censor.rate == 1){
          Censor.time = rexp(N,rate = 1/4202) # NEW
        }else if(censor.rate == 2){
          Censor.time = rexp(N,rate = 1/1332) # NEW
        } else{
          stop("Wrong censoring type.")
        }
      } else {
        stop("Wrong model type.")
      }
    } else { # snrhigh == FALSE
      if (model == "linear"){
        if(censor.rate == 0){
          Censor.time <- rep(Inf,N)
        }else if(censor.rate == 1){
          Censor.time = rexp(N,rate = 1/4202)
        }else if(censor.rate == 2){
          Censor.time = rexp(N,rate = 1/1252) # new
        } else{
          stop("Wrong censoring type")
        }
      }else if (model == "nonlinear"){
        if(censor.rate == 0){
          Censor.time <- rep(Inf,N)
        }else if(censor.rate == 1){
          Censor.time = rexp(N,rate = 1/2202)
        }else if(censor.rate == 2){
          Censor.time = rexp(N,rate = 1/682) # new
        } else{
          stop("Wrong censoring type.")
        }
      } else if (model == "interaction"){
        if(censor.rate == 0){
          Censor.time <- rep(Inf,N)
        }else if(censor.rate == 1){
          Censor.time = rexp(N,rate = 1/4202)
        }else if(censor.rate == 2){
          Censor.time = rexp(N,rate = 1/1322) # new
        } else{
          stop("Wrong censoring type.")
        }
      } else {
        stop("Wrong model type.")
      }
    }
    
  } else { # Basic DGP
    if (snrhigh) {
      if (model == "linear"){
        if(censor.rate == 0){
          Censor.time <- rep(Inf,N)
        }else if(censor.rate == 1){
          Censor.time = rexp(N,rate = 1/4202)
        }else if(censor.rate == 2){
          Censor.time = rexp(N,rate = 1/1270) # NEW
        }else{
          stop("Wrong censoring type")
        }
      }else if (model == "nonlinear"){
        if(censor.rate == 0){
          Censor.time <- rep(Inf,N)
        }else if(censor.rate == 1){
          Censor.time = rexp(N,rate = 1/2202)
        }else if(censor.rate == 2){
          Censor.time = rexp(N,rate = 1/652) # NEW
        }else{
          stop("Wrong censoring type.")
        }
      } else if (model == "interaction"){
        if(censor.rate == 0){
          Censor.time <- rep(Inf,N)
        }else if(censor.rate == 1){
          Censor.time = rexp(N,rate = 1/4202)
        }else if(censor.rate == 2){
          Censor.time = rexp(N,rate = 1/1302) # NEW
        }else{
          stop("Wrong censoring type.")
        }
      } else {
        stop("Wrong model type.")
      }
    } else { # snrhigh == FALSE
      if (model == "linear"){
        if(censor.rate == 0){
          Censor.time <- rep(Inf,N)
        }else if(censor.rate == 1){
          Censor.time = rexp(N,rate = 1/4202) 
        }else if(censor.rate == 2){
          Censor.time = rexp(N,rate = 1/1302) 
        }else{
          stop("Wrong censoring type")
        }
      }else if (model == "nonlinear"){
        if(censor.rate == 0){
          Censor.time <- rep(Inf,N)
        }else if(censor.rate == 1){
          Censor.time = rexp(N,rate = 1/2202) 
        }else if(censor.rate == 2){
          Censor.time = rexp(N,rate = 1/702) 
        }else{
          stop("Wrong censoring type.")
        }
      } else if (model == "interaction"){
        if(censor.rate == 0){
          Censor.time <- rep(Inf,N)
        }else if(censor.rate == 1){
          Censor.time = rexp(N,rate = 1/4202) 
        }else if(censor.rate == 2){
          Censor.time = rexp(N,rate = 1/1302) 
        }else{
          stop("Wrong censoring type.")
        }
      } else {
        stop("Wrong model type.")
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
  
  n = nrow(Data)
  ## noise variables
  Data$X14 <- sample(c(1:5),n,replace = TRUE)
  Data$X15 <- runif(n,0,1)
  Data$X17 <- runif(n,0,1)
  Data$X19 <- sample(c(0,1),n,replace=TRUE)
  Data$X7 = 0
  Data$X8 = 0
  Data$X9 = 0
  Data$X10 = 0
  Data$X11 = 0
  Data$X12 = 0
  Data$X13 = 0
  Data$X16 = 0
  Data$X18 = 0
  Data$X20 = 0
  for (Count in unique(Data$ID)){
    Data[Data$ID==Count, ]$X7 <- runif(1,0,1)
    Data[Data$ID==Count, ]$X8 <- runif(1,1,2)
    Data[Data$ID==Count, ]$X9 <- sample(c(1:5),1,replace=TRUE)
    Data[Data$ID==Count, ]$X10 <- runif(1,0,1)
    Data[Data$ID==Count, ]$X11 <- sample(c(0,1),1,replace=TRUE)
    Data[Data$ID==Count, ]$X12 <- sample(c(0,1,2),1,replace=TRUE)
    
    n0 = nrow(Data[Data$ID==Count, ])
    x13r <- sample(npseu-1,1)
    Data[Data$ID==Count,]$X13 <- c(numeric(x13r),rep(1,npseu-x13r))[1:n0]

    x16r <- sample(npseu-1,1)
    x161 <- sample(c(0,1),1)
    Data[Data$ID==Count,]$X16 <- c(rep(x161,x16r),rep(x161+1,npseu-x16r))[1:n0]
    
    x18r <- sort(sample(npseu-1,2))
    Data[Data$ID==Count,]$X18 <- c(rep(0,x18r[1]),rep(1,x18r[2]-x18r[1]),rep(2,npseu-x18r[2]))[1:n0]
    k = runif(2)
    Data[Data$ID==Count,]$X20 <- k[1] * Data[Data$ID==Count,]$Start + k[2]
  }
  
  info = list(Rtype = model, Coeff = RT$Coeff, Set="nonPH-tvary")
  RET$fullData = Data
  RET$Info = info
  
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
      ## at least 30% of the times are not observable
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

