#####################################################################################################
## Last updated April 11th: add evalData
## Last updated March 8th: add option partial=TRUE
## Last updated March 5th: add option partial=c("None","Partial","Baseline")
#####################################################################################################
itct_term <- function(beta,x1,x2,x3,x4,x5,x6){
  R1 = x1*x2-log(x3+x4)-x6/x5
  R2 = beta[1]*x1+beta[2]*x2+beta[3]*x3+beta[4]*x4+beta[5]*x5+beta[6]*x6
  R3 = cos(pi*(x1+x5))+sqrt(x6+x2)-x3
  R4 = -beta[1]*x1+beta[2]*x2-beta[3]*x3+beta[4]*x4-beta[5]*x5+beta[6]*x6
  R0 = (x4>=0.7)*((x5==2)+(x5==5))*R1+(x4>=0.7)*(1-(x5==2)-(x5==5))*R2+
    (x4<0.7)*((x5==2)+(x5==5))*R3 + (x4<0.7)*(1-(x5==2)-(x5==5))*R4
  return(R0)
}
#####################################################################################################
## Range_T function returns simulated survival time 
Range_T_itct <- function(TALL, DIST, X, U){  # TALL = c(0,TS)
  u = U
  tlen = length(TALL)
  Beta <- c(1,1,-1,-1,0.5,-0.5)
  R0 = sapply(1:tlen, function(i) itct_term(beta=Beta,
                                            x1=X$X1[i],x2=X$X2[i],
                                            x3=X$X3[i],x4=X$X4[i],
                                            x5=X$X5[i],x6=X$X6[i]))
  R0 = exp(R0)
  if(DIST == "Exp"){
    Lambda = 0.05
    Alpha = 0
    V = 0
    
    R = Lambda*R0[-tlen]*(TALL[-1]-TALL[-tlen])
    VEC = c(0,cumsum(R),Inf)
    R.ID <- findInterval(-log(u), VEC)
    TT = -log(u) - VEC[R.ID] 
    TT = TT/Lambda/R0[R.ID] + TALL[R.ID]
  }else if(DIST == "WD"){
    Lambda = 0.14
    V = 0.5
    Alpha = 0
    R = Lambda*R0[-tlen]*(TALL[-1]^V-TALL[-tlen]^V)
    
    VEC = c(0,cumsum(R),Inf)
    
    R.ID <- findInterval(-log(u), VEC)
    TT = -log(u) - VEC[R.ID] 
    TT = (TT/Lambda/R0[R.ID] + TALL[R.ID]^V)^(1/V)
    
  }else if(DIST == "WI"){
    Lambda = 0.0001
    #V = 2
    V = 1.8
    Alpha = 0
    R = Lambda*R0[-tlen]*(TALL[-1]^V-TALL[-tlen]^V)
    
    VEC = c(0,cumsum(R),Inf)
    R.ID <- findInterval(-log(u), VEC)
    TT = -log(u) - VEC[R.ID] 
    TT = (TT/Lambda/R0[R.ID]+TALL[R.ID]^V)^(1/V)
    
  }else if(DIST == "Gtz"){
    Alpha = 0.1
    Lambda = 0.015
    V = 0
    R = Lambda*R0[-tlen]/Alpha*(exp(Alpha*TALL[-1])-exp(Alpha*TALL[-tlen]))
    
    VEC = c(0,cumsum(R),Inf)
    R.ID <- findInterval(-log(u), VEC)
    
    TT = Alpha*(-log(u)-VEC[R.ID])/Lambda/R0[R.ID]
    TT = log(TT + exp(Alpha*TALL[R.ID]))/Alpha
  }
  
  result = list(Time=TT, Row=R.ID, Lambda=Lambda, Beta=Beta, Alpha=Alpha, V=V, Xi = R0)
  return(result)
}

#####################======== Large number of pseudo-subjects ===============#####################
Timevarying_itct_gnrt <- function(N = 200, Distribution = "Exp", censor.rate = 1, npseu = 11,
                                  partial = TRUE){
  Data <- as.data.frame(matrix(NA,npseu*N,87))
  names(Data)<-c("I","ID","X1","X2","X3","X4","X5","X6","X7","X8","X9","X10", 
                 "X11","X12","X13","X14","X15","X16","X17","X18","X19","X20", 
                 "X21","X22","X23","X24","X25","X26","X27","X28","X29","X30", 
                 "X31","X32","X33","X34","X35","X36","X37","X38","X39","X40",
                 "X41","X42","X43","X44","X45","X46","X47","X48","X49","X50",
                 "X51","X52","X53","X54","X55","X56","X57","X58","X59","X60",
                 "X61","X62","X63","X64","X65","X66","X67","X68","X69","X70",
                 "X71","X72","X73","X74","X75","X76","X77","X78","X79","X80", 
                 "Start","Stop","C","Event","Xi")
  Data$ID <- rep(1:N,each=npseu)
  ## time-invariant
  Data$X1 <- rep(sample(c(0,1),N,replace=TRUE),each=npseu)
  Data$X2 <- rep(runif(N,0,1),each=npseu)
  Data$X3 <- rbinom(npseu*N,size = 1,prob = 0.5)
  Data$X4 <- runif(npseu*N,0,1)
  Data$X5 <- sample(c(1:5),npseu*N,replace = TRUE)
  Data$X6 <- 2
  Data$X27 <- 2
  Data$X40 <- 2
  Data$X52 <- 2
  Data$X65 <- 2
  Data$X80 <- 2
  Data$X7 <- rep(runif(N,0,1),each=npseu)
  Data$X8 <- rep(runif(N,1,2),each=npseu)
  Data$X31 <- rep(runif(N,1,2),each=npseu)
  Data$X32 <- rep(runif(N,1,2),each=npseu)
  Data$X42 <- rep(runif(N,1,2),each=npseu)
  Data$X55 <- rep(runif(N,1,2),each=npseu)
  Data$X72 <- rep(runif(N,1,2),each=npseu)
  Data$X9 <- rep(sample(c(1:5),N,replace=TRUE),each=npseu)
  Data$X33 <- rep(sample(c(1:5),N,replace=TRUE),each=npseu)
  Data$X43 <- rep(sample(c(1:5),N,replace=TRUE),each=npseu)
  Data$X56 <- rep(sample(c(1:5),N,replace=TRUE),each=npseu)
  Data$X73 <- rep(sample(c(1:5),N,replace=TRUE),each=npseu)
  Data$X10 <- rep(runif(N,0,1),each=npseu)
  Data$X11 <- rep(sample(c(0,1),N,replace=TRUE),each=npseu)
  Data$X44 <- rep(sample(c(0,1),N,replace=TRUE),each=npseu)
  Data$X57 <- rep(sample(c(0,1),N,replace=TRUE),each=npseu)
  Data$X74 <- rep(sample(c(0,1),N,replace=TRUE),each=npseu)
  Data$X12 <- rep(sample(c(0,1,2),N,replace=TRUE),each=npseu) 
  Data$X45 <- rep(sample(c(0,1,2),N,replace=TRUE),each=npseu) 
  Data$X58 <- rep(sample(c(0,1,2),N,replace=TRUE),each=npseu) 
  Data$X75 <- rep(sample(c(0,1,2),N,replace=TRUE),each=npseu) 
  Data$X21 <- rep(runif(N,0,1),each=npseu)
  Data$X41 <- rep(runif(N,0,1),each=npseu)
  Data$X54 <- rep(runif(N,0,1),each=npseu)
  Data$X71 <- rep(runif(N,0,1),each=npseu)
  Data$X22 <- rep(sample(c(0,1),N,replace=TRUE),each=npseu)
  Data$X34 <- rep(sample(c(0,1),N,replace=TRUE),each=npseu)
  Data$X26 <- rep(sample(c(0,1,2),N,replace=TRUE),each=npseu)
  Data$X35 <- rep(sample(c(0,1,2),N,replace=TRUE),each=npseu)
  ## time-varying 
  Data$X14 <- sample(c(1:5),npseu*N,replace = TRUE)
  Data$X47 <- sample(c(1:5),npseu*N,replace = TRUE)
  Data$X60 <- sample(c(1:5),npseu*N,replace = TRUE)
  Data$X76 <- sample(c(1:5),npseu*N,replace = TRUE)
  Data$X15 <- runif(npseu*N,0,1)
  Data$X17 <- runif(npseu*N,0,1)
  Data$X23 <- runif(npseu*N,0,1)
  Data$X48 <- runif(npseu*N,0,1)
  Data$X61 <- runif(npseu*N,0,1)
  Data$X67 <- runif(npseu*N,0,1)
  Data$X77 <- runif(npseu*N,0,1)
  Data$X68 <- runif(npseu*N,1,2)
  Data$X69 <- runif(npseu*N,1,2)
  Data$X70 <- runif(npseu*N,1,2)
  Data$X79 <- runif(npseu*N,1,2)
  Data$X19 <- sample(c(0,1),npseu*N,replace=TRUE) 
  Data$X51 <- sample(c(0,1),npseu*N,replace=TRUE) 
  Data$X64 <- sample(c(0,1),npseu*N,replace=TRUE) 
  Data$X78 <- sample(c(0,1),npseu*N,replace=TRUE) 
  Data$X24 <- sample(c(0,1),npseu*N,replace=TRUE) 
  Data$X25 <- sample(c(1:5),npseu*N,replace = TRUE)
  Data$X39 <- sample(c(0,1),npseu*N,replace=TRUE) 
  
  Count = 1
  while(Count <= N){
    TS <- rep(0, npseu-1)
    if (Distribution == "Exp"){
      while(any(diff(sort(TS))<0.4)){
        #TS <- sort(rtrunc(npseu-1, spec="beta", a=0.0001, shape1 = 0.001, shape2 = 10))*400
        # TS <- sort(runif(npseu-1, 0.6,50))
        TS <- sort(rtrunc(npseu-1, spec="beta", a=0.0001, b = 9, shape1 = 0.001, shape2 = 20))*800
      }
    } else if (Distribution == "WD"){
      while(any(diff(sort(TS))<0.4)){
        TS <- sort(rtrunc(npseu-1, spec="beta", a=0.0001, shape1 = 0.001, shape2 = 2))*200
      }
    } else if (Distribution == "WI"){
      while(any(diff(sort(TS))<0.4)){
        TS <- sort(rtrunc(npseu-1, spec="beta",a=0.0005, b = 4, shape1 = 0.5, shape2 = 6))*1000
      }
    } else if (Distribution == "Gtz"){
      while(any(diff(sort(TS))<0.4)){
        TS <- sort(rtrunc(npseu-1, spec="beta", a=0.001, b = 4, shape1 = 0.001, shape2 = 6))*200
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
    
    x271 <- sample(c(0,1,2),1)
    x272 <- sample(npseu,1)-1
    Data[Data$ID==Count,]$X27[1:(1+x272)] <- rep(x271*(x271!=2) + 2*(x271==2), x272+1)
    if (x271 == 0 && x272 != (npseu - 1)){
      x273 <- sample(npseu-(1+x272),1)
      Data[Data$ID==Count,]$X27[(1+x272+1):(1+x272+x273)] <- rep(1,x273)
    }
    
    x401 <- sample(c(0,1,2),1)
    x402 <- sample(npseu,1)-1
    Data[Data$ID==Count,]$X40[1:(1+x402)] <- rep(x401*(x401!=2) + 2*(x401==2), x402+1)
    if (x401 == 0 && x402 != (npseu - 1)){
      x403 <- sample(npseu-(1+x402),1)
      Data[Data$ID==Count,]$X40[(1+x402+1):(1+x402+x403)] <- rep(1,x403)
    }
    
    x521 <- sample(c(0,1,2),1)
    x522 <- sample(npseu,1)-1
    Data[Data$ID==Count,]$X52[1:(1+x522)] <- rep(x521*(x521!=2) + 2*(x521==2), x522+1)
    if (x521 == 0 && x522 != (npseu - 1)){
      x523 <- sample(npseu-(1+x522),1)
      Data[Data$ID==Count,]$X52[(1+x522+1):(1+x522+x523)] <- rep(1,x523)
    }
    
    x651 <- sample(c(0,1,2),1)
    x652 <- sample(npseu,1)-1
    Data[Data$ID==Count,]$X65[1:(1+x652)] <- rep(x651*(x651!=2) + 2*(x651==2), x652+1)
    if (x651 == 0 && x652 != (npseu - 1)){
      x653 <- sample(npseu-(1+x652),1)
      Data[Data$ID==Count,]$X65[(1+x652+1):(1+x652+x653)] <- rep(1,x653)
    }
    
    x801 <- sample(c(0,1,2),1)
    x802 <- sample(npseu,1)-1
    Data[Data$ID==Count,]$X80[1:(1+x802)] <- rep(x801*(x801!=2) + 2*(x801==2), x802+1)
    if (x801 == 0 && x802 != (npseu - 1)){
      x803 <- sample(npseu-(1+x802),1)
      Data[Data$ID==Count,]$X80[(1+x802+1):(1+x802+x803)] <- rep(1,x803)
    }
    
    x13r <- sample(npseu-1,1)
    Data[Data$ID==Count,]$X13 <- c(numeric(x13r),rep(1,npseu-x13r))
    x28r <- sample(npseu-1,1)
    Data[Data$ID==Count,]$X28 <- c(numeric(x28r),rep(1,npseu-x28r))
    x46r <- sample(npseu-1,1)
    Data[Data$ID==Count,]$X46 <- c(numeric(x46r),rep(1,npseu-x46r))
    x59r <- sample(npseu-1,1)
    Data[Data$ID==Count,]$X59 <- c(numeric(x59r),rep(1,npseu-x59r))
    
    x16r <- sample(npseu-1,1)
    x161 <- sample(c(0,1),1)
    Data[Data$ID==Count,]$X16 <- c(rep(x161,x16r),rep(x161+1,npseu-x16r))
    x29r <- sample(npseu-1,1)
    x291 <- sample(c(1,2),1)
    Data[Data$ID==Count,]$X29 <- c(rep(x291,x29r),rep(x291+1,npseu-x29r))
    x36r <- sample(npseu-1,1)
    x361 <- sample(c(1,2),1)
    Data[Data$ID==Count,]$X36 <- c(rep(x361,x36r),rep(x361+1,npseu-x36r))
    x49r <- sample(npseu-1,1)
    x491 <- sample(c(1,2),1)
    Data[Data$ID==Count,]$X49 <- c(rep(x491,x49r),rep(x491+1,npseu-x49r))
    x62r <- sample(npseu-1,1)
    x621 <- sample(c(1,2),1)
    Data[Data$ID==Count,]$X62 <- c(rep(x621,x62r),rep(x621+1,npseu-x62r))
    
    
    x18r <- sort(sample(npseu-1,2))
    Data[Data$ID==Count,]$X18 <- c(rep(0,x18r[1]),rep(1,x18r[2]-x18r[1]),rep(2,npseu-x18r[2]))
    x37r <- sort(sample(npseu-1,2))
    Data[Data$ID==Count,]$X37 <- c(rep(0,x37r[1]),rep(1,x37r[2]-x37r[1]),rep(2,npseu-x37r[2]))
    x38r <- sort(sample(npseu-1,2))
    Data[Data$ID==Count,]$X38 <- c(rep(0,x38r[1]),rep(1,x38r[2]-x38r[1]),rep(2,npseu-x38r[2]))
    x50r <- sort(sample(npseu-1,2))
    Data[Data$ID==Count,]$X50 <- c(rep(0,x50r[1]),rep(1,x50r[2]-x50r[1]),rep(2,npseu-x50r[2]))
    x63r <- sort(sample(npseu-1,2))
    Data[Data$ID==Count,]$X63 <- c(rep(0,x63r[1]),rep(1,x63r[2]-x63r[1]),rep(2,npseu-x63r[2]))
    
    
    Data[Data$ID==Count,]$X20 <- k[1] * c(0,TS) + k[2]
    Data[Data$ID==Count,]$X30 <- k[2] * c(0,TS) + k[1]
    
    k2 <- runif(2,-2,-1)
    Data[Data$ID==Count,]$X53 <- k2[1] * c(0,TS) + k2[2]
    Data[Data$ID==Count,]$X66 <- k2[2] * c(0,TS) + k2[1]
    
    Data[Data$ID==Count,]$Start <- c(0,TS)
    Data[Data$ID==Count,]$Stop <- c(TS,NA)
    
    RT <- Range_T_itct(TALL=c(0,TS),DIST=Distribution,X=Data[Data$ID==Count,],U=u)
    t = RT$Time
    rID = RT$Row
    Data[Data$ID==Count,]$Xi <- RT$Xi
    
    if(rID==1){
      Data[Data$ID==Count,][1,]$Stop = t
      Data[Data$ID==Count,][1,]$Event = 1
    }else{
      Data[Data$ID==Count,][1:(rID-1),]$Event=0
      Data[Data$ID==Count,][rID,]$Event = 1
      Data[Data$ID==Count,][rID,]$Stop = t
    }
    # print(Count)
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
  RET$evalData = DATA[,c("ID","X1","X2","X3","X4","X5","X6","Xi","Start","Stop","Event")]
  
  if (npseu == 11){ ### NOT YET ADJUSTED FOR NPSEU = 11
    if(censor.rate == 0){
      Censor.time <- rep(Inf,N)
    }else if(censor.rate == 1){
      if(Distribution == "Exp"){
        Censor.time = rexp(N,rate = 1/85)
      }else if(Distribution == "WD"){
        Censor.time = rexp(N,rate = 1/310)
      }else if(Distribution == "WI"){
        Censor.time = rexp(N,rate = 1/700)
      }else if(Distribution == "Gtz"){
        Censor.time = rexp(N,rate = 1/88)
      }
    }else if(censor.rate == 2){
      if(Distribution == "Exp"){
        Censor.time = rexp(N,rate = 1/18)
      }else if(Distribution == "WD"){
        Censor.time = rexp(N,rate = 1/65)
      }else if(Distribution == "WI"){
        Censor.time = rexp(N,rate = 1/73)
      }else if(Distribution == "Gtz"){
        Censor.time = rexp(N,rate = 1/15)
      }
    }else{
      stop("Wrong censoring type")
    }
  } else if (npseu == 5){
    if(censor.rate == 0){
      Censor.time <- rep(Inf,N)
    }else if(censor.rate == 1){
      if(Distribution == "Exp"){
        Censor.time = rexp(N,rate = 1/130)
      }else if(Distribution == "WD"){
        Censor.time = rexp(N,rate = 1/800)
      }else if(Distribution == "WI"){
        Censor.time = rexp(N,rate = 1/531)
      }else if(Distribution == "Gtz"){
        Censor.time = rexp(N,rate = 1/90)
      }
    }else if(censor.rate == 2){ ## NOT YET ADJUSTED FOR CENSOR.RATE = 2
      if(Distribution == "Exp"){
        Censor.time = rexp(N,rate = 1/20)
      }else if(Distribution == "WD"){
        Censor.time = rexp(N,rate = 1/135)
      }else if(Distribution == "WI"){
        Censor.time = rexp(N,rate = 1/50)
      }else if(Distribution == "Gtz"){
        Censor.time = rexp(N,rate = 1/13)
      }
    }else{
      stop("Wrong censoring type.")
    }
  } else {
    stop("Wrong npseu is given.")
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
    partialInfo = data.frame(matrix(0,nrow = N,ncol = 7))
    names(partialInfo) = c("n","n_unobv","percT_unobv","wXiabsDiff","n_mto","percT_unobv_mto","wXiabsDiff_mto")
    for (Count in IDrev) {
      pN = sum(Data$ID == Count)
      I_count = Data[Data$ID==Count,]$I
      ## at least 30% of the times are not observable
      n_unobv = floor(pN*0.5)
      partialInfo[Count,1]=pN # number of pseudo-subject in the full dataset
      partialInfo[Count,2]=n_unobv # number of unobserved in the partial dataset
      if (pN > 1){
        unobv = sort(sample(1:(pN-1),n_unobv))+1
        tLast = Data[Data$ID==Count,]$Stop[pN]
        EventLast = Data[Data$ID==Count,]$Event[pN]
        
        T_unobv = Data[I_count[unobv],]$Stop
        dt_unobv = Data[I_count[unobv],]$Stop - Data[I_count[unobv],]$Start
        Xi_unobv = Data[I_count[unobv],]$Xi
        partialInfo[Count,3]= sum(dt_unobv)/tLast
        
        Data = Data[-I_count[unobv],]
        if (pN > 2){
          Data[Data$ID==Count,]$Stop[1:(pN-n_unobv-1)] = Data[Data$ID==Count,]$Start[2:(pN-n_unobv)]
        } 
        Data[Data$ID==Count,]$Stop[pN-n_unobv] = tLast
        Data[Data$ID==Count,]$Event[pN-n_unobv] = EventLast
        
        Xi_new_rID = findInterval(T_unobv, c(0,Data[Data$ID==Count,]$Stop),
                                  left.open = TRUE, rightmost.closed = TRUE)
        partialInfo[Count,4]=sum(abs(Data[Data$ID==Count,]$Xi[Xi_new_rID]-Xi_unobv)/Xi_unobv *dt_unobv)/sum(dt_unobv)
      }    
    }
    Data$I = 1:nrow(Data)
    RET$partialData = Data
    
    RET$partialInfo = round(colMeans(partialInfo),digits = 3)
    RET$partialInfo[5] = sum(partialInfo$n_unobv==0)
    RET$partialInfo[6] = round(mean(partialInfo$percT_unobv[partialInfo$percT_unobv!=0]),digits = 3)
    RET$partialInfo[7] = round(mean(partialInfo$wXiabsDiff[partialInfo$wXiabsDiff!=0]),digits = 3)
  }
  
  rm(Data)
  gc()
  return(RET)
  
}
