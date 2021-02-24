#########################################################################
Surv_funct_nonPH <- function(X,info,t){
  TALL = X$Start
  Xi = X$Xi
  R.ID = findInterval(t,TALL)
  
  Lambda = info$Coeff
  if (R.ID == 1){
    R = (Lambda^Xi[1])*(t^Xi[1])
  } else{
    R = Lambda^Xi[1:R.ID]*(c(TALL[2:R.ID]^Xi[1:(R.ID-1)],t^Xi[R.ID])-TALL[1:R.ID]^Xi[1:R.ID])
  }
  
  S = exp(-sum(R))
  return(S)
}
#########################################################################
Surv_funct_PH <- function(X, info, t){
  TALL = X$Start
  Distribution = info$Dist
  Xi = X$Xi
  Lambda = info$Coeff$Lambda
  Alpha = info$Coeff$Alpha
  V = info$Coeff$V
  Beta = info$Coeff$Beta
  R.ID = findInterval(t,TALL)

  if (Distribution == "WD") {
    if (R.ID == 1){
      TD = Lambda * (t^V - TALL[1]^V)
    } else {
      TD = Lambda * (c(TALL[2:R.ID]^V,t^V) - TALL[1:R.ID]^V)
    }
  } else if (Distribution == "WI") {
    if (R.ID == 1){
      TD = Lambda * (t^V - TALL[1]^V)
    } else {
      TD = Lambda * (c(TALL[2:R.ID]^V,t^V) - TALL[1:R.ID]^V)
    }
  } else {
    stop("Wrong distribution is given")
  } ## end loop for distribution
  
  R = Xi[1:R.ID]*TD
  S = exp(-sum(R))
  return(S)
}
### this one does have Xi 
#########################################################################
s_funct <- function(Ni, fulldata, info, tpnt){
  ## This function is to compute the estimated survival probability of the Ni-th subject
  id_uniq <- unique(fulldata$ID)
  ## the i-th data 
  testData = fulldata[fulldata$ID == id_uniq[Ni],]
  
  sfun <- switch(info$Set,
                 "PH" = Surv_funct_PH,
                 "nonPH-tvary" = Surv_funct_nonPH)
  ## Compute the true survival probability of the Ni-th subject with all observations
  surv_t <- function(t){
    sfun(X = testData, info = info, t)
  }
  St <- sapply(tpnt, surv_t)
  return(St)
}

#########################################################################
getSurv <- function(obj, times) {
  # get the survival probability for times from KM curve `obj'
  lt <- length(times)
  nsurv <- times
  
  if (inherits(obj, "survfit")){
    objTime = obj$time
    objSurv = obj$surv
  } else {
    objTime = obj$Time
    objSurv = obj$Survival
  }

  # if the times are the same, return the km-curve
  if(length(times) == length(objTime)) {
    if (all(times == objTime)) return(objSurv)
  }
  
  # otherwise get the km-value for every element of times separatly
  inside <- times %in% objTime
  for (i in (1:lt)) {
    if (inside[i])
      nsurv[i] <- objSurv[objTime == times[i]]
    else  {
      less <- objTime[objTime < times[i]]
      if (length(less) == 0) 
        nsurv[i] <- 1
      else 
        nsurv[i] <- objSurv[objTime == max(less)]
    }
  }
  
  
  nsurv
}

