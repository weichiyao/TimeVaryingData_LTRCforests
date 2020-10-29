#####################################################################################################
## Last updated Sept 28th: it can also deal with indepLTRC
## Last updated May 18th: to match Timevarying_linear_gnrt_0707_partial.R, Xi added 
##  == Also: T.pnt = T.pnt[T.pnt<=Testdata$Stop]
#####################################################################################################
Integrate <- function(f, time.pnt){
  f.value <- sapply(time.pnt, f)
  result <- diff(time.pnt)%*%(f.value[-length(f.value)] + f.value[-1])/2
  return(result)
}
#########################################################################
getSurv <- function(obj, times)
{
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

#####################################################################################################
Surv_funct_tindep <- function(X, t, info){
  Xi = X$Xi
  if (info$Set == "PH"){
    Distribution = info$Dist
    if (Distribution == "Exp"){
      TD = info$Coeff$Lambda * t
    } else if (Distribution == "WD"){
      TD = info$Coeff$Lambda * t^info$Coeff$V
    } else if (Distribution == "WI"){
      TD = info$Coeff$Lambda * t^info$Coeff$V
    } else if (Distribution == "Gtz"){
      TD = info$Coeff$Lambda * exp(info$Coeff$Alpha*t)/info$Coeff$Alpha
    } else {
      stop("Wrong distribution is given when continuous X5 is presented")
    } ## end loop for distribution
    R = TD*Xi
  } else if (info$Set == "nonPH"){
    R = (info$Coeff$Lambda^Xi)*(t^Xi)
  }
  S = exp(-R)
  return(S)
}

#####################################################################################################
Loss_funct_tindep <- function(KM, Testdata, Info, T.pnt){
  ## Compute the true survival probability of the i-th data
  Surv_t <- function(t){
    Surv_funct_tindep(Testdata, t, Info)
  }
  
  if ("survfit" %in% class(KM)){
    integrand <- function(t){
      ( ipred::getsurv(KM,t) - Surv_t(t))^2
    }
    
    Int.value <- Integrate(f = integrand, time.pnt = T.pnt)
    return(Int.value/diff(range(T.pnt)))
    
  } else if ("icfit" %in% class(KM)){
    
    integrand <- function(t){
      (interval::getsurv(t,KM,nonUMLE.method = "right")[[1]]$S - Surv_t(t))^2
    }
    
    Int.value <- Integrate(f = integrand, time.pnt = T.pnt)
    return(Int.value/diff(range(T.pnt)))
    
  } else if (class(KM)=="sp_curves"||"rfsrc" %in% class(KM)||class(KM) == "ranger"){
    integrand <- function(t){
      (sapply(t, Interpolate, Curve = KM ) - Surv_t(t))^2
    }
    
    Int.value <- Integrate(f = integrand, time.pnt = T.pnt)
    return(Int.value/diff(range(T.pnt)))
    
  } else if (class(KM) == "list"){
    integrand <- function(t){
      (sapply(t, Interpolate, Curve = KM ) - Surv_t(t))^2
    }
    
    Int.value <- Integrate(f = integrand, time.pnt = T.pnt)
    return(Int.value/diff(range(T.pnt)))    
  } else if (class(KM)[1] == "ic_ph"||class(KM)[1] == "ic_np"){
    integrand <- function(t){
      (1-icenReg::getFitEsts(fit = KM, q = t) - Surv_t(t))^2
    }
    
    Int.value <- Integrate(f = integrand, time.pnt = T.pnt)
    return(Int.value/diff(range(T.pnt)))   
  } else if (class(KM) == "numeric"){
    KM = KM[1:length(T.pnt)]
    if (length(KM) != length(T.pnt)){
      stop("Please check whether the time of insterest are the same!")
    }
    f.value <- (KM - Surv_t(T.pnt))^2
    Int.value <- diff(T.pnt)%*%(f.value[-length(f.value)] + f.value[-1])/2
    return(Int.value/diff(range(T.pnt)))
  } else {
    stop("wrong class of NPMLE object specified")
  }
}## end of function

L2_tindep <- function(KM, Data, Info, T.pnt){
  n_uniq = nrow(Data)
  
  L2res <- sapply(1:n_uniq, function(Ni){
    # up to the last observed time
    tpnt <- T.pnt[T.pnt <= Data$Stop[Ni]]
    if (class(KM)[1] == "survfitcox"){
      pred <- KM[Ni]
    } else if (class(KM)[1] == "matrix"){
      pred <- KM[,Ni]
    }  else {
      pred <- KM[[Ni]]
    }
    Loss_funct_tindep(pred, Data[Ni, ], Info, tpnt)
  })
  ret <- mean(L2res)
  return(ret)
}

