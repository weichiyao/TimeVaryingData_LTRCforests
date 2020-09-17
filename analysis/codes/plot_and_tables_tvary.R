format_number <- function(x) {
  sprintf('{\\tablenum[table-format = 1.2]{%1.2f}}', x)
}

format_row <- function(xrow) {
  paste(xrow, collapse=" & ")
}


format_all <- function(xall) {
  cat(paste(xall, collapse=" \\\\\n "))
}
######################### ================== mtry =========================== ######################
########### ======= The following code is to reproduce mtry analysis ===== ######################
ddist = c("Exp","WD","WI","Gtz")
modelname <- c("Tree", "Linear", "Nonlinear","Interaction")
cratename <- c("Censoring: 0", "Censoring: 20%", "Censoring: 50%")
ndataname <- c("N = 50", "N = 100", "N = 300", "N = 500")
distname <- c("Exponential","Weibull-D","Weibull-I","Gompertz")
frname = c("CF","RRF","TSF")
mmodel <- 1:4
ccrate = c(0,1,2)
nndata = c(50,100,300,500)
nnpseu = c(11,5)

varN = 20

pp = 1
cc = 2
varN = 20

setting = "nonPH"
DD = 2:4
for (ff in 1:3){
  if (setting == "nonPH") DD=3
  for (dd in DD){
    pdf_file_name <- sprintf("./mtry%s_%s_%1.0fvar_%s_npseu%1.0f_c%1.0f.pdf",
                             frname[ff],setting,varN,ddist[dd],nnpseu[pp],ccrate[cc])
    pdf(pdf_file_name, width=15, height = 18)
    par(mfrow=c(4,3))
    for (nn in 1:4){
      for (mm in 2:4){ 
        aa = NULL
        filename <- sprintf('./L2_%s_%1.0fvar_N%1.0f_%s_m%1.0f_c%1.0f_p%1.0f.rds',
                            setting,varN,nndata[nn],ddist[dd],mmodel[mm],ccrate[cc],nnpseu[pp])
        
        resall <- readRDS(filename)
        Q = matrix(0,ncol = 8,nrow = 500)
        
        if (ff == 1){
          ## L2 errors for LTRCCF with different mtry
          L2mtry = matrix(unlist(resall$L2mtry_cf),ncol = 7,byrow = FALSE)
        } else if (ff == 2){
          ## L2 errors for LTRCRRF with different mtry
          L2mtry = matrix(unlist(resall$L2mtry_rrf),ncol = 7,byrow = FALSE)
        } else {
          ## L2 errors for TSF with different mtry
          L2mtry = matrix(unlist(resall$L2mtry_tsf),ncol = 7,byrow = FALSE)
        }
        
        L2opt = matrix(unlist(resall$L2opt),ncol = 3,byrow = FALSE)
        ## L2 results with mtry = 1,2,3,5,20 and best mtry
        Q[,1:7] = L2mtry[,c(6,5,4,3,2,1,7)]
        
        ## L2 results with mtry tuned by OOB
        Q[,8] = L2[,ff]
        
        ## Get rid of the invalid results    
        Q = Q[rowSums(Q==Inf)==0,]
        Q = Q[rowSums(Q==0)==0,]
        
        boxplot.matrix(Q, main=c(modelname[mm],"",ndataname[nn]), cex.main=2, cex.lab=1.5,
                       xaxt="n", 
                       ylim = aa,
                       ylab = "Integrated L2")
        
        xtick <- c(1:8)
        text(x=xtick,  par("usr")[3],
             labels = c("1","2","3","5","10",
                        "20","Opt","Tuned"), 
             pos = 1, xpd = TRUE, cex = 1.6)
        
      }
    }
    dev.off()
  }
}


######################### ========== Default or Proposed? =========== ######################
########### ======= The following code is to reproduce Table of Default vs Proposed ======= ##########

ddist = c("Exp","WD","WI","Gtz")
modelname <- c("Tree", "Linear", "Nonlinear","Interaction")
cratename <- c("Censoring: 0", "Censoring: 20%", "Censoring: 50%")
ndataname <- c("N = 50", "N = 100", "N = 300", "N = 500")
distname <- c("Bathtub","Log-normal","Weibull-D","Weibull-I")
mmodel <- 1:4
ccrate = c(0,1,2)
nndata = c(50,100,300,500)
nnpseu = c(11,5)

varN = 20

pp = 1
cc = 2
varN = 20
Nfold = 10

setting = "nonPH"
dd = 3
mm = 3
xall = rep(0,8)
for (nn in 1:4){
  filename <- sprintf('./L2_%s_%1.0fvar_N%1.0f_%s_m%1.0f_c%1.0f_p%1.0f.rds',
                         setting,varN,nndata[nn],ddist[dd],mmodel[mm],ccrate[cc],nnpseu[pp])
  resall <- matrix(unlist(readRDS(filename)),nrow=500,byrow=FALSE)
  
  ## L2 results for KM fit
  L2KM = matrix(unlist(resall$caseI$L2KM),ncol = 1)
  
  ## CaseI : L2 results for Cox, cfD, cfP, rrfD, rrfP, tsfD, tsfp
  Ec[,1:7] = matrix(unlist(resall$caseI$L2), ncol = 7, byrow = FALSE)
  ## CaseII : L2 results for Cox, cfD, cfP, rrfD, rrfP, tsfD, tsfp
  Er[,1:7] = matrix(unlist(resall$caseII$L2), ncol = 7, byrow = FALSE)
  
  Ec = matrix(0,ncol = 7,nrow = 500)
  Er = matrix(0,ncol = 7,nrow = 500)
  
  Eall = cbind(L2KM,Ec,Er)
  
  ## Get rid of invarid results
  Eall = Eall[rowSums(Eall==0)==0,]
  Eall = Eall[rowSums(Eall==Inf)==0,]
  print(sprintf("nn=%1.0f",nn))
  Iall = matrix(0,nrow = nrow(Eall), ncol = ncol(Eall))
  for (ll in 1:nrow(Eall)){
    Iall[ll,]=(Eall[ll,]-Eall[ll,1])/Eall[ll,1]
  }
  
  IALL = colMeans(Iall)
  
  xrow = format_number(IALL[-1])
  xall[nn] = sprintf("%1.0f & %s",nndata[nn],format_row(xrow[1:7]))
  xall[nn+4] = sprintf("%1.0f & %s",nndata[nn],format_row(xrow[8:14]))
}

format_all(xall)
######################### ========== Boxplots of L2 =========== ######################
########### ======= The following code is to reproduce boxplots of L2 ======= ##########
ddist = c("Exp","WD","WI","Gtz")
modelname <- c("Tree", "Linear", "Nonlinear","Interaction")
cratename <- c("Censoring: 0", "Censoring: 20%", "Censoring: 50%")
ndataname <- c("N = 50", "N = 100", "N = 300", "N = 500")
distname <- c("Exponential","Weibull-D","Weibull-I","Gompertz")
mmodel <- 1:4
ccrate = c(0,1,2)
nndata = c(50,100,300,500)
nnpseu = c(11,5)

varN = 20

pp = 1
cc = 2

varN = 20
setting = "nonPH"
DD = 2:4
if (setting == "nonPH") DD=3
for (dd in DD){
  pdf_file_name <- sprintf("./L2_%s_%1.0fvar_%s_npseu%1.0f_c%1.0f_ROC.pdf",
                           setting,varN,ddist[dd],nnpseu[pp],ccrate[cc])
  pdf(pdf_file_name, width=17, height = 18)
  par(mfrow=c(4,3))
  for (nn in 1:4){
    for (mm in 2:4){ 
      aa = NULL
      
      filename <- sprintf('./L2_%s_%1.0fvar_N%1.0f_%s_m%1.0f_c%1.0f_p%1.0f.rds',
                          setting,varN,nndata[nn],ddist[dd],mmodel[mm],ccrate[cc],nnpseu[pp])
      resall <- matrix(unlist(readRDS(filename)),nrow=500,byrow=FALSE)
      
      ## L2 results for KM fit
      L2KM = matrix(unlist(resall$caseI$L2KM),ncol = 1)
      
      ## CaseI : L2 results for Cox, cfD, cfP, rrfD, rrfP, tsfD, tsfp
      Ec[,1:7] = matrix(unlist(resall$caseI$L2), ncol = 7, byrow = FALSE)
      ## CaseII : L2 results for Cox, cfD, cfP, rrfD, rrfP, tsfD, tsfp
      Er[,1:7] = matrix(unlist(resall$caseII$L2), ncol = 7, byrow = FALSE)
      
      
      E = cbind(L2KM, Ec, Er)
      
      E = E[rowSums(E==Inf)==0,]
      E = E[rowSums(E==0)==0,]
      ########## ======== presented way 1 ========= #######
      boxplot.matrix(E, main=c(modelname[mm],"",ndataname[nn]), cex.main=2, cex.lab=1.5,
                     xaxt="n", 
                     ylim = aa,
                     ylab = "Integrated L2")
      
      abline(v = 4.5, lty = 2, col = 4, lwd = 3)
      abline(h = E[, 1], lty = 2, col = 2, lwd = 3)
      
      xtick <- c(1:8)
      text(x=xtick,  par("usr")[3],
           labels = c("Cox","cfP","rrfP","tsfP",
                      "Cox","cfP","rrfP","tsfP"),
           pos = 1, xpd = TRUE, cex = 1.6)
      
    }
  }
  dev.off()
}
######################### ========== IBS based CV rule =========== ######################
########### ======= The following code is to reproduce Table of Summary of IBS-based CV rules ======= ############
ddist = c("Exp","WD","WI","Gtz")
modelname <- c("Tree", "Linear", "Nonlinear","Interaction")
cratename <- c("Censoring: 0", "Censoring: 20%", "Censoring: 50%")
ndataname <- c("Subject size N = 50", "Subject size N = 100", "Subject size N = 300", "Subject size N = 500")
distname <- c("Exponential","Weibull-D","Weibull-I","Gompertz")
mmodel <- 1:4
ccrate = c(0,1,2)
nndata = c(50,100,300,500)
nnpseu = c(11,5)

varN = 20

pp = 1
cc = 2

varN = 20
Nfold = 10
nn = 3

setting = "PH"
DD = 3
if (setting == "PH") DD = c(2:4)
for (dd in DD){
  for (mm in c(2,3,4)){
    Qc = matrix(0, ncol = 7, nrow = 500)
    Qr = matrix(0, ncol = 7, nrow = 500)
    filename <- sprintf('./L2_%s_%1.0fvar_N%1.0f_%s_m%1.0f_c%1.0f.rds',
                        setting,varN,nndata[nn],ddist[dd],mmodel[mm],ccrate[cc])
    
    resall <- readRDS(filename)
    L2c = matrix(unlist(resall$caseI$L2),nrow = 500,byrow = FALSE)
    L2r = matrix(unlist(resall$caseII$L2),nrow = 500,byrow = FALSE)
    
    ## Get the L2 errors of Cox, cfP, rrfP, tsfP for CaseI
    Qc[, 1:4] = L2c[, c(1,3,5,7)]
    ## Get the L2 errors of Cox, cfP, rrfP, tsfP for CaseII
    Qr[, 1:4] = L2r[, c(1,3,5,7)]
    
    IILc = matrix(unlist(resall$caseI$ibsCVerr), ncol = 4, byrow = FALSE)
    IILr = matrix(unlist(resall$caseII$ibsCVerr), ncol = 4, byrow = FALSE)
    
    idx_iLc = apply(IILc, 1, which.min)
    idx_iLr = apply(IILr, 1, which.min)
    
    nz_ic = which(as.numeric(rowSums(IILc != 0) < 4) > 0)
    nz_ir = which(as.numeric(rowSums(IILr != 0) < 4) > 0)
    
    for (ll in 1:500){
      ## L2 error of the method chosen by OOB for CaseI
      Qc[ll,5] = Qc[ll,idx_iLc[ll]]
      
      ## L2 error of the method chosen by OOB for CaseII
      Qr[ll,5] = Qr[ll,idx_iLr[ll]]
    }
    
  
    Qc[nz_ic, 5] = 0 
    Qr[nz_ir, 5] = 0 
    
    Qc = Qc[rowSums(Qc[,1:5]==0)==0,]
    Qr = Qr[rowSums(Qr[,1:5]==0)==0,]
    Qc = Qc[rowSums(Qc==Inf)==0,]
    Qr = Qr[rowSums(Qr==Inf)==0,]
    
    k_r = 0
    for (ll in 1:nrow(Qr)){
      L2min = min(Qr[ll, 1:4])
      L2max = max(Qr[ll, 1:4])
      Qr[ll,6] = abs(Qr[ll,5] - L2min)/L2min
      Qr[ll,7] = abs(Qr[ll,5] - L2max)/L2max
      if (Qr[ll,5]==L2min){
        k_r = k_r+1
      }
    }
    p_r = k_r/nrow(Er)
    
    k_c = 0
    for (ll in 1:nrow(Qc)){
      L2min = min(Qc[ll,1:4])
      L2max = max(Qc[ll,1:4])
      Qc[ll,6] = abs(Qc[ll,5]-L2min)/L2min
      Qc[ll,7] = abs(Qc[ll,5]-L2max)/L2max
      if (Qc[ll,5]==L2min){
        k_c = k_c+1
      }
    }
    p_c = k_c/nrow(Ec)
    
    
    Ire = matrix(0,nrow = 1,ncol = 6)
    Ire[1,1] =  p_c
    Ire[1,2:3] =  colMeans(Qc)[6:7]
    Ire[1,4] =  p_r
    Ire[1,5:6] =  colMeans(Qr)[6:7]
    Ire = format(round(Ire,digits=2), nsmall = 2, scientific=FALSE)
    Ire = format_row(Ire)
    if (mm == 2){
      cat(sprintf("\\multirow{9}{*}{%s}& \\multirow{3}{*}{%s}& %s & %s\\\\\n",setting,distname[dd],modelname[mm],Ire))
    } else if (mm == 3){
      cat(sprintf("&& %s & %s\\\\\n",modelname[mm],Ire))
    } else {
      cat(sprintf("&& %s & %s\\\\\n",modelname[mm],Ire))
    }
  }
  
}


################## ========= Plot Pseudo-subject vs Subject (side by side, only the forests) ========= ##################################################
########### ======= The following code is to reproduce the plot of Pseudo-subject vs Subject ======= ######################
ddist = c("Exp","WD","WI","Gtz")
modelname <- c("Tree", "Linear", "Nonlinear","Interaction")
cratename <- c("Censoring: 0", "Censoring: 20%", "Censoring: 50%")
ndataname <- c("N = 50", "N = 100", "N = 300", "N = 500")
distname <- c("Exponential","Weibull-D","Weibull-I","Gompertz")
mmodel <- 1:4
ccrate = c(0,1,2)
nndata = c(50,100,300,500)
nnpseu = c(11,5)
pp = 1
cc = 2

varN = 20
Nfold = 10
dd=3
for (setting in c("PH","nonPH")){
  pdf_file_name <- sprintf("./subseu_%s_%1.0fvar_%s_npseu%1.0f_c%1.0f.pdf",
                           setting,varN,ddist[dd],nnpseu[pp],ccrate[cc])
  pdf(pdf_file_name, width=15, height = 15)
  par(mfrow=c(4,3))
  aa = NULL
  for (nn in 1:4){
    for (mm in 2:4){ 
      filename <- sprintf('./L2_%s_%1.0fvar_N%1.0f_%s_m%1.0f_c%1.0f_p%1.0f.rds',
                          setting,varN,nndata[nn],ddist[dd],mmodel[mm],ccrate[cc],nnpseu[pp])
      resall_cf <- matrix(unlist(readRDS(filename)$caseI$L2mtry$cf), nrow=500, byrow=FALSE)[,3]
      resall_rrf <- matrix(unlist(readRDS(filename)$caseI$L2mtry$rrf), nrow=500, byrow=FALSE)[,3]
      resall_tsf <- matrix(unlist(readRDS(filename)$caseI$L2mtry$tsf), nrow=500, byrow=FALSE)[,3]
      resall_seu <- matrix(unlist(readRDS(filename)$caseI$L2seu), nrow=500, byrow = FALSE)
      E = cbind(resall_seu[,1],resall_cf,resall_seu[,2],resall_rrf,resall_seu[,3],resall_tsf)
      
      E = E[rowSums(E==0)==0,]
      #################
      boxplot.matrix(E, main=c(modelname[mm],"",ndataname[nn]), cex.main=2, cex.lab=1.5,
                     xaxt="n", 
                     ylim = aa,
                     ylab = "Integrated L2")
      
      abline(v=2.5,lty = 2, col = 4, lwd = 2)
      abline(v=4.5,lty = 2, col = 4, lwd = 2)
      
      xtick <- c(1:6)
      text(x=xtick,  par("usr")[3],
           labels = c("1","2","3","4","5","6"), pos = 1, xpd = TRUE, cex = 1.4)
    }
  }
  
  dev.off()
}




############################################################################################