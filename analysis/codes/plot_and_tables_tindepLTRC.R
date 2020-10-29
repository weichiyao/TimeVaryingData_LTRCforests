######################### ========== mtry =========================== ######################
############### ===== The following code is to reproduce mtry analysis===== ######################
modelname <- c("Linear", "Nonlinear","Interaction")
cratename <- c("Censoring: 20%", "Censoring: 50%")
ndataname <- c("N = 100", "N = 300", "N = 500", "N = 1000")
frname = c("LTRCCIF","LTRCRRF","LTRCTSF")
L2name <- c("old","new")
mmodel <- 1:3
ccrate = c(1,2)
nndata = c(100,300,500,1000)


cc = 1

setting = "PH"

L2ll = 2 ###### Decided to go with L2 up to the last observed time point
for (ff in 1:3){
    pdf_file_name <- sprintf("./LTRCmtry%s_%s_c%1.0f_%s.pdf",
                             frname[ff],setting,ccrate[cc],L2name[L2ll])
    pdf(pdf_file_name, width=15, height = 18)
    par(mfrow=c(4,3))
    for (nn in 1:4){
      for (mm in 1:3){ 
        aa = NULL
        filename <- sprintf('./L2LTRC_%s_N%1.0f_m%1.0f_c%1.0f.rds',
                            setting,nndata[nn],mmodel[mm],ccrate[cc])
        
        Q = matrix(0,ncol = 8,nrow = 500)
        
        if (ff == 1){
          ## L2 errors for cforest with different mtry
          L2mtryM = matrix(unlist(resall$L2mtry_cf),ncol = 7,byrow = FALSE)
        } else if (ff == 2){
          ## L2 errors forr RRF with different mtry
          L2mtry = matrix(unlist(resall$L2mtry_rrf),ncol = 7,byrow = FALSE)
        } else {
          ## L2 errors for TSF with different mtry
          L2mtry = matrix(unlist(resall$L2mtry_tsf),ncol = 7,byrow = FALSE)
        }
        
        ## L2 results with mtry = 1,2,3,5,20 and best mtry
        Q[,1:7] = L2mtry[, c(6,5,4,3,2,1,7)]
        
        ## L2 results with mtry tuned by OOB
        Q[,8] = L2[, ff*2+1]
        
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

######################### ========== Default or Proposed? =========== ######################
########### ==== The following code is to reproduce Table of Default vs Proposed======= ##########
modelname <- c("Linear", "Nonlinear","Interaction")
cratename <- c("Censoring: 20%", "Censoring: 50%")
ndataname <- c("N = 100", "N = 300", "N = 500", "N = 1000")
mmodel <- 1:3
ccrate = c(1,2)
nndata = c(100,300,500,1000)

cc = 1
Nfold = 10

setting = "PH"
mm = 4
xall = rep(0,8)
for (nn in 1:4){
  filename <- sprintf('./L2LTRC_%s_N%1.0f_m%1.0f_c%1.0f.rds',
                      setting,nndata[nn],mmodel[mm],ccrate[cc])
  
  resall <- readRDS(filename)
  Q = matrix(0, ncol = 7, nrow = 500)
  ## L2 results for Cox, cfD, cfP, rrfD, rrfP, tsfD, tsfp
  Q[,1:7] = matrix(unlist(resall$L2), ncol = 7, byrow = FALSE)
  
  Q = Q[rowSums(Q==0)==0,]
  Q = Q[rowSums(Q==Inf)==0,]
  print(sprintf("nn=%1.0f",nn))
  Iall = matrix(0,nrow = nrow(Q), ncol = ncol(Q))
  for (ll in 1:nrow(Q)){
    Iall[ll,1:7]=(Q[ll,1:7]-Q[ll,1])/Q[ll,1]
  }
  
  for (ll in 1:nrow(Q)){
    Iall[ll,8:14]=(Q[ll,8:14]-Q[ll,8])/Q[ll,8]
  }
  
  IALL = colMeans(Iall)
  
  xrow = format_number(IALL[-c(1,8)])
  xall[nn] = sprintf("%1.0f & %s",nndata[nn],format_row(xrow[1:6]))
  xall[nn+4] = sprintf("%1.0f & %s",nndata[nn],format_row(xrow[7:12]))
}

format_all(xall)


######################### ========== Boxplots of L2 =========== ##################################
##### ==== The following code is to reproduce the plots of performance comparison, including IBSCV==== #######
modelname <- c("Linear", "Nonlinear","Interaction")
cratename <- c("Censoring: 20%", "Censoring: 50%")
ndataname <- c("N = 100", "N = 300", "N = 500", "N = 1000")


mmodel <- 1:3
ccrate = c(1,2)
nndata = c(100,300,500,1000)

cc = 1

setting = "nonPH"

pdf_file_name <- sprintf("./IBSCV_%s_c%1.0f.pdf",
                         setting,ccrate[cc])
pdf(pdf_file_name, width=15, height = 18)
par(mfrow=c(4,3))
for (nn in 1:4){
  for (mm in 1:3){
    aa = NULL
    Q = matrix(0,ncol = 6, nrow = 500)
    filename <- sprintf('./L2LTRC_%s_N%1.0f_m%1.0f_c%1.0f.rds',
                        setting,nndata[nn],mmodel[mm],ccrate[cc])
    
    resall <- readRDS(filename)
    L2 = matrix(unlist(resall$L2),nrow = 500,byrow = FALSE)
    
    ## Get the L2 errors of Cox, cifP, rrfP, tsfP
    Q[,1:4]=L2[,c(1,3,5,7)]
    
    IIL = matrix(unlist(resall$ibsCVerr), ncol = 4, byrow = FALSE)
    
    idx_iL = apply(IIL, 1, which.min)
    
    nz_i = which(as.numeric(rowSums(IIL!=0)<4)>0)
    
    for (ll in 1:500){
      ## L2 error of the best method 
      Q[ll,5] = min(Q[ll,1:4])
      ## L2 error of the method chosen by OOB
      Q[ll,6] = Q[ll,idx_iL[ll]]
    }
    Q[nz_i,6] = 0 
    ks[1,nz_i] = 0
    Q = Q[rowSums(Q==0)==0,]
    
    boxplot.matrix(Q, main=c(modelname[mm],"",ndataname[nn]), cex.main=2, cex.lab=1.5,
                   xaxt="n", 
                   ylim = aa,
                   ylab = "Integrated L2")
    
    xtick <- c(1:6)
    text(x=xtick,  par("usr")[3],
         labels = c("1","2","3","4","Opt","CV"),#
         pos = 1, xpd = TRUE, cex = 1.6)
  }
}

dev.off()









########################################################################################

