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
modelname <- c("Linear", "Nonlinear","Interaction")
cratename <- c("Censoring: 20%", "Censoring: 50%")
ndataname <- c("N = 100", "N = 300", "N = 500")
frname = c("LTRCCIF","LTRCRRF","LTRCTSF")
mmodel <- 1:3
ccrate = c(1,2)
nndata = c(100,300,500)

cc = 2

setting = "nonPH"
for (ff in 1:3){
  
    pdf_file_name <- sprintf("./mtry%s_%s_c%1.0f.pdf",
                             frname[ff],setting,ccrate[cc])
    pdf(pdf_file_name, width=15, height = 14)
    par(mfrow=c(3,3))
    for (nn in 1:3){
      for (mm in 1:3){ 
        aa = NULL
        filename <- sprintf('./L2_%s_N%1.0f_m%1.0f_c%1.0f.rds',
                            setting,nndata[nn],mmodel[mm],ccrate[cc])
        
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


######################### ========== Default or Proposed? =========== ######################
########### ======= The following code is to reproduce Table of Default vs Proposed ======= ##########
modelname <- c("Linear", "Nonlinear","Interaction")
cratename <- c("Censoring: 20%", "Censoring: 50%")
ndataname <- c("N = 100", "N = 300", "N = 500")
mmodel <- 1:3
ccrate = c(1,2)
nndata = c(100,300,500)

cc = 1
setting = "nonPH"
mm = 3
xall = rep(0,6)
for (nn in 1:3){
  filename <- sprintf('./L2_%s_N%1.0%_m%1.0f_c%1.0f.rds',
                         setting,nndata[nn],mmodel[mm],ccrate[cc])
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
modelname <- c("Linear", "Nonlinear","Interaction")
cratename <- c("Censoring: 20%", "Censoring: 50%")
ndataname <- c("N = 100", "N = 300", "N = 500")
mmodel <- 1:3
ccrate = c(1,2)
nndata = c(100,300,500)

cc = 1

setting = "nonPH"

pdf_file_name <- sprintf("./L2_%s_c%1.0f_ROC.pdf",
                         setting,ccrate[cc])
pdf(pdf_file_name, width=17, height = 18)
par(mfrow=c(3,3))
for (nn in 1:3){
  for (mm in 1:3){ 
    aa = NULL
    
    filename <- sprintf('./L2_%s_N%1.0f_m%1.0f_c%1.0f.rds',
                        setting,nndata[nn],mmodel[mm],ccrate[cc])
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
         labels = c("Cox","LTRC CIF(P)","LTRC RRF(P)","LTRC TSF(P)",
                    "Cox","LTRC CIF(P)","LTRC RRF(P)","LTRC TSF(P)"),
         pos = 1, xpd = TRUE, cex = 1.6)
    
  }
}
dev.off()

######################### ========== IBS based CV rule =========== ######################
########### ======= The following code is to reproduce Table of Summary of IBS-based CV rules ======= ############
modelname <- c("Linear", "Nonlinear","Interaction")
cratename <- c("Censoring: 20%", "Censoring: 50%")
ndataname <- c("N = 100", "N = 300", "N = 500")
mmodel <- 1:3
ccrate = c(1,2)
nndata = c(100,300,500)


cc = 1
nn = 3
setting = "PH"

for (mm in 1:3){
  Qc = matrix(0, ncol = 7, nrow = 500)
  Qr = matrix(0, ncol = 7, nrow = 500)
  filename <- sprintf('./L2_%s_N%1.0f_m%1.0f_c%1.0f.rds',
                      setting,nndata[nn],mmodel[mm],ccrate[cc])
  
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
  if (mm == 1){
    cat(sprintf("\\multirow{9}{*}{%s}& \\multirow{3}{*}{%s}& %s & %s\\\\\n",setting,distname[dd],modelname[mm],Ire))
  } else if (mm == 2){
    cat(sprintf("&& %s & %s\\\\\n",modelname[mm],Ire))
  } else {
    cat(sprintf("&& %s & %s\\\\\n",modelname[mm],Ire))
  }
}




################## ========= Plot Pseudo-subject vs Subject (side by side, only the forests) ========= ##################################################
########### ======= The following code is to reproduce the plot of Pseudo-subject vs Subject ======= ######################
modelname <- c("Linear", "Nonlinear","Interaction")
cratename <- c("Censoring: 20%", "Censoring: 50%")
ndataname <- c("N = 100", "N = 300", "N = 500")
mmodel <- 1:3
ccrate = c(1,2)
nndata = c(100,300,500)
cc = 1

for (setting in c("PH","nonPH")){
  pdf_file_name <- sprintf("./subseu_%s_c%1.0f.pdf",
                           setting,ccrate[cc])
  pdf(pdf_file_name, width=15, height = 15)
  par(mfrow=c(3,3))
  aa = NULL
  for (nn in 1:3){
    for (mm in 1:3){ 
      filename <- sprintf('./L2_%s_N%1.0f_m%1.0f_c%1.0f.rds',
                          setting,nndata[nn],mmodel[mm],ccrate[cc])
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




###################### =========== Main effects plots ========= #################################
modelname <- c("Linear", "Nonlinear","Interaction")
cratename <- c("20%", "50%")
ndataname <- c("N = 100", "N = 300", "N = 500")
ssetting <- c("PH", "nonPH")
ssettingname <- c("PH", "non-PH")
SNRname <- c("Low", "High")
knowledgename <- c("Full", "Half")
scenname <- c("2TI + 1TV", "2TI + 4TV")
mmodel <- 1:3
ccrate = c(1,2)
nndata = c(100,300,500)

resaggre <- 
  lapply(1:2, function(ee){ #scenario
    if (ee == 1){ # "2TI + 1TV"
      varfile = "_vT"
    } else {
      varfile = ""
    }
    
    lapply(1:2, function(rr){
      
      if (rr == 1){ # Low
        snrfile = ""
      } else {
        snrfile =  "_sH"
      }
      lapply(1:2, function(ss){
        lapply(1:2, function(cc){
          lapply(1:3, function(nn){
            lapply(1:3, function(mm){
              filename <- sprintf('./L2%s%s_%s_N%1.0f_m%1.0f_c%1.0f.rds',
                                  snrfile, varfile, ssetting[ss], nndata[nn], mmodel[mm], ccrate[cc])
              
              resall <- readRDS(filename)
              ## KM, CIFPI, RRFPI, CIFPII, RRFPII
              Q <- as.matrix(resall$L2[, c(1,3,4,9,10)])
              Q = Q[rowSums(Q==Inf)==0, ]
              Q = Q[rowSums(Q==0)==0, ]
              ########## ======== presented way 1 ========= #######
              xrow <- colMeans((Q[, c(2,3,4,5)] - Q[, 1]) / Q[, 1])
              res = data.frame(matrix(0, ncol = 9, nrow = 2))
              names(res) <- c("LTRCCIF", "LTRCRRF",
                              "scenario","SNR","setting","model","crate","ntrain","knowledge")
              res$model <- modelname[mm]
              res$setting <- ssettingname[ss]
              res$ntrain <- nndata[nn]
              res$crate <- cratename[cc]
              res$scenario <- scenname[ee]
              res$SNR <- SNRname[rr]
              res$knowledge <- c("Full", "Half")
              res$LTRCCIF <- xrow[c(1,3)]
              res$LTRCRRF <- xrow[c(2,4)]
              print(sprintf("ee=%1.0f, rr = %1.0f, nn = %1.0f, ss = %1.0f, mm = %1.0f, nz = %1.0f", 
                            ee, rr, nn, ss, mm, sum(rowSums(Q==0)==0)))
              res
              
            }) %>% bind_rows()
          }) %>% bind_rows()
        }) %>% bind_rows()
      }) %>% bind_rows()
    }) %>% bind_rows()
  }) %>% bind_rows()


resallcombined <- resaggre %>% as_tibble() 
resallcombined$ntrain = factor(resallcombined$ntrain, levels = c("100", "300", "500"))
# resallcombined$scenario = factor(resallcombined$scenario, levels = c("2TI + 1TV", "2TI + 4TV"))
resallcombined$SNR = factor(resallcombined$SNR, levels = c("Low", "High"))
# resallcombined$knowledge = factor(resallcombined$knowledge, levels = c("Full", "Half"))
resallcombined$crate = factor(resallcombined$crate, levels = c("20%", "50%"))
resallcombined$setting = factor(resallcombined$setting, levels = c("PH", "non-PH"))
resallcombined$model = factor(resallcombined$model, levels = c("Linear", "Nonlinear","Interaction"))
######################### ========== main effect plot ========== #########################
create_tblplot <- function(dblfinal, method){
  
  tblplot <- lapply(method, function(met){
    # Formula <- as.formula(met ~ scenario + SNR + setting + model + crate + ntrain + knowledge)
    Formula <- as.formula(met ~ scenario + SNR + setting + crate + ntrain + knowledge)
    Formula <- as.formula(paste(c(met, "~", Formula[[3]]), collapse = " "))
    mod  <- lm(Formula, data = dblfinal)
    
    # mean of all
    MEall <- mean(dblfinal[[met]])
    
    MEdf <- as.data.frame(effects::allEffects(mod))
    namMEdf <- names(MEdf)
    # namMEdfnew <- c("Scenario", "SNR", "Setting", "Relationship", "Censoring rate", "Sample size", "Knowledge")
    namMEdfnew <- c("Scenario", "SNR", "Setting", "Censoring rate", "Sample size", "Knowledge")
    
    lapply(1:length(MEdf), function(mm){
      MEdf[[mm]][, 1:2] %>%
        mutate(cat = namMEdfnew[mm]) %>%
        rename("subcat" = namMEdf[mm],
               "Mean" = fit)
    }) %>%
      bind_rows() %>%
      select(cat, subcat, Mean) %>%
      as_tibble() %>%
      mutate(method = met) %>%
      mutate(meanall = MEall)
  }) %>%
    bind_rows()  %>%
    mutate(method = recode(method,
                           "LTRCCIF" = "LTRC CIF(P)",
                           "LTRCRRF" = "LTRC RRF(P)"))
  
  p1 <- tblplot %>%
    # mutate(cat = factor(cat, levels = c("Scenario", "Relationship", "Distribution", "SNR", 
    #                                     "Autocorrelation", "Censor rate", "Sample size"))) %>%
    mutate(subcat = factor(subcat, levels = unique(subcat))) %>% # this is important to made ntrain = 200, 1000, 5000 in order
    # arrange(match(subcat, c("200", "1000", "5000"))) %>%
    ggplot(aes(subcat, Mean)) +
    geom_point(group = 1, color = "steelblue") +
    geom_line(group = 1, color = "steelblue") +
    facet_grid(cols = vars(cat), 
               rows = vars(method), 
               scales = "free_x") +
    geom_hline(aes(yintercept = meanall), linetype = "dashed", color = "#999999") + 
    geom_hline(aes(yintercept = 0), linetype = "solid", color = "gray40") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x = element_blank(),
          strip.background = element_rect(fill="white")) 
  
  print(p1) # has to print, or do nothing 
  # dev.off()
}

dblfinal <- resallcombined %>% filter(model == "Interaction")
create_tblplot(dblfinal = dblfinal, 
               method = c("LTRCCIF", "LTRCRRF"))



