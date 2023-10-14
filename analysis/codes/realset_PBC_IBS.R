# ###############################################
# setwd("./utils/randomForestSRC")
# devtools::load_all()
# 
# setwd("./utils/transformation/mlt_mod")
# devtools::load_all()
# 
# setwd("./utils/transformation/tram_mod")
# devtools::load_all()
# 
# setwd("./utils/transformation/trtf_mod")
# devtools::load_all()
#  
# ###############################################
# library(ipred) # for kfoldcv
# library(survival)
# library(LTRCforests)
# 
# source("./data/pbc_complete.R")
# source("./utils/tsf_tvary_funct.R")
# source("./utils/Loss_funct_tvary.R")
###############################################
# comp_scores <- function(jobname, seed, Nfold, ifold, ntree=100L){
#   filename <- sprintf('./%s_seed%1.0f_Nfold%1.0f_%1.0f.rds',
#                       jobname, seed, Nfold, ifold)
#   print(filename)
###############################################
## grab the array id value from the environment variable passed from sbatch

slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
#slurm_ncpus <- Sys.getenv("SLURM_CPUS_PER_TASK")
#typeof(slurm_arrayid )
#cat(slurm_arrayid)

# coerce the value to an integer
jj <- strtoi(slurm_arrayid)
setwd("/home/wy635/time-varyingRSF/wy_code/randomForestSRC")
# setwd("/Users/wyao/Dropbox/RESEARCH/time-varyingRSF/wy_code/randomForestSRC")
# create_package("rpart")
# clone the source code from the package repo, add your modifications, load the package with
# devtools::load_all() (must be inside the package folder) and test it
devtools::load_all()

setwd("/home/wy635/time-varyingRSF/wy_code/transformation/mlt_changed")
# setwd("/Users/wyao/Dropbox/RESEARCH/time-varyingRSF/transformation/mlt_changed")
# clone the source code from the package repo, add your modifications, load the package with
# devtools::load_all() (must be inside the package folder) and test it
devtools::load_all()

setwd("/home/wy635/time-varyingRSF/wy_code/transformation/tram_changed")
# setwd("/Users/wyao/Dropbox/RESEARCH/time-varyingRSF/transformation/tram_changed")
# clone the source code from the package repo, add your modifications, load the package with
# devtools::load_all() (must be inside the package folder) and test it
devtools::load_all()

setwd("/home/wy635/time-varyingRSF/wy_code/transformation/trtf_changed")
# setwd("/Users/wyao/Dropbox/RESEARCH/time-varyingRSF/transformation/trtf_changed")
# create_package("rpart")
# clone the source code from the package repo, add your modifications, load the package with
# devtools::load_all() (must be inside the package folder) and test it
devtools::load_all()



###############################################
library(ipred) # for kfoldcv
library(survival)
library(LTRCforests)

source("/home/wy635/time-varyingRSF/realset_PBC/pbc_complete.R")
source("/home/wy635/time-varyingRSF/realset_PBC/tsf_tvary_funct.R")
source("/home/wy635/time-varyingRSF/realset_PBC/Loss_funct_tvary.R")
###############################################
comp_scores <- function(seed, Nfold, ifold, ntree=100L){
  filename <- sprintf('/home/wy635/time-varyingRSF/realset_PBC/data-out/%s_seed%1.0f_Nfold%1.0f_%1.0f.rds',
                      jobname, seed, Nfold, ifold)
  print(filename)
  
  RES <- try(readRDS(filename), silent = T)
  
  if (class(RES) == "try-error"){
    ibsCVerr <- data.frame(matrix(0, nrow = 1, ncol = 4))
    names(ibsCVerr) <- c("cox", "cif", "rrf", "tsf")
    RES <- NULL
    RES$IBS <- ibsCVerr
    RES$BS <- NULL
  }  
  
  if (sum(RES$IBS[1,]==0)>0){
    ## Make the PBC dataset 
    fullDATA <- make_realset(0)
    
    #### Time points to be evaluated  
    Tpnt_max = c(0, sort(unique(fullDATA$Stop)))
    if (jobname == "ti"){
      fullDATA <- make_timeinvariant(fullDATA)
      Formula_TD = Surv(Stop, Event ) ~ age+edema+alk.phos+chol+ast+platelet+spiders+hepato+ascites+albumin+bili+protime
 
    } else{
      Formula_TD = Surv(Start, Stop, Event, type = "counting") ~ age+edema+alk.phos+chol+ast+platelet+spiders+hepato+ascites+albumin+bili+protime
 
    }
    #####################################################################################################################
    Formula = Surv(Start, Stop, Event) ~ age+edema+alk.phos+chol+ast+platelet+spiders+hepato+ascites+albumin+bili+protime
    
    
    set.seed(seed)
    # unique ID of the subject
    id_uniq = unique(fullDATA$ID)
    # Number of subjects
    N = length(id_uniq)
    
    kfold = kfoldcv(k = Nfold, N = N)
    
    dataCV = rep(list(0), Nfold)
    id_uniq_spd = sample(id_uniq, N)
    # Create the lists, each contains the data ID for the corresponding fold 
    for (kk in 1:Nfold){
      dataCV[[kk]] = id_uniq_spd[((kk - 1) * kfold[kk] + 1):(kk * kfold[kk])]
    }
    
    
    idxCV = which(fullDATA$ID %in% dataCV[[ifold]])
    # estimate on the fold jj
    newdata_b = fullDATA[idxCV, ]
    # use the rest of the folds to train
    data_b = fullDATA[-idxCV, ]
    
    Sobj_b = Surv(newdata_b$Start, newdata_b$Stop, newdata_b$Event)
    
    id_uniq_new = unique(newdata_b$ID)
    n_uniq_new = length(id_uniq_new)
    
     
    tau_b_max = sapply(1:n_uniq_new, function(i){max(Tpnt_max)})
    ## tsf
    if (RES$IBS[1,4]==0){
      set.seed(seed)
      
      modelT <- tsf_wy(formula = Formula_TD, data = data_b,
                       ntree = ntree, stepFactor = 1.5)
      predT <- predict_tsf_wy(object = modelT, newdata = newdata_b)
      rm(modelT)
      
      Shat = shat(newdata_b, pred = predT, tpnt = Tpnt) 
      predTnew_max = list(survival.probs = Shat, survival.times = Tpnt_max, survival.tau = tau_b_max)
      rm(Shat) 
      RES$IBS[1,4] <- sbrier_ltrc(obj = Sobj_b, pred = predTnew_max, id = newdata_b$ID)
      RES$BS$tsf   <- sbrier_ltrc(obj = Sobj_b, pred = predTnew_max, id = newdata_b$ID, type = "BS")
      
      saveRDS(RES,filename)
      
      rm(predTnew_max)
      gc()
    }
    print(sprintf("Fold%1.0f - TSF ibs=%1.3f", ifold, RES$IBS[1,4]))
    
    ## Cox model
    if (RES$IBS[1,1]==0){
      set.seed(seed)
      #### Training
      Coxfit <- coxph(formula = Formula, data = data_b)
      
      #### Prediction
      predT <- survival::survfit(Coxfit, newdata = newdata_b)
      rm(Coxfit)
      gc()
      Shat = shat(newdata_b, pred = predT, tpnt = Tpnt) 
      rm(predT)
      predTnew_max = list(survival.probs = Shat, survival.times = Tpnt_max, survival.tau = tau_b_max)
      rm(Shat)
      gc()
      
      RES$IBS[1,1] <- sbrier_ltrc(obj = Sobj_b, pred = predTnew_max, id = newdata_b$ID)
      RES$BS$cox   <- sbrier_ltrc(obj = Sobj_b, pred = predTnew_max, id = newdata_b$ID, type = "BS")
      
      saveRDS(RES,filename)
      rm(predTnew_max)
    }
    print(sprintf("Fold%1.0f - Cox ibs=%1.3f", ifold, RES$IBS[1,1]))
    
    ## cif
    if (RES$IBS[1,2]==0){
      set.seed(seed)
      modelT <- ltrccif(formula = Formula, data = data_b, id = ID, 
                        ntree = ntree, stepFactor = 1.5)
      
      predT_max <- predictProb(object = modelT, newdata = newdata_b, newdata.id = ID, 
                               time.eval = Tpnt_max, 
                               time.tau = tau_b_max)
      rm(modelT)
      
      RES$IBS[1,2] <- sbrier_ltrc(obj = Sobj_b, pred = predT_max, id = newdata_b$ID) 
      RES$BS$cif   <- sbrier_ltrc(obj = Sobj_b, pred = predT_max, id = newdata_b$ID, type = "BS")
      saveRDS(RES, filename)
      
      
      rm(predT_max)
      gc()
    }
    print(sprintf("Fold%1.0f - CIF ibs=%1.3f", ifold, RES$IBS[1,2]))
    
    ## rrf
    if (RES$IBS[1,3]==0){
      set.seed(seed)
      modelT <- ltrcrrf(formula = Formula, data = data_b, id = ID, 
                        ntree = ntree, stepFactor = 1.5)
      
      predT_max <- predictProb(object = modelT, newdata = newdata_b, newdata.id = ID, 
                               time.eval = Tpnt_max, 
                               time.tau = tau_b_max)
      rm(modelT)
      
      RES$IBS[1,3] <- sbrier_ltrc(obj = Sobj_b, pred = predT_max, id = newdata_b$ID) 
      RES$BS$rrf   <- sbrier_ltrc(obj = Sobj_b, pred = predT_max, id = newdata_b$ID, type = "BS")
      saveRDS(RES,filename)
      
      rm(predT_max)
      gc()
    }
    print(sprintf("Fold%1.0f - RRF ibs=%1.3f", ifold, RES$IBS[1,3]))
    
  }
  
  
}


Nfold=10
ntree=100L
seed=1021

comp_scores(jobname="ti", seed=seed, Nfold=Nfold, ifold=jj, ntree=ntree)

###########################################################################
library(ggplot2)
library(reshape2)
library(patchwork) # To display 2 charts together
library(tidyverse)
library(dplyr)
 

 

##############################################################################
filename = "/Users/wyao/Dropbox/RESEARCH/time-varyingRSF/report12_Denis/realset/tvti_all.rds"
RES <- readRDS(filename)
tlen = length(RES$bsCVerr_tv$t)
dataGP = data.frame(matrix(0, nrow = tlen, ncol = 10))
names(dataGP) = c("Time", "CIF-TV", "RRF-TV", "TSF-TV", "Extended Cox",
                  "CIF", "RRF", "TSF", "Cox", "propAtRisk")
dataGP$Time = RES$bsCVerr_tv$t

dataGP$`CIF-TV` = RES$bsCVerr_tv$cif 
dataGP$`RRF-TV` = RES$bsCVerr_tv$rrf
dataGP$`TSF-TV` = RES$bsCVerr_tv$tsf
dataGP$`Extended Cox` = RES$bsCVerr_tv$cox
dataGP$`CIF` = RES$bsCVerr_ti$cif 
dataGP$`RRF` = RES$bsCVerr_ti$rrf
dataGP$`TSF` = RES$bsCVerr_ti$tsf
dataGP$`Cox` = RES$bsCVerr_ti$cox

##########################################################################
#### ==== With the proportion at risk curve  ===== ####


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
n = 4
cols = gg_color_hue(n)
# Value used to transform the data
coeff <- 0.22
fullDATA <- make_realset(0)
dataGP$`propAtRisk` = compute_propAtRisk(DATA=fullDATA, Tpnt=dataGP$Time)
dataGP %>% 
  pivot_longer(-c(Time,propAtRisk), names_to = "Model") %>%
  mutate(
    Covariate = if_else(
      .data$Model %in% c("CIF-TV", "RRF-TV", "TSF-TV", "Extended Cox"), 'Time-varying', 'Time-invariant'
    )
  ) %>%
  mutate(
    Method = case_when(
      .data$Model %in% c("CIF-TV", "CIF") ~ "Conditional Inference",
     .data$Model %in% c("RRF-TV", "RRF") ~ "Relative Risk",
     .data$Model %in% c("TSF-TV", "TSF") ~ "Transformation",
     .data$Model %in% c("Extended Cox", "Cox") ~ "Cox Proportional Hazards"
    )
  ) %>%
  ggplot(aes(x=Time)) + 
  geom_step(aes(y = propAtRisk), color = "midnightblue", linetype = "dotdash", size = 0.7) +
  geom_line(aes(y = value / coeff, color = Method, linetype = Covariate), 
            size = 0.65) +  
  scale_color_manual(values=cols)+
  scale_linetype_manual(values=c("dotdash", "solid"))+
  scale_y_continuous(
    # Features of the first axis
    name = "Proportion at risk (%)",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*coeff, name = "Brier score-based cross-validation error")
  ) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))+
  theme_bw() + 
  theme(legend.title = element_text(size = 11, face = "bold"))+
  theme(legend.text = element_text(size = 9)) 

# 7.5 x 4.5
#### ==== Without the proportion at risk curve  ===== ####


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
n = 4
cols = gg_color_hue(n)
 
dataGP %>% 
  pivot_longer(-c(Time,propAtRisk), names_to = "Model") %>%
  mutate(
    Covariate = if_else(
      .data$Model %in% c("CIF-TV", "RRF-TV", "TSF-TV", "Extended Cox"), 'Time-varying', 'Time-invariant'
    )
  ) %>%
  mutate(
    Method = case_when(
      .data$Model %in% c("CIF-TV", "CIF") ~ "Conditional Inference",
      .data$Model %in% c("RRF-TV", "RRF") ~ "Relative Risk",
      .data$Model %in% c("TSF-TV", "TSF") ~ "Transformation",
      .data$Model %in% c("Extended Cox", "Cox") ~ "Cox Proportional Hazards"
    )
  ) %>%
  ggplot(aes(x=Time)) + 
  geom_line(aes(y = value  , color = Method, linetype = Covariate), 
            size = 0.65) +  
  scale_color_manual(values=cols)+
  scale_linetype_manual(values=c("dotdash", "solid"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))+
  theme_bw() + 
  theme(legend.title = element_text(size = 11, face = "bold"))+
  theme(legend.text = element_text(size = 9)) +
  ylab("Brier score based cross validation error") + 
  xlim(0,4800)
# 7 x 4
#### ==== Without the proportion at risk curve ===== ####
cols <- c("gray", "black")
dataGP %>% 
  pivot_longer(-c(Time,propAtRisk), names_to = "Model") %>%
  mutate(
    Covariate = if_else(
      .data$Model %in% c("CIF-TV", "RRF-TV", "TSF-TV", "Extended Cox"), 'Time-varying', 'Time-invariant'
    )
  ) %>%
  mutate(
    Method = case_when(
      .data$Model %in% c("CIF-TV", "CIF") ~ "Conditional Inference",
      .data$Model %in% c("RRF-TV", "RRF") ~ "Relative Risk",
      .data$Model %in% c("TSF-TV", "TSF") ~ "Transformation",
      .data$Model %in% c("Extended Cox", "Cox") ~ "Cox Proportional Hazards"
    )
  ) %>%
  ggplot(aes(x=Time)) + 
  geom_line(aes(y = value, color = Covariate, linetype = Method), 
            size = 0.65) +  
  scale_color_manual(values=cols)+
  scale_linetype_manual(values=c("dotted", "solid", "twodash", "dashed")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))+
  theme_bw() + 
  theme(legend.title = element_text(size = 11, face = "bold"))+
  theme(legend.text = element_text(size = 9)) +
  theme(legend.position = "top") + 
  # scale_x_continuous(limits = c(0, 4800)) +
  ylab("Brier score based cross validation error") + 
  xlim(0, 4800) 
  

# 11 x 6 