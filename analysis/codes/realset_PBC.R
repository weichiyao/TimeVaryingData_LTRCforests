## Make the PBC dataset 
DATA = make_realset(0)

#####################################################################################################################
Formula = Surv(Start, Stop, Event) ~ age+edema+alk.phos+chol+ast+platelet+spiders+hepato+ascites+albumin+bili+protime
Formula_TD = Surv(Start, Stop, Event, type = "counting") ~ age+edema+alk.phos+chol+ast+platelet+spiders+hepato+ascites+albumin+bili+protime

Nfold = 10

# time points of interest
Tpnt = c(0, sort(unique(DATA$Stop)))

BSt = NULL
#####################################################################################################################
Survobj = Surv(DATA$Start, DATA$Stop, DATA$Event)
## LTRCRRF
modelT <- ltrcrrf(formula = Formula, data = DATA, id = ID, stepFactor = 1.5)
predT <- predictProb(modelT, time.eval = Tpnt)
BSt$rrf <- sbrier_ltrc(obj = Survobj, id = DATA$ID, pred = predT, type = "BS")

## LTRCCIF
modelT <- ltrccif(formula = Formula, data = DATA, id = ID, stepFactor = 1.5)
predT <- predictProb(modelT, time.eval = Tpnt)
BSt$cf <- sbrier_ltrc(obj = Survobj, id = DATA$ID, pred = predT, type = "BS")

## Cox
modelT = coxph(formula = Formula, data = DATA)
predT = survfit(modelT, newdata = DATA)
Shat = shat(data = DATA, pred = predT, tpnt = Tpnt)
predTnew = list(survival.probs = Shat, 
                survival.times = Tpnt, 
                survival.tau = rep(max(Tpnt), length(unique(DATA$ID))))
BSt$cx = sbrier_ltrc(obj = Survobj, id = DATA$ID, pred = predTnew, type = "BS")
  
## LTRCTSF
modelT <- tsf_wy(formula = Formula_TD, data = DATA, stepFactor = 1.5)
predT <- predict_tsf_wy(object = modelT)

Shat = shat(DATA, pred = predT, tpnt = Tpnt)
predTnew = list(survival.probs = Shat, 
                survival.times = Tpnt, 
                survival.tau = rep(max(Tpnt), length(unique(DATA$ID))))
rm(Shat)
BSt$tsf = sbrier_ltrc(obj = Survobj, id = DATA$ID, pred = predTnew, type = "BS")

################ ======== Plot -- Using ggplot2 ========== ###############
library(ggplot2)
library(dplyr)
library(reshape2)
library(patchwork) # To display 2 charts together
##########################################################################
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
n = 4
cols = gg_color_hue(n)
cols = cols[c(2,1,3,4)]

##############################################################################
N = length(unique(DATA$ID))
id_uniq = sort(unique(DATA$ID))
data_sbrier <- data.frame(matrix(0, nrow = N, ncol = 3))
names(data_sbrier) <- c("ID", "times", "cens")
data_sbrier$ID <- id_uniq
data = DATA
for (ii in 1:N){
  data_sbrier[ii, ]$times <- max(data[data$ID==id_uniq[ii], ]$Stop)
  data_sbrier[ii, ]$cens <- sum(data[data$ID==id_uniq[ii], ]$Event)
}


tlen = length(Tpnt)
propAtRisk = rep(0, tlen)

for(i in 1:tlen){
  propAtRisk[i] = sum(data_sbrier$times > Tpnt[i]) / length(data_sbrier$times)
}


dataGP = data.frame(matrix(0, nrow = tlen, ncol = 6))
names(dataGP) = c("Time", "LTRC CIF", "LTRC RRF", "LTRC TSF", "Cox", "propAtRisk")
dataGP$Time = Tpnt
dataGP$`LTRC CF` = BSt$cf
dataGP$`LTRC RRF` = BSt$rrf
dataGP$`TSF` = BSt$tsf
dataGP$`Cox` = BSt$cx
dataGP$`propAtRisk` = propAtRisk


# Reshape for the ggplot2
melted = melt(dataGP, id.vars=c("Time","propAtRisk"),
              measure.vars = c("LTRC CIF", "LTRC RRF", "LTRC TSF", "Cox"))
names(melted)[3] = "Model"

# Value used to transform the data
coeff <- 0.22
# A few constants
ggplot(melted, aes(x=Time)) + 
  geom_step(aes(y = propAtRisk), color = "midnightblue", linetype = "solid", size = 0.7)
  geom_line(aes(y = value / coeff, color = Model, linetype = Model, group=Model), 
            size = 0.65) + 
  scale_linetype_manual(values=c("solid", "twodash", "F1", "dotdash"))+
  scale_color_manual(values=cols)+
  scale_y_continuous(
    # Features of the first axis
    name = "Proportion at risk (%)",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*coeff, name = "Brier score at time t")
  ) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))+
  #theme(legend.position = c(0.85, 0.85)) +
  theme(legend.title = element_text(size = 11, face = "bold"))+
  theme(legend.text = element_text(size = 9))



