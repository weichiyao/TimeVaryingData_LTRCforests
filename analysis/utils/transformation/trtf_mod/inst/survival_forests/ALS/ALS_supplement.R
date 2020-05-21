## ----libraries ,message=FALSE, size='scriptsize', cache=TRUE, cache.lazy=FALSE----
## load ALS data
## The data is available to registered users of https://nctu.partners.org/ProACT
## Register, download the data, and run
## R> library("TH.data"); demo("PROACT");
load("./ALS/ALSsurvdata.rda")
## attach packages
library("survival")
library("tram")
library("party")
library("trtf")
## to make results reproducible
set.seed(290875)

## ----data_preprocessing, size='scriptsize', cache=TRUE, cache.lazy=FALSE----
## remove items with missing values
ALScc <- ALSsurvdata[complete.cases(ALSsurvdata),]
## convert time variable from days to years
ALScc$time_onset_treatment <- ALScc$time_onset_treatment / 365
## too few levels for Coxph and Survreg
levels(ALScc$race) <- c("Caucasian", rep("Other", 3))

## stratification control for predictive models
grps <- with(ALScc, interaction(cens, Riluzole))

## output as Surv object
ALScc$y <- with(ALScc, Surv(survival.time, cens))
ALScc$cens <- ALScc$survival.time <- NULL
dfprog <- dfpred <- ALScc
dfprog$Riluzole <- NULL
## transform Riluzole variable to 0-1 variable
dfpred$Riluzole <- as.numeric(dfpred$Riluzole == 'Yes')

## ----settings, size='scriptsize', cahce=TRUE, cache.lazy=FALSE-----------
## number of randomly preselected variables for splitting
## set equal to the square root of the total number of variables
mtry <- ceiling(sqrt(ncol(ALScc) - 2))

### joint parameters for all forest procecdures
## number of trees
ntree <-250
## maximal depth of trees
nodedepth <- 10
## number of obs in a terminal node
minbucket <- 20

### Bernstein approximation parameters
## order of polynomials
ord <- 5
## indicator of whether log-transformation is applied
## to the time variable before the approximation
logf <- TRUE
## the largest time values are not used for approximation
## to keep it stable
sup <- c(sqrt(.Machine$double.eps), quantile(ALScc$y[,1], .95))

## set to integer value if what to learn forests with
## parallel computation
NCORES <- NULL
## Warning: when non-null value is used for NCORES it may
## offen cause issues. Non-parallel learning is stable.

### parameters for party/partykit
## with testype = "Univ" and mincriterion = 0 test-splitting
## procedure is free from early stops caused by p-values
ctrl_partykit <- partykit:::ctree_control(teststat = "quad",
                                          testtype = "Univ",
                                          mincriterion = 0,
                                          minbucket = minbucket,
                                          maxdepth = nodedepth)

## time points at which survival function is estimated
stime <- seq(from = 0, to = max(ALScc$y[,1]), length.out = 100)

## ----prognostic_cox_weibull_models, cache=TRUE, cache.lazy=FALSE, size='scriptsize'----
method <- c("Cox", "Weibull_alpha", "Cox_alpha", 
            "W_alpha", "Bs_alpha", "W_theta", "Bs_theta")
method <- factor(method, levels = method, labels = method)
distr <- expand.grid(time = stime, method = method)

## Cox PH model with constant PH
coxconst <- Coxph(y ~ 1, dfprog, order = ord, log_first = logf, support = sup)
logLik(coxconst)
distr$surv[distr$method == "Cox"] <- predict(coxconst, 
                                             newdata = data.frame(1), 
                                             q = stime, 
                                             type = "survivor")
## Cox PH model with prognostic PH
coxph <- Coxph(y ~ ., dfprog, order = ord, log_first = logf, support = sup)
logLik(coxph)
distr$surv[distr$method == "Cox_alpha"] <- predict(coxph, 
                                                   newdata = dfprog[1,], 
                                                   q = stime, 
                                                   type = "survivor")
## Weibull model
srg <- Survreg(y ~ ., dfprog, order = ord)
logLik(srg)
distr$surv[distr$method == "Weibull_alpha"] <- predict(srg, 
                                                       newdata = dfprog[1,], 
                                                       q = stime, 
                                                       type = "survivor")

## ----Bs_theta, cache=TRUE, cache.lazy=FALSE,size='scriptsize'------------
tf_B  <- traforest(as.mlt(Coxph(y ~ 1, 
                                data = dfprog, 
                                order = ord, 
                                log_first = logf, 
                                sup = sup)),
                   formula = y ~ ., 
                   data = dfprog,
                   ntree = ntree, 
                   mtry = mtry,
                   control = ctrl_partykit, 
                   cores = NCORES)
logLik(tf_B, OOB = TRUE)

distr$surv[distr$method == "Bs_theta"] <- predict(tf_B, 
                                                  dfprog[1,], 
                                                  type = 'survivor', 
                                                  q = stime)[[1]]
## The weights of the learned forest are used for learning other
## prognostic forests
initw <- do.call("cbind", tf_B$weights)
storage.mode(initw) <- "double"

## ----W_theta, cache=TRUE, cache.lazy=FALSE,size='scriptsize'-------------
tf_W <- traforest(as.mlt(Survreg(y ~ 1, data = dfprog)),
                  formula = y ~ ., 
                  data = dfprog,
                  ntree = ntree, 
                  mtry = mtry,
                  control = ctrl_partykit, 
                  weights = initw,
                  cores = NCORES)
logLik(tf_W, OOB = TRUE)

distr$surv[distr$method == "W_theta"] <- predict(tf_W, 
                                                 dfprog[1,], 
                                                 type = 'survivor', 
                                                 q = stime)[[1]]

## ----W_alpha, cache=TRUE, cache.lazy=FALSE,size='scriptsize'-------------
tf_W_alpha <- traforest(as.mlt(Survreg(y ~ 1, data = dfprog)),
                        formula = y ~ ., 
                        data = dfprog, 
                        parm = "(Intercept)",
                        ntree = ntree, 
                        mtry = mtry,
                        control = ctrl_partykit, 
                        weights = initw,
                        cores = NCORES)
logLik(tf_W_alpha, OOB = TRUE)

distr$surv[distr$method == "W_alpha"] <- predict(tf_W_alpha, 
                                                 dfprog[1,], 
                                                 type = 'survivor', 
                                                 q = stime)[[1]]

## ----Bs_alpha, cache=TRUE, cache.lazy=FALSE,size='scriptsize'------------
dfprog$alpha <- 1
tf_B_alpha <- traforest(m <- as.mlt(Coxph(y ~ alpha, 
                                          data = dfprog, 
                                          order = ord, 
                                          log_first = logf,
                                          support = sup, 
                                          fixed = c("alpha" = 0))),
                        formula = y | alpha ~ ., 
                        data = dfprog,
                        parm = "alpha",
                        mltargs = list(fixed = c("alpha" = 0)),
                        ntree = ntree, 
                        mtry = mtry,
                        control = ctrl_partykit,
                        weights = initw,
                        cores = NCORES)
logLik(tf_B_alpha, OOB = TRUE)

## a hack to make predictions
cf1 <- predict(tf_B_alpha, dfprog[1,], type = 'coef')[[1]]
cf <- coef(m, fixed = TRUE)
cf[names(cf1)] <- cf1
coef(m) <- cf
distr$surv[distr$method == "Bs_alpha"] <- predict(m, 
                                                  newdata = data.frame(alpha = 0), 
                                                  type = 'survivor', 
                                                  q = stime)

## ----prognostic_plot, fig.width=8, fig.height=6, size='scriptsize', cache=TRUE, cache.lazy=FALSE----
library("lattice")
xyplot(surv ~ time, data = distr, group = method, type = "l", auto.key = TRUE)

## ----predictive_cox_weibull_models, cache=TRUE, cache.lazy=FALSE, size='scriptsize'----
## a formula that states time-independent part of log-hazard as a linear combination
## of variables, treatmenent and treatment-variable interactionsn
fm <- names(dfpred)
fm <- paste(fm[!(fm %in% c("y", "Riluzole"))], collapse = "+")
fmy <- as.formula(paste("y ~ Riluzole +", fm, " + Riluzole * (", fm, ")"))

## grid methods vs time points
method <- c("Cox_alpha_beta", "Weibull_alpha_beta", 
            "Bs_theta_theta", "W_theta_theta", 
            "Bs_theta_beta", "W_theta_beta",
            "Bs_alpha_beta", "W_alpha_beta")
method <- factor(method, levels = method, labels = method)
distr <- expand.grid(time = stime, Riluzole = 0:1, method = method)

## for the first patient from ALS dataset prediction is made in
## the abcence and in the presence of treatment by Riluzole
nd <- dfpred[c(1, 1),]
nd$Riluzole <- 0:1

## predictive Cox PH model
tcoxph <- Coxph(fmy , dfpred, order = ord, log_first = logf, support = sup)
logLik(tcoxph)
distr$surv[distr$method == "Cox_alpha_beta"] <-  c(predict(tcoxph,
                                                           newdata = nd, 
                                                           q = stime, 
                                                           type = "survivor"))

## predictive Weibull PH model
tweib <- Survreg(fmy, dfpred)
logLik(tweib)
distr$surv[distr$method == "Weibull_alpha_beta"] <- c(predict(tweib,
                                                              newdata = nd, 
                                                              q = stime, 
                                                              type = "survivor"))

## ----Bs_theta_theta, cache=TRUE, cache.lazy=FALSE, size='scriptsize'-----
Bs_theta_theta <- traforest(as.mlt(Coxph(y | Riluzole ~ 1, 
                                         data = dfpred, 
                                         order = ord,
                                         log_first = logf, 
                                         support = sup)),
                            formula = y | Riluzole ~ ., 
                            data = dfpred,
                            mtry = mtry, 
                            ntree = ntree,
                            strata = grps, 
                            control = ctrl_partykit,
                            cores = NCORES)
logLik(Bs_theta_theta, OOB = TRUE)
distr$surv[distr$method == "Bs_theta_theta"] <- c(predict(Bs_theta_theta,
    newdata = dfpred[1,], mnewdata = data.frame(Riluzole = 0:1),
    type = 'survivor', q = stime)[[1]])

initw <- do.call("cbind", Bs_theta_theta$weights)
storage.mode(initw) <- "double"

## ----W_theta_theta, cache=TRUE, cache.lazy=FALSE, size='scriptsize'------
W_theta_theta <- traforest(as.mlt(Survreg(y | Riluzole ~ 1, data = dfpred)),
                           formula = y | Riluzole ~ ., 
                           data = dfpred,
                           mtry = mtry, 
                           ntree = ntree,
                           strata = grps, 
                           control = ctrl_partykit,
                           weights = initw, 
                           cores = NCORES)
logLik(W_theta_theta, OOB = TRUE)
distr$surv[distr$method == "W_theta_theta"] <- c(predict(W_theta_theta,
                                                         newdata = dfpred[1,], 
                                                         mnewdata = data.frame(Riluzole = 0:1),
                                                         type = 'survivor', 
                                                         q = stime)[[1]])

## ----Bs_theta_beta, cache=TRUE, cache.lazy=FALSE, size='scriptsize'------
Bs_theta_beta <- traforest(as.mlt(Coxph(y ~ Riluzole, 
                                        data = dfpred, 
                                        order = ord,
                                        log_first = logf, 
                                        support = sup)),
                           formula = y | Riluzole ~ ., 
                           data = dfpred,
                           mtry = mtry, 
                           ntree = ntree,
                           strata = grps, 
                           control = ctrl_partykit,
                           weights = initw, 
                           cores = NCORES)
logLik(Bs_theta_beta, OOB = TRUE)
distr$surv[distr$method == "Bs_theta_beta"] <- c(predict(Bs_theta_beta,
                                                         newdata = dfpred[1,], 
                                                         mnewdata = data.frame(Riluzole = 0:1),
                                                         type = 'survivor', 
                                                         q = stime)[[1]])

## ----W_theta_beta, cache=TRUE, cache.lazy=FALSE, size='scriptsize'-------
W_theta_beta <- traforest(as.mlt(Survreg(y ~ Riluzole, data = dfpred)),
                          formula = y | Riluzole ~ ., 
                          data = dfpred,
                          mtry = mtry, 
                          ntree = ntree,
                          strata = grps, 
                          control = ctrl_partykit,
                          weights = initw, 
                          cores = NCORES)
logLik(W_theta_beta, OOB = TRUE)
distr$surv[distr$method == "W_theta_beta"] <- c(predict(W_theta_beta,
                                                        newdata = dfpred[1,], 
                                                        mnewdata = data.frame(Riluzole = 0:1),
                                                        type = 'survivor', 
                                                        q = stime)[[1]])

## ----W_alpha_beta, cache=TRUE, cache.lazy=FALSE, size='scriptsize'-------
W_alpha_beta <- traforest(as.mlt(Survreg(y ~ Riluzole, data = dfpred)),
                          formula = y | Riluzole ~ ., 
                          data = dfpred,
                          parm = c("(Intercept)", "Riluzole"),
                          mtry = mtry, 
                          ntree = ntree,
                          strata = grps, 
                          control = ctrl_partykit,
                          weights = initw, 
                          cores = NCORES)
logLik(W_alpha_beta, OOB = TRUE)
distr$surv[distr$method == "W_alpha_beta"] <- c(predict(W_alpha_beta,
                                                        newdata = dfpred[1,], 
                                                        mnewdata = data.frame(Riluzole = 0:1),
                                                        type = 'survivor', 
                                                        q = stime)[[1]])

## ----Bs_alpha_beta, cache=TRUE, cache.lazy=FALSE, size='scriptsize'------
dfpred$int <- 1
Bs_alpha_beta <- traforest(m <- as.mlt(Coxph(y ~ int + Riluzole, 
                                             data = dfpred, 
                                             fixed = c("int" = 0), 
                                             order = ord,
                                             log_first = logf, 
                                             support = sup)),
                           formula = y | int + Riluzole ~ .,
                           data = dfpred, 
                           parm = c("int", "Riluzole"),
                           mltargs = list(fixed = c("int" = 0)), 
                           mtry = mtry,
                           ntree = ntree, 
                           strata = grps,
                           control = ctrl_partykit, 
                           weigths = initw,
                           cores = NCORES)
logLik(Bs_alpha_beta, OOB = TRUE)
cf1 <- predict(Bs_alpha_beta, dfprog[1,], type = 'coef')[[1]]
cf <- coef(m, fixed = TRUE)
cf[names(cf1)] <- cf1
coef(m) <- cf
distr$surv[distr$method == "Bs_alpha_beta"] <- c(predict(m, 
                                                         newdata = data.frame(int = 0, 
                                                                              Riluzole = 0:1), 
                                                         type = 'survivor', 
                                                         q = stime))

## ----predictive_plot, fig.width = 8, fig.height = 6, size='scriptsize', cache=TRUE, cache.lazy=FALSE----
library("lattice")
xyplot(surv ~ time | method, data = distr, group = Riluzole, type = "l", auto.key = TRUE)

