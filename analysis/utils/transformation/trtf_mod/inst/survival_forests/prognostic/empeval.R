start_time <- Sys.time()
library("parallel")
ncores <- 30

set.seed(290875)

source("competitors.R")
load("simulated_data.rda")

### number of repetitions
nsim <- length(simdat[[1]])

ret <- args[rep(1:nrow(args), each = nsim),]
simdat <- do.call("c", simdat)
ret$tf_W <- ret$tf_B <- ret$tf_W_alpha <- ret$tf_B_alpha <- 0
ret$cf <- ret$rsf <- ret$L1 <- ret$rg <- 0
ret$loglik <- unlist(do.call("c", loglik))

### filename to save estimated likelihoods
ret_filename <- 'results_empeval.rda'

### TSF Bs(theta)
trtf_fun <- function(i) {
    ## number of randomly preseleted variables to perform a split
    ## for low-dimensional data it is set to the number of variables
    ## for high-dimensional data it is set to the square root of the number of variables
    mtry <- ifelse(ret[i, "p"] <= 5, 
                   ncol(simdat[[i]]$train) - 1, 
                   ceiling(sqrt(ncol(simdat[[i]]$train) - 1)))
    
    ## learning and validation
    tf <- mytraforest(simdat[[i]]$train, simdat[[i]]$test, mtry = mtry)
    ## the learned forest is saved to be used in other methods
    save(tf, file=paste0(tf_dir, '/tf', i, '.rda'))
    unclass(tf$logLik_NN)
}
### the folder to save forest models
tf_folder <- 'results_tf'
dir.create(tf_folder)
tf_dir <- paste(getwd(), tf_folder, sep='/')

### learn and validate for each repetition of the simulated data
tf <- mclapply(1:length(simdat), trtf_fun, mc.cores = ncores)
ret$tf_B <- unlist(tf)
print('Traforest B is done!')
save(ret, file = ret_filename)

### DSF W(theta)
trtf_fun_Seibold <- function(i) {
    ## the corresponding TSF BS(theta) forest is loaded
    load(paste0(tf_dir, '/tf', i, '.rda'))
    ## learning and validation
    tf <- mytraforest_Seibold(tf, simdat[[i]]$train, simdat[[i]]$test)
    unclass(tf$logLik_NN)                        
}

### learn and validate for each repetition of the simulated data
tf <- mclapply(1:length(simdat), trtf_fun_Seibold, mc.cores = ncores)
ret$tf_W <- unlist(tf)
print('Traforest W is done!')
save(ret, file = ret_filename)

### TSF Bs(alpha)
trtf_fun_B <- function(i) {
    ## the corresponding TSF BS(theta) forest is loaded
    load(paste0(tf_dir, '/tf', i, '.rda'))
    ## learning and validation
    tf <- mytraforest_B(tf, simdat[[i]]$train, simdat[[i]]$test)
    unclass(tf$logLik_NN)
}

### learn and validate for each repetition of the simulated data
tf <- mclapply(1:length(simdat), trtf_fun_B, mc.cores = ncores)
ret$tf_B_alpha <- unlist(tf)
print('Traforest B_alpha is done!')
save(ret, file = ret_filename)

### DSF W(alpha)
trtf_fun_W <- function(i) {
    ## the corresponding TSF Bs(theta) forest is loaded
    load(paste0(tf_dir, '/tf', i, '.rda'))
    ## learning and validation
    tf <- mytraforest_W(tf, simdat[[i]]$train, simdat[[i]]$test)
    unclass(tf$logLik_NN)
}

### learn and validate for each repetition of the simulated data
tf <- mclapply(1:length(simdat), trtf_fun_W, mc.cores = ncores)
ret$tf_W_alpha <- unlist(tf)
print('Traforest W_alpha is done!')
save(ret, file = ret_filename)

### Ranger
library("ranger")
rg_fun <- function(i) {
    ## the corresponding TSF Bs(theta) forest is loaded
    load(paste0(tf_dir, '/tf', i, '.rda'))
    ## learning and validation
    rg <- myranger(tf, simdat[[i]]$train, simdat[[i]]$test)
    unclass(rg$logLik_NN)
}

### learn and validate for each repetition of the simulated data
rg <- mclapply(1:length(simdat), rg_fun, mc.cores = ncores)
ret$rg <- unlist(rg)
print('Ranger is done!')
save(ret, file = ret_filename)

### CForest
cf_fun <- function(i) {
    ## the corresponding TSF Bs(theta) forest is loaded
    load(paste0(tf_dir, '/tf', i, '.rda'))
    ## learning and validation
    cf <- mycforest(tf, simdat[[i]]$train, simdat[[i]]$test)
    unclass(cf$logLik_NN)
}

### learn and validate for each repetition of the simulated data
cf <- mclapply(1:length(simdat), cf_fun, mc.cores = ncores)
ret$cf <- unlist(cf)
print('Cforest is done!')
save(ret, file = ret_filename)

### RSF
library("randomForestSRC")
rsf_fun <- function(i) {
    ## the corresponding TSF Bs(theta) forest is loaded
    load(paste0(tf_dir, '/tf', i, '.rda'))
    ## learning and validation
    rsf <- myrfsrc(tf, simdat[[i]]$train, simdat[[i]]$test)
    unclass(rsf$logLik_NN)
}

### learn and validate for each repetition of the simulated data
rsf <- mclapply(1:length(simdat), rsf_fun, mc.cores = ncores)
ret$rsf <- unlist(rsf)
print('RSF is done!')

### L1 is a modified randomForestSRC, therefore we
### detach it to prevent the conflicts
detach(package:randomForestSRC)
save(ret, file = ret_filename)

### L1-splitting forest
### the library is not freely available
if {FALSE}
library("L1SRC")

L1_fun <- function(i) {
    ## the corresponding TSF Bs(theta) forest is loaded
    load(paste0(tf_dir, '/tf', i, '.rda'))
    ## learning and validation
    rsf <- myL1(tf, simdat[[i]]$train, simdat[[i]]$test)
    unclass(rsf$logLik_NN)
}

### learn and validate for each repetition of the simulated data
L1 <- mclapply(1:length(simdat), L1_fun, mc.cores = ncores)
ret$L1 <- unlist(L1)
print('L1 is done!')
}

### for each repetition of the simulated data
### the parameters of the corresponding model,
### the true model likelihood and
### the estimated likelihoods of all forest models are saved
save(ret, file = ret_filename)

end_time <- Sys.time()
print(end_time - start_time)

sessionInfo()
