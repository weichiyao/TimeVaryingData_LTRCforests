
options(warn = -1L)

###################################################
### code chunk number 1: setup
###################################################

library("trtf")
library("quantregForest")

nmax <- Inf

pdf("applications.pdf", width = 12, height = 8)

data("abalone", package = "AppliedPredictiveModeling")
response <- "Rings"
abalone[[response]] <- as.numeric(abalone[[response]])

ns <- 100
fm <- Rings ~ Type + LongestShell + Diameter + Height + WholeWeight + ShuckedWeight + 
        VisceraWeight + ShellWeight
mtry <- ceiling(length(all.vars(fm[[3]])) / 3)
var_m <- numeric_var("Rings", support = quantile(abalone[[response]], prob = c(.1, .9)), 
                     add = range(abalone[[response]]) - quantile(abalone[[response]], prob = c(.1, .9)), 
                     bounds = c(0, Inf))

B_m <- Bernstein_basis(var_m, order = 4, ui = "increasing")
uc_ctm_AB <- ctm(B_m, data = abalone, todistr = "Normal")
uc_mlt_AB <- mlt(uc_ctm_AB, data = abalone, scale = FALSE)

c_ctm_AB <- ctm(B_m, data = abalone, todistr = "Normal", shift = fm[c(1, 3)])
c_mlt_AB <- mlt(c_ctm_AB, data = abalone, scale = TRUE, maxit = 2500)

tt_AB <- trafotree(uc_ctm_AB, formula = fm, data = abalone, 
                control = ctree_control(mincriterion = .95, minsplit = ns*2, minbucket = ns, nmax = nmax),
                mltargs = list(maxit = 10000, scale = TRUE, gtol = 1e-3, trace = FALSE))

if (FALSE) {
tf_AB <- traforest(uc_ctm_AB, formula = fm, data = abalone, ntree = 100, trace = TRUE, 
                  mltargs = list(maxit = 10000, scale = TRUE, gtol = 1e-3, trace = FALSE),
                  control = ctree_control(mincriterion = 0, minsplit = ns*2, minbucket = ns, nmax = nmax),
                                          mtry = mtry)

qrf_AB <- quantregForest(x = abalone[, all.vars(fm[[3]])], y = abalone[, all.vars(fm[[2]])], nodesize = ns,
                                       mtry = mtry, ntree = 100, keep.inbag = TRUE)
}

data("BostonHousing2", package = "mlbench")
response <- "cmedv"
BostonHousing2[[response]] <- as.numeric(BostonHousing2[[response]])

ns <- 40
fm <- cmedv ~ crim + zn + indus + chas + nox + rm + age + dis + rad + tax + ptratio + b + lstat
mtry <- ceiling(length(all.vars(fm[[3]])) / 3)
var_m <- numeric_var("cmedv", support = quantile(BostonHousing2[[response]], prob = c(.1, .9)), 
                     add = range(BostonHousing2[[response]]) - quantile(BostonHousing2[[response]], prob = c(.1, .9)), 
                     bounds = c(0, Inf))

B_m <- Bernstein_basis(var_m, order = 4, ui = "increasing")
uc_ctm_BH <- ctm(B_m, data = BostonHousing2, todistr = "Normal")
uc_mlt_BH <- mlt(uc_ctm_BH, data = BostonHousing2, scale = FALSE, trace = FALSE)

c_ctm_BH <- ctm(B_m, data = BostonHousing2, todistr = "Normal", shift = fm[c(1, 3)])
c_mlt_BH <- mlt(c_ctm_BH, data = BostonHousing2, scale = TRUE, trace = FALSE)

tt_BH <- trafotree(uc_ctm_BH, formula = fm, data = BostonHousing2, 
                control = ctree_control(mincriterion = .95, minsplit = 2*ns, minbucket = ns, nmax = nmax),
                mltargs = list(maxit = 10000, scale = TRUE, gtol = 1e-3, trace = FALSE))

if (FALSE) {
tf_BH <- traforest(uc_ctm_BH, formula = fm, data = BostonHousing2, ntree = 100, 
                  control = ctree_control(mincriterion = 0, minsplit = 2*ns, minbucket = ns,
                                          nmax = nmax), mtry = mtry, 
                  trace = TRUE, mltargs = list(maxit = 10000, scale = TRUE, gtol = 1e-3, trace = FALSE))

qrf_BH <- quantregForest(x = BostonHousing2[, all.vars(fm[[3]])], y = BostonHousing2[, all.vars(fm[[2]])], 
                                       nodesize = ns, mtry = mtry, ntree = 100, keep.inbag = TRUE)
}

data("BigMac2003", package = "alr3")
response <- "BigMac"
BigMac2003[[response]] <- as.numeric(BigMac2003[[response]])

ns <- 20
fm <- BigMac ~ Bread + Rice + FoodIndex + Bus + Apt + TeachGI + 
        TeachNI + TaxRate + TeachHours
mtry <- ceiling(length(all.vars(fm[[3]])) / 3)
var_m <- numeric_var("BigMac", support = quantile(BigMac2003[[response]], prob = c(.1, .9)),
                     add = range(BigMac2003[[response]]) - quantile(BigMac2003[[response]], prob = c(.1, .9)), 
                     bounds = c(0, Inf))

B_m <- Bernstein_basis(var_m, order = 4, ui = "increasing")
uc_ctm_BM <- ctm(B_m, data = BigMac2003, todistr = "Normal")
uc_mlt_BM <- mlt(uc_ctm_BM, data = BigMac2003, scale = FALSE, trace = FALSE)

c_ctm_BM <- ctm(B_m, data = BigMac2003, todistr = "Normal", shift = fm[c(1, 3)])
c_mlt_BM <- mlt(c_ctm_BM, data = BigMac2003, scale = TRUE, trace = FALSE)

tt_BM <- trafotree(uc_ctm_BM, formula = fm, data = BigMac2003, 
                control = ctree_control(mincriterion = .95, minsplit = 2 * ns, minbucket = ns, nmax = nmax),
                mltargs = list(maxit = 10000, scale = TRUE, gtol = 1e-3, trace = FALSE))

if (FALSE) {
tf_BM <- traforest(uc_ctm_BM, formula = fm, data = BigMac2003, ntree = 100, 
                  control = ctree_control(mincriterion = 0, minsplit = 2*ns, minbucket = ns, 
                  nmax = nmax), mltargs = list(maxit = 10000, scale = TRUE, gtol = 1e-3, trace = FALSE),
                  trace = TRUE, mtry = mtry)

qrf_BM <- quantregForest(x = BigMac2003[, all.vars(fm[[3]])], 
    y = BigMac2003[, all.vars(fm[[2]])], nodesize = ns, mtry = mtry, ntree = 100, keep.inbag = TRUE)
}

data("Ozone", package = "mlbench")
Ozone <- subset(Ozone, complete.cases(Ozone))
Ozone <- as.data.frame(lapply(Ozone, function(x) {
    x <- x[, drop = TRUE]
    if (is.factor(x)) return(as.ordered(x))
    x
}))
response <- "V4"
Ozone[[response]] <- as.numeric(Ozone[[response]])

ns <- 20
fm <- V4 ~ V1 + V2 + V3 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 + V13
mtry <- ceiling(length(all.vars(fm[[3]])) / 3)
var_m <- numeric_var("V4", support = quantile(Ozone[[response]], prob = c(.1, .9)), 
                     add = range(Ozone[[response]]) - quantile(Ozone[[response]], prob = c(.1, .9)), 
                     bounds = c(0, Inf))

B_m <- Bernstein_basis(var_m, order = 4, ui = "increasing")
uc_ctm_Ozone <- ctm(B_m, data = Ozone, todistr = "Normal")
uc_mlt_Ozone <- mlt(uc_ctm_Ozone, data = Ozone, scale = FALSE)

c_ctm_Ozone <- ctm(B_m, data = Ozone, todistr = "Normal", shift = fm[c(1, 3)])
c_mlt_Ozone <- mlt(c_ctm_Ozone, data = Ozone, scale = TRUE, maxit = 10000)

tt_Ozone <- trafotree(uc_ctm_Ozone, formula = fm, data = Ozone, 
                control = ctree_control(mincriterion = .95, minsplit = 2*ns, minbucket = ns, nmax = nmax),
                mltargs = list(maxit = 10000, scale = TRUE, gtol = 1e-3, trace = FALSE))

tf_Ozone <- traforest(uc_ctm_Ozone, formula = fm, data = Ozone, ntree = 100, 
                  control = ctree_control(mincriterion = 0, minsplit = 2*ns,
                  minbucket = ns, nmax = nmax), trace = TRUE, mtry = mtry, 
                  mltargs = list(maxit = 10000, scale = TRUE, gtol = 1e-3, trace = FALSE))

qrf_Ozone <- quantregForest(x = Ozone[, all.vars(fm[[3]])], y = Ozone[, all.vars(fm[[2]])], 
                                          nodesize = ns, mtry = mtry, ntree = 100, keep.inbag = TRUE)


library("lattice")
lwd <- 1.5
col <- rgb(.1, .1, .1, .05)
colR <- rgb(.75, 0, 0, .8)
colB <- rgb(0, 0, 0.5, .8)
trellis.par.set(list(plot.symbol = list(col="blue",pch=18, cex=0.5),
                     box.rectangle = list(col=1),
                     box.umbrella = list(lty=1, col=1),
                     strip.background = list(col = "white")))
ltheme <- canonical.theme(color = FALSE)     ## in-built B&W theme
ltheme$strip.background$col <- "transparent" ## change strip bg
lattice.options(default.theme = ltheme)
bpch <- 18
bfill <- FALSE

###################################################
### code chunk number 4: Ozone-uctm
###################################################
q <- mkgrid(uc_mlt_Ozone, n = 100)[[1]]
d <- predict(uc_mlt_Ozone, newdata = Ozone, q = q, type="density")
plot(d ~ q, type = "l", ylab = "Density", xlab = "Daily Maximum One-hour-average Ozone Reading")
rug(Ozone$V4, col = rgb(.1, .1, .1, .1), lwd = 2)


###################################################
### code chunk number 5: Ozone-ltm
###################################################
q <- mkgrid(uc_mlt_Ozone, n = 100)[[1]]

d <- predict(c_mlt_Ozone, newdata = Ozone, q = q, type="distribution")
lp <- c(predict(c_mlt_Ozone, newdata = Ozone, q = 0, terms = "bshift"))
nd <- expand.grid(q = q, lp = -lp)
nd$d <- c(d)
pfun <- function(x, y, z, subscripts, at, ...) {
    panel.contourplot(x, y, z, subscripts, at = 1:9/10, ...)
    panel.xyplot(x = -lp, y = Ozone$V4, pch = 19, cex = 1.5,
                 col = rgb(.1, .1, .1, .1), ...)
}
print(contourplot(d ~ lp + q, data = nd, panel = pfun, xlab = "Linear Predictor", 
                  ylab = "Observed", ylim = c(0, 40)))


###################################################
### code chunk number 6: Ozone-trtree
###################################################
plot(tt_Ozone, tp_args = list(type = "density", id = FALSE, ylines = 0, K = 100))


###################################################
### code chunk number 8: Ozone-forest
###################################################
nx <- lapply(Ozone, function(x) {
    if (is.factor(x))
        return(factor(levels(x)[which.max(table(x))], levels = levels(x), ordered = is.ordered(x)))
    return(median(x))
})

nx$V8 <- 20:100
nx <- expand.grid(nx)

q1 <- predict(qrf_Ozone, newdata = nx, what = c(.2, .5, .8))

q2 <- predict(tf_Ozone, newdata = nx, K = 100, type = "quantile", p = c(.2, .5,.8))
q2 <- do.call("rbind", q2)

nx <- rbind(cbind(nx["V8"], q = q1[,1], p = "0.2", model = "QRF", stringsAsFactors = FALSE),
            cbind(nx["V8"], q = q1[,2], p = "0.5", model = "QRF", stringsAsFactors = FALSE),
            cbind(nx["V8"], q = q1[,3], p = "0.8", model = "QRF", stringsAsFactors = FALSE),
            cbind(nx["V8"], q = q2[,1], p = "0.2", model = "Forest", stringsAsFactors = FALSE),
            cbind(nx["V8"], q = q2[,2], p = "0.5", model = "Forest", stringsAsFactors = FALSE),
            cbind(nx["V8"], q = q2[,3], p = "0.8", model = "Forest", stringsAsFactors = FALSE))
nx$model <- factor(as.character(nx$model), levels = c("QRF", "Forest"), labels = c("QRF", "Forest"))

print(xyplot(q ~ V8 | model, data = nx, groups = p, type = "l", lty = 1:3, col = "black",
             xlab = "Sandburg, CA, Temperature in Fahrenheit",
             ylab = "Daily Maximum One-hour-average Ozone Reading"))

varimp(tf_Ozone)

###################################################
### code chunk number 2: abalone-ltm
###################################################
q <- mkgrid(uc_mlt_AB, n = 100)[[1]]

d <- predict(c_mlt_AB, newdata = abalone, q = q, type="distribution")
lp <- c(predict(c_mlt_AB, newdata = abalone, q = 0, terms = "bshift"))
nd <- expand.grid(q = q, lp = -lp)
nd$d <- c(d)
pfun <- function(x, y, z, subscripts, at, ...) {
    panel.contourplot(x, y, z, subscripts, at = 1:9/10, ...)
    panel.xyplot(x = -lp, y = abalone$Rings, pch = 19, cex = 1.5,
                 col = rgb(.1, .1, .1, .1), ...)
}
print(contourplot(d ~ lp + q, data = nd, panel = pfun, xlab = "Linear Predictor",
                  ylab = "Observed"))


###################################################
### code chunk number 3: abalone-trtree
###################################################
plot(tt_AB, tp_args = list(type = "density", id = FALSE, ylines = 0, K = 100))


###################################################
### code chunk number 6: BigMac-ltm
###################################################
q <- mkgrid(uc_mlt_BM, n = 100)[[1]]

d <- predict(c_mlt_BM, newdata = BigMac2003, q = q, type="distribution")
lp <- c(predict(c_mlt_BM, newdata = BigMac2003, q = 0, terms = "bshift"))
nd <- expand.grid(q = q, lp = -lp)
nd$d <- c(d)
pfun <- function(x, y, z, subscripts, at, ...) {
    panel.contourplot(x, y, z, subscripts, at = 1:9/10, ...)
    panel.xyplot(x = -lp, y = BigMac2003$BigMac, pch = 19, cex = 1.5,
                 col = rgb(.1, .1, .1, .1), ...)
}
print(contourplot(d ~ lp + q, data = nd, panel = pfun, xlab = "Linear Predictor",
                  ylab = "Observed"))


###################################################
### code chunk number 7: BigMac-trtree
###################################################
plot(tt_BM, tp_args = list(type = "density", id = FALSE, ylines = 0, K = 100))


###################################################
### code chunk number 10: BH-ltm
###################################################
q <- mkgrid(uc_mlt_BH, n = 100)[[1]]

d <- predict(c_mlt_BH, newdata = BostonHousing2, q = q, type="distribution")
lp <- c(predict(c_mlt_BH, newdata = BostonHousing2, q = 0, terms = "bshift"))
nd <- expand.grid(q = q, lp = -lp)
nd$d <- c(d)
pfun <- function(x, y, z, subscripts, at, ...) {
    panel.contourplot(x, y, z, subscripts, at = 1:9/10, ...)
    panel.xyplot(x = -lp, y = BostonHousing2$cmedv, pch = 19, cex = 1.5,
                 col = rgb(.1, .1, .1, .1), ...)
}
print(contourplot(d ~ lp + q, data = nd, panel = pfun, xlab = "Linear Predictor",
                  ylab = "Observed"))


###################################################
### code chunk number 11: BH-trtree
###################################################
plot(tt_BH, tp_args = list(type = "density", id = FALSE, ylines = 0, K = 100))

dev.off()

sessionInfo()
