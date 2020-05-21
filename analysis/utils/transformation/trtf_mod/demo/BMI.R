## ----setup, echo = FALSE, results = "hide", message = FALSE, warning = FALSE----
### from CRAN
library("lattice")
library("latticeExtra")
library("multcomp")
library("memisc")
library("Matrix")
options(Matrix.warn = FALSE)
library("colorspace")
library("grid")

library("libcoin")
library("inum")
library("partykit")
library("ATR")
library("trtf")
library("mlt")

pdf("BMI.pdf", width = 12, height = 8)

col <- diverge_hcl(2, h = c(246, 40), c = 120, l = c(65, 90), alpha = .75)
trellis.par.set(list(plot.symbol = list(col="blue",pch=18, cex=0.5),
                     box.rectangle = list(col=1),
                     box.umbrella = list(lty=1, col=1),
                     strip.background = list(col = "white")))
ltheme <- canonical.theme(color = FALSE)     ## in-built B&W theme
ltheme$strip.background$col <- "transparent" ## change strip bg
lattice.options(default.theme = ltheme)

set.seed(290875)

### original SGB data may NOT be distributed.
### We use a simulated dataset. BMI values were simulated
### from the transformation forest. All explanatory
### variables (except sex and smoking) were permuted.
### Results will differ from what is printed in the manuscript.
frda <- system.file("BMI_artificial.rda", package = "trtf")
load(frda)
xvars <- c("age40", "edu", "swiss", "frveg", "physcat", "agramtag", "regl2")

# learn$age <- as.double(learn$age)
vn <- c(age = "Age", 
        agramtag = "Alcohol intake", 
        frveg = "Fruit and vegetables", 
        physcat = "Physical activity", 
        edu = "Education", 
        swiss = "Nationality", 
        regl2 = "Region",
        sex = "Sex",
        smoking = "Smoking")
learn$wght <- with(learn, wght / sum(wght) * nrow(learn))

order <- 5
d2c <- function(x) formatC(round(x, 2), digits = 2, format = "f")

## ----ss-ecdf, fig = "pdf", echo = FALSE, fig.height = 5------------------
myprepanel <- function (x, y, f.value = NULL, ...) 
{
    ans <- prepanel.default.qqmath(x, f.value = f.value, distribution = qunif)
    with(ans, list(xlim = ylim, ylim = c(0, 1), dx = dy, dy = dx))
}


mypanel <- function (x, y, f.value = NULL, type = "s", groups = NULL, qtype = 7, 
    ref = TRUE, ...) 
{
    if (ref) {
        reference.line <- trellis.par.get("reference.line")
        do.call(panel.abline, c(list(h = c(0, 1)), reference.line))
    }
    x <- as.numeric(x)
    distribution <- qunif
    nobs <- sum(!is.na(x))
    if (!is.null(groups)) {
        panel.superpose(x, y = NULL, f.value = f.value, type = type, 
            distribution = distribution, qtype = qtype, groups = groups, 
            panel.groups = panel.ecdfplot, ...)
    }
    else if (nobs) {
        if (is.null(f.value)) {
            panel.xyplot(x = sort(x), y = cumsum(y[order(x)]) / sum(y),
                type = type, ...)
        }
        else {
            p <- if (is.numeric(f.value)) 
                f.value
            else f.value(nobs)
            panel.xyplot(x = quantile(x, p, names = FALSE, type = qtype, 
                na.rm = TRUE), y = distribution(p), type = type, 
                ...)
        }
    }
}

#library("latticeExtra")
xyplot(wght~bmi | smoking + sex, data = learn, col = col[1], lwd = 3, xlab = "BMI", 
       panel = mypanel, prepanel = myprepanel, ylab = "Empirical CDF")
# ecdfplot(~ bmi | smoking + sex, data = learn, col = col[1], lwd = 3, xlab = "BMI")

## ----ss1, echo = FALSE---------------------------------------------------
### lm; cell-means model
lm_ss1 <- lm(bmi ~ smoking:sex - 1, data = learn, weights = learn$wght)
ll_ss1 <- logLik(lm_ss1)
sd_ss1 <- summary(lm_ss1)$sigma

## ----ss1-tab, echo = FALSE, results = "asis"-----------------------------
ci <- d2c(confint(lm_ss1))
tb <- paste(d2c(coef(lm_ss1)), paste(" (", ci[,1], "--",  ci[,2], sep = ""), ")", sep = "")
tb <- matrix(tb, ncol = 2)
attr(tb, "row.vars") <- list(Smoking = levels(learn$smoking))
attr(tb, "col.vars") <- list(Sex = levels(learn$sex))
class(tb) <- "ftable"
toLatex(tb, digits = 3, useDcolum = FALSE)

## ----ss1-plot, fig = "pdf", echo = FALSE, fig.height = 5-----------------
### set-up grid
b <- 150:400 / 10
sm <- sort(unique(learn$smoking))
se <- sort(unique(learn$sex))
nd <- expand.grid(bmi = b, smoking = sm, sex = se)

### compute normal probabilities
nd$mean <- predict(lm_ss1, newdata = nd)
nd$sd <- summary(lm_ss1)$sigma
nd$prob1 <- pnorm(nd$bmi, mean = nd$mean, sd = nd$sd)


plotfun <- function(prob, data, ...) {
    fm <- as.formula(paste(prob, "~ bmi | smoking + sex"))
    xyplot(fm, data = data, type = "l", 
        panel = function(x, y, subscripts, ...) {
            tmp <- subset(learn, sex == unique(nd[subscripts, "sex"]) & 
                                     smoking == unique(nd[subscripts, "smoking"]))
            mypanel(tmp$bmi, tmp$wght, lwd = 3, col = col[1])
            panel.xyplot(x, y, ...)
    }, col = col[2],  xlab = "BMI", ylab = "CDF", lwd = 2, ...)
}
plotfun("prob1", nd)#, main = "lm constant var")

## ----ss2, echo = FALSE---------------------------------------------------
### varying variances
sd_ss2 <- c(tapply(resid(lm_ss1), iss <- with(learn, interaction(smoking, sex)), sd))
ll_ss2 <- sum(learn$wght * dnorm(learn$bmi, mean = coef(lm_ss1)[iss], sd = sd_ss2[iss], log = TRUE))

## ----ss2-plot, fig = "pdf", echo = FALSE, fig.height = 5-----------------
nd$sd2 <- rep(sd_ss2, each = length(b))
nd$prob2 <- pnorm(nd$bmi, mean = nd$mean, sd = nd$sd2)
plotfun("prob2", nd)# , main = "lm varying var")

## ----ss3, echo = FALSE, cache = TRUE, message = FALSE, warning = FALSE, results = "hide"----
bSMK <- as.basis(~ smoking - 1, data = learn)
bSEX <- as.basis(~ sex - 1, data = learn)
vBMI <- numeric_var("bmi", bounds = c(0, Inf),
    support = quantile(learn$bmi, prob = c(.01, .99), na.rm = TRUE),
    add = c(-5, 5))
bBMI <- Bernstein_basis(vBMI, order = order, ui = "increasing")
ctm_ss3 <- ctm(bBMI, interacting = b(sm = bSMK, sex = bSEX),
               data = learn, todistr = "Normal")
mlt_ss3 <- mlt(ctm_ss3, data = learn, scale = TRUE, weights = learn$wght)
ll_ss3 <- logLik(mlt_ss3)

## ----ss3-plot, fig = "pdf", echo = FALSE, message = FALSE, warning = FALSE, results = "hide", fig.height = 5----
nd$prob3 <- predict(mlt_ss3, newdata = nd, type = "distribution")
plotfun("prob3", nd)#, main = "nonlin mlt | smoking + sex")

## ----ss3-trafo-plot, fig = "pdf", echo = FALSE, message = FALSE, warning = FALSE, results = "hide", fig.height = 5----
nd3 <- nd
nd3$trafo3 <- predict(mlt_ss3, newdata = nd, type = "trafo")
bBMIlin <- as.basis(~ bmi, data = learn)
ctm_ss3a <- ctm(bBMIlin, interacting = b(sm = bSMK, sex = bSEX),
               data = learn, todistr = "Normal")
mlt_ss3a <- mlt(ctm_ss3a, data = learn, scale = TRUE, weights = learn$wght)
ll_ss3a <- logLik(mlt_ss3a)
nd3a <- nd
nd3a$trafo3 <- predict(mlt_ss3a, newdata = nd, type = "trafo")
nd3$type <- "Non-normal"
nd3a$type <- "Normal"
ndd <- rbind(nd3, nd3a)
xyplot(trafo3 ~ bmi | smoking + sex, data = ndd, type = "l", groups = type,
       col = col, lwd = 3, xlab = "BMI", ylab = expression(paste("Transformation Function ", h)))
# plotfun("trafo3", nd)#, main = "nonlin mlt | smoking + sex")

## ----ss3-plot-density, fig = "pdf", echo = FALSE, message = FALSE, warning = FALSE, results = "hide", fig.height = 5----
nd$density <- predict(mlt_ss3, newdata = nd, type = "density")
xyplot(density ~ bmi | smoking + sex, data = nd, type = "l", 
       col = col[2],  xlab = "BMI", ylab = "Density", lwd = 2)

## ----ss4, echo = FALSE, cache = TRUE, message = FALSE, warning = FALSE, results = "hide"----
bSEXtrt <- as.basis(~ sex, data = learn)
bSMKtrt <- as.basis(~ smoking, data = learn, remove_intercept = TRUE, negative = TRUE)
### non-normal transformation model, conditional on sex with sex-specific smoking shift
ctm_ss4 <- ctm(bBMI, interacting = bSEX,
               shifting = b(smktrt = bSMKtrt, sextrt = bSEXtrt), todistr = "Normal")
mlt_ss4 <- mlt(ctm_ss4, data = learn, scale = TRUE, weights = learn$wght)
ll_ss4 <- logLik(mlt_ss4)

## ----ss4-plot, fig = "pdf", echo = FALSE, message = FALSE, warning = FALSE, results = "hide", fig.height = 5----
nd$prob4 <- predict(mlt_ss4, newdata = nd, type = "distribution")
plotfun("prob4", nd)#, main = "ninlin mlt | sex - smoking(sex)")

## ----ss5, echo = FALSE, cache = TRUE, message = FALSE, warning = FALSE, results = "hide"----
### non-normal transformation model, conditional on sex with sex-specific smoking shift
ctm_ss5 <- ctm(bBMI, interacting = bSEX,
               shifting = b(smktrt = bSMKtrt, sextrt = bSEXtrt), todistr = "Logistic")
mlt_ss5 <- mlt(ctm_ss5, data = learn, scale = TRUE, weights = learn$wght)
ll_ss5 <- logLik(mlt_ss5)

## ----ss5-tab, echo = FALSE, results = "asis", warning = FALSE------------
K <- glht(lm(bmi ~ smoking, data = learn), linfct = mcp(smoking = "Dunnett"))$linfct[, -1]
X <- model.matrix(b(smktrt = bSMKtrt, sextrt = bSEXtrt), data = expand.grid(smoking = sm, sex = se))
### bSMKtrt has negative = TRUE, but we need -X %*% beta here!
X <- X * (-1)
K <- contrMat(table(learn$smoking), "Dunnett")
class(K) <- "matrix"
K <- as(K, "Matrix")
K2 <- as(bdiag(K, K), "matrix")
rownames(K2) <- c(paste("Female", rownames(K), sep = ":"), paste("Male", rownames(K), sep = ":"))
cf <- coef(mlt_ss5)
ci <- confint(glht(mlt_ss5, vcov. = vcov, linfct = K2 %*% cbind(matrix(0, ncol = length(cf) - ncol(X), nrow = nrow(X)), X)), 
              calpha = univariate_calpha())
tb <- paste(d2c(exp(ci$confint[,1])), paste(" (", d2c(exp(ci$confint[,2])), "--",
d2c(exp(ci$confint[,3])), sep = ""), ")", sep = "")
tb <- rbind(c(1, 1), matrix(tb, ncol = 2))
attr(tb, "row.vars") <- list(Smoking = levels(learn$smoking))
attr(tb, "col.vars") <- list(Sex = levels(learn$sex))
class(tb) <- "ftable"
toLatex(tb, digits = 3, useDcolum = FALSE)

## ----ss6, echo = FALSE, cache = TRUE, message = FALSE, warning = FALSE, results = "hide"----
### non-normal transformation model, conditional on sex with smoking shift
ctm_ss6 <- ctm(bBMI, interacting = bSEX,
               shifting = bSMKtrt, todistr = "Logistic")
mlt_ss6 <- mlt(ctm_ss6, data = learn, scale = TRUE, weights = learn$wght)
ll_ss6 <- logLik(mlt_ss6)

## ----cmpx_tree, echo = FALSE, cache = TRUE-------------------------------
mod0 <- ctm(bBMI, todistr = "Logistic")
mtmp <- mlt(mod0, data = learn)
xvars2 <- xvars
xvars2[xvars2 == "age40"] <- "age"
fm <- as.formula(paste("bmi ~ sex + smoking + ", 
                       paste(xvars2, collapse = "+")))
ctrl <- ctree_control(maxdepth = 4, nmax = 100, 
                      alpha = 1 / 10000)
tlearn <- learn
levels(tlearn$edu) <- c("I", "II", "III")
levels(tlearn$smoking) <- c("N", "F", "L", "M", "H")
tlearn$smoking <- ordered(as.character(tlearn$smoking), 
                          levels = c("N", "F", "L", "M", "H"),
                          labels = c("N", "F", "L", "M", "H"))
levels(tlearn$physcat) <- c(">2", "1-2", "0")
tlearn$physcat <- ordered(as.character(tlearn$physcat),
                          levels = c("0", "1-2", ">2"),
                          labels = c("0", "1-2", ">2"))
iw <- tlearn$wght
iw[iw < .5] <- 1
iw <- as.integer(round(iw))
trtr <- trafotree(mod0, formula = fm, data = tlearn,
                  control = ctrl, weights = iw)
### needs newdata to compute logLik for unbinned data
ll_cmpx_tree <- logLik(trtr, newdata = tlearn)

## ----cmpx-tree-plot, echo = FALSE, fig = "pdf", fig.width = 13, fig.height = 8----
### plot(trtr, tp_args = list(type = "density", id = FALSE, ylines = 0, K = 100, fill = col[1]))

node_BMI <- function(obj, col = "black", bg = "white", fill = "transparent",
                     ylines = 2, id = TRUE, mainlab = NULL, gp = gpar(), K = 20,
                     type = c("trafo", "distribution", "survivor",  
                              "density", "logdensity", "hazard",    
                              "loghazard", "cumhazard", "quantile"),
                     flip = FALSE, axes = TRUE, xaxis = NULL, ...)
{
    mod <- obj$model
    q <- mkgrid(mod, n = K)[[mod$response]]
    type <- match.arg(type)

    if (type %in% c("distribution", "survivor")) {
        yscale <- c(0, 1)
    } else {
        pr <- trtf:::predict.trafotree(obj, q = q, type = type, ...)
        yscale <- range(pr) * c(1, 1.05)
    }
    xscale <- range(q)
    axes <- rep_len(axes, 2)

    ### panel function for ecdf in nodes
    rval <- function(node) {

        nid <- id_node(node)
        dat <- data_party(obj, nid)
        wn <- dat[["(weights)"]]   

        cf <- obj$coef[as.character(nid),]
        coef(mod) <- cf
        y <- as.vector(predict(mod, newdata = data.frame(1), q = q, type = type))

        ## set up plot
        q <- q - xscale[1]
        q <- q / diff(xscale)
        y <- y - yscale[1]
        y <- y / diff(yscale)

        top_vp <- viewport(layout = grid.layout(nrow = 1, ncol = 2,
                           widths = unit(c(1, ylines), c("null", "lines")),
                           heights = unit(1, "null")),
                           width = unit(1, "npc"),
                           height = unit(1, "npc"), # - unit(2, "lines"),
                           name = paste("node_mlt", nid, sep = ""), gp = gp)

        pushViewport(top_vp)
        grid.rect(gp = gpar(fill = bg, col = 0))

        plot <- viewport(layout.pos.col=1, layout.pos.row=1,
                         xscale = if(flip) yscale else xscale,
                         yscale = if(flip) xscale else yscale,
                         name = paste0("node_mlt", nid, "plot"),
                         clip = FALSE)

        pushViewport(plot)
        if(axes[1]) {
            if (is.null(xaxis)) {
                 grid.xaxis()
            } else {
                 if (nid == xaxis)
                     grid.xaxis()
            }
        }
        if(axes[2]) grid.yaxis(main = FALSE)
        grid.rect(gp = gpar(fill = "transparent"))
        grid.clip()
        if(flip) {
          if(fill != "transparent") {
            grid.polygon(c(min(y), y, min(y)), c(q[1], q, q[K]), gp = gpar(col = col, fill = fill))
          } else {
            grid.lines(y, q, gp = gpar(col = col))
          }
        } else {
          if(fill != "transparent") {
            grid.polygon(c(q[1], q, q[K]), c(min(y), y, min(y)), gp = gpar(col = col, fill = fill))
          } else {
            grid.lines(q, y, gp = gpar(col = col))
          }
        }
        upViewport(2)
    }

    return(rval)
}
class(node_BMI) <- "grapcon_generator"

nid <- min(nodeids(trtr, terminal = TRUE))

### weights are caseweights for ctree but not for this application
tmp <- unclass(trtr)
tmp$data[["(weights)"]] <- 1L
tmp$fitted[["(weights)"]] <- 1L
class(tmp) <- class(trtr)
plot(rotate(tmp), terminal_panel = node_BMI, nobs.loc='top',
     tp_args = list(type = "density", id = FALSE, ylines = 3, K = 100, fill = col[1], xaxis = nid))


## ----cmpx_forest, echo = FALSE, cache = TRUE, results = "hide", message = FALSE, warning = FALSE----
trtf <- traforest(mod0, formula = fm, data = learn, weights = learn$wght,
                  control = ctree_control(nmax = 100, mincriterion = 0, minsplit = 500), 
                  ntree = 100, mtry = 9, trace = TRUE)

## ----cmpx_vi, echo = FALSE, cache = TRUE---------------------------------
vi <- varimp(trtf)

## ----cmpx_ll, echo = FALSE, cache = TRUE---------------------------------
set.seed(290875)
n <- nrow(learn)
cf <- predict(trtf, newdata = learn, type = "coef")
i <- which(sapply(cf, function(x) any(abs(x) > 10)))
if (length(i) > 0) 
    cf[i] <- predict(trtf, newdata = learn[i,], type = "coef")
ll_cmpx_forest <- logLik(trtf, coef = cf)

## ----cmpx_predict, echo = FALSE, cache = TRUE----------------------------
### prepare figure
nd2 <- lapply(learn[,xvars2], function(x) {
    if (is.factor(x)) return(sort(x)[1])
    if (is.data.frame(x)) return(NULL)
    median(x)
})
nd2$sex <- sort(unique(learn$sex))
nd2$smoking <- sort(unique(learn$smoking))

nd3 <- nd2
nd3$age <- 18:70
nd4 <- expand.grid(nd3)
bmi <- 15:40

p <- predict(trtf, newdata = nd4, type = "distribution", q = bmi)

## ----cmpx-varimp-plot, echo = FALSE, fig = "pdf", fig.height = 5---------
vi <- vi[c("sex", "smoking", xvars2)]
vi[is.na(vi)] <- 0
# vi <- vi / sum(vi)
names(vi) <- vn[names(vi)]
vi <- vi[order(vi)]
barchart(~ vi, xlim = c(0, max(vi) * 1.1), xlab = "Mean Decrease
Log-Likelihood", col = col[1])

## ----cmpx-age-plot, echo = FALSE, fig = "pdf", fig.height = 5------------
tmp <- expand.grid(c(list(bmi = bmi), nd3))
tmp$p <- unlist(p)

contourplot(p ~ age + bmi | smoking + sex, data = tmp, 
            panel = function(x, y, z, ...) {
                panel.contourplot(x, y, z, ...)
#                 panel.rug(x, col = rgb(.1, .1, .1, .01))
            },
            at = 1:9/10, ylim = c(18, 34), xlab = "Age", ylab = "BMI",
            labels = list(cex = .5))

## ----cmpx-ctm, echo = FALSE, cache = TRUE, message = FALSE, warning = FALSE, results = "hide"----
bAGE <- Bernstein_basis(numeric_var("age", support = c(20, 70)), order = order)
xfm <- as.formula(paste("~", paste(xvars[xvars != "age40"], collapse = "+")))
ctm_ctm <- ctm(bBMI, interacting = b(sex = bSEX, age = bAGE),
               shifting = c(sstrt = b(smktrt = bSMKtrt, sextrt = bSEX),
                            lp = as.basis(xfm, data = learn, remove_intercept = TRUE, negative = TRUE)),
             todistr = "Logistic")
mlt_ctm <- mlt(ctm_ctm, data = learn, scale = TRUE, weights = learn$wght)
ll_ctm <- logLik(mlt_ctm)

## ----cmpx-ctm-plot, echo = FALSE, fig = "pdf", fig.height = 5------------
nd3 <- nd2
nd3$age <- 18:70
nd3$bmi <- 15:35
nd3 <- nd3[c("bmi", names(nd3)[names(nd3) != "bmi"])]

p <- predict(mlt_ctm, newdata = nd3, type = "distribution")
tmp <- do.call("expand.grid", nd3)
tmp$p <- c(p)

contourplot(p ~ age + bmi | smoking + sex, data = tmp,
            at = 1:9/10, ylim = c(18, 34), xlab = "Age", ylab = "BMI",
            labels = list(cex = .5))

## ----dr, echo = FALSE, cache = TRUE, message = FALSE, warning = FALSE, results = "hide"----
bAGE <- as.basis(~ age, data = learn)

### distribution regression | age 
ctm_dr <- ctm(bBMI, interacting = b(sex = bSEX, age = bAGE),
            shifting = c(sstrt = b(smktrt = bSMKtrt, sextrt = bSEX),
                            lp = as.basis(xfm, data = learn, remove_intercept = TRUE, negative = TRUE)),
            todistr = "Logistic")
mlt_dr <- mlt(ctm_dr, data = learn, scale = TRUE, weights = learn$wght)
ll_dr <- logLik(mlt_dr)

## ----dr-plot, echo = FALSE, fig = "pdf", fig.height = 5------------------
nd3 <- nd2
nd3$age <- 18:70
nd3$bmi <- 15:35
nd3 <- nd3[c("bmi", names(nd3)[names(nd3) != "bmi"])]

p <- predict(mlt_dr, newdata = nd3, type = "distribution")
tmp <- do.call("expand.grid", nd3)
tmp$p <- c(p)

contourplot(p ~ age + bmi | smoking + sex, data = tmp,
            at = 1:9/10, ylim = c(18, 34), xlab = "Age", ylab = "BMI",
            labels = list(cex = .5))

## ----smpl, echo = FALSE, cache = TRUE, message = FALSE, warning = FALSE, results = "hide"----
bAGE <- as.basis(~ age - 1, data = learn, negative = TRUE)
ctm_smpl <- ctm(bBMI, interacting = bSEX,
            shifting = c(sstrt = b(smktrt = bSMKtrt, sextrt = bSEX),
                         sexage = b(sex = bSEX, age = bAGE),
                         lp = as.basis(xfm, data = learn, remove_intercept = TRUE, negative = TRUE)),
            todistr = "Logistic")
mlt_smpl <- mlt(ctm_smpl, data = learn, scale = TRUE, weights = learn$wght)
ll_smpl <- logLik(mlt_smpl)

ciMale <- confint(mlt_smpl)["sexMale:age",]
ciFemale <- confint(mlt_smpl)["sexFemale:age",]

ageMale <- paste("$", d2c(exp(coef(mlt_smpl))["sexMale:age"]), "$ ($", d2c(exp(ciMale)[1]),
  "-", d2c(exp(ciMale)[2]), "$)", sep = "")

ageFemale <- paste("$", d2c(exp(coef(mlt_smpl))["sexFemale:age"]), "$ ($", d2c(exp(ciFemale)[1]),
  "-", d2c(exp(ciFemale)[2]), "$)", sep = "")



## ----summary-tab, echo = FALSE, results = "asis"-------------------------
cf_ctm <- coef(mlt_ctm)
cf_dr <- coef(mlt_dr)
cf_smpl <- coef(mlt_smpl)

xparm <- names(cf_ctm)[-grep("^Bs", names(cf_ctm))]

ci_ctm <- exp(confint(mlt_ctm)[xparm,])
ci_dr <- exp(confint(mlt_dr)[xparm,])
ci_smpl <- exp(confint(mlt_smpl)[xparm,])

tab <- cbind(paste(d2c(exp(cf_ctm[xparm])), " (", d2c(ci_ctm[,1]), "--", d2c(ci_ctm[,2]), ")", sep = ""),
paste(d2c(exp(cf_dr[xparm])), " (", d2c(ci_dr[,1]), "--", d2c(ci_dr[,2]), ")", sep = ""),
paste(d2c(exp(cf_smpl[xparm])), " (", d2c(ci_smpl[,1]), "--", d2c(ci_smpl[,2]), ")", sep = ""))

tab2 <- rbind(" ", rbind("$1$", tab))
rownames(tab2) <- 1:nrow(tab2)
rownames(tab2)[1:2] <- c("empty", "one")
rownames(tab2)[-(1:2)] <- xparm

tab2 <- tab2[c("empty", "one", "smokingFormer:sexFemale",
               "smokingLight:sexFemale",
               "smokingMedium:sexFemale",
               "smokingHeavy:sexFemale",
               "empty", "one", "smokingFormer:sexMale",
               "smokingLight:sexMale",
               "smokingMedium:sexMale",
               "smokingHeavy:sexMale",
               "agramtag", 
               "empty", "one", "frvegLow", 
               "empty", "one", "physcat1-2 d/w", "physcatnone",
               "empty", "one", "eduSecondary II", "eduTertiary",
               "empty", "one", "swissForeign", 
               "empty", "one", "regl2French", "regl2Italian"),,drop = FALSE]
rownames(tab2) <- c("Smoking (Females)", "\\hspace{.4cm}Never", "\\hspace{.4cm}Former", "\\hspace{.4cm}Light", "\\hspace{.4cm}Medium", "\\hspace{.4cm}Heavy",
                    "Smoking (Males)", "\\hspace{.4cm}Never", "\\hspace{.4cm}Former", "\\hspace{.4cm}Light", "\\hspace{.4cm}Medium", "\\hspace{.4cm}Heavy",
                    "Alcohol intake (g/d)",
                    "Fruit and vegetables", "\\hspace{.4cm}High", "\\hspace{.4cm}Low",
                    "Physical activity", "\\hspace{.4cm}High", "\\hspace{.4cm}Moderate", "\\hspace{.4cm}Low",
                    "Education", "\\hspace{.4cm}Mandatory (I)", "\\hspace{.4cm}Secondary (II)", "\\hspace{.4cm}Tertiary (III)",
                    "Nationality", "\\hspace{.4cm}Swiss", "\\hspace{.4cm}Foreign",
                    "Region", "\\hspace{.4cm}German speaking", "\\hspace{.4cm}French speaking", "\\hspace{.4cm}Italian speaking")
colnames(tab2) <- c("Model (\\ref{mod:ctm})", "Model (\\ref{mod:dr})", "Model (\\ref{mod:smpl})")
toLatex(tab2, digits = 3, useDcolum = FALSE)

## ----ss3-high, echo = FALSE, results = "hide", cache = TRUE--------------
bBMIh <- Bernstein_basis(vBMI, order = 2 * order, ui = "increasing")
ctm_ss3h <- ctm(bBMIh, interacting = b(sm = bSMK, sex = bSEX),
               data = learn, todistr = "Normal")
mlt_ss3h <- mlt(ctm_ss3h, data = learn, scale = TRUE, weights = learn$wght)
ll_ss3h <- logLik(mlt_ss3h)

dev.off()
