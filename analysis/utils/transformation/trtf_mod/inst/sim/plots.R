
library("lattice")
library("grid")
lwd <- 1.5
col <- rgb(.1, .1, .1, .05)
colR <- rgb(.75, 0, 0, .8)
colRlight <- rgb(.75, 0, 0, .1)
colB <- rgb(0, 0, 0.75, .8)
trellis.par.set(list(plot.symbol = list(col="black",pch=18, cex=0.75),
                     box.rectangle = list(col=1),
                     box.umbrella = list(lty=1, col=1),
                     strip.background = list(col = "white")))
ltheme <- canonical.theme(color = FALSE)     ## in-built B&W theme
ltheme$strip.background$col <- "transparent" ## change strip bg
lattice.options(default.theme = ltheme)

fn <- function(n) formatC(n, big.mark = "{,}")

plotfun <- function(file, models = c("gt", "tt"), labels = models,
                    wplot = c("ll", "L.1", "L.5", "L.9"), 
                    hypo = c("H2", "H2", "H2", "H0"), main = file) {

  wplot <- match.arg(wplot, several.ok = TRUE)

  load(file)

  df0 <- subset(df, tau < 0.1 & model %in% models)
  df0$model <- df0$model[,drop = TRUE]
  df0$model <- factor(as.character(df0$model), levels = models, labels = labels)

  df0$p <- factor(df0$p, levels = sort(unique(df0$p)),
                  labels = c("Low dimensional", "High dimensional"))
  df0$prod_mu <- factor(df0$prod_mu, levels = 0:1, labels = c("Constant mean", "Varying mean"))
  df0$prod_mu <- factor(as.character(df0$prod_mu), 
                        levels = rev(levels(df0$prod_mu)), 
                        labels = rev(levels(df0$prod_mu)))
  df0$prod_sigma <- factor(df0$prod_sigma, levels = 0:1, 
                           labels = c("Constant variance", "Varying variance"))
  df0$prod_sigma <- factor(as.character(df0$prod_sigma), 
                           levels = rev(levels(df0$prod_sigma)), 
                           labels = rev(levels(df0$prod_sigma)))

  df0gt <- subset(df0, model == labels[models == "gt"])
  nm <- names(df0gt)
  nm[match(c("ll", "L.1", "L.5", "L.9"), nm)] <- paste0("gt", c("ll", "L.1", "L.5", "L.9"))
  names(df0gt) <- nm
  df0gt$model <- NULL
  df0 <- merge(df0, df0gt, by = c("id", "p", "prod_mu", "prod_sigma"))
  df0 <- subset(df0, model != labels[models == "gt"])
  df0$model <- df0$model[, drop = TRUE]
  df0$rightlabel <- with(df0, interaction(prod_mu, prod_sigma))
  levels(df0$rightlabel) <- hypo

  mypanel <- function(x, y, groups, subscripts, ...) {
    lim <- current.panel.limits(unit = "native")
    at <- mean(ylim <- lim$ylim)
    lrect(lim$xlim[1], lim$ylim[1], lim$xlim[2], lim$ylim[2], col = "white")
    panel.bwplot(x = x, y = y, ...)
#    panel.abline(h = median(y[x == "Truth"]), lty = 3)
    panel.abline(h = 0, lty = 3)
    lab <- as.character(unique(df0$rightlabel[subscripts]))
    if (unique(df0$p[subscripts]) == "High dimensional")
        panel.axis("right", at = at, labels = lab, 
                   ticks = FALSE, outside = TRUE)
    tapply(1:length(y), groups[subscripts], function(i) {
        llines(x = 1:nlevels(x), y = y[i][order(x[i])], 
               col = rgb(0.1, 0.1, 0.1, 0.1))
    })
    lrect(lim$xlim[1], lim$ylim[2], lim$xlim[2], lim$ylim[2] * 100, 
          col = "white", border = "white")
  }

  myprepanel <- function(df, response) {
      df$response <- response
      function(x, y, groups, subscripts, ...) {
          ps <- df[subscripts, "prod_sigma"]
          pm <- df[subscripts, "prod_mu"]
          tdf <- subset(df, prod_sigma == ps & prod_mu == pm)
          bp <- boxplot(split(tdf[["response"]], tdf$model), plot = FALSE)
          ylim <- range(bp$stats) * c(1, 1.2)
#          ylim <- pmin(500, range(bp$stats) * c(.995, 1))
          if (ylim[1] > 0) ylim[1] <- 0
          list(ylim = ylim) 
      }
  }

### we need clipping off for outside panel.axis but this
### draws outside the box. 
ps <-  list(clip = list(panel = "off", strip = "on"),
                    layout.widths = list(right.padding = 5))
#                    layout.heights = list(top.padding = 0,
#                        main = 0, main.key.padding = 0, key.top = 0))

sc <- function(which.given, ..., factor.levels) {
    if (which.given == 1) {
        strip.default(which.given = which.given, ..., factor.levels)
    } else {
        if (which.given == 2) {
            strip.default(which.given = 2, ..., 
               factor.levels = c(expression(paste("Varying ", mu)), expression(paste("Constant ", mu))))
        } else {
            strip.default(which.given = 3, ..., 
               factor.levels = c(expression(paste("Varying ", sigma)), expression(paste("Constant ", sigma))))
        }
   }
}

  if ("ll" %in% wplot)
  print(bwplot(I(- ll + gtll) ~ model | p + prod_mu + prod_sigma, data = df0, 
        panel = mypanel, groups = id, do.out = FALSE,
        prepanel = myprepanel(df0, with(df0, - ll + gtll)),
        par.settings = ps, strip = sc,
        page = function(...) grid.text(main, x = .5, y = .980, 
                                       gp = gpar(cex = 1.25)),
        ylab = "Negative Log-likelihood difference", scales = list(y =
        "free", x = list(rot = 45)), layout = c(2, 4)))

  if ("L.1" %in% wplot)
  print(bwplot(L.1 - gtL.1 ~ model | p + prod_mu + prod_sigma, data = df0, 
        panel = mypanel, groups = id, do.out = FALSE,
        prepanel = myprepanel(df0, with(df0, L.1 - gtL.1)),
        par.settings = ps, strip = sc,
        page = function(...) grid.text(main, x = .5, y = .975, 
                                       gp = gpar(cex = 1.25)),
        ylab = "10% quantile risk difference", scales = list(y =
        "free", x = list(rot = 45)), layout = c(2, 4)))

  if ("L.5" %in% wplot)
  print(bwplot(L.5 - gtL.5 ~ model | p + prod_mu + prod_sigma, data = df0, 
        panel = mypanel, groups = id, do.out = FALSE,
        prepanel = myprepanel(df0, with(df0, L.5 - gtL.5)),
        par.settings = ps, strip = sc,
        page = function(...) grid.text(main, x = .5, y = .975, 
                                       gp = gpar(cex = 1.25)),
        ylab = "Absolute error difference", scales = list(y =
        "free", x = list(rot = 45)), layout = c(2, 4)))

  if ("L.9" %in% wplot)
  print(bwplot(L.9 - gtL.9 ~ model | p + prod_mu + prod_sigma, data = df0, 
        panel = mypanel, groups = id, do.out = FALSE,
        prepanel = myprepanel(df0, with(df0, L.9 - gtL.9)),
        par.settings = ps, strip = sc,
        page = function(...) grid.text(main, x = .5, y = .975, 
                                       gp = gpar(cex = 1.25)),
        ylab = "90% quantile risk difference", scales = list(y =
        "free", x = list(rot = 45)), layout = c(2, 4)))
}

forests <- c("gt", "tfbagg", "tfrf", "tfbaggBern", "tfrfBern", "cfbagg", "cfrf", "bagg", "rf")
trees <- c("gt", "ct", "tt", "ttExSplit", "ttM", "ttBern", "ttBernExSplit")
models <- c("gt", "ct", "tt", "ttBern", # "tfbagg", 
            "cfrf", 
            "rf", "rfBern",
            "tfrf", # "tfbaggBern", 
            "tfrfBern") # "cfbagg", 
            # "bagg", 
labels <- c("Truth", "CTree (P = 2)", "TTree (P = 2)", "TTree (P = 6)", 
            #"TBagging (lin)", 
            "CForest (P = 2)", 
            #"Bagging", 
            "RForest (P = 2)", "RForest (P = 6)",
            "TForest (P = 2)", 
            #"TBagging (nonlin)", 
            "TForest (P = 6)") 
            # "CBagging", 
logmodels <- c(
            "gt", 
            #"ct", 
            # "tt", 
            "ttBern", # "tfbagg", 
            "rfBern",
            # "tfrf", # "tfbaggBern", 
            "tfrfBern") # "cfbagg", 
            # "cfrf", 
            # "bagg", 
            #"rf", 
loglabels <- c(
            "Truth", 
            #"CTree (P = 2)", "TTree (P = 2)", 
            "TTree (P = 6)", 
            #"TBagging (lin)", 
            #"TForest (P = 2)", 
            #"TBagging (nonlin)", 
            "RForest (P = 6)",

            "TForest (P = 6)") 
            # "CBagging", 
            #"CForest (P = 2)", 
            #"Bagging", 
            #"RForest (P = 2)", 



### R code from vignette source '/home/thothorn/transforest/paper/empeval.Rnw'

###################################################
### code chunk number 2: 2d-ll
###################################################
plotfun("2d.rda_out.rda", models = models, labels = labels, wplot = "ll",
        main = "Tree-Structured Normal (H1a)", hypo = c("H2c", "H2c", "H2b", "H2a"))


###################################################
### code chunk number 3: ln-2d-ll
###################################################
plotfun("lognormal_2d.rda_out.rda", models = logmodels, labels = loglabels, wplot = "ll",
        main = "Tree-Structured Log-Normal (H1a)", hypo = c("H2d", "H2d", "H2b", "H2a"))


###################################################
### code chunk number 4: fm-ll
###################################################
plotfun("friedman.rda_out.rda", models = models, labels = labels, wplot = "ll",
        main = "Non-Linear Normal (H1b)", hypo = c("H2c", "H2c", "H2b", "H2a"))


###################################################
### code chunk number 5: ln-fm-ll
###################################################
plotfun("lognormal_friedman.rda_out.rda", models = logmodels, labels = loglabels, wplot = "ll",
         main = "Non-Linear Log-Normal (H1b)", hypo = c("H2d", "H2d", "H2b", "H2a"))


###################################################
### code chunk number 2: 2d-L1
###################################################
plotfun("2d.rda_out.rda", models = models, labels = labels, wplot = "L.1",
        main = "Tree-Structured Normal (H1a)", hypo = c("H2c", "H2c", "H2b", "H2a"))


###################################################
### code chunk number 3: 2d-L5
###################################################
plotfun("2d.rda_out.rda", models = models, labels = labels, wplot = "L.5",
        main = "Tree-Structured Normal (H1a)", hypo = c("H2c", "H2c", "H2b", "H2a"))


###################################################
### code chunk number 4: 2d-L9
###################################################
plotfun("2d.rda_out.rda", models = models, labels = labels, wplot = "L.9",
        main = "Tree-Structured Normal (H1a)", hypo = c("H2c", "H2c", "H2b", "H2a"))


###################################################
### code chunk number 5: ln-2d-L1
###################################################
plotfun("lognormal_2d.rda_out.rda", models = logmodels, labels = loglabels, wplot = "L.1",
        main = "Tree-Structured Log-Normal (H1a)", hypo = c("H2d", "H2d", "H2b", "H2a"))


###################################################
### code chunk number 6: ln-2d-L5
###################################################
plotfun("lognormal_2d.rda_out.rda", models = logmodels, labels = loglabels, wplot = "L.5",
        main = "Tree-Structured Log-Normal (H1a)", hypo = c("H2d", "H2d", "H2b", "H2a"))


###################################################
### code chunk number 7: ln-2d-L9
###################################################
plotfun("lognormal_2d.rda_out.rda", models = logmodels, labels = loglabels, wplot = "L.9",
        main = "Tree-Structured Log-Normal (H1a)", hypo = c("H2d", "H2d", "H2b", "H2a"))


###################################################
### code chunk number 8: fm-L1
###################################################
plotfun("friedman.rda_out.rda", models = models, labels = labels, wplot = "L.1",
        main = "Non-Linear Normal (H1b)", hypo = c("H2c", "H2c", "H2b", "H2a"))


###################################################
### code chunk number 9: fm-L5
###################################################
plotfun("friedman.rda_out.rda", models = models, labels = labels, wplot = "L.5",
        main = "Non-Linear Normal (H1b)", hypo = c("H2c", "H2c", "H2b", "H2a"))


###################################################
### code chunk number 10: fm-L9
###################################################
plotfun("friedman.rda_out.rda", models = models, labels = labels, wplot = "L.9",
        main = "Non-Linear Normal (H1b)", hypo = c("H2c", "H2c", "H2b", "H2a"))


###################################################
### code chunk number 11: ln-fm-L1
###################################################
plotfun("lognormal_friedman.rda_out.rda", models = logmodels, labels = loglabels, wplot = "L.1",
         main = "Non-Linear Log-Normal (H1b)", hypo = c("H2d", "H2d", "H2b", "H2a"))


###################################################
### code chunk number 12: ln-fm-L5
###################################################
plotfun("lognormal_friedman.rda_out.rda", models = logmodels, labels = loglabels, wplot = "L.5",
         main = "Non-Linear Log-Normal (H1b)", hypo = c("H2d", "H2d", "H2b", "H2a"))


###################################################
### code chunk number 13: ln-fm-L9
###################################################
plotfun("lognormal_friedman.rda_out.rda", models = logmodels, labels = loglabels, wplot = "L.9",
         main = "Non-Linear Log-Normal (H1b)", hypo = c("H2d", "H2d", "H2b", "H2a"))

