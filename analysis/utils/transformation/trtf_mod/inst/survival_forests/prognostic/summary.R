
load("results_empeval.rda")

i <- 1:nrow(ret)
nm <- c("p", "prod_shift", "prod_scale", "loglik")
m <- colnames(ret)
m <- m[!(m %in% nm)]

meth <- c('rsf', 'rg', 'cf', 'tf_W_alpha', 'tf_B_alpha', 'tf_W', 'tf_B', 'L1')
lab <- c('RSF NP(\U0001D6FC)', 
         'Ranger NP(\U0001D6FC)', 
         'CForest NP(\U0001D6FC)', 
         'TF W(\U0001D6FC)', 
         'TF Bs(\U0001D6FC)', 
         'TF W(\U0001D6DD)',
         'TF Bs(\U0001D6DD)',
         'L1') 


x <- expand.grid(i = i, m = m)
r <- ret[x$i,nm]
r$m <- factor(x$m, levels = meth, label = lab)
r$fll <- do.call("c", ret[,m])
r$id <- factor(i)
for (n in nm[1:3]) r[[n]] <- factor(r[[n]])

r$setup <- with(r, interaction(prod_shift, prod_scale))
levels(r$setup) <- lev <- c("No", "PH", "Non-PH", "Combined")

r$setup <- factor(as.character(r$setup), levels = rev(lev), labels = rev(lev))
levels(r$p) <- c("Low", "High")

r$fll[r$fll == 0] <- NA

library("lattice")

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

mypanel <- function(x, y, groups, subscripts, ...) {
    panel.abline(h = 0, lty = 3)
    panel.bwplot(x = x, y = y, ...)
    tapply(1:length(y), groups[subscripts], function(i) {
        llines(x = 1:nlevels(x), y = y[i][order(x[i])], 
               col = rgb(0.1, 0.1, 0.1, 0.1))
    })
}

bwplot(I(fll - loglik) ~ m | p + setup, data = r, groups = id,
       scales = list(x = list(rot = 45)), layout = c(2, 4), 
       panel = mypanel, ylim = c(-200, 25), 
       ylab = "Log-likelihood Difference")

