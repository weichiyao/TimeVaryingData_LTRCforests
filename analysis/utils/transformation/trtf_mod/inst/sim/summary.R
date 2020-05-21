
library("lattice")

plotfun <- function(file, doplot = TRUE) {

load(file)

parm <- lapply(out, function(x) x$parm)

df <- NULL

for (i in 1:length(parm)) {

    tmp <- c()

    chkRsk <- out[[i]][["checkRisk"]]
    out[[i]][["L.1"]] <- chkRsk[,,1]
    out[[i]][["L.5"]] <- chkRsk[,,2]
    out[[i]][["L.9"]] <- chkRsk[,,3]
    nm <- names(out[[i]])
    out[[i]] <- out[[i]][nm[nm != "checkRisk"]]

    for (j in 2:length(out[[i]])) {
        dummy <- as.data.frame(out[[i]][[j]])
        vars <- 1:ncol(dummy)
        vnames <- colnames(dummy)[vars]
        dummy$id <- factor(1:nrow(dummy))

        dummy <- reshape(dummy, varying = list(vars), direction = "long", 
                         idvar = "id", timevar = "model", 
                         v.names = names(out[[i]])[j])
        dummy$model <- factor(dummy$model, levels = vars, labels = vnames)

        if (j == 2) { 
            tmp <- dummy
        } else {
            tmp <- merge(tmp, dummy, by = c("id", "model"))
        }
    }
    tmp <- cbind(tmp, as.data.frame(parm[i])[rep(1, nrow(tmp)),,drop = FALSE])
    df <- rbind(df, tmp)
}

save(df, file = paste(file, "_out.rda", sep = ""))

if (doplot) {
df <- subset(df, model != "ttBern")
df <- subset(df, model != "ttBernExSplit")

pdf(paste(file, ".pdf", sep = ""))
print(bwplot(ll ~ model | p + tau + prod_mu + prod_sigma, data = df, scales = list(y =
"free", x = list(rot = 45)), main = file))#, ylim = c(-800, -350)))
print(bwplot(time ~ model | p + tau + prod_mu + prod_sigma, data = df, scales = list(y =
"free", x = list(rot = 45)), main = file))#, ylim = c(-800, -350)))
print(bwplot(L.1 ~ model | p + tau + prod_mu + prod_sigma,, data = df,
scales = list(y = "free", x = list(rot = 45)), main = file))
print(bwplot(L.5 ~ model | p + tau + prod_mu + prod_sigma,, data = df,
scales = list(y = "free", x = list(rot = 45)), main = file))
print(bwplot(L.9 ~ model | p + tau + prod_mu + prod_sigma,, data = df,
scales = list(y = "free", x = list(rot = 45)), main = file))
dev.off()

}

}

plotfun("2d.rda", doplot = FALSE)

plotfun("lognormal_2d.rda", doplot = FALSE)

plotfun("friedman.rda", doplot = FALSE)

plotfun("lognormal_friedman.rda", doplot = FALSE)

plotfun("timings.rda", doplot = FALSE)

