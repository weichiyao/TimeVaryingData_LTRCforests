load('results_empeval.rda')

for (name in names(ret)){
    ret[, name] <- as.numeric(ret[, name])
}

meth <- c('rsf', 'rg', 'cf', 'tf_W_alpha', 'tf_B_alpha', 'tf_W', 'tf_B', 'L1')
tick_meth <- c('RSF NP(\U0001D6FC)', 
               'Ranger NP(\U0001D6FC)', 
               'CForest NP(\U0001D6FC)', 
               'TF W(\U0001D6FC)', 
               'TF Bs(\U0001D6FC)', 
               'TF W(\U0001D6DD)',
               'TF Bs(\U0001D6DD)',
               'L1') 

### expression(W(bold(vartheta)))


args <- unique(ret[, c('p', 'prod_shift', 'prod_scale')])

LogLikDiff <- ret[, meth] - ret$loglik


cairo_pdf(file = 'empeval_boxplot.pdf')
par(mfrow=c(4,2), oma = c(10, 6, 4, 8),
    mar=c(0, 0, 0, 0))


for (i in 1:nrow(args)){
  arg <- args[i,]
  idx <- (ret$prod_shift == arg$prod_shift) & (ret$prod_scale == arg$prod_scale) & (ret$p == arg$p)
  boxplot(LogLikDiff[idx, ], las=2, axes=FALSE, frame.plot=TRUE, ylim=c(-320, 15), at = c(1:3, 4.5, 5.5, 7, 8, 9.5))
  if (i > 6) {
    axis(1, col = "grey40", col.axis = "grey20", at = c(1:3, 4.5, 5.5, 7, 8, 9.5), labels=tick_meth, las = 2, cex.axis = 1.5)
  }
  if (i %in% c(1, 3, 5, 7)){
    axis(2, col = "grey40", col.axis = "grey20", at=c(-300, -250, -200, -150, -100, -50, 0))
  }
  if (i == 1){ mtext("low", side=3, line=1, cex = 1.2)}
  if (i == 2) {
    mtext('high', side=3, line=1, cex = 1.2)
    mtext('No', side=4, line=1, las=2, cex = 1.2)}
  if (i == 4) {
    mtext('PH', side=4, line=1, las=2, cex = 1.2)
  }
  if (i == 6){
    mtext('Non-PH', side=4, line=1, las=2, cex = 1.2)
  }
  if(i == 8){
    mtext('Combined', side=4, line=1, las=2, cex = 1.2)
  }
  if(i == 3){
    mtext('log-likelihood difference', outer=TRUE, side=2, line=3, cex = 1.2)
  }
}
dev.off()

