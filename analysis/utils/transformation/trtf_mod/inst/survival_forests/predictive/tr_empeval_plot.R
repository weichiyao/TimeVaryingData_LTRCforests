load('tr_results_empeval.rda')

meth <- c('w_alpha_beta', 'bs_alpha_beta', 'w_theta_beta', 'bs_theta_beta', 'w_theta_theta', 'bs_theta_theta')

tick_meth <- c('TF W(\U0001D6FC, \U0001D6FD)', 'TF Bs(\U0001D6FC, \U0001D6FD)', 
               'TF W(\U0001D6DD, \U0001D6FD)', 'TF Bs(\U0001D6DD, \U0001D6FD)',
               expression('TF W(\U0001D6DD, \U0001D6DD'['tr']*')'), expression('TF Bs(\U0001D6DD, \U0001D6DD'['tr']*')'))

#args <- data.frame(p = rep(sort(unique(ret$p)), times = 3),
#                   prod_shift_int = rep(c(1, 0, 1), each = 2),
#                   prod_scale_int = rep(c(0, 1, 1), each = 2))

#args$prod_shift_trt <- args$prod_shift_int
#args$prod_scale_trt <- args$prod_scale_int

args <- unique(ret[, c('p', 'prod_shift_int', 'prod_shift_trt', 'prod_scale_int', 'prod_scale_trt')])

LogLikDiff <- ret[, meth] -ret$loglik


cairo_pdf(file = 'tr_boxplots.pdf')
par(mfrow=c(3,2), oma = c(8, 6, 4, 8),
    mar=c(0, 0, 0, 0))


for (i in 1:nrow(args)){
  arg <- args[i,]
  idx <- (ret$prod_shift_int == arg$prod_shift_int) & 
         (ret$prod_shift_trt == arg$prod_shift_trt) & 
         (ret$prod_scale_int == arg$prod_scale_int) &
         (ret$prod_scale_trt == arg$prod_scale_trt) &
         (ret$p == arg$p)
  
  boxplot(LogLikDiff[idx, ], axes=FALSE, frame.plot=TRUE, ylim = c(-420, -20), at = c(1, 2, 3.5, 4.5, 6, 7))
  if (i > 4) {
    axis(1, col = "grey40", col.axis = "grey20", at = c(1, 2, 3.5, 4.5, 6, 7), labels=tick_meth, las=2, cex.axis = 1.5)
  }
  if (i %% 2 == 1){
    axis(2, col = "grey40", col.axis = "grey20", at=c(-50, -100, -150, -200, -250, -300, -350, -400))
  }
  if (i == 1){ mtext("low", side=3, line=1, cex = 1.2)}
  if (i == 2) {
    mtext('high', side=3, line=1, cex = 1.2)
    mtext('PH', side=4, line=1, las=2, cex = 1.2)}
  if (i == 4) {
    mtext('Non-PH', side=4, line=1, las=2, cex = 1.2)
  }
  if (i == 6){
    mtext('Combined', side=4, line=1, las=2, cex = 1.2)
  }
  if(i == 3){
    mtext('log-likelihood difference', outer=TRUE, side=2, line=3, cex = 1.2)
  }
}
dev.off()
