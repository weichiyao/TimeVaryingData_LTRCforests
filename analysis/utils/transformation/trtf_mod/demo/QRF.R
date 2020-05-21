
library("quantregForest")
set.seed(290875)

pdf("QRF.pdf", width = 12, height = 8)

ntree <- 100
n <- 10000
p <- 10
x <- runif(n)
y <- rnorm(n, sd = 1 + (x > .5))
X <- matrix(runif(n * p), ncol = p)

qrf <- quantregForest(y = y, x = data.frame(x = x, X), mtry = p + 1, 
                      nodesize = 25, ntree = ntree, replace = FALSE)

xn <- 1:50 / 51
nX <- matrix(.5, nrow = length(xn), ncol = ncol(X))
nd <- data.frame(x = xn, nX)
qQRF <- predict(qrf, newdata = nd, what = c(.1, .9))

library("trtf")
var_y <- numeric_var("y", support = c(-5, 5))
B_y <- Bernstein_basis(var_y, order = 5, ui = "increasing")
m_y <- ctm(B_y)
trf <- traforest(m_y, formula = y ~ ., data = data.frame(y = y, x = x, X), 
                 ntree = ntree,
                 control = ctree_control(mincriterion = 0, 
                     minsplit = 25, minbucket = 10), 
                 mtry = p + 1, trace = TRUE)

trt <- trafotree(m_y, formula = y ~ ., data = data.frame(y = y, x = x, X))

qTRT <- predict(trt, newdata = nd, type = "quantile", prob = c(.1, .9))
qTRT <- t(qTRT)

qTRF <- list(p1 = predict(trf, newdata = nd, type = "quantile", prob = .1),
             p2 = predict(trf, newdata = nd, type = "quantile", prob = .9))

lwd <- 1.5
col <- rgb(.1, .1, .1, .05)
colR <- rgb(.75, 0, 0, .8)
colRlight <- rgb(.75, 0, 0, .1)
colB <- rgb(0, 0, 0.75, .8)

q1 <- qnorm(.1, sd = 1 + (x > .5))
q9 <- qnorm(.9, sd = 1 + (x > .5))
plot(y ~ x, pch = 19, col = c(col, colRlight)[(y < q1 | y > q9) + 1], cex = .5)
lines(xn, qQRF[,1], lty = 1, lwd = lwd * 1.5, col = "black")
lines(xn, qQRF[,2], lty = 1, lwd = lwd * 1.5, col = "black")
lines(xn, qTRT[,1], lty = 2, lwd = lwd * 1.5, col = colB)
lines(xn, qTRT[,2], lty = 2, lwd = lwd * 1.5, col = colB)
lines(xn, unlist(qTRF[[1]]), lty = 1, lwd = lwd * 1.5, col = colB)
lines(xn, unlist(qTRF[[2]]), lty = 1, lwd = lwd * 1.5, col = colB)
legend("topleft", lty = c(1, 2, 1), lwd = c(lwd * 1.5, lwd * 1.5, lwd * 1.5),
       col = c("black", colB, colB),
       legend = c("Quantile Regression Forest", "Transformation Tree", "Transformation Forest"),
       bty = "n")

dev.off()
