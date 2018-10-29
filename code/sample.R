rm(list = ls(all = TRUE))
source("model.R")
data(Boston)
y = Boston$medv
x = Boston[, c(1, 8)]
fit = regression(y = y, x = x)
for(j in 1:(length(x) + 1)){
    cat("beta", j - 1, ": ", median(fit$beta[j, ]), "\n")
}
