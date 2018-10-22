rm(list = ls(all = TRUE))
source("model.R")
data(Boston)
y = Boston$medv
x = Boston[, c(1, 8)]
fit = regression(y = y, x = x)
