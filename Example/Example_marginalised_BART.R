# my implementation in R of this marginalisation (quite slow!)
install_github("ebprado/ExtensionsBART/myBART")
library(myBART)

library(dbarts)
library(dplyr)

# Simulate data -------------------------------------------
f <- function(x) {
  10 * sin(pi * x[,1] * x[,2]) + 20 * (x[,3] - 0.5)^2 +
    10 * x[,4] + 5 * x[,5]
}

set.seed(99)
sigma <- 1.0
n     <- 200

x  <- matrix(runif(n * 10), n, 10)
Ey <- f(x)
y  <- rnorm(n, Ey, sigma)

# run BART --------------------------------------------------------
set.seed(99)
bartFit <- bart(x, y, ntree = 10, nskip = 100, ndpost = 13, keeptrees = TRUE)

# Get trees -------------------------------------------------------
all_trees = bartFit$fit$getTrees()

# Select only trees that use x4 (let's say we're interested in x4)
var_interest = 4

aux = all_trees %>% filter(var != -1) %>% group_by(sample, tree, var) %>% summarise(n=n())
aux = aux %>% group_by(sample,tree) %>% mutate(n=n(), keep=ifelse(var==var_interest & n==1, TRUE, FALSE))
aux = aux %>% select(sample, tree, keep)

trees_with_x4_only = left_join(all_trees,aux,by = c('sample','tree')) %>% filter(keep==TRUE)

trees_with_x4_only %>% group_by(sample,tree) %>% summarise(n=n())

# -----------------------------------------------------------------
## GAM ------------------------------------------------------------

library(mgcv)
dat = data.frame(y,x)
b <- gam(y~s(X1,X2)+s(X3)+s(X4)+s(X5),data=dat)

plot(b,pages=1,residuals=TRUE)  ## show partial residuals
# zoom in on x4
plot(b,pages=0,residuals=TRUE, select=3, main='Friedman data: GAM fit - marginal effect of x4')

# -----------------------------------------------------------------
## My implementation of this marginalisation (quite slow)
# -----------------------------------------------------------------
mybart = myBART::bart(x,y,ntrees = 200, nburn = 100, npost=100)

#----------
VAR_AUX = as.character(var_interest)
data_test = dat[,-1]
gam_bart_overall_prediction = myBART::predict_mybart(mybart, data_test, type = 'mean')
gam_bart_posterior_samples = myBART::marginal_predict_mybart(mybart, VAR_AUX, data_test, type = 'all')
gam_bart_pred = apply(gam_bart_posterior_samples$prediction,2,mean)
gam_bart_pred_sd = apply(gam_bart_posterior_samples$prediction,2,sd)


diff = mean(apply(mybart$y_hat, 2,mean)) - mean(gam_bart_pred)
gam_bart_pred = gam_bart_pred + diff

# residuals
y_scaled = (dat[,'y'] - mean(dat[,'y']))
x_aux = data_test[,var_interest]
y_hat_scaled = gam_bart_pred - mean(gam_bart_pred)
db = data.frame(y = dat[,'y'],
                y_scaled = y_scaled,
                x = x_aux,
                y_hat_scaled_marginalised = y_hat_scaled,
                res = y_scaled - gam_bart_overall_prediction,
                y_hat_overall = gam_bart_overall_prediction)
db = db %>% arrange(x)

#### plot
points(db$x, db$y_hat_scaled_marginalised, col=4)
auxx = (db$y - db$y_hat_overall)
points(db$x, db$y_hat_scaled_marginalised + auxx, col=4, cex = 0.5)
legend(0.1, 5, legend=c("GAM", "BART"), pch=c(NA,1), lty=c(1,0),
       col=c("black", "blue"), cex=1)

# bartFit$fit$state
# bartFit$fit$getTrees()
# bartFit$fit$printTrees()
# bartFit$fit$plotTree(1,1)
# bartFit$fit$setControl()
