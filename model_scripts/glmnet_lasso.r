## script to run lasso (ridge) models for regression (continuous target) with glmnet
## first stub: to be further developed

# Load libraries
library("glmnet")
library("ggplot2")
library("tidyverse")
library("doParallel")

options(expressions = 5e5)

##############
## PARAMETERS
##############
basefolder = "/home/filippo/Documents/tania/PI"
# fname = "PI/data.RData" ## input data
# outfile = "PI/model.RData"
fname = "dati/subset.RData" ## input data
outfile = "results/model.RData"
false_reps = 0 ## for subsetting: the higher this number, the fewer the SNPs retained in the analysis (set this to zero to use the whole dataset)
ncpus = 4
nfolds = 5 ## n. of folds to fine-tune lambda for lasso-penalised regression
alpha = 1 ## alpha = 1 ##(1 − α)/2 * ||β||^2 + α||β||_1): alpha = 1 --> lasso; alpha = 0 --> ridge
###################
###################

## Import datasets
writeLines(" - loading data")
inputname = file.path(basefolder, fname)
load(inputname)

writeLines(" - subsetting the data")
## subsetting the data (columnwise)
y = dataset$PIresINT
temp <- dataset[,-1]
vec <- sample(c(TRUE,rep(FALSE,false_reps)), size = ncol(temp), replace = TRUE)
subset <- bind_cols("Y" = y, temp[,vec, with = FALSE])
rm(temp, dataset)

## data cleaning (uncomment if needed)
writeLines(" - data cleaning")
res = colMeans(subset[,-1])  ## remove SNPs with average close to 0 or 2 (likely monomorphic or close to be)
res = c(TRUE, as.vector(!(res < 0.05 | res > 1.95))) 
subset = subset[,res, with = FALSE]

print(paste("N. of SNPs remaining after subsetting and cleaning", ncol(subset)-1))

###############################
## LASSO MODEL
##############################
print("Using glmnet")

## training / test split
writeLines(" - split data")

n_train = round(nrow(subset)*0.8, digits = 0)
vec = sample(1:nrow(subset), n_train)
train = subset[vec,]
test = subset[-vec,]

rm(subset)
gc()

# Setup parallel backend to use n. processors
cl <- makeCluster(ncpus)
registerDoParallel(cl)

## Cross validation (fine-tuning of lambda)
writeLines(" - fine-tuning lambda through k-fold CV")
print(paste("N. of folds used for CV:", nfolds))

y <- train |> pull(Y)
X <- train |>
  select(-Y) |>
  as.matrix()

tX <- test |>
  select(-Y) |>
  as.matrix()

cvfit <- cv.glmnet(X, y, alpha = 1, type.measure = "mse", nfolds = nfolds, family = "gaussian")
plot(cvfit)

lambdamin = cvfit$lambda.min
print(paste("Best lambda is", lambdamin))


#### fit the final model
writeLines(" - fit the final model")

#lambda can be provided, but is typically not and the program constructs a sequence. 
#Supplying a value of lambda overrides this. WARNING: use with care. 
#Avoid supplying a single value for lambda (for predictions after CV use predict() instead). 
#Supply instead a decreasing sequence of lambda values. 
#glmnet relies on its warms starts for speed, and its often faster to fit a whole path than compute a single fit.

lambdaseq = seq(lambdamin-5*(0.1*lambdamin), lambdamin+0.1*lambdamin, 0.1*lambdamin/5)
fit <- glmnet(X, y, lambda = lambdaseq, family = "gaussian", alpha = alpha)
plot(fit, label = TRUE)

print(fit)
lambda =  lambdamin
snps = coef(fit, s = lambda) |> row.names()
coeffs = as.vector(coef(fit, s = lambda)[,"s1"])
coeffs = data.frame("snp"=snps, "coef" = coeffs)
coeffs = filter(coeffs, abs(coef) > 0)

y_hat = predict(fit, newx = tX, s = lambda)
print(paste("correlation between y and y_hat:",cor(test$Y, y_hat[,1])))

y_hat2 = predict(cvfit, newx = tX, s = lambdamin)
print(paste("correlation between y and y_hat2:",cor(test$Y, y_hat2[,1])))

print("DONE!")
