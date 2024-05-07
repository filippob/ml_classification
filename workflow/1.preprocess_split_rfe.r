
library("caret")
library("tidyverse")
library("data.table")
library("doParallel")
library("randomForest")

## PARAMETERS
args = commandArgs(trailingOnly=TRUE)
if (length(args) >= 1) {
  
  #loading the parameters
  source(args[1])
  # source("Analysis/hrr/config.R")
  
} else {
  #this is the default configuration, used for development and debug
  writeLines('Using default config')
  
  #this dataframe should be always present in config files, and declared
  #as follows
  config = NULL
  config = rbind(config, data.frame(
    basefolder = "Documents/tania/kmers",
    datafile = "data/genome_characteristics_vertebrate_invertebrate_real.csv",
    resfolder = "Analysis",
    suffix = "1.preprocess_split_rfe",
    split = 0.8, ## proportion of data for training
    folds = 10,
    nrepeates = 1,
    ncpu = 4, ## n. of cores to use for parallel computations (if any)
    force_overwrite = FALSE
  ))
}

HOME <- Sys.getenv("HOME")
basefolder = file.path(HOME, config$basefolder)
outdir = file.path(basefolder,config$resfolder)

writeLines(" - saving configuration parameters to log file")
fname = paste(config$suffix, ".config.r", sep="")
fname = file.path(outdir, fname)
fwrite(x = config, file = fname)

print(paste("saved configuration to file", fname))

# Import dataset
writeLines(" - importing data")

fname = file.path(basefolder, config$datafile)
print(paste("reading data from file", fname))

kmers <- fread(fname)
str(kmers) |> print()

### preprocessing
writeLines(" - data preprocessing and environment setting")

temp <- kmers %>%
  mutate(Outcome = as.factor(Outcome)) |>
  select(-ID)

temp %>%
  count(Outcome) |>
  print()


if (config$ncpu > 1) {
  
  print(paste("parallelising on", config$ncpu, "cores"))
  cl <- makeCluster(config$ncpu)
  registerDoParallel(cl)
}

writeLines(" - splitting data in training and test sets")
partition_index = createDataPartition(temp$Outcome, p=config$split, list=FALSE)

# Select 20% of the data for validation
validation <- temp[-partition_index,]

# Use the remaining 80% of data to train and test the models
training <- temp[partition_index,]

print("distribution of target variable in the TRAINING and TEST sets")

training %>%
  count(Outcome) |>
  print()

validation %>%
  count(Outcome) |> 
  print()

to_save <- list("validation"=validation, "training"=training)


# Feature selection (rfe): automatically selecting a subset of the most predictive
# features on train/test dataset
writeLines(" - feature selection")

# 10-fold cross validation (cv) repeated 1000 times
# Use different seed for cv in rfe and train
rfectrl_rf <- rfeControl(functions=rfFuncs, method="repeatedcv", number=config$folds, repeats=config$nrepeates,
                         verbose = FALSE, returnResamp = "final", allowParallel = TRUE)


nfeatures = ncol(training)-1
results_rfe <- rfe(x = training[,-1], y = training$Outcome, sizes=c(1:nfeatures), rfeControl=rfectrl_rf)
print(results_rfe$fit)

to_save[["rfe"]] = results_rfe

# Plot the results
writeLines(" - plotting results of RFE")
p <- plot(results_rfe, type=c("g", "o"))

to_save[["plot_rfe"]] = rfe

fname = "RFE.png"
dir.create(file.path(outdir, "figures"), showWarnings = FALSE, recursive = TRUE)
fname = file.path(outdir, "figures", fname)

png(filename=fname)
plot(results_rfe, type=c("g", "o"))
dev.off()

# List the chosen features
writeLines(" - selected features from RFE")
predictors(results_rfe) |> print()
to_save[["rfe_predictors"]] = predictors(results_rfe)

## save results to R object
writeLines(" - saving results and data to R object")
fname = paste(config$suffix, ".RData", sep="")
fname = file.path(outdir, fname)
save(to_save, file = fname)

print(paste("saved results to file", fname))

print("DONE!!")


