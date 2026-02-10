
library("dplyr")
library("tidyverse")
library("tidymodels")
library("data.table")


###################
# script parameters
###################
base_folder <- "~/Documents/tania/vertrebrates_invertebrates/"
infile = "cleaned_kmer_frequencies.csv"
paramfile = "best_hyperparameters_tuning_lasso.csv"
outdir = "results"
method = "lasso" ## "rf", "lasso"
class_weights = TRUE

## preprocessing
collinearity_threshold = 0.95
nproc <- 4

### write out chosen parameters
print("###################################")
print(paste("Parameter file is:", paramfile))
print(paste("Method is :", method))
temp = ifelse(class_weights, "Yes", "No")
print(paste("Class weights:", temp))
print("###################################")

# create output folder
fname = file.path(base_folder, outdir)
if (!dir.exists(fname)) {
  dir.create(fname, recursive = TRUE, showWarnings = FALSE)
}

## input data
writeLines(" - reading input data")

fname = file.path(base_folder, infile)
inpdata <- fread(fname)

## model parameters from fine-tuning
fname = file.path(base_folder, outdir, paramfile)
params = fread(fname)
print("Loaded model parameters (from fine-tuning):")
print(params)

## 
writeLines(" - residual preprocessing ... ")

preprocess_recipe <- recipe(Outcome ~ ., data = inpdata) %>%
  update_role(Species, new_role = "Species") %>%  # Species will not be used as predictor
  step_zv(all_numeric(), -all_outcomes()) %>%
  step_corr(all_predictors(), threshold = collinearity_threshold)

# Calculate weights
if (class_weights) {
  
  class_freq <- table(inpdata$Outcome)
  class_weights <- sum(class_freq) / (length(class_freq) * class_freq)
  weight_vector = c(as.numeric(class_weights["invertebrate"]), as.numeric(class_weights["vertebrate"]))
  
} else weight_vector =c(1,1)

print("class weights")
print(weight_vector)

if (method == "lasso") {
  
  mod_spec <- logistic_reg(
    penalty = params$lambda, 
    mixture = 1 ## 1: lasso; 0: ridge
  ) |>
    set_mode("classification") |>
    set_engine("glmnet",
               classwt = weight_vector)
} else if (method == "rf") {
  
  mod_spec <- rand_forest(
    mtry = params$mtry,
    trees = params$ntrees,
    min_n = params$min_n
  ) %>%
    set_mode("classification") %>%
    set_engine("randomForest", 
               classwt = weight_vector)
}

print("model specification")
print(mod_spec)

model_wf <- workflow() %>%
  add_recipe(preprocess_recipe) %>%   
  add_model(mod_spec)

print("workflow")
print(model_wf)

writeLines(" - fitting the model to the training data")
model_fit <- fit(model_wf, inpdata)

writeLines(" - saving results to file")
fname = file.path(base_folder, outdir, "trained_model.RData")
save(model_fit, file = fname)

print("DONE!!")
