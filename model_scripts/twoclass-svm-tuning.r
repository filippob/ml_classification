library("repr")
library("dials")
library("readxl")
library("themis")
library("lobstr")
library("butcher")
library("ggplot2")
library("kernlab")
library("tidyverse")
library("tidymodels")
library("data.table")

## Parameters
basefolder <- "/home/filippo/Documents/tania/probiotics"
input_file <- "splits/train_set.csv"
outdir = "results"
# input_file <- "data/subset.csv"
nproc <- 4
split_ratio <- 0.75
method_cv <- "repeatedcv"
k_folds <- 5
nrepeats_cv <- 5
nlevels = 10 ## levels of hyperparameters to test
target_var = "Label"
id_vars = c("Organism", "Taxon", "Definition")

## Import dataset
writeLines(" - reading the data ...")
fname <- file.path(basefolder, input_file)
ext = gsub("^.*\\.","",fname)
if (ext %in% c("xlsx")) {
  
  dataset <- readxl::read_xlsx(fname)
} else dataset <- fread(fname)

print(paste("Data size (n. records):",nrow(dataset)))

## Data preparation
dataset <- dataset %>% 
  mutate({{target_var}} := as.factor(.data[[target_var]]))

## Setup parallel backend
writeLines(" - parallelise on the number of declared CPUs")
cl <- parallel::makeCluster(nproc)
doParallel::registerDoParallel(cl)

## training/test split
writeLines(" - splitting the data in training/test stes")
dt_split <- initial_split(dataset, strata = !!target_var, prop = split_ratio)
dt_train <- training(dt_split)
dt_test <- testing(dt_split)

dt_train$Label |> table()
dt_test$Label |> table()

#### Preprocessing

# We use tidymodels to build a recipe for data preprocessing:
#   
# -   remove correlated variables
# -   remove non informative variables (zero variance)
# -   standardize all variables
# -   impute missing data (Random Forest does not handle missing data)

writeLines(" - preprocessing ... ")
svm_recipe <- dt_train %>%
  recipe(reformulate(".", response = target_var)) %>%
  update_role(all_of(id_vars), new_role = "ID") |>
  step_corr(all_predictors(), threshold = 0.95) %>%
  step_zv(all_predictors(), -all_outcomes()) %>%
  step_nzv(all_predictors()) |>
  step_normalize(all_numeric(), -all_outcomes()) |>
  step_novel(all_nominal_predictors()) |>
  step_unknown(all_nominal_predictors()) |>
  step_upsample(all_outcomes(), over_ratio = 1)

prep_dt <- prep(svm_recipe)
print(prep_dt)

training_set <- juice(prep_dt)
table(training_set$Label)

#### Model building

# We now specify the structure of our model:
#   
# -   the problem at hand (classification)
# -   the engine (R package)
# 
# Then we put this in a workflow together with the preprocessing recipe
writeLines(" - fine-tuning of the hyperparameters ...")

svm_spec <-
  svm_rbf(cost = tune(), rbf_sigma = tune()) %>%
  set_mode("classification") %>%
  set_engine("kernlab")

tune_wf <- workflow() %>%
  add_recipe(svm_recipe) |>
  add_model(svm_spec)

#### Tuning the hyperparameters

# We use k-fold cross-validation to tune the hyperparameters in the training set
dt_folds <- vfold_cv(dt_train, v = k_folds, repeats = nrepeats_cv, strata = !!target_var)

## grid of hyperparameter values
svm_hyperparams <- parameters(
  cost(range = c(-4, 4)),
  rbf_sigma(range = c(-5, -3))
)
svm_grid <- grid_regular(svm_hyperparams, levels = nlevels)

head(svm_grid)
nrow(svm_grid)

fine_tune_res <- tune_grid(
  tune_wf,
  metrics = metric_set(roc_auc, accuracy, mcc),
  resamples = dt_folds,
  grid = svm_grid
)

tuning_res <- fine_tune_res |>
  collect_metrics()

print(tuning_res)

## from "repr"
options(repr.plot.width=14, repr.plot.height=8)

g <- fine_tune_res %>%
  collect_metrics() %>%
  filter(.metric == "mcc") %>%
  mutate(rbf_sigma = factor(rbf_sigma)) %>%
  ggplot(aes(cost, mean, color = rbf_sigma)) +
  geom_line(alpha = 0.5, size = 1.5) +
  geom_point() +
  labs(y = "MCC")

print(g)

## from "butcher"
obj_size(fine_tune_res)
weigh(fine_tune_res)

cleaned_finetune <- butcher(fine_tune_res, verbose = TRUE)
obj_size(cleaned_finetune)
weigh(cleaned_finetune)

writeLines(" - saving results")
## tuned model
fname <- file.path(basefolder, outdir, "twoclass_tuned_model.RData")
to_save = list(cleaned_finetune)
to_save[[2]] = svm_spec
to_save[[3]] = svm_recipe
to_save[[4]] = dt_split
save(to_save, file = fname)

## tuning results
fname <- file.path(basefolder, outdir, "twoclass_tuning_results.csv")
fwrite(x = tuning_res, file = fname, sep = ",")

## tuning plot
fname <- file.path(basefolder, outdir, "twoclass_tuning.png")
ggsave(filename = fname, plot = g, device = "png", width = 7, height = 6)

print("DONE!!")

