library("vip")
library("ggplot2")
library("tidyverse")
library("tidymodels")
library("data.table")
library("randomForest")

## Parameters
basefolder <- "/home/filippo/Documents/tania/kmers"
input_file <- "data/dataset_v_i.csv"
nproc <- 4
split_ratio <- 0.80
method_cv <- "repeatedcv"
k_folds <- 10
# nrepeats_rfe <- 10
nrepeats_cv <- 1
sampling_method <- "up"
positive_class <- "vertebrate"

## Import dataset
fname <- file.path(basefolder, input_file)
dataset <- fread(fname)
head(dataset)

## Data preparation
dataset <- dataset %>% 
  mutate(Outcome = as.factor(Outcome), ID = as.factor(ID))

if(sum(is.na(dataset)) == 0) print("No missing data in the dataset: OK!")

## Setup parallel backend
cl <- parallel::makeCluster(nproc)
doParallel::registerDoParallel(cl)

kmer_dt <- select(dataset, -c(ID))
kmer_split <- initial_split(kmer_dt, strata = Outcome, prop = split_ratio)
kmer_train <- training(kmer_split)
kmer_test <- testing(kmer_split)

#### Preprocessing

# We use tidymodels to build a recipe for data preprocessing:
#   
# -   remove correlated variables
# -   remove non informative variables (zero variance)
# -   standardize all variables
# -   impute missing data (Random Forest does not handle missing data)

library("corrplot")
V <- cor(kmer_train[,-1])
corrplot(V)

kmer_recipe <- kmer_train %>%
  recipe(Outcome ~ .) %>%
  # step_corr(all_predictors(), threshold = 0.99) %>%
  step_zv(all_numeric(), -all_outcomes()) %>%
  step_normalize(all_numeric(), -all_outcomes())

prep_kmer <- prep(kmer_recipe)
print(prep_kmer)

training_set <- juice(prep_kmer)
head(training_set)

#### Model building

# We now specify the structure of our model:
#   
# -   hyperparameters to tune: `mtry` (number of features to sample for each tree) and `min_n` (minimum number of data points in a node to allow further splitting)
# -   number of trees in the forest
# -   the problem at hand (classification)
# -   the engine (R package)
# 
# Then we put this in a workflow together with the preprocessing recipe

tune_spec <- rand_forest(
  mtry = tune(),
  trees = 500,
  min_n = tune()
) %>%
  set_mode("classification") %>%
  set_engine("randomForest")

tune_wf <- workflow() %>%
  add_formula(Outcome ~ .) %>%
  add_model(tune_spec)

#### Tuning the hyperparameters

# We use k-fold cross-validation to tune the hyperparameters in the training set

trees_folds <- vfold_cv(training_set, v = k_folds, repeats = nrepeats_cv)

tune_res <- tune_grid(
  tune_wf,
  resamples = trees_folds,
  grid = 10 ## n. of tuning combinations
)

print(tune_res)

library("repr")
options(repr.plot.width=14, repr.plot.height=8)

tune_res %>%
  collect_metrics() %>%
  filter(.metric == "roc_auc") %>%
  select(mean, min_n, mtry) %>%
  pivot_longer(min_n:mtry,
               values_to = "value",
               names_to = "parameter"
  ) %>%
  ggplot(aes(value, mean, color = parameter)) +
  geom_point(show.legend = FALSE) +
  facet_wrap(~parameter, scales = "free_x") +
  labs(x = NULL, y = "AUC")


# We now try to start from $\sqrt{p}$ (classification problem)
m <- round(sqrt(ncol(training_set)-1),0)
print(m)
rf_grid <- grid_regular(
  mtry(range = c(m-3, m+6)),
  min_n(range = c(3, 9)),
  levels = c(8,4)
)

print(rf_grid)

regular_res <- tune_grid(
  tune_wf,
  resamples = trees_folds,
  grid = rf_grid
)


regular_res %>%
  collect_metrics() %>%
  filter(.metric == "roc_auc") %>%
  mutate(min_n = factor(min_n)) %>%
  ggplot(aes(mtry, mean, color = min_n)) +
  geom_line(alpha = 0.5, size = 1.5) +
  geom_point() +
  labs(y = "AUC")



