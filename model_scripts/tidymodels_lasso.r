library("vip")
library("glmnet")
library("ggplot2")
library("tidyverse")
library("tidymodels")
library("data.table")

## Parameters
basefolder <- "/home/filippo/Documents/tania/kmers"
input_file <- "data/dataset_v_i.csv"
nproc <- 4
split_ratio <- 0.80
method_cv <- "repeatedcv"
k_folds <- 10
nrepeats_cv <- 5
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
fname <- file.path(basefolder, "correlation_plot.png")
png(fname)
corrplot(V)
dev.off()

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

tune_spec <- logistic_reg(
  mode = "classification",
  penalty = tune(),
  mixture = 1      # α = 1 (Lasso)
) %>%
  set_engine("glmnet")

tune_wf <- workflow() %>%
  add_formula(Outcome ~ .) %>%
  add_model(tune_spec)

#### Tuning the hyperparameters

# We use k-fold cross-validation to tune the hyperparameters in the training set

trees_folds <- vfold_cv(training_set, v = k_folds, repeats = nrepeats_cv)

lasso_grid <- grid_regular(
  penalty(range = c(-7, -0.5)),
  levels = 10
)

head(lasso_grid)
nrow(lasso_grid)

regular_res <- tune_grid(
  tune_wf,
  metrics = metric_set(roc_auc, accuracy, mcc),
  resamples = trees_folds,
  grid = lasso_grid
)

regular_res |>
  collect_metrics() |>
  filter(.metric == "mcc") |>
  print(n = 20)

library("repr")
options(repr.plot.width=14, repr.plot.height=8)

regular_res %>%
  collect_metrics() %>%
  # filter(.metric == "mcc") %>%
  ggplot(aes(penalty, mean, color = .metric)) +
  geom_line(alpha = 0.5, size = 1.5) +
  geom_point() +
  labs(y = "MCC")

best_auc <- select_best(x = regular_res, metric = "mcc")
show_best(regular_res, metric = "mcc")

# 2.  finalise the model:
final_lasso <- finalize_model(
  tune_spec,
  best_auc
)

print(final_lasso)

# 3.  finalise the workflow and fit it to the initial split (training and test data):
final_wf <- workflow() %>%
  add_recipe(kmer_recipe) %>%
  add_model(final_lasso)

final_res <- final_wf %>%
  last_fit(kmer_split, metrics = metric_set(roc_auc, accuracy, mcc, brier_class))

# 4.  evaluate the fine-tuned RF model:
print(final_res)
final_res %>%
  collect_metrics()


# 5.  get variable importance:
final_res %>% 
  pluck(".workflow", 1) %>%   
  extract_fit_parsnip() %>% 
  #vip(num_features = 20, geom = "point")
  vip()

#### Predictions

# We collect the predictions on the test set: for each test observations we get the probabilities of belonging to each of the four classes.
final_res %>%
  collect_predictions()

cm <- final_res %>%
  collect_predictions() %>%
  conf_mat(Outcome, .pred_class)

print(cm)

fname <- file.path(basefolder, "confusion_matrix.png")
png(fname)
autoplot(cm, type = "heatmap")
dev.off()

print("DONE!!")




