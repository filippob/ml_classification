library("vip")
library("ggplot2")
library("kernlab")
library("tidyverse")
library("tidymodels")
library("data.table")

## Parameters
basefolder <- "/home/filippo/Documents/tania/kmers"
input_file <- "data/dataset_v_i.csv"
# input_file <- "data/subset.csv"
nproc <- 4
split_ratio <- 0.80
method_cv <- "repeatedcv"
k_folds <- 5
nrepeats_cv <- 1
positive_class <- "vertebrate"

## Import dataset
writeLines(" - reading the data ...")
fname <- file.path(basefolder, input_file)
dataset <- fread(fname)
head(dataset)

## Data preparation
dataset <- dataset %>% 
  mutate(Outcome = as.factor(Outcome), ID = as.factor(ID))

if(sum(is.na(dataset)) == 0) print("No missing data in the dataset: OK!")

## Setup parallel backend
writeLines(" - parallelise on the number of declared CPUs")
cl <- parallel::makeCluster(nproc)
doParallel::registerDoParallel(cl)

## training/test split
writeLines(" - splitting the data in training/test stes")
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

writeLines(" - preprocessing ... ")
kmer_recipe <- kmer_train %>%
  recipe(Outcome ~ .) %>%
  step_corr(all_predictors(), threshold = 0.99) %>%
  step_zv(all_numeric(), -all_outcomes()) %>%
  step_normalize(all_numeric(), -all_outcomes())

prep_kmer <- prep(kmer_recipe)
print(prep_kmer)

training_set <- juice(prep_kmer)
head(training_set)

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
  add_formula(Outcome ~ .) %>%
  add_model(svm_spec)

#### Tuning the hyperparameters

# We use k-fold cross-validation to tune the hyperparameters in the training set
trees_folds <- vfold_cv(training_set, v = k_folds, repeats = nrepeats_cv)
# trees_folds <- vfold_cv(kmer_train, v = k_folds, repeats = nrepeats_cv)

# We now try to start from $\sqrt{p}$ (classification problem)

svm_grid <- grid_regular(
  cost(range = c(-6, 1)),
  rbf_sigma(range = c(-6, -2)),
  levels = c(5,5)
)

head(svm_grid)
nrow(svm_grid)

fine_tune_res <- tune_grid(
  tune_wf,
  metrics = metric_set(roc_auc, accuracy, mcc),
  resamples = trees_folds,
  grid = svm_grid
)

fine_tune_res |>
  collect_metrics() |>
  print()

library("repr")
options(repr.plot.width=14, repr.plot.height=8)

fine_tune_res %>%
  collect_metrics() %>%
  filter(.metric == "mcc") %>%
  mutate(rbf_sigma = factor(rbf_sigma)) %>%
  ggplot(aes(cost, mean, color = rbf_sigma)) +
  geom_line(alpha = 0.5, size = 1.5) +
  geom_point() +
  labs(y = "MCC")


best_auc <- select_best(x = fine_tune_res, metric = "mcc")
show_best(fine_tune_res, metric = "mcc")

# 2.  finalise the model:
final_svm <- finalize_model(
  svm_spec,
  best_auc
)

print(final_svm)

# 3.  finalise the workflow and fit it to the initial split (training and test data):
## final workflow for model accuracy
final_wf <- workflow() %>%
  add_recipe(kmer_recipe) %>%
  add_model(final_svm)

final_res <- final_wf %>%
  last_fit(kmer_split, metrics = metric_set(roc_auc, accuracy, mcc, brier_class))

## final workflow for variable importance (right interface for DALEX::explain_tidymodels)
base_wf <- workflow() %>%
  add_formula(Outcome ~ .)

svm_fitted <- base_wf %>%
  add_model(final_svm) %>%
  fit(training_set)

# svm_fitted <- base_wf %>%
#   add_model(final_svm) %>%
#   fit(kmer_train)

# 4.  evaluate the fine-tuned RF model:
print(final_res)
final_res %>%
  collect_metrics()

# 5.  get variable importance:
library("DALEX")
library("DALEXtra")

## variable importance from training set (default)
temp <- training_set
temp$Outcome = as.numeric(temp$Outcome) - 1 
tmp <- select(temp, -Outcome)

## variable importance from entire dataset (uncomment to use)
# temp <- kmer_dt
# temp$Outcome = as.numeric(temp$Outcome) - 1 
# tmp <- select(temp, -Outcome)

# explainer_svm <- explain_tidymodels(svm_fitted, data = temp[,-17], y = temp$Outcome)
explainer_svm <- explain_tidymodels(svm_fitted, data = tmp, y = temp$Outcome)

vip <- model_parts(explainer = explainer_svm)

## plot
fname <- file.path(basefolder, "variable_importance.png")
png(fname)
plot(vip)
dev.off()

## varimp results
fname <- file.path(basefolder, "variable_importance.csv")
fwrite(x = vip, file = fname)

#####################################
## from: https://www.tmwr.org/explain
#####################################

# shap_val <- 
#   predict_parts(
#     explainer = explainer_svm, 
#     new_observation = kmer_test, 
#     type = "shap",
#     B = 20
#   )

# library("forcats")
# shap_val %>%
#   group_by(variable) %>%
#   mutate(mean_val = mean(contribution)) %>%
#   ungroup() %>%
#   mutate(variable = fct_reorder(variable, abs(mean_val))) %>%
#   ggplot(aes(contribution, variable, fill = mean_val > 0)) +
#   geom_col(data = ~distinct(., variable, mean_val), 
#            aes(mean_val, variable), 
#            alpha = 0.5) +
#   geom_boxplot(width = 0.5) +
#   theme(legend.position = "none") +
#   scale_fill_viridis_d() +
#   labs(y = NULL)


#### Predictions

# We collect the predictions on the test set: for each test observations we get the probabilities of belonging to each of the four classes.
final_res %>%
  collect_predictions() |>
  print()

cm <- final_res %>%
  collect_predictions() %>%
  conf_mat(Outcome, .pred_class)

print(cm)

fname <- file.path(basefolder, "confusion_matrix.png")
png(fname)
autoplot(cm, type = "heatmap")
dev.off()

print("DONE!!")

