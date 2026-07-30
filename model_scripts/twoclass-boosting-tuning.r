library("repr")
library("dials")
library("glmnet")
library("readxl")
library("themis")
library("lobstr")
library("GGally")
library("xgboost")
library("butcher")
library("ggplot2")
library("kernlab")
library("tidyverse")
library("tidymodels")
library("data.table")

## Parameters
tmstmp = as.integer(Sys.time())
basefolder <- "/home/filippo/Documents/tania/probiotics"
input_file <- "splits/train_set.csv"
outdir = "results/twoclass/boosting"
nproc <- 4
method_cv <- "repeatedcv"
k_folds <- 5
nrepeats_cv <- 3
nlevels = 25 ## levels of hyperparameters to test
target_var = "Label"
id_vars = c("Organism", "Taxon", "Definition")


## log file
fname = paste(tmstmp, ".twoclass-boosting-tuning.r.log", sep="")
logn = file.path(basefolder, "log", fname)
fileConn <- file(logn)
data_list = list("timestamp"=tmstmp, "basefolder"=basefolder, "input_file"=input_file, "outdir"=outdir, 
                 "nproc"=nproc, "method_cv"=method_cv, "k_folds"=k_folds, "nrepeats_cv"=nrepeats_cv, 
                 "nlevels"=nlevels, "target_var"=target_var, "id_vars"=id_vars)
lines <- paste0(names(data_list), ": ", unlist(data_list))
writeLines(lines, fileConn)
close(fileConn)

## make timestamp directory
dir.create(file.path(basefolder, outdir, tmstmp), showWarnings = FALSE, recursive = TRUE)

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
# writeLines(" - splitting the data in training/test stes")
# dt_split <- initial_split(dataset, strata = !!target_var, prop = split_ratio)
# dt_train <- training(dt_split)
# dt_test <- testing(dt_split)
# 
# dt_train$Label |> table()
# dt_test$Label |> table()

dt_train <- dataset
rm(dataset)

#### Preprocessing

# We use tidymodels to build a recipe for data preprocessing:
#   
# -   remove correlated variables
# -   remove non informative variables (zero variance)
# -   standardize all variables
# -   impute missing data (Random Forest does not handle missing data)

writeLines(" - preprocessing ... ")
model_recipe <- dt_train %>%
  recipe(reformulate(".", response = target_var)) %>%
  update_role(all_of(id_vars), new_role = "ID") |>
  step_zv(all_predictors(), -all_outcomes()) %>%
  step_nzv(all_predictors()) |>
  step_corr(all_predictors(), threshold = 0.95) %>%
  step_normalize(all_numeric(), -all_outcomes()) |>
  step_novel(all_nominal_predictors()) |>
  step_unknown(all_nominal_predictors()) |>
  step_upsample(all_outcomes(), over_ratio = 1)

prep_dt <- prep(model_recipe)
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

## n. of variables (features) in the dataset
m = (ncol(training_set) - 4)

tune_spec <- boost_tree(
  trees = tune(),
  tree_depth = tune(), 
  min_n = tune(),
  loss_reduction = 0,                     ## first three: model complexity
  sample_size = tune(), 
  mtry = floor(sqrt(m)),         ## randomness
  learn_rate = tune()                          ## step size
) %>%
  set_engine("xgboost") %>%
  set_mode("classification")

tune_wf <- workflow() %>%
  add_recipe(model_recipe) |>
  add_model(tune_spec)

#### Tuning the hyperparameters

# We use k-fold cross-validation to tune the hyperparameters in the training set
dt_folds <- vfold_cv(dt_train, v = k_folds, repeats = nrepeats_cv, strata = !!target_var)

## grid of hyperparameter values
hyper_grid <- grid_space_filling(
  trees(range = c(500,1000)),
  tree_depth(),
  min_n(),
  # loss_reduction = 0,
  # loss_reduction(),
  sample_size = sample_prop(range = c(0.25,0.75)),
  learn_rate(),
  size = nlevels
)

head(hyper_grid)
nrow(hyper_grid)

fine_tune_res <- tune_grid(
  tune_wf,
  metrics = metric_set(roc_auc, accuracy, mcc),
  resamples = dt_folds,
  grid = hyper_grid
)

tuning_res <- fine_tune_res |>
  collect_metrics() |>
  filter(!is.na(mean))

print(tuning_res)

collect_metrics(fine_tune_res) %>%
  filter(.metric == "roc_auc") %>%
  ggplot(
    aes(log10(learn_rate), mean)
  ) +
  geom_point() +
  geom_smooth(se = FALSE)


## from "repr"
options(repr.plot.width=14, repr.plot.height=8)

dd <- fine_tune_res %>%
  collect_metrics() %>%
  filter(.metric == "mcc", !is.na(mean)) %>%
  select(c(trees, min_n, tree_depth, learn_rate, sample_size, mean)) |>
  mutate(learn_rate = log10(learn_rate)) |>
  gather(key = "hyperparameter", value = "value", -mean) |>
  rename(mcc = mean) 


g <- ggplot(dd, aes(x = value, y = mcc)) +
  geom_point(aes(color = hyperparameter)) + 
  facet_wrap(~hyperparameter, scales = "free_x") + 
  xlab(NULL)

print(g)

## from "butcher"
obj_size(fine_tune_res)
weigh(fine_tune_res)

cleaned_finetune <- butcher(fine_tune_res, verbose = TRUE)
obj_size(cleaned_finetune)
weigh(cleaned_finetune)

writeLines(" - saving results")
## tuned model
fname <- file.path(basefolder, outdir, tmstmp, "twoclass_tuned_model.RData")
to_save = list(cleaned_finetune)
to_save[[2]] = tune_spec
to_save[[3]] = model_recipe
to_save[[4]] = tmstmp
to_save[[5]] = m
save(to_save, file = fname)

## tuning results
fname <- file.path(basefolder, outdir, tmstmp, "twoclass_tuning_results.csv")
fwrite(x = tuning_res, file = fname, sep = ",")

## tuning plot
fname <- file.path(basefolder, outdir, tmstmp, "twoclass_tuning.png")
ggsave(filename = fname, plot = g, device = "png", width = 7, height = 6)

print("DONE!!")

