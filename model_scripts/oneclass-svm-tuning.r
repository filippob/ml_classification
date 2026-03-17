library("repr")
library("dplyr")
library("dials")
library("purrr")
library("readxl")
library("themis")
library("lobstr")
library("butcher")
library("ggplot2")
library("kernlab")
library("discrim")
library("rsample")
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
normal_class <- "Probiotic"
target_var = "Label"
to_remove = c("Organism", "Taxon", "Definition")

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

if(sum(is.na(dataset)) == 0) {
  print("No missing data in the dataset: OK!")
  } else writeLines("###########################################\n 
                    !! Warning: there are missing data in your dataframe: check if this may be a problem !! 
                    \n ###########################################")

##################################
## CHECK WHAT IS MISSING
# vec <- colSums(is.na(dataset)) > 0
# dataset[,vec]
# 
# vec <- rowSums(is.na(dataset)) > 0
# x <- dataset[vec,]
##################################

dataset |>
  group_by(Label) |>
  summarise(N = n())

## Setup parallel backend
writeLines(" - parallelise on the number of declared CPUs")
cl <- parallel::makeCluster(nproc)
doParallel::registerDoParallel(cl)

## training/test split
writeLines(" - splitting the data in training/test stes")
# dt_split <- initial_split(dataset, strata = !!target_var, prop = split_ratio)
# dt_train <- training(dt_split)
# dt_test <- testing(dt_split)

dt_train_normal <- subset(dataset, Label == normal_class)
train_normal_class_label <- dataset$Label

dt_test <- subset(dataset, Label != normal_class)
# dt_test <- dt_test |> bind_rows(subset(dt_train, Label != normal_class))
test_label <- dt_test$Label

dt_train_normal$Label <- NULL
dt_test$Label <- NULL

dt_train_normal <- dt_train_normal |> select(!all_of(to_remove))
dt_test_reduced <- dt_test |> select(!all_of(to_remove))

#### Preprocessing

# We use tidymodels to build a recipe for data preprocessing:
#   
# -   remove correlated variables
# -   remove non informative variables (zero variance)
# -   standardize all variables
# -   impute missing data (Random Forest does not handle missing data)

writeLines(" - preprocessing ... ")
svm_recipe <- dt_train_normal %>%
  recipe(~ ., data = dt_train_normal) %>%
  # update_role(all_of(id_vars), new_role = "ID") |>
  step_corr(all_predictors(), threshold = 0.95) %>%
  # step_zv(all_predictors(), -all_outcomes()) %>%
  step_nzv(all_predictors()) |>
  step_normalize(all_numeric(), -all_outcomes()) |>
  step_novel(all_nominal_predictors()) |>
  step_unknown(all_nominal_predictors())

prep_dt <- prep(svm_recipe)
print(prep_dt)

training_set <- juice(prep_dt)

#### Model building

# We now specify the structure of our model:
#   
# -   the problem at hand (classification)
# -   the engine (R package)
# 
# Then we put this in a workflow together with the preprocessing recipe
writeLines(" - fine-tuning of the hyperparameters ...")

## 1) grid of values for the hyperparameters to tune
nu_vals <- seq(0.005, 0.2, by = 0.02)
sigma_vals <- 2^seq(-12, -2, length.out = 10)

grid <- expand.grid(
  nu = nu_vals,
  sigma = sigma_vals
)

## 2) cv folds
folds <- vfold_cv(training_set, v = k_folds, repeats = nrepeats_cv)

## 3) scoring function (one-class only in the training data)
score_model <- function(train_data, test_data, nu, sigma) {
  
  tryCatch({
    
    model <- ksvm(
      x = as.matrix(train_data),
      type = "one-svc",
      kernel = rbfdot(sigma = sigma),
      nu = nu,
      scaled = TRUE
    )
    
    pred <- predict(model, as.matrix(test_data))
    
    mean(pred)
    
  }, error = function(e) {
    NA_real_
  })
  
}

## 4) actual fine-tuning
results <- grid %>%
  mutate(score = map2_dbl(nu, sigma, function(nu, sigma) {
    
    fold_scores <- map_dbl(folds$splits, function(split) {
      
      train_data <- analysis(split)
      test_data  <- assessment(split)
      
      score_model(train_data, test_data, nu, sigma)
      
    })
    
    print(paste("Total n. of folds (folds x repeates):",  length(fold_scores)))
    mean(fold_scores, na.rm = TRUE)
    
  }))


## 5) check results
g <- ggplot(results, aes(x = sigma, y = score)) + geom_point() + facet_wrap(~nu)
print(g)

## tuning plot
fname <- file.path(basefolder, outdir, "oneclass_tuning.png")
ggsave(filename = fname, plot = g, device = "png", width = 7, height = 6)

## 6) select best hyperparameters
best_params <- results %>%
  arrange(desc(score)) %>%
  slice(1)

print(best_params)

writeLines(" - saving results")
## tuned model
fname <- file.path(basefolder, outdir, "oneclass_tuned_model.RData")
to_save = list(results)
to_save[[2]] = svm_recipe
save(to_save, file = fname)

print("DONE!!")




