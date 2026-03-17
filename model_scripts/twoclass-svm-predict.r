
library("themis")
library("yardstick")
library("tidymodels")
library("data.table")

## Parameters
basefolder <- "/home/filippo/Documents/tania/probiotics"
tuned_model <- "results/twoclass_tuned_model.RData"
train_set = "splits/train_set.csv"
test_set = "splits/test_set.csv"
outdir = "results"
nproc <- 4
id_vars = c("Organism", "Taxon", "Definition")
positive_class <- "Probiotic"
target_var = "Label"
flag_manual = TRUE ## for manual explicit workflow

## load tuned model
writeLines(" - loading tuned model")
fname = file.path(basefolder, tuned_model)
load(fname)

fine_tune_res <- to_save[[1]]
svm_spec = to_save[[2]]
svm_recipe = to_save[[3]]
dt_split = to_save[[4]]

## read training and test data
writeLines(" - read training data")
fname = file.path(basefolder, train_set)
train <- fread(fname)

writeLines(" - read test data")
fname = file.path(basefolder, test_set)
test <- fread(fname)

writeLines(" - show best model")
best_model <- select_best(x = fine_tune_res, metric = "mcc")
show_best(fine_tune_res, metric = "mcc")

# 2.  finalise the model:
writeLines(" - finalise model")

## Data preparation
train <- train %>% 
  mutate({{target_var}} := factor(.data[[target_var]]))

prepped_rec <- prep(svm_recipe)

final_svm <- finalize_model(
  svm_spec,
  best_model
)

# 3.  finalise the workflow and fit it to the initial split (training and test data):
## final workflow for model accuracy

if (flag_manual) {
  
  ## with no data augmentation
  training_set <- bake(prepped_rec, new_data = train)
  training_set <- training_set |> select(!all_of(id_vars))
  # table(training_set$Label)
  
  final_wf <- workflow() %>%
    add_formula(reformulate(".", response = target_var)) %>%
    add_model(final_svm)
  
  print(final_wf)
  
  final_res <- final_wf |> fit(data = training_set)
  
} else {
  
  final_wf <- workflow() %>%
    add_recipe(svm_recipe) %>%
    add_model(final_svm)
  
  # training_set <- juice(prepped_rec)
  # table(training_set$Label)
  
  final_res <- final_wf |> fit(data = train)
}

print(final_res)

final_res %>%
  extract_fit_parsnip()

# final_res <- final_wf %>%
#   last_fit(dt_split, metrics = metric_set(roc_auc, accuracy, mcc, brier_class))


## evaluate on test
test <- test %>% 
  mutate({{target_var}} := factor(.data[[target_var]]))


if (flag_manual) {
  
  test_set <- bake(prepped_rec, new_data = test)
  test_set <- test_set |> select(!all_of(id_vars))
  
  preds = predict(final_res, test_set, type="prob")
  temp <- test_set |> select(!!target_var) |> pull()
  preds$.pred_class = colnames(preds)[max.col(preds)]
  preds$.pred_class = gsub(".pred_","",preds$.pred_class)
  preds$.pred_class = as.factor(preds$.pred_class)
  preds <- preds |> mutate({{target_var}} := temp)
  
  auc = roc_auc(preds, truth = Label, .pred_Nonprobiotic)
  acc = accuracy(preds, truth = Label, estimate = .pred_class)
  mcc = mcc(preds, truth = Label, estimate = .pred_class)
  brier_score = brier_class(preds, truth = Label, .pred_Probiotic)
  
} else {
  
  preds = predict(final_res, test, type="prob")
  temp <- test |> select(!!target_var) |> pull()
  preds <- preds |> mutate({{target_var}} := temp)
  preds <- preds |> bind_cols(predict(final_res, test, type="class"))
  
  auc = roc_auc(preds, truth = Label, .pred_Nonprobiotic)
  acc = accuracy(preds, truth = Label, estimate = .pred_class)
  mcc = mcc(preds, truth = Label, estimate = .pred_class)
  brier_score = brier_class(preds, truth = Label, .pred_Probiotic)
}

df_metrics <- bind_rows(acc,auc,mcc,brier_score)

print("Performance metrics")
print(df_metrics)

fname <- file.path(basefolder, outdir, "twoclass-metrics.csv")
fwrite(x = df_metrics, file = fname)

print("Confusion Matrix")
cm <- preds |>
  conf_mat(!!target_var, .pred_class)

print(cm)

fname <- file.path(basefolder, outdir, "twoclass-confusion_matrix.png")
g <- autoplot(cm, type = "heatmap") + scale_fill_gradientn(colours = c("lightyellow", "yellow", "orange", "red"))
print(g)

ggsave(filename = fname, plot = g, device = "png", width = 5, height = 4.5)

writeLines(" - error analysis")

preds <- test |>
  select(Label, Organism, Taxon, Definition) |>
  rename(Label_orig = Label) |>
  bind_cols(preds)

print(paste("N. of mismatches between labels from the test set and from collected predictions", sum(preds$Label != preds$.pred_class)))

errors <- preds |>
  filter(Label != .pred_class)

fname = file.path(basefolder, outdir, "twoclass-all_predictions.csv")
fwrite(x = preds, file = fname, sep = "\t")

fname = file.path(basefolder, outdir, "twoclass-errors.csv")
fwrite(x = errors, file = fname, sep = "\t")

print("DONE!!")
