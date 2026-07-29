
library("dplyr")
library("ggpubr")
library("ggplot2")
library("data.table")

basefolder <- "/home/filippo/Documents/tania/probiotics"
problem = "oneclass"
method = "svm"
outdir = file.path("results", problem, method)
positive_class = "Probiotic"
negative_class = "Nonprobiotic"

runs <- list.files(file.path(basefolder, outdir))

### model fine-tuning
writeLines(" - parsing results from fine-tuning")
bestmod = data.frame(NULL)
for (d in runs) {
  
  temp = paste(problem, "_tuning_results.csv", sep="")
  fname = file.path(outdir, d, temp)
  print(fname)
  if (file.exists(fname)) {
    temp = fread(fname)
    if (problem == "twoclass") {
      temp = filter(temp, .metric == "mcc")
      temp = filter(temp, mean == max(mean))
    }
    temp$tmstmp = as.character(d)
    
    bestmod = rbind.data.frame(bestmod, temp) 
  }
}

### accuracy
writeLines(" - parsing performance metrics")
metrics = data.frame(NULL)
for (d in runs) {
  
  temp = paste(problem, "-metrics.csv", sep="")
  fname = file.path(outdir, d, temp)
  print(fname)
  temp = fread(fname)
  metrics = rbind.data.frame(metrics, temp)
}

### errors
writeLines(" - parsing predictions and errors")

res = data.frame(NULL)
errors = data.frame(NULL)
for (d in runs) {
  
  temp = paste(problem, "-all_predictions-test_set.csv", sep="")
  fname = file.path(outdir, d, temp)
  print(fname)
  preds = fread(fname)
  if (problem == "oneclass") preds <- rename(preds, .pred_class = predicted_label)
  
  ## positive labels
  temp <- filter(preds, Label == positive_class)
  fn = filter(temp, Label != .pred_class) |> nrow()
  tp = filter(temp, Label == .pred_class) |> nrow()
  
  ## negative labels
  temp <- filter(preds, Label == negative_class)
  fp = filter(temp, Label != .pred_class) |> nrow()
  tn = filter(temp, Label == .pred_class) |> nrow()
  
  temp = data.frame("TP"=tp, "TN"=tn, "FP"=fp, "FN"=fn)
  temp$tmstmp = as.integer(d)
  res = rbind.data.frame(res, temp)
  
  temp <- filter(preds, Label != .pred_class)
  temp$tmstmp = as.integer(d)
  errors = rbind.data.frame(errors, temp)
}

### discovery
writeLines(" - parsing new probiotics")

newprobs = data.frame(NULL)
for (d in runs) {
  
  temp = paste(problem, "-newprobiotics-filtered_merged_bits26_testset.csv", sep="")
  fname = file.path(outdir, d, temp)
  print(fname)
  probs = fread(fname)
  probs$tmstmp = as.integer(d)
  
  newprobs = rbind.data.frame(newprobs, probs)
}

### WRITE OUT ###
writeLines(" - writing out results")
## make timestamp directory
resfolder = file.path(basefolder, "combined_results")
dir.create(file.path(resfolder, problem, method), showWarnings = FALSE, recursive = TRUE)

fname = file.path(resfolder, problem, method, "bestmodel.csv")
print(fname)
fwrite(x = bestmod, file = fname)

fname = file.path(resfolder, problem, method, "metrics.csv")
print(fname)
fwrite(x = metrics, file = fname)

fname = file.path(resfolder, problem, method, "prediction_results.csv")
print(fname)
fwrite(x = res, file = fname)

fname = file.path(resfolder, problem, method, "errors.csv")
print(fname)
fwrite(x = errors, file = fname)

fname = file.path(resfolder, problem, method, "new_probiotics.csv")
print(fname)
fwrite(x = newprobs, file = fname)
