
######################################################
## script to make figures and tables for PROB-AI paper
######################################################
library("dplyr")
library("ggpubr")
library("ggplot2")
library("data.table")
library("tidymodels")

basefolder <- "/home/filippo/Documents/tania/probiotics"
resdir = "combined_results"
twoclass_results = "twoclass/svm"
oneclass_results = "oneclass/svm"
outdir = "figures_tables"

dir.create(file.path(basefolder, outdir), showWarnings = FALSE)

###################
## MODEL TUNING ###
###################
fname = file.path(basefolder, resdir, twoclass_results, "bestmodel.csv")
bestmodel = fread(fname)
bestmodel <- bestmodel |>
  select(cost, rbf_sigma, mean) |>
  gather(key = "hyperparameter", value = "value", -mean) |>
  mutate("problem" = "twoclass") |>
  rename(accuracy = mean)

fname = file.path(basefolder, resdir, oneclass_results, "bestmodel.csv")
temp = fread(fname)
temp <- temp |>
  select(-tmstmp) |>
  gather(key = "hyperparameter", value = "value", -score) |>
  mutate("problem" = "oneclass") |>
  rename(accuracy = score)

bestmodel <- rbind.data.frame(bestmodel, temp)

p <- ggplot(bestmodel, aes(x = value, y = accuracy)) + geom_jitter(aes(color=hyperparameter))
p <- p + facet_wrap(~problem, scales = "free")
p

fname = file.path(outdir, "hyperparameters_accuracy.png")
ggsave(filename = fname, plot = p, device = "png")

D <- bestmodel |>
  group_by(problem, hyperparameter, value) |>
  summarise(N = n()) |>
  arrange(value) |>
  mutate(value = as.character(round(value, 3)))

p <- ggplot(filter(D, problem == "oneclass"), aes(x = value, y = N))
p <- p + geom_bar(aes(fill=hyperparameter), stat="identity")
p <- p + facet_wrap(~hyperparameter, scales = "free")
# p

q <- ggplot(filter(D, problem == "twoclass"), aes(x = value, y = N))
q <- q + geom_bar(aes(fill=hyperparameter), stat="identity")
q <- q + facet_wrap(~hyperparameter, scales = "free")
# q

gg <- ggarrange(p,q, labels = c("one","two"), ncol = 1)

fname = file.path(outdir, "hyperparameters_frequency.png")
ggsave(filename = fname, plot = gg, device = "png")

##################################
## MODEL PERFORMANCE          ####
##################################
fname = file.path(basefolder, resdir, twoclass_results, "metrics.csv")
metrics = fread(fname)
metrics$problem = "twoclass"

fname = file.path(basefolder, resdir, oneclass_results, "metrics.csv")
temp = fread(fname)
temp$problem = "oneclass"

metrics <- rbind.data.frame(metrics, temp)

p <- ggplot(metrics, aes(x = .metric, y = .estimate))
p <- p + geom_jitter(aes(color=.metric))
p <- p + facet_wrap(~problem, scales="free")
p

fname = file.path(outdir, "model_performance.png")
ggsave(filename = fname, plot = p, device = "png")

##################################
## CONFUSION MATRIX           ####
##################################
fname = file.path(basefolder, resdir, twoclass_results, "prediction_results.csv")
preds = fread(fname)
preds$problem = "twoclass"

fname = file.path(basefolder, resdir, oneclass_results, "prediction_results.csv")
temp = fread(fname)
temp$problem = "oneclass"

preds <- rbind.data.frame(preds,temp)

preds <- preds |>
  group_by(problem) |>
  summarise(FP = mean(FP), FN = mean(FN), TP = mean(TP), TN = mean(TN))

## round floats
preds <- preds %>%
  mutate(across(where(is.numeric),
      ~ if_else(abs(.x - round(.x)) < .Machine$double.eps^0.5,
        .x,
        round(.x, 2)
)))

temp <- filter(preds, problem == "twoclass") 
temp <- temp |>
  gather(key = "pred", value = "count", -problem)

temp$label = c("other","probiotic","probiotic","other")
temp$pred_class = c("probiotic","other","probiotic","other")

temp$pred_class = factor(temp$pred_class, levels = rev(c("other", "probiotic")))
temp$label = factor(temp$label, levels = c("other", "probiotic"))

p <- ggplot(temp, aes(x = label, y = pred_class, fill = count)) +
  geom_tile(color = "white") +
  geom_text(aes(label = ifelse(count %% 1 == 0, as.character(count), sprintf("%.2f", count)), size = 6)) +
  scale_fill_gradient(low = "lightyellow", high = "red") +
  labs(
    x = "Truth",
    y = "Prediction",
    fill = "Count"
  ) + theme_minimal()

temp <- filter(preds, problem == "oneclass") 
temp <- temp |>
  gather(key = "pred", value = "count", -problem)

temp$label = c("other","probiotic","probiotic","other")
temp$pred_class = c("probiotic","other","probiotic","other")

temp$pred_class = factor(temp$pred_class, levels = rev(c("other", "probiotic")))
temp$label = factor(temp$label, levels = c("other", "probiotic"))

q <- ggplot(temp, aes(x = label, y = pred_class, fill = count)) +
  geom_tile(color = "white") +
  geom_text(aes(label = ifelse(count %% 1 == 0, as.character(count), sprintf("%.2f", count)), size = 6)) +
  scale_fill_gradient(low = "lightyellow", high = "red") +
  labs(
    x = "Truth",
    y = "Prediction",
    fill = "Count"
  ) + theme_minimal()

gg <- ggarrange(q,p, ncol = 2, labels = c("one", "two"), common.legend = TRUE)

fname = file.path(outdir, "confusion_matrix.png")
ggsave(filename = fname, plot = gg, device = "png")
