
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

D <- bestmodel |>
  group_by(problem, hyperparameter, value) |>
  summarise(N = n()) |>
  arrange(value) |>
  mutate(value = as.character(round(value, 3)))

p <- ggplot(filter(D, problem == "oneclass"), aes(x = value, y = N))
p <- p + geom_bar(aes(fill=hyperparameter), stat="identity")
p <- p + facet_wrap(~hyperparameter, scales = "free")
p

q <- ggplot(filter(D, problem == "twoclass"), aes(x = value, y = N))
q <- q + geom_bar(aes(fill=hyperparameter), stat="identity")
q <- q + facet_wrap(~hyperparameter, scales = "free")
q

gg <- ggarrange(p,q, labels = c("one","two"), ncol = 1)
