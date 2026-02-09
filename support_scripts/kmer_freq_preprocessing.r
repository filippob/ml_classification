## script to clean and preprocess data for classification models

library("dplyr")
library("GGally")
library("corrplot")
library("tidyverse")
library("tidymodels")
library("data.table")

# script parameters
n_repliche <- 25
base_folder <- "~/Documents/tania/vertrebrates_invertebrates/"
infile = "genome_with_taxonomy_from_taxize.csv"
outdir = "results"


# create output folder
fname = file.path(base_folder, outdir)
if (!dir.exists(fname)) {
  dir.create(fname, recursive = TRUE, showWarnings = FALSE)
}

## input data
writeLines(" - reading input data")

fname = file.path(base_folder, infile)
inpdata <- fread(fname)

## FEATURE ENGINEERING
writeLines(" - creating new features")
# Create new average features
inpdata <- inpdata %>%
  mutate(AA_TT_avg = (`AA (%)` + `TT (%)`) / 2) %>%
  mutate(CC_GG_avg = (`CC (%)` + `GG (%)`) / 2) %>%
  mutate(TG_CA_avg = (`TG (%)` + `CA (%)`) / 2) %>%
  mutate(CT_AG_avg = (`CT (%)` + `AG (%)`) / 2) %>%
  mutate(AC_GT_avg = (`AC (%)` + `GT (%)`) / 2) %>%
  mutate(GA_TC_avg = (`GA (%)` + `TC (%)`) / 2) %>%
  select(-c(`AA (%)`, `TT (%)`, `CC (%)`,`GG (%)`, `TG (%)`, `CA (%)`, `CT (%)`, `AG (%)`,
            `AC (%)`, `GT (%)`, `GA (%)`, `TC (%)`)) %>%
  select(`Species`, `Outcome`, everything()) %>%
  mutate(Outcome = as.factor(Outcome)) |>
  select(-c(Superkingdom, Kingdom, Phylum, Class, Order, Family))

## DESCRIPTION
## average value of numerical variables
print("average values of features")
inpdata |>
  summarise_if(is.numeric, mean, na.rm = TRUE)


temp <- inpdata |> select(is.numeric)
M <- cor(temp)

fname = file.path(base_folder, outdir, "correlation_plot.png")
png(fname, width = 8, height = 8, units = "in", res = 200)
corrplot(M)
dev.off()

g <- ggpairs(temp) 
ggsave(filename = file.path(base_folder, outdir, "data_distribution.png"), plot = g, device = "png")

# ============================================================================
# 4. PREPROCESSING RECIPE
# ============================================================================
writeLines(" - preprocessing ... ")

preprocess_recipe <- recipe(Outcome ~ ., data = inpdata) %>%
  update_role(Species, new_role = "Species") %>%  # Species will not be used as predictor
  step_zv(all_numeric(), -all_outcomes()) %>%
  step_corr(all_predictors(), threshold = 0.9) %>%
  step_normalize(all_numeric(), -all_outcomes()) 

# Apply preprocessing
prep_recipe <- prep(preprocess_recipe)
print(prep_recipe)
inpdata_preprocessed <- bake(prep_recipe, new_data = NULL)

fname = file.path(base_folder, "cleaned_kmer_frequencies.csv")
fwrite(x = inpdata_preprocessed, file = fname, sep = ",", col.names = TRUE)

print("DONE!")
