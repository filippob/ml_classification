
## script to make a test set

library("dplyr")
library("readxl")
library("data.table")

## Parameters
basefolder <- "/home/filippo/Documents/tania/probiotics"
input_file <- "data/filtered_merged_bits26_testset.csv"
outdir = "splits"
target_var = NULL
nclass_minority = NULL
nclass_majority = NULL
to_remove = c("Assembly", "Class", "Locus_ID", "Version") ## columns to remove


## Import dataset
writeLines(" - reading the data ...")
fname <- file.path(basefolder, input_file)
ext = gsub("^.*\\.","",fname)
if (ext %in% c("xlsx")) {
  
  dataset <- readxl::read_xlsx(fname)
} else dataset <- fread(fname)

print(paste("Data size (n. records):",nrow(dataset)))

## Data preparation
## making target variable a factor

if (!is.null(target_var)) {
  
  print("making the target variable a factor")
  dataset <- dataset %>% 
    mutate({{target_var}} := as.factor(.data[[target_var]]))
}


## removing unnecessary columns
dataset <- dataset |> select(!all_of(to_remove))

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

## remove missing
dataset <- na.omit(dataset) ## since we have two probiotics with missing data (potential duplicates): TO  CHECK !!!
if (!is.null(target_var)) {
  dataset |>
    group_by(Label) |>
    summarise(N = n())
}

## training/test split
writeLines(" - splitting the data in training and test sets")

if (!is.null(nclass_minority) & !is.null(nclass_majority)) {
  
  print(paste("N. of minority class records to be sampled for test:", nclass_minority))
  print(paste("N. of majority class records to be sampled for test:", nclass_majority))
  
  probiotics <- subset(dataset, Label == "Probiotic")
  nonprobiotics <- subset(dataset, Label == "Nonprobiotic")
  
  vec <- sample(1:nrow(probiotics), nclass_minority)
  test_probs = probiotics[vec,]
  train_probs = probiotics[-vec,]
  
  vec <- sample(1:nrow(nonprobiotics), nclass_majority)
  test_nonprobs = nonprobiotics[vec,]
  train_nonprobs = nonprobiotics[-vec,]
  
  test <- test_probs |> bind_rows(test_nonprobs)
  train <- train_probs |> bind_rows(train_nonprobs)
  
  writeLines(" - writing out test and training sets")
  fname = file.path(basefolder, outdir, "test_set.csv")
  fwrite(x = test, file = fname)
  
  fname = file.path(basefolder, outdir, "train_set.csv")
  fwrite(x = train, file = fname)
} else {
  
  writeLines(" - writing out independent test set")
  fname = file.path(basefolder, outdir, basename(input_file))
  fwrite(x = dataset, file = fname)
}


print("DONE!!")
