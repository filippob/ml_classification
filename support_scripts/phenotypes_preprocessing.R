library("data.table")
library("dplyr")

###################################
## read arguments from command line
###################################
allowed_parameters = c(
  'basefolder',
  'fname', # input file name
  'cols' # column names to be selected
)

args <- commandArgs(trailingOnly = TRUE)

print(args)
for (p in args){
  pieces = strsplit(p, '=')[[1]]
  #sanity check for something=somethingElse
  if (length(pieces) != 2){
    stop(paste('badly formatted parameter:', p))
  }
  if (pieces[1] %in% allowed_parameters)  {
    assign(pieces[1], pieces[2])
    next
  }
  
  #if we get here, is an unknown parameter
  stop(paste('bad parameter:', pieces[1]))
}

basefolder = "/home/tbobbo/Colombo"
fname = "data/PI_clinvar.csv"
cols = 'count,PIresINT'

print(paste("basefolder:", basefolder))
print(paste("input file name:", fname))

## read input file
input = fread(file.path(basefolder, fname))

## select columns


## filter data


## write output file
