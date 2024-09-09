#!/bin/sh

## GLOBAL PARAMS FROM EXTERNAL FILE
source config.sh

## LOCAL PARAMS
inputfile='Analysis/3.imputing/pi_imputed'
outdir='Analysis/4.covariance_matrix'

echo "path to Plink is $plink"

echo " - creating output directory"
mkdir -p $outdir

## make kinship matrix
echo " - calculating the kinship matrix ... "
$plink --bfile $inputfile --maf 0.01 --make-rel square gz --out=$outdir/pi_kinship

echo "DONE!"
