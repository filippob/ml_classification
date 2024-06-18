#!/bin/sh

## GLOBAL PARAMS FROM EXTERNAL FILE
source config.sh

## LOCAL PARAMS
inputfile='imputed'
outdir='Analysis'

echo "path to Plink is $plink"

echo " - creating output directory"
mkdir -p $outdir

## make kinship matrix
echo " - calculating the kinship matrix ... "
$plink --bfile imputed --make-rel square gz out=$outdir/kinship

echo "DONE!"
