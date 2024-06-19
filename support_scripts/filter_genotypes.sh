#!/bin/sh

## GLOBAL PARAMS (FROM CONFIG FILE)
source config.sh

## LOCAL PARAMS
inputfile='PHENOMENA_28122022'
outdir='Analysis'
keepf="keep.fam"

MAC=100
MIND=0.05

echo "path to Plink is: $plink"

echo " - creating output directory"
mkdir -p $outdir

echo " - filtering genotype data ... "
$plink --file $inputfile --allow-extra-chr --chr 1-23 --mac $MAC --mind $MIND --geno 0.01 --bp-space 1 --snps-only 'just-acgt' --keep $keepf --make-bed --out ${outdir}/filtered

echo "DONE!"
