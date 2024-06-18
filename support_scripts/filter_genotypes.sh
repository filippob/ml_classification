#!/bin/sh

## SOFTWARE & FILES
plink=/home/filippo/Downloads/plink
inputfile='PHENOMENA_28122022'
outdir='Analysis'

## PARAMETERS
MAC=100
MIND=0.05

echo " - creating output directory"
mkdir -p $outdir

echo " - filtering genotype data ... "
$plink --file $inputfile --allow-extra-chr --chr 1-23 --mac $MAC --mind $MIND --geno 0.01 --bp-space 1 --recode --out ${outdir}/filtered

echo "DONE!"
