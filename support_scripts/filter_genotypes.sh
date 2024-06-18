#!/bin/sh

## GLOBAL PARAMS (FROM CONFIG FILE)
source config.sh

## SOFTWARE & FILES
inputfile='PHENOMENA_28122022'
outdir='Analysis'

## PARAMETERS
MAC=100
MIND=0.05

echo $plink
echo $beagle
echo $basefolder

echo " - creating output directory"
mkdir -p $outdir

echo " - filtering genotype data ... "
#$plink --file $inputfile --allow-extra-chr --chr 1-23 --mac $MAC --mind $MIND --geno 0.01 --bp-space 1 --make-bed --out ${outdir}/filtered

echo "DONE!"
