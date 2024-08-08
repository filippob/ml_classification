#!/bin/sh

## GLOBAL PARAMS (FROM CONFIG FILE)
source config.sh ## basefolder and paths to software from here

## LOCAL PARAMS
inputfile='Analysis/1.cleaning/pi_cleaned'
outdir='Analysis/2.filtering'
keepf="data/keep.fam"

MAF=0.01
GENO=0.01
HWE=1e-06

echo "path to Plink is: $plink"

echo " - creating output directory"
mkdir -p $basefolder/$outdir

echo " - filtering genotype data ... "
$plink --file $basefolder/$inputfile --allow-extra-chr --chr 1-22 --geno $GENO --maf $MAF --hwe $HWE --bp-space 1 --snps-only 'just-acgt' --keep $keepf --make-bed --out ${outdir}/pi_filtered

echo "DONE!"
