#!/bin/sh

## GLOBAL PARAMS FROM EXTERNAL FILE
source config.sh ##basefolder and paths to software

## LOCAL PARAMS
inputfile='Analysis/2.filtering/pi_filtered'
outdir='Analysis/3.imputing'

echo "path to Plink is $plink"
echo "path to Beagle is $beagle"


echo " - creating output directory"
mkdir -p $basefolder/$outdir

## prepare data
echo " - converting to VCF"
$plink --bfile $basefolder/$inputfile --recode vcf-iid --out $outdir/pi_filtered

## impute
echo " - imputing missing genotype data ... "
java -Xmx8g -jar $beagle gt=$basefolder/$outdir/pi_filtered.vcf out=$basefolder/$outdir/pi_imputed

## post-processing
echo " - convert back to plink (binary)"
$plink --vcf $basefolder/$outdir/pi_imputed.vcf.gz --make-bed --out $outdir/pi_imputed

echo "DONE!"
