#!/bin/sh

## GLOBAL PARAMS FROM EXTERNAL FILE
source config.sh

## LOCAL PARAMS
inputfile='filtered'
outdir='Analysis'

echo "path to Plink is $plink"
echo "path to Beagle is $beagle"


echo " - creating output directory"
mkdir -p $outdir

## prepare data
echo " - converting to VCF"
$plink --bfile $inputfile --recode vcf-iid --out $outdir/filtered

## impute
echo " - imputing missing genotype data ... "
java -Xmx8g -jar $beagle gt=$outdir/filtered.vcf out=$outdir/imputed

## post-processing
echo " - convert back to plink (binary)"
$plink --vcf $outdir/imputed.vcf.gz --make-bed --out $outdir/imputed

echo "DONE!"
