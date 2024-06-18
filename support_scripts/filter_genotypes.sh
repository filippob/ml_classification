#!/bin/sh

plink=/home/filippo/Downloads/plink
inputfile='PHENOMENA_28122022'
MAC=100
MIND=0.05

echo " - filtering genotype data ... "
$plink --file $inputfile --allow-extra-chr --chr 1-23 --mac $MAC --mind $MIND --geno 0.01 --bp-space 0 --recode --out filtered

echo "DONE!"
