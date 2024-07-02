#!/bin/sh

## GLOBAL PARAMS (FROM CONFIG FILE)
source config.sh

## LOCAL PARAMS
prjfolder="$HOME/Colombo"
inputfile='data/genotyped_nofilters/PHENOMENA_28122022.ped'
outdir='Analysis/1.cleaning'


## Rimuovi i controlli (HG00264 & HG00097) con record duplicati nei genotipi (assenti nel dataset dei fenotipi) per poter usare plink: 5 record con family ID 0 e stesso individual ID (HG00264 o HG00097) e 24 record ciascuno HG00XXX_n. (da 1 a 24) con individual e family ID uguali

fname=$prjfolder/$inputfile

echo "N. of samples in the original file"
wc -l $fname

#grep 'HG00264' genotyped_nofilters/PHENOMENA_28122022.ped | less -S
grep -v 'HG00264' $fname > temp.ped
#grep 'HG00097' temp.ped | less -S

outfile="$prjfolder/$outdir/pi_cleaned.ped"
grep -v 'HG00097' temp.ped > $outfile

echo "N. of samples after cleaning"
wc -l $outfile

echo "copying and renaming the .map file"
basename="${inputfile%.*}"
mapfile="${prjfolder}/{$basename}.map"

cp $mapfile "$prjfolder/$outdir/pi_cleaned.map"

echo "DONE!"
