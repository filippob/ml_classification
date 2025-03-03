This code can download the gbff files from refseq if you know the GCF number of the organism you want to download.
1. Create a .txt file with the GCF numbers of all the genomes you want to download. The numbers should not contain spaces and each GCF number in in a new line.
2. run the tool in the command line with the code,
```
   ./gfc_downloader.py -in gcf_numbers.txt -out path/to/download/folder
```
