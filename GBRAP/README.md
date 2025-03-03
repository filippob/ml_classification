This tool can download the gbff files from refseq if you know the GCF number of the organism you want to download.
1. Create a .txt file with the GCF numbers of all the genomes you want to download. The numbers should not contain spaces and each GCF number should be in a new line.
2. Go to the directory containing the  tool with the command,
   
```
cd /path/to/the/folder
```
*If necessary, first set the correct permissions to excute the script using the command, 

```
chmod +x gcf_downloader.py
```

Run the tool in the command line with the code,

```
./gcf_downloader.py -in gcf_numbers.txt -out path/to/download/folder
```
