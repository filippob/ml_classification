## GFC Downloader
   
The gcf_downloader tool can download the gbff files from refseq if you know the GCF number of the organism you want to download. If a gbff file for the GCF number already exists in the output folder, the tool will not download it again. 

### Usage instructions

1. Create a .txt file with the GCF numbers of all the genomes you want to download. The numbers should not contain spaces and each GCF number should be in a new line.
2. Go to the directory containing the  tool using the command,
   
```
cd /path/to/the/folder
```
*If necessary, set the correct permissions to excute the script using the command, 

```
chmod +x gcf_downloader.py
```

3. Run the tool in the command line using the command,

```
./gcf_downloader.py -in gcf_numbers.txt -out path/to/download/folder
```

## NC to GCF

The nc_2_gcf tool can retrieve the NCBI assembly accession (GCF) corresponding to a list of NC accession numbers. It queries the NCBI Entrez API and outputs the results into a text file.

### Usage Instructions



1. Create a CSV file containing NC accession numbers, where each number is on a new line and without spaces.
2. Go to the directory containing the tool using the command:
```
cd /path/to/the/folder
```
*If necessary, set the correct permissions to execute the script using the command:
```
chmod +x ncbi_accession_tool.py
```
3. Run the tool in the command line using the command:
```
./nc_2_gcf.py -in input_nc_numbers.csv -out gcf_numbers.txt -email your_email@example.com
```
Input your email adress to comply with NCBI's API requirements.

#### Note
The tool includes a 0.5-second delay between each request to NCBI to avoid rate limits. Please increase the delay in the script if required.
If no corresponding GCF accession is found for an NC accession, the script will notify and skip it in the output list.
