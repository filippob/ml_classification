#!/usr/bin/env python3

"""
Created on Mon Mar  3 17:10:54 2025

@author: sach_s__macbook
"""

'''This tool downloads the gbff files from refseq when a .txt file with the GCF numbers of the genome is provided'''

import os
import re
import requests
import gzip
import shutil
import argparse
from pathlib import Path
from urllib.parse import urljoin

def extract_triplets(gcf_number):
    match = re.search(r'GCF_(\d+)\.', gcf_number)
    if not match:
        print(f"Invalid GCF number format: {gcf_number}")
        return None
    digits = match.group(1)
    triplets = [digits[i:i+3] for i in range(0, len(digits), 3)]
    return '/'.join(triplets)

def download_genome(gcf_number, output_dir="downloads"):
    triplet_path = extract_triplets(gcf_number)
    if not triplet_path:
        return
    
    base_url = f"https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/{triplet_path}/{gcf_number}/"
    genome_file = f"{gcf_number}_genomic.gbff"
    genome_gz_file = genome_file + ".gz"

    
    output_gbff_path = Path(output_dir) / genome_file  # Uncompressed file path
    output_gz_path = Path(output_dir) / genome_gz_file  # Compressed file path
    
    # Check if the file already exists
    if output_gz_path.exists() or output_gbff_path.exists():
        print(f"gbff file already exists for {gcf_number} in the output directory, skipping download")
        return
    
    
    genome_url = urljoin(base_url, genome_gz_file)

    response = requests.head(genome_url)
    if response.status_code != 200:
        print(f"Genome file not found for {gcf_number}")
        return
    
    os.makedirs(output_dir, exist_ok=True)
   
    
    with requests.get(genome_url, stream=True) as r:
        r.raise_for_status()
        with open(output_gz_path, 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
    
    print(f"Downloaded: {output_gz_path}")
    
    # Extract the .gz file
    with gzip.open(output_gz_path, 'rb') as f_in:
        with open(output_gbff_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    
    os.remove(output_gz_path)  # Remove the .gz file after extraction
    print(f"Extracted: {output_gbff_path}")

def main():
    # Set up command-line argument parsing
    parser = argparse.ArgumentParser(description="Download and extract genomes based on GCF numbers.")
    parser.add_argument("-in", "--input_file", required=True, help="Path to the .txt file containing the GCF numbers (one per line).")
    parser.add_argument("-out", "--output_dir", required=True, help="Directory where the downloaded genomes will be stored.")
    
    args = parser.parse_args()
    
    # Read the GCF numbers from the input file
    with open(args.input_file, 'r') as file:
        gcf_list = [line.strip() for line in file.readlines() if line.strip()]
    
    # Download genomes for each GCF number in the list
    for gcf in gcf_list:
        download_genome(gcf, args.output_dir)

if __name__ == "__main__":
    main()
