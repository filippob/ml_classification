#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 11 16:32:34 2025

@author: sach_s__macbook
"""

import csv
import time
import argparse
from Bio import Entrez

def get_assembly_accession(nc_accession):
    """Retrieve the NCBI assembly accession (GCF) corresponding to an NC accession number."""
    try:
        # First, find the nucleotide record to ensure it exists
        handle = Entrez.esearch(db="nucleotide", term=nc_accession, retmax=1)
        record = Entrez.read(handle)
        handle.close()

        if not record["IdList"]:
            return None
        
        # Then, search for the corresponding assembly
        handle = Entrez.esearch(db="assembly", term=nc_accession, retmax=1)
        record = Entrez.read(handle)
        handle.close()

        if not record["IdList"]:
            print(f"No assembly found for {nc_accession}.")
            return None

        assembly_id = record["IdList"][0]
        
        # Retrieve assembly summary
        handle = Entrez.esummary(db="assembly", id=assembly_id)
        summary = Entrez.read(handle)
        handle.close()

        #print (summary) #un-comment to get a summary of the info inside Entrez API
        
        # Extract the Assembly Acession number + the name from the FTP path
        # (The code retrives the ftp path directly instead of the Accession number because there are different accession versions 
        #it is better to get the aceesion version included in the ftp path.)
        
        if "DocumentSummarySet" in summary and "DocumentSummary" in summary["DocumentSummarySet"]:
            doc_summary = summary["DocumentSummarySet"]["DocumentSummary"]
            if doc_summary and "FtpPath_RefSeq" in doc_summary[0]:
                ftp_path = doc_summary[0]["FtpPath_RefSeq"]
                # Split by '/' and return the last part which contains the Assembly accession_number+name
                return ftp_path.split("/")[-1]

        print(f"Assembly accession not found in summary for {nc_accession}.")
        return None

    except Exception as e:
        print(f"Error retrieving {nc_accession}: {e}")
        return None

def process_nc_list(input_file, output_file):
    """Process a CSV file containing NC accession numbers and write the corresponding GCF assembly accessions."""
    with open(input_file, newline='') as csvfile:
        reader = csv.reader(csvfile)
        
        with open(output_file, 'w') as textfile:
            for row in reader:
                nc_accession = row[0].strip()  # Remove extra spaces or invisible characters
                assembly_accession = get_assembly_accession(nc_accession)
                
                if assembly_accession:
                    textfile.write(f"{assembly_accession}\n")
                    print(f"{nc_accession} -> {assembly_accession}")
                else:
                    print(f"No GCF found for {nc_accession}. Skipping entry.")

                time.sleep(0.5)  # Delay to avoid NCBI API rate limits

def main():
    parser = argparse.ArgumentParser(description="Retrieve NCBI assembly accession numbers from NC accession numbers.")
    parser.add_argument("-in", "--input", required=True, help="Input CSV file containing NC accession numbers.")
    parser.add_argument("-out", "--output", required=True, help="Output text file to save GCF accession numbers.")
    parser.add_argument("-email", required=True, help="Email address for NCBI Entrez API.")
    
    args = parser.parse_args()
    
    # Set email for Entrez API
    Entrez.email = args.email
    
    process_nc_list(args.input, args.output)

if __name__ == "__main__":
    main()