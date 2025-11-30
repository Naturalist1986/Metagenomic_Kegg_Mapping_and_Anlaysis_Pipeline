#!/usr/bin/env python3
"""
Calculate mean single-copy gene counts from KEGG-annotated TSV files

This script:
1. Reads summed KEGG hit TSV files from a directory
2. Filters for single-copy KEGG genes
3. Calculates the mean count for each file (representing estimated genome count)
4. Outputs results to an Excel file
"""

import os
import argparse
import pandas as pd
import numpy as np
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm  # For progress indication

# Single-copy KEGG genes of interest
kegg_numbers = [
    "K09748", "K03687", "K00962", "K02864", "K02994", "K02996", "K03438", "K02835",
    "K02836", "K02968", "K07042", "K11749", "K02879", "K02888", "K03110", "K03531",
    "K02834", "K11753", "K03550", "K03664", "K15429", "K09710", "K08316",
    "K03545", "K02357", "K03501", "K02939", "K02990", "K02520", "K03218", "K03703",
    "K07447", "K01937", "K01872", "K02313", "K03544", "K01870", "K01869", "K01874",
    "K01875", "K04485", "K03177", "K01883", "K03595", "K01892", "K01000",
    "K01887", "K01876", "K00604", "K01889", "K01890", "K02519", "K02838", "K02495",
    "K03723", "K06187", "K03702", "K03631", "K03551", "K03655", "K02338",
    "K02945", "K02528", "K03075", "K02601", "K01756", "K03106", "K03070", "K03073",
    "K03076", "K02988", "K02992", "K02887", "K02112",
    "K02890", "K02470", "K02469", "K02871", "K02876", "K02895", "K01924", "K01925",
    "K02340", "K02878", "K02863", "K02886", "K00088", "K02316", "K03596",
    "K06207", "K03625", "K02600", "K03553", "K03043", "K03040",
    "K03685", "K02860", "K03046", "K02343", "K04075", "K03979",
    "K00942", "K03977", "K02906", "K02948",  "K25706", "K14742", "K02926"
]

def process_file(file_path):
    """For each file, calculate the mean of the summed KEGG gene counts distribution."""
    try:
        # Extract run_accession from filename (remove the suffix)
        filename = os.path.basename(file_path)
        run_accession = filename.replace('_kegg_hits_summed.tsv', '')
        
        # 1. Skip empty files
        if os.stat(file_path).st_size == 0:
            print(f"Skipping empty file: {file_path}")
            return {
                'run_accession': run_accession,
                'num_genomes': None
            }
        
        # 2. Read the file
        data = pd.read_csv(file_path, sep='\t', header=None, names=["KO", "Count"])
        
        # 3. Ensure "Count" column is numeric and drop invalid values
        data["Count"] = pd.to_numeric(data["Count"], errors="coerce")
        data = data.dropna()  # Remove rows with NaN values
        
        # 4. Filter rows with single-copy KEGG numbers
        filtered_data = data[data["KO"].isin(kegg_numbers)]
        if filtered_data.empty:
            print(f"No matching single-copy KEGGs in file: {file_path}")
            return {
                'run_accession': run_accession,
                'num_genomes': None
            }
        
        # 5. Compute sum of hits per KEGG gene
        sum_hits = filtered_data.groupby("KO")["Count"].sum()
        sum_hits = pd.to_numeric(sum_hits, errors="coerce").dropna()
        
        # 6. Calculate the mean of the distribution
        mean_value = sum_hits.mean()
        
        # 7. Return results
        return {
            'run_accession': run_accession,
            'num_genomes': mean_value
        }
    
    except Exception as e:
        print(f"Error processing file {file_path}: {e}")
        filename = os.path.basename(file_path)
        run_accession = filename.replace('_kegg_hits_summed.tsv', '')
        return {
            'run_accession': run_accession,
            'num_genomes': None
        }

def main():
    """Main function to process command-line arguments and run the pipeline"""
    parser = argparse.ArgumentParser(
        description='Calculate mean single-copy gene counts from KEGG-annotated TSV files'
    )
    parser.add_argument('--input-dir', required=True,
                       help='Directory containing _kegg_hits_summed.tsv files')
    parser.add_argument('--output', required=True,
                       help='Output Excel file path (e.g., kegg_stats.xlsx)')

    args = parser.parse_args()

    input_dir = args.input_dir
    output_file = args.output

    # Validate input directory
    if not os.path.isdir(input_dir):
        print(f"ERROR: Input directory does not exist: {input_dir}")
        return 1

    # Collect files
    files = [
        os.path.join(input_dir, f)
        for f in os.listdir(input_dir)
        if f.endswith('_kegg_hits_summed.tsv')
    ]

    if not files:
        print(f"WARNING: No files matching pattern '*_kegg_hits_summed.tsv' found in {input_dir}")
        return 1

    print(f"Found {len(files)} files to process")

    # Process in parallel
    results = []
    with ProcessPoolExecutor() as executor:
        for result in tqdm(
            executor.map(process_file, files),
            total=len(files),
            desc="Processing Files"
        ):
            results.append(result)

    # Convert to a DataFrame and save
    results_df = pd.DataFrame(results)
    results_df.to_excel(output_file, index=False)
    print(f"\nProcessed {len(results_df)} files")
    print(f"Results saved to: {output_file}")

    return 0

if __name__ == "__main__":
    exit(main())
