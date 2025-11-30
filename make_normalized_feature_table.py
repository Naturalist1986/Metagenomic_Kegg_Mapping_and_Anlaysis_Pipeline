#!/usr/bin/env python3
"""
Normalize KEGG hits by genome counts

This script:
1. Reads KEGG stats file (containing run_accession and num_genomes)
2. Processes TSV files with KEGG hits
3. Normalizes hits by dividing by num_genomes
4. Outputs a normalized feature table
"""

import os
import argparse
import pandas as pd
from multiprocessing import Pool, Manager
from tqdm import tqdm

def process_tsv_file(args):
    file_name, tsv_folder, kegg_stats_data, progress_queue = args
    run_accession = file_name.replace('_kegg_hits_summed.tsv', '')  # Extract run_accession from file name

    try:
        # Check if run_accession exists in KEGG stats
        relevant_stats = kegg_stats_data[kegg_stats_data['run_accession'] == run_accession]
        if not relevant_stats.empty:
            num_genomes = relevant_stats.iloc[0]['num_genomes']

            # Read the TSV file
            tsv_path = os.path.join(tsv_folder, file_name)
            if os.stat(tsv_path).st_size == 0:  # Skip empty files
                progress_queue.put(1)
                print(f"Skipping empty file: {file_name}")
                return run_accession, {}

            tsv_data = pd.read_csv(tsv_path, sep='\t')

            # Normalize the hits
            tsv_data['normalized_hits'] = tsv_data['sum_num_hits'] / num_genomes

            # Collect normalized values
            normalized_results = {}
            for _, row in tsv_data.iterrows():
                kegg_number = row['kegg_number']
                normalized_value = row['normalized_hits']

                normalized_results[kegg_number] = normalized_value

            progress_queue.put(1)
            return run_accession, normalized_results
    except Exception as e:
        progress_queue.put(1)
        print(f"Error processing file {file_name}: {e}")
    return run_accession, {}

def normalize_kegg_hits(tsv_folder, kegg_stats_path, output_file):
    # Read KEGG stats Excel file
    kegg_stats_data = pd.read_excel(kegg_stats_path)

    # Clean column names to avoid KeyError due to mismatched column names
    kegg_stats_data.columns = kegg_stats_data.columns.str.strip().str.lower()

    # Ensure 'run_accession' and 'num_genomes' are present
    if 'run_accession' not in kegg_stats_data.columns or 'num_genomes' not in kegg_stats_data.columns:
        raise KeyError("The KEGG stats file must contain 'run_accession' and 'num_genomes' columns.")

    # Validate that output_file is not a directory
    if os.path.isdir(output_file):
        output_file = os.path.join(output_file, "normalized_kegg_results.tsv")  # Create default file name in the directory

    # List all TSV files in the folder
    tsv_files = [file_name for file_name in os.listdir(tsv_folder) if file_name.endswith(".tsv")]

    # Use multiprocessing to process files in parallel with progress indication
    with Manager() as manager:
        progress_queue = manager.Queue()
        with Pool(processes=4) as pool:
            # Start a progress bar in the main process
            with tqdm(total=len(tsv_files), desc="Processing TSV files", unit="file") as pbar:
                def update_progress():
                    while True:
                        progress_queue.get()
                        pbar.update()
                        if pbar.n == len(tsv_files):
                            break

                results = pool.map_async(
                    process_tsv_file,
                    [(file_name, tsv_folder, kegg_stats_data, progress_queue) for file_name in tsv_files]
                )

                update_progress()
                results = results.get()

    # Combine results from all processes
    combined_results = {}
    all_kegg_numbers = set()
    for run_accession, result in results:
        if not result:  # Skip empty results
            continue
        for kegg_number, normalized_value in result.items():
            if kegg_number not in combined_results:
                combined_results[kegg_number] = {}
            combined_results[kegg_number][run_accession] = normalized_value
            all_kegg_numbers.add(kegg_number)

    # Ensure all KEGG numbers are accounted for in all run_accessions
    all_run_accessions = [file_name.replace('_kegg_hits_summed.tsv', '') for file_name in tsv_files]
    for kegg_number in all_kegg_numbers:
        for run_accession in all_run_accessions:
            if run_accession not in combined_results[kegg_number]:
                combined_results[kegg_number][run_accession] = 0

    # Create a DataFrame with KEGG numbers as rows and run_accessions as columns
    result_df = pd.DataFrame.from_dict(combined_results, orient='index')
    result_df.index.name = 'kegg_number'
    result_df.reset_index(inplace=True)

    # Ensure the columns are sorted by run_accession for consistency
    columns_order = ['kegg_number'] + sorted(all_run_accessions)
    result_df = result_df[columns_order]

    # Save the result to the output file
    result_df.to_csv(output_file, sep='\t', index=False, float_format='%.6f')

def main():
    """Main function to process command-line arguments and run normalization"""
    parser = argparse.ArgumentParser(
        description='Normalize KEGG hits by genome counts from stats file'
    )
    parser.add_argument('--input-dir', required=True,
                       help='Directory containing _kegg_hits_summed.tsv files')
    parser.add_argument('--kegg-stats', required=True,
                       help='Excel file with run_accession and num_genomes columns')
    parser.add_argument('--output', required=True,
                       help='Output TSV file for normalized results')

    args = parser.parse_args()

    # Validate input directory
    if not os.path.isdir(args.input_dir):
        print(f"ERROR: Input directory does not exist: {args.input_dir}")
        return 1

    # Validate KEGG stats file
    if not os.path.isfile(args.kegg_stats):
        print(f"ERROR: KEGG stats file does not exist: {args.kegg_stats}")
        return 1

    # Call the normalization function
    try:
        normalize_kegg_hits(args.input_dir, args.kegg_stats, args.output)
        print(f"\nNormalization complete! Results saved to: {args.output}")
        return 0
    except Exception as e:
        print(f"ERROR: {e}")
        return 1

if __name__ == "__main__":
    exit(main())
