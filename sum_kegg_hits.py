#!/usr/bin/env python

import argparse
import pandas as pd

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Sum KEGG hits from a Diamond KEGG hits file')
    parser.add_argument('--input', required=True, help='Input Diamond KEGG hits file')
    parser.add_argument('--output', required=True, help='Output file for summed KEGG hits')
    
    args = parser.parse_args()
    
    try:
        # Read the input file with column names
        # Expected format: qseqid, kegg_number, num_hits
        # Set low_memory=False to avoid DtypeWarning
        data = pd.read_csv(args.input, sep='\t', names=['qseqid', 'kegg_number', 'num_hits'], 
                          low_memory=False)
        
        # Make sure num_hits is converted to float before summing
        data['num_hits'] = pd.to_numeric(data['num_hits'], errors='coerce')
        
        # Drop any rows where conversion failed (NaN values)
        data = data.dropna(subset=['num_hits'])
        
        # Sum the num_hits values for each kegg_number
        ko_sums = data.groupby('kegg_number')['num_hits'].sum().reset_index()
        
        # Rename columns for output to match expected format
        ko_sums.columns = ['kegg_number', 'sum_num_hits']
        
        # Save to output file WITH HEADERS
        ko_sums.to_csv(args.output, sep='\t', index=False, header=True)
        print(f"Summed KEGG hits saved to {args.output}")
        
    except Exception as e:
        print(f"Error processing file: {e}")
        return 1
    
    return 0

if __name__ == "__main__":
    exit(main())