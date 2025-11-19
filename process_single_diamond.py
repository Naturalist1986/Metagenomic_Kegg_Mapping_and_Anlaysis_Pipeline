#!/usr/bin/env python3
"""
Process Diamond BLASTX output against KEGG database

This script:
1. Reads Diamond BLASTX output (tsv.gz format)
2. Extracts KEGG orthology (KO) numbers from subject IDs
3. Parses KEGG database file to map gene IDs to KO numbers
4. Outputs qseqid, kegg_number, num_hits for each hit
"""

import argparse
import gzip
import sys
import re
from collections import defaultdict

def parse_kegg_dat(kegg_file):
    """
    Parse KEGG prokaryotes.dat file to create mapping of gene ID to KO number

    Expected format:
    gene_id  KO:K##### (plus other annotations)

    Returns:
        dict: {gene_id: 'K#####'}
    """
    gene_to_ko = {}

    print(f"Parsing KEGG database file: {kegg_file}", file=sys.stderr)

    try:
        with open(kegg_file, 'r') as f:
            for i, line in enumerate(f):
                if i % 1000000 == 0 and i > 0:
                    print(f"  Processed {i:,} lines, found {len(gene_to_ko):,} gene-KO mappings",
                          file=sys.stderr)

                parts = line.strip().split('\t')
                if len(parts) < 2:
                    continue

                gene_id = parts[0]
                annotation = parts[1]

                # Extract KO number from annotation (format: KO:K##### or just K#####)
                ko_match = re.search(r'K\d{5}', annotation)
                if ko_match:
                    ko_number = ko_match.group(0)
                    gene_to_ko[gene_id] = ko_number

        print(f"Completed parsing: {len(gene_to_ko):,} gene-KO mappings loaded", file=sys.stderr)

    except FileNotFoundError:
        print(f"ERROR: KEGG database file not found: {kegg_file}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"ERROR parsing KEGG database: {e}", file=sys.stderr)
        sys.exit(1)

    return gene_to_ko

def process_diamond_output(input_file, output_file, gene_to_ko):
    """
    Process Diamond BLASTX output to extract KEGG hits

    Args:
        input_file: Diamond output (tsv.gz)
        output_file: Output TSV (qseqid, kegg_number, num_hits)
        gene_to_ko: Dictionary mapping gene IDs to KO numbers
    """
    print(f"Processing Diamond output: {input_file}", file=sys.stderr)

    hit_counts = defaultdict(lambda: defaultdict(int))
    total_lines = 0
    skipped_lines = 0
    no_ko_mapping = 0

    try:
        # Open input file (handles .gz compression)
        if input_file.endswith('.gz'):
            f_in = gzip.open(input_file, 'rt')
        else:
            f_in = open(input_file, 'r')

        with f_in:
            for i, line in enumerate(f_in):
                total_lines += 1

                if i % 1000000 == 0 and i > 0:
                    print(f"  Processed {i:,} lines", file=sys.stderr)

                # Skip header line
                if line.startswith('qseqid') or line.startswith('#'):
                    skipped_lines += 1
                    continue

                parts = line.strip().split('\t')
                if len(parts) < 2:
                    skipped_lines += 1
                    continue

                qseqid = parts[0]  # Query sequence ID
                sseqid = parts[1]  # Subject (KEGG gene) ID

                # Extract gene ID from subject ID (may have prefixes like 'gnl|...')
                gene_id = sseqid.split('|')[-1] if '|' in sseqid else sseqid

                # Map gene to KO number
                if gene_id in gene_to_ko:
                    ko_number = gene_to_ko[gene_id]
                    hit_counts[qseqid][ko_number] += 1
                else:
                    # Try alternative formats
                    # Sometimes the ID might have version or other suffixes
                    gene_id_base = gene_id.split('.')[0]
                    if gene_id_base in gene_to_ko:
                        ko_number = gene_to_ko[gene_id_base]
                        hit_counts[qseqid][ko_number] += 1
                    else:
                        no_ko_mapping += 1

        print(f"Completed processing Diamond output:", file=sys.stderr)
        print(f"  Total lines: {total_lines:,}", file=sys.stderr)
        print(f"  Skipped lines: {skipped_lines:,}", file=sys.stderr)
        print(f"  Lines without KO mapping: {no_ko_mapping:,}", file=sys.stderr)
        print(f"  Unique query sequences: {len(hit_counts):,}", file=sys.stderr)

    except Exception as e:
        print(f"ERROR processing Diamond output: {e}", file=sys.stderr)
        sys.exit(1)

    # Write output
    print(f"Writing output to: {output_file}", file=sys.stderr)

    try:
        with open(output_file, 'w') as f_out:
            # No header - sum_kegg_hits.py expects no header
            total_hits = 0
            for qseqid in sorted(hit_counts.keys()):
                for ko_number in sorted(hit_counts[qseqid].keys()):
                    num_hits = hit_counts[qseqid][ko_number]
                    f_out.write(f"{qseqid}\t{ko_number}\t{num_hits}\n")
                    total_hits += 1

            print(f"Wrote {total_hits:,} KEGG hits to output file", file=sys.stderr)

    except Exception as e:
        print(f"ERROR writing output file: {e}", file=sys.stderr)
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(
        description='Process Diamond BLASTX output to extract KEGG hits'
    )
    parser.add_argument('--input', required=True,
                       help='Input Diamond BLASTX output file (tsv or tsv.gz)')
    parser.add_argument('--output', required=True,
                       help='Output file (qseqid, kegg_number, num_hits)')
    parser.add_argument('--kegg', required=True,
                       help='KEGG database file (prokaryotes.dat)')

    args = parser.parse_args()

    print("="*60, file=sys.stderr)
    print("Diamond KEGG Hit Processor", file=sys.stderr)
    print("="*60, file=sys.stderr)

    # Step 1: Parse KEGG database
    gene_to_ko = parse_kegg_dat(args.kegg)

    if not gene_to_ko:
        print("ERROR: No gene-KO mappings found in KEGG database", file=sys.stderr)
        return 1

    # Step 2: Process Diamond output
    process_diamond_output(args.input, args.output, gene_to_ko)

    print("="*60, file=sys.stderr)
    print("Processing complete!", file=sys.stderr)
    print("="*60, file=sys.stderr)

    return 0

if __name__ == '__main__':
    sys.exit(main())
