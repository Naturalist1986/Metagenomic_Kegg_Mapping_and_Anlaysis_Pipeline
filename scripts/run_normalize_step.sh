#!/bin/bash
#SBATCH --job-name=normalize_kegg
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=2:00:00
##SBATCH --account=ACCOUNT_NAME

################################################################################
# Normalization Step - Create Normalized Feature Table
#
# This script normalizes KEGG hits by genome equivalents to create
# the final feature table for downstream analysis
################################################################################

set -euo pipefail

# Function for logging
log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

# Load configuration
if [ $# -eq 0 ]; then
    log_message "ERROR: No configuration file provided"
    exit 1
fi

CONFIG_FILE=$1
source "$CONFIG_FILE"

# Set up directories
SUM_KEGG_DIR="${OUTPUT_BASE_DIR}/kegg_summed"
FINAL_OUTPUT_DIR="${OUTPUT_BASE_DIR}/final_results"

mkdir -p "$FINAL_OUTPUT_DIR"

log_message "=========================================="
log_message "Creating Normalized Feature Table"
log_message "=========================================="

# Check if input directory exists
if [ ! -d "$SUM_KEGG_DIR" ]; then
    log_message "ERROR: Input directory not found: $SUM_KEGG_DIR"
    exit 1
fi

# Check if kegg_stats.xlsx exists
KEGG_STATS_FILE="${SUM_KEGG_DIR}/kegg_stats.xlsx"
if [ ! -f "$KEGG_STATS_FILE" ]; then
    log_message "ERROR: KEGG stats file not found: $KEGG_STATS_FILE"
    log_message "Please run the mean single-copy step first"
    exit 1
fi

OUTPUT_FILE="${FINAL_OUTPUT_DIR}/normalized_kegg_feature_table.tsv"

# Check if output already exists
if [ -f "$OUTPUT_FILE" ]; then
    log_message "Output file already exists: $OUTPUT_FILE"
    log_message "Removing and regenerating..."
    rm -f "$OUTPUT_FILE"
fi

# Locate the normalization script
NORMALIZE_SCRIPT="${PIPELINE_DIR}/make_normalized_feature_table.py"

if [ ! -f "$NORMALIZE_SCRIPT" ]; then
    log_message "ERROR: Python script not found: $NORMALIZE_SCRIPT"
    exit 1
fi

log_message "Input directory: $SUM_KEGG_DIR"
log_message "KEGG stats file: $KEGG_STATS_FILE"
log_message "Output file: $OUTPUT_FILE"
log_message ""
log_message "Running normalization..."

# Run the script with command-line arguments
python3 "$NORMALIZE_SCRIPT" \
    --input-dir "$SUM_KEGG_DIR" \
    --kegg-stats "$KEGG_STATS_FILE" \
    --output "$OUTPUT_FILE"

if [ $? -ne 0 ]; then
    log_message "ERROR: Normalization failed"
    exit 1
fi

# Verify output exists
if [ -f "$OUTPUT_FILE" ] && [ -s "$OUTPUT_FILE" ]; then
    OUTPUT_SIZE=$(du -h "$OUTPUT_FILE" | cut -f1)
    NUM_ROWS=$(wc -l < "$OUTPUT_FILE")
    NUM_KEGG=$(( NUM_ROWS - 1 ))  # Subtract header
    NUM_SAMPLES=$(head -1 "$OUTPUT_FILE" | awk -F'\t' '{print NF-1}')  # Subtract KEGG column

    log_message ""
    log_message "=========================================="
    log_message "SUCCESS!"
    log_message "=========================================="
    log_message "Output file: $OUTPUT_FILE"
    log_message "File size: $OUTPUT_SIZE"
    log_message "Number of KEGG features: $NUM_KEGG"
    log_message "Number of samples: $NUM_SAMPLES"
    log_message "=========================================="
else
    log_message "ERROR: Output file was not created or is empty"
    exit 1
fi

log_message ""
log_message "=========================================="
log_message "PIPELINE COMPLETE!"
log_message "=========================================="
log_message "Final results are in: $FINAL_OUTPUT_DIR"
log_message ""
log_message "Next steps:"
log_message "  1. Review the normalized feature table: $OUTPUT_FILE"
log_message "  2. Check SingleM results: ${OUTPUT_BASE_DIR}/singlem/"
log_message "  3. Review treatment group summaries: ${OUTPUT_BASE_DIR}/singlem/combined_by_treatment/"
log_message "=========================================="

exit 0
