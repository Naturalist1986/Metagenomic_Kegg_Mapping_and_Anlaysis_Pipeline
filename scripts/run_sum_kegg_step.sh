#!/bin/bash
#SBATCH --job-name=sum_kegg
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=1:00:00

################################################################################
# Sum KEGG Hits Step
#
# This script sums KEGG hits for each sample
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

# Calculate array index
OFFSET=${ARRAY_OFFSET:-0}
LINE_NUM=$((SLURM_ARRAY_TASK_ID + OFFSET))

# Set up directories
POST_DIAMOND_DIR="${OUTPUT_BASE_DIR}/post_diamond_processing"
SUM_KEGG_DIR="${OUTPUT_BASE_DIR}/kegg_summed"

mkdir -p "$SUM_KEGG_DIR"

# Get list of files to process
shopt -s nullglob
FILES=( ${POST_DIAMOND_DIR}/*_kegg_hits.tsv )
shopt -u nullglob

if [ ${#FILES[@]} -eq 0 ]; then
    log_message "ERROR: No KEGG hits files found in $POST_DIAMOND_DIR"
    exit 1
fi

if [ $LINE_NUM -ge ${#FILES[@]} ]; then
    log_message "ERROR: Array index $LINE_NUM exceeds number of files (${#FILES[@]})"
    exit 1
fi

# Get file to process
FILE_TO_PROCESS="${FILES[$LINE_NUM]}"
BASENAME=$(basename "$FILE_TO_PROCESS" _kegg_hits.tsv)
OUTPUT_FILE="${SUM_KEGG_DIR}/${BASENAME}_kegg_hits_summed.tsv"

log_message "=========================================="
log_message "Processing: $FILE_TO_PROCESS"
log_message "Output: $OUTPUT_FILE"
log_message "=========================================="

# Skip if output already exists
if [ -f "$OUTPUT_FILE" ] && [ -s "$OUTPUT_FILE" ]; then
    log_message "Output file already exists. Skipping."
    exit 0
fi

# Run the Python script
SUM_SCRIPT="${PIPELINE_DIR}/sum_kegg_hits.py"

if [ ! -f "$SUM_SCRIPT" ]; then
    log_message "ERROR: Python script not found: $SUM_SCRIPT"
    exit 1
fi

python "$SUM_SCRIPT" --input "$FILE_TO_PROCESS" --output "$OUTPUT_FILE"

if [ $? -ne 0 ]; then
    log_message "ERROR: Processing file $FILE_TO_PROCESS failed"
    exit 1
fi

log_message "Successfully summed KEGG hits for $BASENAME"

exit 0
