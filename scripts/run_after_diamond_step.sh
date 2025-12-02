#!/bin/bash
#SBATCH --job-name=after_diamond
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=1-00:00:00
##SBATCH --account=ACCOUNT_NAME

################################################################################
# Post-Diamond KEGG Hit Processing Step
#
# This script processes Diamond BLASTX output to extract KEGG hits
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
DIAMOND_OUTDIR="${OUTPUT_BASE_DIR}/diamond_kegg"
POST_DIAMOND_DIR="${OUTPUT_BASE_DIR}/post_diamond_processing"

mkdir -p "$POST_DIAMOND_DIR"

# Get list of files to process
shopt -s nullglob
FILES=( ${DIAMOND_OUTDIR}/*.tsv.gz )
shopt -u nullglob

if [ ${#FILES[@]} -eq 0 ]; then
    log_message "ERROR: No .tsv.gz files found in $DIAMOND_OUTDIR"
    exit 1
fi

if [ $LINE_NUM -ge ${#FILES[@]} ]; then
    log_message "ERROR: Array index $LINE_NUM exceeds number of files (${#FILES[@]})"
    exit 1
fi

# Get file to process
FILE_TO_PROCESS="${FILES[$LINE_NUM]}"
BASENAME=$(basename "$FILE_TO_PROCESS" .tsv.gz)
OUTPUT_FILE="${POST_DIAMOND_DIR}/${BASENAME}_kegg_hits.tsv"
TEMP_OUTPUT="${OUTPUT_FILE}.tmp"
LOCK_FILE="${OUTPUT_FILE}.lock"

log_message "=========================================="
log_message "Processing: $FILE_TO_PROCESS"
log_message "Output: $OUTPUT_FILE"
log_message "=========================================="

# Check if output already exists
if [ -f "$OUTPUT_FILE" ] && [ -s "$OUTPUT_FILE" ]; then
    log_message "Output file already exists and has content. Skipping."
    exit 0
fi

# Use file locking to prevent concurrent processing
(
    flock -n 200 || {
        log_message "Another job is already processing this file. Exiting."
        exit 0
    }

    # Clean up any previous temporary files
    rm -f "$TEMP_OUTPUT"

    log_message "Memory available: $(free -h | grep Mem | awk '{print $7}')"
    log_message "Starting processing with $SLURM_CPUS_PER_TASK CPUs"

    # Check if process_single_diamond.py exists in the pipeline directory
    PROCESS_SCRIPT="${PIPELINE_DIR}/process_single_diamond.py"

    if [ ! -f "$PROCESS_SCRIPT" ]; then
        log_message "ERROR: Python script not found: $PROCESS_SCRIPT"
        exit 1
    fi

    # Run the processing script
    if python "$PROCESS_SCRIPT" \
        --input "$FILE_TO_PROCESS" \
        --output "$TEMP_OUTPUT" \
        --kegg "$KEGG_DAT_FILE"; then

        # Move temporary file to final location only if successful
        mv "$TEMP_OUTPUT" "$OUTPUT_FILE"
        log_message "Successfully processed $BASENAME"

        # Validate output file
        if [ ! -s "$OUTPUT_FILE" ]; then
            log_message "WARNING: Output file is empty after processing"
            exit 1
        fi

        # Log output file size
        OUTPUT_SIZE=$(du -h "$OUTPUT_FILE" | cut -f1)
        log_message "Output file size: $OUTPUT_SIZE"
    else
        EXIT_CODE=$?
        log_message "ERROR: Processing failed with exit code $EXIT_CODE"
        rm -f "$TEMP_OUTPUT"
        exit $EXIT_CODE
    fi

) 200>"$LOCK_FILE"

# Clean up lock file
rm -f "$LOCK_FILE"

log_message "Job completed successfully"

exit 0
